//#include <stdlib.h>
#include "Define.h"
#include "BasisFunction.h"
#include "PiecewisePol.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "matio.h"
#include "Zmatrices.h"
#include <omp.h>
#ifdef OS_WIN
#include <time.h>
void start_timing(clock_t& t)
{
	t = clock();
}
void finish_timing(clock_t& t)
{
	t = clock() - t;
	printf("\n\nComplete. The elapsed time is %f seconds\n\n",
		((float)t) / CLOCKS_PER_SEC);
}
#else
#include <sys/time.h>	// this isn't in MSVC
void start_timing(struct timeval *start_time)
{
	gettimeofday(start_time,NULL);
}
void finish_timing(struct timeval *start_time)
{
	struct timeval end_time;
	gettimeofday(&end_time,NULL);
	printf("\n\nComplete. The elapsed time is %f seconds\n\n",
		end_time.tv_sec-start_time->tv_sec+(end_time.tv_usec-start_time->tv_usec)/1e6);
}
#endif

void CreateStruct(mat_t **file_ptr, matvar_t **struct_in,
	const char *filename, const char *mat_name,
	const char **fieldnames, unsigned int nfields)
{
	mat_t *file_p = Mat_CreateVer(filename, NULL, MAT_FT_DEFAULT);
	*file_ptr = file_p;
	if (NULL == *file_ptr) {
		fprintf(stderr, "Error opening MAT file %s\n", filename);
		return;
	}

	size_t struct_dims[2] = { 1, 1 };		// default size for a scalar Matlab struct

	*struct_in = Mat_VarCreateStruct(mat_name, 2, struct_dims, fieldnames, nfields);
	if (NULL == *struct_in) {
		fprintf(stderr, "Error creating %s struct\n", mat_name);
		Mat_Close(*file_ptr);
		return;
	}
}

void InsertCubeIntoStruct(matvar_t **struct_in, const char *name, cube &D)
{
	//write cube D
	size_t dims[3] = { D.n_rows, D.n_cols, D.n_slices };
	int tot = D.n_elem;
	double *tmp = new double[tot];
	for (int i = 0; i<tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp, 0);
	Mat_VarSetStructFieldByName(*struct_in, name, 0, field);
	delete tmp;
}

void InsertMatrixIntoStruct(matvar_t **struct_in, const char *name, MATRIX &D)
{
	
	//size_t dims[2] = { D.n_cols, D.n_rows };		// row major
	size_t dims[2] = { D.n_rows, D.n_cols };		// column major

	int tot = D.n_elem;
	double *tmp = new double[tot];
	
	for (int i = 0; i<tot; i++)
		tmp[i] = D(i);

	matvar_t *field;
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp, 0);
	Mat_VarSetStructFieldByName(*struct_in, name, 0, field);
	delete tmp;
}

void FinishStruct(mat_t **file_ptr, matvar_t **struct_in)
{
	Mat_VarWrite(*file_ptr, *struct_in, MAT_COMPRESSION_ZLIB);

	Mat_VarFree(*struct_in);
	Mat_Close(*file_ptr);
}

bool ReadGeometry(GEOMETRY& geometry, double& dt, double& c, const char* fname)
{
	mat_t *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening %s file %s\n", fname);
		return false;
	}

	// Get dt
	matvar = Mat_VarRead(matfp, "dt");
	memcpy(&dt,matvar->data,matvar->data_size);

	// Get c
	matvar = Mat_VarRead(matfp, "c");
	memcpy(&c, matvar->data, matvar->data_size);

	// Get geometry
	matvar_t* tmpp[100];				// TODO: arbitrary max number of edges
	EDGE edge;
	int index = 0;
	while (index<(int)matvar->dims[0])
	{
		matvar = Mat_VarRead(matfp, "geometry");
		if (matvar == NULL)
			return false;
		tmp = Mat_VarGetStructFieldByName(matvar, "he_idx", index);
		memcpy(&edge.he_idx, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "a", index);
		memcpy(&edge.a, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "b", index);
		memcpy(&edge.b, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "l", index);
		memcpy(&edge.l, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "t", index);
		memcpy(&edge.t, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "n", index);
		memcpy(&edge.n, tmp->data, tmp->nbytes);
		geometry.push_back(edge);
		tmpp[index] = matvar;
		index++;
	}

	for (int i = 0; i < index; i++){
		Mat_VarFree(tmpp[i]);
	}
	Mat_Close(matfp);
	return true;
}

void TempConvs_example()
{
	printf("\nComputing TempConvs example...");
	double	c = 1;								// speed of light
	double	dt = 0.1 / c;						// timestep
	VECTOR P = linspace<vec>(1e-6, 1, 1e+6);    // distances
	VECTOR P1 = P / c;

	int	k[] = { 0, 1, 2, 6, 9 };		// chosen timesteps to plot
	const UINT degree = 2;					// basis degree
	CTempConvs tempconvs;

	// Create the time basis functions
	CLagrange_interp timeBasis = CLagrange_interp(dt, degree);
	CLagrange_interp timeBasis_D = timeBasis;
	CLagrange_interp timeBasis_Nh = timeBasis;
	timeBasis_Nh.integrate();
	CLagrange_interp timeBasis_Ns = timeBasis;
	timeBasis_Ns.diff();
	// Pad the interpolators so the sizes match (Ns and D to match Nh)
	timeBasis_Ns = timeBasis_Nh.pad_coeffs(timeBasis_Ns);
	timeBasis_D = timeBasis_Nh.pad_coeffs(timeBasis_D);

	// Loop over all chosen timesteps to get graph of convolution for that timestep
	MATRIX Fh(P.size(), sizeof(k) / 4, fill::zeros);
	MATRIX Fs(P.size(), sizeof(k) / 4, fill::zeros);
	MATRIX dF(P.size(), sizeof(k) / 4, fill::zeros);

	// start timing
	#ifdef OS_WIN
	clock_t t;
	#else
	struct timeval * t = new struct timeval;
	#endif
	start_timing(t);

	for (int K = 0; K<sizeof(k) / 4; K++)
	{
		// Shifted time basis - shift and flip the time basis functions to get T(k*dt-t)
		CLagrange_interp	shiftedTB_D = timeBasis_D;
		shiftedTB_D.translate(k[K] * dt, -1);
		CLagrange_interp	shiftedTB_Nh = timeBasis_Nh;
		shiftedTB_Nh.translate(k[K] * dt, -1);
		CLagrange_interp	shiftedTB_Ns = timeBasis_Ns;
		shiftedTB_Ns.translate(k[K] * dt, -1);

		// Compute the temporal convolutions
		VECTOR vFh(P.size(), fill::zeros), vFs(P.size(), fill::zeros), vdF(P.size(), fill::zeros);
		P /= c;
		tempconvs.compute2(P1, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, vFh, vFs, vdF);

		Fh.col(K) = vFh;
		Fs.col(K) = vFs;
		dF.col(K) = vdF;
		dF = dF / c;
	}
	finish_timing(t);
	

	// Output MAT file
	mat_t *matfpT = NULL;
	matvar_t *matvar=NULL;
	const char *fieldnames[3] = { "Fh", "Fs", "dF" };
	unsigned nfields = 3;
	CreateStruct(&matfpT, &matvar, "./results/TempConvs_example.mat", "TempConvs2", fieldnames, nfields);
	
	InsertMatrixIntoStruct(&matvar, "Fh", Fh);
	InsertMatrixIntoStruct(&matvar, "Fs", Fs);
	InsertMatrixIntoStruct(&matvar, "dF", dF);

	FinishStruct(&matfpT, &matvar);
}

bool Zmatrices_example()
{
	// Load the geometry
	GEOMETRY geometry;
	double dt, c;
	if (!ReadGeometry(geometry, dt, c, "./input/cyl_res24.mat"))
	{
		fprintf(stderr, "Error opening geometry. ");
		return false;
	}

	//simulation params
	UINT N_T = 20;
	//c = 1; dt = 0.1;
	UINT outer_points_sp = 250;
	UINT inner_points_sp = 251;
	UINT outer_points = 50;
	UINT inner_points = 51;
	UINT Lagrange_degree = 1;

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	Z_matrices.status_int = 1;		// define when status indicator updates
	//Z_matrices.use_cheat();

	//	 Basis and test functions applied to each segment
	//	 S = transverse plane, Z = z-directed, d = S divergence
	CBasisFunction BasisFunction;
	Z_matrices.basis_function_Z = BasisFunction.createSquare(geometry, true);
	Z_matrices.basis_function_S = BasisFunction.createHat(geometry, true);
	Z_matrices.test_function_Z = BasisFunction.createSquare(geometry, false);
	Z_matrices.test_function_S = BasisFunction.createHat(geometry, false);

	// Lagrange interpolators (temporal basis functions)
	CLagrange_interp timeBasis = CLagrange_interp(dt, Lagrange_degree);
	Z_matrices.timeBasis_D = timeBasis;
	CLagrange_interp timeBasis_Nh = timeBasis;
	timeBasis_Nh.integrate();
	Z_matrices.timeBasis_Nh = timeBasis_Nh;
	CLagrange_interp timeBasis_Ns = timeBasis;
	timeBasis_Ns.diff();
	Z_matrices.timeBasis_Ns = timeBasis_Ns;

	Z_matrices.z_outer_points_sp = outer_points_sp;
	Z_matrices.z_inner_points_sp = inner_points_sp;
	Z_matrices.z_outer_points = outer_points;
	Z_matrices.z_inner_points = inner_points;

	// do main computation and time
	cube S, D, Dp, Nh, Ns;
	printf("\n%s\n\n", "Computing operators...");
	#ifdef OS_WIN
	clock_t t;
	#else
	struct timeval * t = new struct timeval;
	#endif
	start_timing(t);
	Z_matrices.compute(S, D, Dp, Nh, Ns);
	finish_timing(t);

	// Output MAT file as a struct that stores the operators
	mat_t *matfpZ = NULL;
	matvar_t *matvar = NULL;
	const unsigned nfields = 5;
	const char *fieldnames[nfields] = { "S", "D", "Dp", "Nh", "Ns" };
	CreateStruct(&matfpZ, &matvar, "./results/cyl_res24.mat", "Z_Matrices", fieldnames, nfields);
	InsertCubeIntoStruct(&matvar, "S", S);
	InsertCubeIntoStruct(&matvar, "D", D);
	InsertCubeIntoStruct(&matvar, "Dp", Dp);
	InsertCubeIntoStruct(&matvar, "Nh", Nh);
	InsertCubeIntoStruct(&matvar, "Ns", Ns);
	FinishStruct(&matfpZ, &matvar);

	// free memory
	S.clear(), D.clear(), Dp.clear(), Ns.clear(), Nh.clear();

	return true;
}


int main(int argc, char* argv[])
{
	//omp setup
	int nthreads, tid;
#pragma omp parallel private (tid)
	{
		tid = omp_get_thread_num();
		if (tid == 0)
		{
			nthreads = omp_get_num_threads();
			printf("\nTotal threads = %i\n\n", nthreads);
		}
	}
	omp_set_num_threads(nthreads);


	//TempConvs_example();

	if (!Zmatrices_example()){ return -1; }

	return 0;
}
