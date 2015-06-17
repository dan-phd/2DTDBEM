//#include <stdlib.h>
#include "Define.h"
#include "BasisFunction.h"
#include "PiecewisePol.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "matio.h"
#include "windows.h"		// for timing
#include "Zmatrices.h"

typedef unsigned long DWORD;

void CreateStruct(mat_t **file_ptr, matvar_t **struct_in,
	char *filename, const char *mat_name,
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

BOOL ReadGeometry(GEOMETRY& geometry, char* fname)
{
	mat_t *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening %s file %s\n", fname);
		return FALSE;
	}
	EDGE edge;
	int index = 0;
	matvar = Mat_VarReadNext(matfp);
	matvar_t* tmpp[100];				// TODO: arbitrary max number of edges
	while (index<(int)matvar->dims[0])
	{
		matvar = Mat_VarRead(matfp, matvar->name);
		if (matvar == NULL)
			return	FALSE;
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
	return	TRUE;
}

void TempConvs_example()
{
	printf("\nComputing TempConvs example...");
	double	c = 1;								// speed of light
	double	dt = 0.1 / c;						// timestep
	VECTOR P = linspace<vec>(1e-2, 1, 1e+2);    // distances
	VECTOR P1 = P / c;

	int	k[] = { 0, 1, 2, 6, 9 };		// chosen timesteps to plot
	UINT degree = 2;					// basis degree

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
	DWORD t_time;	t_time = GetTickCount();
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
		VECTOR	vFh, vFs, vdF;
		CTempConvs tempconvs;
		P /= c;
		tempconvs.compute(P1, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, vFh, vFs, vdF);

		Fh.col(K) = vFh;
		Fs.col(K) = vFs;
		dF.col(K) = vdF;
		dF = dF / c;
	}
	printf("\nComplete. The elapsed time is %d milliseconds\n", GetTickCount() - t_time);
	

	// Output MAT file
	mat_t *matfpT = NULL;
	matvar_t *matvar=NULL;
	const char *fieldnames[3] = { "Fh", "Fs", "dF" };
	unsigned nfields = 3;
	CreateStruct(&matfpT, &matvar, "./results/TempConvs_example.mat", "TempConvs", fieldnames, nfields);
	
	InsertMatrixIntoStruct(&matvar, "Fh", Fh);
	InsertMatrixIntoStruct(&matvar, "Fs", Fs);
	InsertMatrixIntoStruct(&matvar, "dF", dF);

	FinishStruct(&matfpT, &matvar);
}

BOOL Zmatrices_example()
{
	//simulation params
	double c = 2.997924580105029e8;
	UINT N_T = 15;
	double dt = 2.821273163357656e-12;
	UINT outer_points_sp = 150;
	UINT inner_points_sp = 151;
	UINT outer_points = 10;
	UINT inner_points = 11;
	UINT Lagrange_degree = 1;

	// Load the geometry
	GEOMETRY geometry;
	if (!ReadGeometry(geometry, "./input/geometry_circle.mat"))
	{
		fprintf(stderr, "Error opening geometry. ");
		return FALSE;
	}

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	Z_matrices.status_int = 1;
	//Z_matrices.use_cheat();

	//	 Basis and test functions applied to each segment
	//	 S = transverse plane, Z = z-directed, d = S divergence
	CBasisFunction BasisFunction;
	Z_matrices.basis_function_Z = BasisFunction.createSquare(geometry, true);
	Z_matrices.basis_function_S = BasisFunction.createHat(geometry, false);
	Z_matrices.test_function_Z = BasisFunction.createSquare(geometry, true);
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

	// do main computation
	cube S, D, Dp, Nh, Ns;
	printf("\n%s\n\n", "Computing operators...");
	DWORD t_time;	t_time = GetTickCount();
	Z_matrices.compute(S, D, Dp, Nh, Ns);
	printf("\n\nComplete. The elapsed time is %d milliseconds\n\n", GetTickCount() - t_time);

	// Output MAT file as a struct that stores the operators
	mat_t *matfpZ = NULL;
	matvar_t *matvar = NULL;
	const unsigned nfields = 5;
	const char *fieldnames[nfields] = { "S", "D", "Dp", "Nh", "Ns" };
	CreateStruct(&matfpZ, &matvar, "./results/Zmatrices1.mat", "Z_matrices", fieldnames, nfields);
	InsertCubeIntoStruct(&matvar, "S", S);
	InsertCubeIntoStruct(&matvar, "D", D);
	InsertCubeIntoStruct(&matvar, "Dp", Dp);
	InsertCubeIntoStruct(&matvar, "Nh", Nh);
	InsertCubeIntoStruct(&matvar, "Ns", Ns);
	FinishStruct(&matfpZ, &matvar);

	// free memory
	S.clear(), D.clear(), Dp.clear(), Ns.clear(), Nh.clear();

	return TRUE;
}


int main(int argc, char* argv[])
{
	//TempConvs_example();

	if (!Zmatrices_example()){ return -1; }

	return 0;
}
