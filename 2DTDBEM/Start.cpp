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

mat_t	 *matfpL = NULL, *matfpB = NULL, *matfpT = NULL;
mat_t	 *matfpZ = NULL;

void	MatWriteLagrange(CLagrange_interp& L, const char* name)
{
	matvar_t *matvar, *field;
	size_t    dims[2] = { 1, 1 }, struct_dims[2] = { 1, 1 };
	const char *fieldnames[4] = { "partition", "coeffs", "degree", "dt" };
	unsigned nfields = 4;

	if (matfpL == NULL)
		return;

	matvar = Mat_VarCreateStruct(name, 2, struct_dims, fieldnames, nfields);
	if (NULL == matvar) {
		fprintf(stderr, "Error creating variable for 'a'\n");
		Mat_Close(matfpL);
		return;
	}
	dims[0] = 1;
	dims[1] = L.m_partition.size();
	double *tmp1 = new double[L.m_partition.size()];
	for (int i = 0; i<(int)L.m_partition.size(); i++)
		tmp1[i] = L.m_partition(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "partition", 0, field);
	delete tmp1;
	MATRIX	t = L.m_coeffs.t();
	dims[0] = t.n_cols;
	dims[1] = t.n_rows;
	double *tmp2 = new double[L.m_coeffs.n_rows*L.m_coeffs.n_cols];

	for (int i = 0; i<(int)t.n_rows; i++)
		for (int j = 0; j<(int)t.n_cols; j++)
			tmp2[i*t.n_cols + j] = t(i, j);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp2, 0);
	Mat_VarSetStructFieldByName(matvar, "coeffs", 0, field);
	delete tmp2;
	dims[0] = 1;
	dims[1] = 1;
	field = Mat_VarCreate(NULL, MAT_C_UINT32, MAT_T_UINT32, 2, dims, &L.m_degree, 0);
	Mat_VarSetStructFieldByName(matvar, "degree", 0, field);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, &L.m_dt, 0);
	Mat_VarSetStructFieldByName(matvar, "dt", 0, field);
	Mat_VarWrite(matfpL, matvar, MAT_COMPRESSION_ZLIB);
	Mat_VarFree(matvar);
}
void	MatWriteTempconvs(MATRIX& Fh, MATRIX& Fs, MATRIX& dF)
{
	matvar_t *matvar, *field;
	size_t    dims[2] = { 1, 1 }, struct_dims[2] = { 1, 1 };
	const char *fieldnames[4] = { "Fh", "Fs", "dF" };
	unsigned nfields = 3;

	if (matfpT == NULL)
		return;

	matvar = Mat_VarCreateStruct("TempConvs", 2, struct_dims, fieldnames, nfields);
	if (NULL == matvar) {
		fprintf(stderr, "Error creating variable \n");
		Mat_Close(matfpT);
		return;
	}
	dims[0] = Fh.n_rows;
	dims[1] = Fh.n_cols;
	double *tmp1 = new double[Fh.n_rows*Fh.n_cols];
	for (int i = 0; i<(int)Fh.n_rows*Fh.n_cols; i++)
		tmp1[i] = Fh(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "Fh", 0, field);
	delete tmp1;

	dims[0] = Fs.n_rows;
	dims[1] = Fs.n_cols;
	tmp1 = new double[Fs.n_rows*Fs.n_cols];
	for (int i = 0; i<(int)Fs.n_rows*Fs.n_cols; i++)
		tmp1[i] = Fs(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "Fs", 0, field);
	delete tmp1;

	dims[0] = dF.n_rows;
	dims[1] = dF.n_cols;
	tmp1 = new double[dF.n_rows*dF.n_cols];
	for (int i = 0; i<(int)dF.n_rows*dF.n_cols; i++)
		tmp1[i] = dF(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "dF", 0, field);
	delete tmp1;

	Mat_VarWrite(matfpT, matvar, MAT_COMPRESSION_ZLIB);
	Mat_VarFree(matvar);
}
void	MatWriteZMatrix(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns, char* str)
{
	matfpZ = Mat_CreateVer(str, NULL, MAT_FT_DEFAULT);
	if (NULL == matfpZ) {
		fprintf(stderr, "Error opening MAT file %s\n", str);
		return;
	}
	matvar_t *matvar, *field;
	size_t    dims[3] = { 1, 1, 1 }, struct_dims[3] = { 1, 1, 1 };
	const char *fieldnames[5] = { "S", "D", "Dp", "Nh", "Ns" };
	unsigned nfields = 5;
	if (matfpZ == NULL)
		return;
	matvar = Mat_VarCreateStruct("Z_Matrices", 3, struct_dims, fieldnames, nfields);
	if (NULL == matvar) {
		fprintf(stderr, "Error creating variable for 'a'\n");
		Mat_Close(matfpZ);
		return;
	}
	//write cube D
	dims[0] = D.n_rows; 	dims[1] = D.n_cols; 	dims[2] = D.n_slices;
	double *tmp1 = new double[D.n_rows * D.n_cols * D.n_slices];
	for (int i = 0; i<D.n_rows * D.n_cols * D.n_slices; i++){
		tmp1[i] = D(i);
	}
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "D", 0, field);
	delete tmp1;

	//write cube S
	dims[0] = S.n_rows; 	dims[1] = S.n_cols; 	dims[2] = S.n_slices;
	tmp1 = new double[S.n_rows * S.n_cols * S.n_slices];
	for (int i = 0; i<S.n_rows * S.n_cols * S.n_slices; i++)
		tmp1[i] = S(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "S", 0, field);
	delete tmp1;

	//write cube Dp
	dims[0] = Dp.n_rows; 	dims[1] = Dp.n_cols; 	dims[2] = Dp.n_slices;
	tmp1 = new double[Dp.n_rows * Dp.n_cols * Dp.n_slices];
	for (int i = 0; i<Dp.n_rows * Dp.n_cols * Dp.n_slices; i++)
		tmp1[i] = Dp(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "Dp", 0, field);
	delete tmp1;

	//write cube Nh
	dims[0] = Nh.n_rows; 	dims[1] = Nh.n_cols; 	dims[2] = Nh.n_slices;
	tmp1 = new double[Nh.n_rows * Nh.n_cols * Nh.n_slices];
	for (int i = 0; i<Nh.n_rows * Nh.n_cols * Nh.n_slices; i++)
		tmp1[i] = Nh(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "Nh", 0, field);
	delete tmp1;

	//write cube Ns
	dims[0] = Ns.n_rows; 	dims[1] = Ns.n_cols; 	dims[2] = Ns.n_slices;
	tmp1 = new double[Ns.n_rows * Ns.n_cols * Ns.n_slices];
	for (int i = 0; i<Ns.n_rows * Ns.n_cols * Ns.n_slices; i++)
		tmp1[i] = Ns(i);
	field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 3, dims, tmp1, 0);
	Mat_VarSetStructFieldByName(matvar, "Ns", 0, field);
	delete tmp1;

	Mat_VarWrite(matfpZ, matvar, MAT_COMPRESSION_ZLIB);
	Mat_VarFree(matvar);
	Mat_Close(matfpZ);

}
void	MatWriteBasis(POLYMAT& geo, const char* name)
{
	if (matfpB == NULL)
		return;
	double *tmp2 = NULL;
	for (int Iter = 0; Iter < (int)geo.size(); Iter++)
	{
		char tname[64];
		sprintf(tname, "%s%d", name, Iter);
		matvar_t *matvar, *field;
		size_t    dims[2] = { 1, 1 }, struct_dims[2] = { 1, 1 };
		const char *fieldnames[3] = { "partition", "coeffs", "degree" };
		unsigned nfields = 3;

		if (matfpB == NULL)
			return;

		matvar = Mat_VarCreateStruct(tname, 2, struct_dims, fieldnames, nfields);
		if (NULL == matvar) {
			fprintf(stderr, "Error creating variable for 'a'\n");
			Mat_Close(matfpL);
			return;
		}
		dims[0] = 1;
		dims[1] = geo.at(Iter).m_partition.size();
		double *tmp1 = new double[geo.at(Iter).m_partition.size()];
		for (int i = 0; i<(int)geo.at(Iter).m_partition.size(); i++)
			tmp1[i] = geo.at(Iter).m_partition(i);
		field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp1, 0);
		Mat_VarSetStructFieldByName(matvar, "partition", 0, field);
		delete tmp1;
		MATRIX	t = geo.at(Iter).m_coeffs.t();
		dims[0] = t.n_cols;
		dims[1] = t.n_rows;

		tmp2 = new double[t.n_rows*t.n_cols];

		for (int i = 0; i<(int)t.n_rows; i++)
			for (int j = 0; j<(int)t.n_cols; j++)
				tmp2[i*t.n_cols + j] = t(i, j);
		field = Mat_VarCreate(NULL, MAT_C_DOUBLE, MAT_T_DOUBLE, 2, dims, tmp2, 0);
		Mat_VarSetStructFieldByName(matvar, "coeffs", 0, field);
		dims[0] = 1;
		dims[1] = 1;
		field = Mat_VarCreate(NULL, MAT_C_UINT32, MAT_T_UINT32, 2, dims, &geo.at(Iter).m_degree, 0);
		Mat_VarSetStructFieldByName(matvar, "degree", 0, field);
		Mat_VarWrite(matfpB, matvar, MAT_COMPRESSION_ZLIB);
		Mat_VarFree(matvar);
		delete tmp2;
	}
}

void	Lagrange_example()
{
	double	dt = 0.1;
	UINT	degree = 2;
	CLagrange_interp	L(dt, degree);
	VECTOR	x = L.makevec(-2, 6, 0.01);
	x *= dt;
	VECTOR	y = L.eval(x);
	printf("*****Lagrange interpolator function of degree %i*****\n", degree);
	printf("-----Creating Lagrange_example.mat...-----\n");
	MatWriteLagrange(L, "./results/MyL");

	CLagrange_interp	shiftedL = L;
	shiftedL.translate(3 * dt, -1);	// shift and scale
	y = shiftedL.eval(x);
	MatWriteLagrange(shiftedL, "./results/MyshiftedL");

	CLagrange_interp	dL = L;
	dL.diff();	// differentiated
	y = dL.eval(x);
	MatWriteLagrange(dL, "./results/MydL");

	CLagrange_interp	intL = L;
	intL.integrate();	// differentiated
	y = intL.eval(x);
	MatWriteLagrange(intL, "./results/MyintL");

	CLagrange_interp	d2L = dL;
	d2L.diff();	// double differentiated
	y = d2L.eval(x);
	MatWriteLagrange(d2L, "./results/Myd2L");

	CLagrange_interp	int2L = intL;
	int2L.integrate();	// double differentiated
	y = int2L.eval(x);
	MatWriteLagrange(int2L, "./results/Myint2L");

	CLagrange_interp	dL2 = intL.pad_coeffs(dL);
	y = dL2.eval(x);
	MatWriteLagrange(dL2, "./results/MydL2");

	printf("-----Lagrange interpolator function end-----\n");
}

BOOL	ReadGeometry(GEOMETRY& geometry, char*fname)
{
	mat_t    *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;

	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening %s file %s\n", fname);
		return FALSE;
	}
	EDGE	edge;
	int index = 0;
	matvar = Mat_VarReadNext(matfp);
	while (index<(int)matvar->dims[0])
	{
		matvar = Mat_VarRead(matfp, matvar->name);
		if (matvar == NULL)
			return	FALSE;
		tmp = Mat_VarGetStructFieldByName(matvar, "a", index);
		memcpy(&edge.a, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "b", index);
		memcpy(&edge.b, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "l", index);
		memcpy(&edge.l, tmp->data, tmp->nbytes);
		geometry.push_back(edge);
		index++;
	}
	Mat_VarFree(matvar);
	Mat_Close(matfp);
	return	TRUE;
}

BOOL	ReadGeometry(ZGEOMETRY& zgeometry, char*fname)
{
	mat_t    *matfp = NULL;
	matvar_t *matvar = NULL, *tmp = NULL;
	matfp = Mat_Open(fname, MAT_ACC_RDONLY);

	if (NULL == matfp)
	{
		fprintf(stderr, "Error opening %s file %s\n", fname);
		return FALSE;
	}
	ZEDGE	z_edge;
	int index = 0;
	matvar = Mat_VarReadNext(matfp);
	matvar_t* tmpp[100];
	while (index<(int)matvar->dims[0])
	{
		matvar = Mat_VarRead(matfp, matvar->name);
		if (matvar == NULL)
			return	FALSE;
		tmp = Mat_VarGetStructFieldByName(matvar, "he_idx", index);
		memcpy(&z_edge.he_idx, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "a", index);
		memcpy(&z_edge.a, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "b", index);
		memcpy(&z_edge.b, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "l", index);
		memcpy(&z_edge.l, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "t", index);
		memcpy(&z_edge.t, tmp->data, tmp->nbytes);
		tmp = Mat_VarGetStructFieldByName(matvar, "n", index);
		memcpy(&z_edge.n, tmp->data, tmp->nbytes);
		zgeometry.push_back(z_edge);
		tmpp[index] = matvar;
		index++;
		//		free(matvar);
		//		matvar=Mat_VarRead( matfp, matvar->name );
	}
	for (int i = 0; i < index; i++){
		Mat_VarFree(tmpp[i]);
	}
	Mat_Close(matfp);
	return	TRUE;
}

void	Basis_example(char*fname)
{
	// Create basis functions
	printf("*****Basis function*****\n");
	printf("-----Creating Basis_example.mat...-----\n");
	GEOMETRY	geometry;
	ReadGeometry(geometry, fname);

	CBasisFunction	BasisFunction;
	POLYMAT	hat_function = BasisFunction.createHat(geometry, FALSE);
	MatWriteBasis(hat_function, "./results/hat");
	POLYMAT square_function = BasisFunction.createSquare(geometry, FALSE);
	MatWriteBasis(square_function, "./results/square");
	BasisFunction.divergence(hat_function, geometry);
	MatWriteBasis(hat_function, "./results/div");
	printf("-----Basis function end-----\n");
}

void	TempConvs_example()
{
	printf("*****TempConvs function*****\n");
	printf("-----Creating Tempconvs_example.mat...-----\n");
	double	c = 3e8;         // speed of light
	double	dt = 0.1 / c;         // timestep
	vec P = linspace<vec>(1e-4, 1, 1e+4);     // distances
	vec P1 = P / c;

	int	k[] = { 0, 1, 2, 6, 9 };      // chosen timesteps to plot
	UINT	degree = 2;         // basis degree

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
	mat Fh(P.size(), sizeof(k) / 4, fill::zeros);
	mat Fs(P.size(), sizeof(k) / 4, fill::zeros);
	mat dF(P.size(), sizeof(k) / 4, fill::zeros);
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
		CTempConvs	tempconvs;
		P /= c;
		tempconvs.TempConvs(P1, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, vFh, vFs, vdF);

		Fh.col(K) = vFh;
		Fs.col(K) = vFs;
		dF.col(K) = vdF;
		dF = dF / c;
	}
	MatWriteTempconvs(Fh, Fs, dF);
	printf("-----TempConvs function end-----\n");

}

BOOL Zmatrices_example()
{
	//simulation params
	double c = 2.997924580105029e8;
	UINT N_T = 5;
	double dt = 2.821273163357656e-12;
	UINT outer_points_sp = 500;
	UINT inner_points_sp = 501;
	UINT outer_points = 100;
	UINT inner_points = 101;
	UINT Lagrange_degree = 1;

	// Load the geometry
	ZGEOMETRY	zgeometry;
	if (!ReadGeometry(zgeometry, ".\input\geometry_circle.mat"))
	{
		return false;
	}

	// Using object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, zgeometry, c);

	//	 Basis and test functions applied to each segment
	//	 S = transverse plane, Z = z-directed, d = S divergence
	CBasisFunction	BasisFunction;
	Z_matrices.basis_function_Z = BasisFunction.createSquare(zgeometry, TRUE);
	Z_matrices.basis_function_S = BasisFunction.createHat(zgeometry, FALSE);
	Z_matrices.test_function_Z = BasisFunction.createSquare(zgeometry, TRUE);
	Z_matrices.test_function_S = BasisFunction.createHat(zgeometry, FALSE);

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

	//count the elapsed times per miliseconds	

	// do Zmatrices.compute()
	cube S, D, Dp, Nh, Ns;
	printf("%s\n\n", "Now, computing S1, D1, Nh1, Ns1, Z1 ...");
	DWORD t_time;	t_time = GetTickCount();
	Z_matrices.compute(S, D, Dp, Nh, Ns, Z_matrices);
	printf("%s\n", "Complete computing S1, D1, Nh1, Ns1, Z1 .");
	printf("The elapsed time is %d milliseconds\n\n", GetTickCount() - t_time);
	MatWriteZMatrix(S, D, Dp, Nh, Ns, "./results/Zmatrices1.mat");
	S.clear(), D.clear(), Dp.clear(), Ns.clear(), Nh.clear();

	return true;

}
int main(int argc, char* argv[])
{
	// Zmatrices example
	if (!Zmatrices_example()){ return 0; }


	// Lagrange example
	/*
	matfpL = Mat_CreateVer("./results/Lagrange_example.mat", NULL, MAT_FT_DEFAULT);
	if (NULL == matfpL) {
		fprintf(stderr, "Error opening MAT file %s\n", "Largrange_example.mat");
		return -1;
	}
	Lagrange_example();
	Mat_Close(matfpL);
	//*/


	// TempConvs example
	/*
	matfpT = Mat_CreateVer("./results/TempConvs_example.mat", NULL, MAT_FT_DEFAULT);
	if (NULL == matfpT) {
		fprintf(stderr, "Error opening MAT file %s\n", "TempConvs_example.mat");
		return -1;
	}
	DWORD t_time;	t_time = GetTickCount();
	TempConvs_example();
	printf("The elapsed time is %d milliseconds\n", GetTickCount() - t_time);
	Mat_Close(matfpT);
	//*/


	// Basis example
	/*
	matfpB = Mat_CreateVer("./results/Basis_example.mat", NULL, MAT_FT_DEFAULT);
	if (NULL == matfpB) {
		fprintf(stderr, "Error opening MAT file %s\n", "Basis_example.mat");
		return -1;
	}
	if (argc<2)
	{
		//printf( "cannot find argument! Usage: mat2c_new.exe geometry.mat\n" );
		//return	-1;
		
		Basis_example("./input/geometry.mat");
	}
	else {
		Basis_example(argv[1]);
	}
	Mat_Close(matfpB);
	//*/


	return 0;
}
