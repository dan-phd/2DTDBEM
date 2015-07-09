#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"


bool TempConvs_example()
{
	printf("\nComputing TempConvs example...");
	double	c = 1;								// speed of light
	double	dt = 0.1 / c;						// timestep
	VECTOR P = linspace<vec>(1e-6, 1, (arma::uword)1e+6);    // distances
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
	matvar_t *matvar = NULL;
	const char *fieldnames[3] = { "Fh", "Fs", "dF" };
	unsigned nfields = 3;
	CreateStruct(&matfpT, &matvar, "./results/TempConvs_example.mat", "TempConvs2", fieldnames, nfields);

	InsertMatrixIntoStruct(&matvar, "Fh", Fh);
	InsertMatrixIntoStruct(&matvar, "Fs", Fs);
	InsertMatrixIntoStruct(&matvar, "dF", dF);

	FinishStruct(&matfpT, &matvar);

	return true;
}

bool Zmatrices_example()
{
	// Load the geometry
	GEOMETRY geometry;
	double dt, c;
	if (!ReadGeometry(geometry, dt, c, "./input/cyl_res24_dual.mat"))
	{
		fprintf(stderr, "Error opening geometry. ");
		return false;
	}

	//simulation params
	UINT N_T = 10;
	c = 3e8; dt = 0.1 / c;
	UINT outer_points_sp = 50;
	UINT inner_points_sp = 51;
	UINT outer_points = 4;
	UINT inner_points = 5;
	UINT Lagrange_degree = 1;

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	Z_matrices.status_int = 1;		// define when status indicator updates
	//Z_matrices.use_cheat();

	//	 Basis and test functions applied to each segment
	//	 S = transverse plane, Z = z-directed, d = S divergence
	CBasisFunction BasisFunction;
	/*
	Z_matrices.basis_function_Z = BasisFunction.createSquare(geometry, true);
	Z_matrices.basis_function_S = BasisFunction.createHat(geometry, true);
	Z_matrices.test_function_Z = BasisFunction.createSquare(geometry, false);
	Z_matrices.test_function_S = BasisFunction.createHat(geometry, false);*/

	Z_matrices.basis_function_Z = BasisFunction.createDualSquare(geometry, true);
	Z_matrices.basis_function_S = BasisFunction.createDualHat(geometry, false);
	Z_matrices.test_function_Z = BasisFunction.createDualSquare(geometry, true);
	Z_matrices.test_function_S = BasisFunction.createDualHat(geometry, false);

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


int start_tests(int argc, char* argv[])
{

	setup_omp();

	// re-define the arguments for easier comparisons
	vector<string> args(argv, argv + argc);

	// args[0] is the program filename, argv[1] will be "test", so we start parsing at 2
	for (size_t i = 2; i < args.size(); ++i)
	{
		if (args[i] == "Zmatrices")
		{
			if (!Zmatrices_example()){ return -1; }
		}
		else if (args[i] == "TempConvs")
		{
			if (!TempConvs_example()){ return -1; }
		}
		else {
			printf("Not enough or invalid arguments, please try again.\n");
			return -1;
		}
	}

	return 0;
}
