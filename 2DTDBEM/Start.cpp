#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"
#include "Tests.h"


int main(int argc, char* argv[])
{
	
	// run tests if specified, usage Start.exe test Zmatrices
	if (argc>1) return start_tests(argc, argv);

	long long a = NaN;

	// Load the geometry
	GEOMETRY geometry; double dt, c;
	if (!ReadGeometry(geometry, dt, c, "./input/squ_cyl_res80_dual.mat"))
	{
		fprintf(stderr, "Error opening geometry. "); return false;
	}

	// Output filename
	const char* result_file = "./results/squ_cyl_res80_dual.mat";

	// Simulation parameters
	UINT N_T = 10000;
	UINT outer_points_sp = 500;
	UINT inner_points_sp = 501;
	UINT outer_points = 100;
	UINT inner_points = 101;
	UINT Lagrange_degree = 1;

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	Z_matrices.status_int = 100;		// define when status indicator updates
	//Z_matrices.use_cheat();

	// Basis and test functions applied to each segment
	// S = transverse plane, Z = z-directed, d = S divergence
	// 2nd argument (true/false) defines whether to scale by edge length
	CBasisFunction BasisFunction;
	Z_matrices.basis_function_Z = BasisFunction.createDualSquare(geometry, true);
	Z_matrices.basis_function_S = BasisFunction.createDualHat(geometry, true);
	Z_matrices.test_function_Z = BasisFunction.createDualSquare(geometry, false);
	Z_matrices.test_function_S = BasisFunction.createDualHat(geometry, false);



	// Below this line does not need to be touched
	//--------------------------------------------------------------------------
	setup_omp();

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
	CreateStruct(&matfpZ, &matvar, result_file, "Z_Matrices", fieldnames, nfields);
	InsertCubeIntoStruct(&matvar, "S", S);
	InsertCubeIntoStruct(&matvar, "D", D);
	InsertCubeIntoStruct(&matvar, "Dp", Dp);
	InsertCubeIntoStruct(&matvar, "Nh", Nh);
	InsertCubeIntoStruct(&matvar, "Ns", Ns);
	FinishStruct(&matfpZ, &matvar);

	// free memory
	S.clear(), D.clear(), Dp.clear(), Ns.clear(), Nh.clear();

	

	return 0;
}
