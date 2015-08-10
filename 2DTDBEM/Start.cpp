#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"
#include "Tests.h"

void run_Zmatrices_calculation(Zmatrices Z_matrices, double dt, UINT Lagrange_degree,
	UINT outer_points_sp, UINT inner_points_sp, UINT outer_points, UINT inner_points, const char* result_file)
{
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

	printf("\nOutput file location: %s\n\n", result_file);
}

int main(int argc, char* argv[])
{

	// Defaults
	char filename[100] = "dans_geometry_2.mat";
	UINT N_T = 10000;
	UINT outer_points = 50;

	// parse arguments, usage:
	// Start.exe filename geometry.mat N_T 100 outer_points 5
	// Start.exe test Zmatrices TempConvs
	if (argc > 1)
	{
		// re-define the arguments for easier comparisons
		vector<string> args(argv, argv + argc);

		// args[0] is the program filename, argv[1] will be "test" or filename
		// we start parsing at 2
		if (args[1] == "test")
		{
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
		}
		else {
			
			size_t i = 1;			// start parsing at 1
			while (i < args.size())
			{
				if (args[i] == "filename")
				{
					string str = args[i + 1];
					strcpy(filename, str.c_str());
					i = i + 2;
				}
				else if (args[i] == "outer_points")
				{
					const char* tmp_arg = args[i + 1].c_str();
					outer_points = strtol(tmp_arg, NULL, 10);
					i = i + 2;
				}
				else if (args[i] == "N_T")
				{
					const char* tmp_arg = args[i + 1].c_str();
					N_T = strtol(tmp_arg, NULL, 10);
					i = i + 2;
				}
				else {
					printf("Invalid arguments, please try again with filename, outer_points, or N_T\n");
					return -1;
				}
			}
		}

	}
	else {
		printf("\nUsing default arguments...\n\n");
	}

	
	setup_omp();

	// Load the geometry
	printf("Opening geometry file: %s\n", filename);
	GEOMETRY geometry; double dt, c, num_shapes_;
	string str1 = "./input/"; str1 = str1 + filename; // strcat(str1, filename);
	if (!ReadGeometry(geometry, dt, c, num_shapes_, str1.c_str()) )
	{
		fprintf(stderr, "\nError opening geometry. \n"); return false;
	}
	UINT num_shapes = (UINT)num_shapes_;
	printf("Success!\n");

	// Output filename
	string str2 = "./results/"; //strcat(str2, filename);
	string str3 = str2+filename;
	char result_file[100];
	strcpy(result_file, str3.c_str());
	//const char* result_file = str2;

	// Simulation parameters
	UINT inner_points = outer_points + 1;
	UINT outer_points_sp = 5 * outer_points;
	UINT inner_points_sp = 5 * outer_points + 1;
	UINT Lagrange_degree = 1;
	printf("\nSimulation parameters:"
		"\n\tN_T = %i"
		"\n\touter_points = %i"
		"\n\tinner_points = %i"
		"\n\touter_points_sp = %i"
		"\n\tinner_points_sp = %i"
		"\n\tLagrange_degree = %i\n\n", 
		N_T, outer_points, inner_points, outer_points_sp, inner_points_sp, Lagrange_degree);

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	Z_matrices.status_int = 1;		// define when status indicator updates
	//Z_matrices.use_cheat(); omp_set_num_threads(1); printf("\nUsing cheat!\n");

	// Basis and test functions applied to each segment
	// S = transverse plane, Z = z-directed, d = S divergence
	// 2nd argument (true/false) defines whether to scale by edge length
	CBasisFunction BasisFunction;
	Z_matrices.basis_function_Z = BasisFunction.createDualSquare(geometry, false);
	Z_matrices.basis_function_S = BasisFunction.createDualHat(geometry, false, num_shapes);
	Z_matrices.test_function_Z = BasisFunction.createDualSquare(geometry, true);
	Z_matrices.test_function_S = BasisFunction.createDualHat(geometry, true, num_shapes);


	run_Zmatrices_calculation(Z_matrices, dt, Lagrange_degree,
		outer_points_sp, inner_points_sp, outer_points, inner_points, result_file);

	
	return 0;
}
