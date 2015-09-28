#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"
#include "Tests.h"
#include "optionparser.h"


void run_Zmatrices_calculation(Zmatrices Z_matrices, double dt, UINT Lagrange_degree,
	UINT outer_points_sp, UINT inner_points_sp, UINT outer_points, UINT inner_points, const char* result_file);

// Argument types for option parser
struct Arg : public option::Arg
{
	static void printError(const char* msg1, const option::Option& opt, const char* msg2)
	{
		fprintf(stderr, "%s", msg1);
		fwrite(opt.name, opt.namelen, 1, stderr);
		fprintf(stderr, "%s", msg2);
	}

	static option::ArgStatus Unknown(const option::Option& option, bool msg)
	{
		if (msg) printError("Unknown option '", option, "'\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Required(const option::Option& option, bool msg)
	{
		if (option.arg != 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires an argument\n");
		return option::ARG_ILLEGAL;
	}

	static option::ArgStatus Numeric(const option::Option& option, bool msg)
	{
		char* endptr = 0;
		if (option.arg != 0 && strtol(option.arg, &endptr, 10)){};
		if (endptr != option.arg && *endptr == 0)
			return option::ARG_OK;

		if (msg) printError("Option '", option, "' requires a numeric argument\n");
		return option::ARG_ILLEGAL;
	}
};

// List of accepted arguments to parse, along with usage
enum  optionIndex { UNKNOWN, HELP, FILENAME, NUM_TIMESTEPS, QUAD_POINTS, TESTS, CHEAT };
const option::Descriptor usage[] = {
	{ UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: 2DTDBEM [options]\n\n"
	"Options:" },
	{ HELP, 0, "h", "help", Arg::None, "  -h,  \t--help  \tPrint usage and exit." },
	{ FILENAME, 0, "f", "file", Arg::Required, "  -f <arg>, \t--file=<arg>  \tInput mesh filename (required)." },
	{ NUM_TIMESTEPS, 0, "t", "timesteps", Arg::Numeric, "  -t <num>, \t--timesteps=<num>  \tNumber of timesteps. [10000]" },
	{ QUAD_POINTS, 0, "q", "quadrature_points", Arg::Numeric, "  -q <num>, \t--quadrature_points=<num>"
		" \tNumber of Gaussian quadrature points to use for integrating. [25]" },
	{ CHEAT, 0, "c", "cheat", Arg::None, "  -c,  \t--cheat  \tUse cheat for faster computation"
		" (only applicable for symmetric cylinder since matrices are SPD)." },
	{ TESTS, 0, "T", "test", Arg::Required, "  -T <args>, \t--test <args>  \tPerform tests." },
	{ UNKNOWN, 0, "", "", Arg::None,
	"\nExamples:\n"
	"  ./bin/2DTDBEM --file mesh.mat \n"
	"  ./bin/2DTDBEM --file=mesh.mat --timesteps=5000 --quadrature_points=10 \n"
	"  ./bin/2DTDBEM -fmesh.mat -t5000 -q10 \n"
	"  ./bin/2DTDBEM --test TempConvs \n"
	"\nThe input file is a specific Matlab type file which contains boundary edges, dt, c, and number of shapes.\n"
	},
	{ 0, 0, 0, 0, 0, 0 } };

int main(int argc, char* argv[])
{
	char filename[100];
	UINT N_T = 10000;
	UINT outer_points = 25;
	bool cheat = false;

	// Parse arguments
	argc -= (argc > 0); argv += (argc > 0); // skip program name argv[0] if present
	option::Stats stats(usage, argc, argv);
	std::vector<option::Option> options(stats.options_max);
	std::vector<option::Option> buffer(stats.buffer_max);
	option::Parser parse(usage, argc, argv, &options[0], &buffer[0]);

	if (parse.error())
		return 1;

	if (options[HELP] || argc == 0)
	{
		int columns = getenv("COLUMNS") ? atoi(getenv("COLUMNS")) : 80;
		option::printUsage(fwrite, stdout, usage, columns);
		return 0;
	}

	string temp_str;
	for (int i = 0; i < parse.optionsCount(); ++i)
	{
		option::Option& opt = buffer[i];
		switch (opt.index())
		{
		case HELP:
			// not possible, because handled further above and exits the program
		case CHEAT:
			cheat = true;
		case NUM_TIMESTEPS:
			if (opt.arg)
				N_T = strtol(opt.arg, NULL, 10);
			break;
		case QUAD_POINTS:
			if (opt.arg)
				outer_points = strtol(opt.arg, NULL, 10);
			break;
		case FILENAME:
			temp_str = opt.arg;
			strcpy(filename, temp_str.c_str());
			break;
		case UNKNOWN:
			// not possible because Arg::Unknown returns ARG_ILLEGAL
			// which aborts the parse with an error
			break;
		}
	}

	for (option::Option* o = options[TESTS]; o; o = o->next())
	{
		string this_arg = o->arg;

		if (this_arg == "TempConvs")
		{
			if (!TempConvs_example()){ return -1; }
		}
		else {
			printf("Not enough or invalid arguments, please try again.\n");
			return 1;
		}
		// tests completed
		return 0;
	}

	
	// Find out how many computing nodes we can use, then use them all
	setup_omp();

	// Load the geometry
	printf("Opening geometry file: %s\n", filename);
	GEOMETRY geometry; double dt, c, num_shapes_;
	string str1 = "./input/"; str1 = str1 + filename;
	if (!ReadGeometry(geometry, dt, c, num_shapes_, str1.c_str()) )
	{
		fprintf(stderr, "\nError opening geometry. \n"); return false;
	}
	UINT num_shapes = (UINT)num_shapes_;
	printf("Success! Number of distinct shapes = %i\n", num_shapes);

	// Output filename
	string str2 = "./results/";
	string str3 = str2+filename;
	char result_file[100];
	strcpy(result_file, str3.c_str());

	// Simulation parameters
	UINT inner_points = outer_points + 1;
	UINT outer_points_sp = 5 * outer_points;
	UINT inner_points_sp = 5 * outer_points + 1;
	UINT Lagrange_degree = 1;
	printf("\nSimulation parameters:"
		"\n\tdt = %e"
		"\n\tN_T = %i"
		"\n\touter_points = %i"
		"\n\tinner_points = %i"
		"\n\touter_points_sp = %i"
		"\n\tinner_points_sp = %i"
		"\n\tLagrange_degree = %i\n\n", 
		dt, N_T, outer_points, inner_points, outer_points_sp, inner_points_sp, Lagrange_degree);

	// Create Z_matrices object
	Zmatrices Z_matrices = Zmatrices(N_T, dt, geometry, c);
	if (cheat == true)
	{
		Z_matrices.use_cheat(); omp_set_num_threads(1); printf("\nUsing cheat!\n");
	}

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
	printf("Compressing matrices into Matlab form...\n");
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

	printf("\nDone. Output file location: %s\n\n", result_file);
}
