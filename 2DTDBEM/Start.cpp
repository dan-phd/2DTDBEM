#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"
#include "Tests.h"
#include "optionparser.h"

void run_Zmatrices_calculation(Zmatrices& Z_matrices, UINT Lagrange_degree,
	UINT outer_points, const char* result_file);

void run_Zmatrices_calculation_scattered_field(Zmatrices& Z_matrices, UINT Lagrange_degree,
	GRID& rho, const char* result_file);

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
enum  optionIndex { UNKNOWN, HELP, FILENAME, NUM_TIMESTEPS, QUAD_POINTS, TESTS, DEGREE, CHEAT, SUFFIX, SCATTERED };
const option::Descriptor usage[] = {
	{ UNKNOWN, 0, "", "", Arg::Unknown, "USAGE: ./bin/2DTDBEM [options]\n\n"
	"Options:" },
	{ HELP, 0, "h", "help", Arg::None, "  -h,  \t--help  \tPrint usage and exit." },
	{ FILENAME, 0, "f", "file", Arg::Required, "  -f <arg>, \t--file=<arg>  \tInput mesh filename, without extension (required)." },
	{ NUM_TIMESTEPS, 0, "t", "timesteps", Arg::Numeric, "  -t <num>, \t--timesteps=<num>  \tNumber of timesteps. [1000]" },
	{ QUAD_POINTS, 0, "q", "quadrature_points", Arg::Numeric, "  -q <num>, \t--quadrature_points=<num>"
	" \tNumber of Gaussian quadrature points used on the outer integral. [25]" },
	{ DEGREE, 0, "d", "degree", Arg::Numeric, "  -d <num>, \t--degree=<num>"
	" \tLagrange interpolator degree to use for the temporal convolutions. [1]" },
	{ SUFFIX, 0, "s", "suffix", Arg::Optional, "  -s <arg>, \t--suffix=<arg>  \tTo attach to end of result filename." },
	{ CHEAT, 0, "c", "cheat", Arg::None, "  -c,  \t--cheat  \tUse cheat for faster computation"
	" (only applicable for symmetric cylinder since matrices are SPD)." },
	{ SCATTERED, 0, "S", "scattered", Arg::None, "  -S <args>, \t--scattered <args>  \tCompute scattered field." },
	{ TESTS, 0, "T", "test", Arg::Required, "  -T <args>, \t--test <args>  \tPerform tests." },
	{ UNKNOWN, 0, "", "", Arg::None,
	"\nExamples:\n"
	"  ./bin/2DTDBEM --file cyl_res21 \n"
	"  ./bin/2DTDBEM --file=cyl_res21 --timesteps=300 --quadrature_points=4 --cheat \n"
	"  ./bin/2DTDBEM -fcyl_res21 -t300 -q4 -c \n"
	"  ./bin/2DTDBEM --scattered -fscattered_mesh -t5000  \n"
	"  ./bin/2DTDBEM --test computeConvolutions \n"
	"\nThe input file is a specific Matlab type file which contains boundary edges, dt, c, number of shapes, and an option to decide whther or not to use dual basis functions."
	"\nThe results folder contains the output files, which have the same name as the input, plus an optional suffix.\n"
	},
	{ 0, 0, 0, 0, 0, 0 } };

int main(int argc, char* argv[])
{
	char filename[100], suffix[20];
	UINT N_T = 1000;
	UINT outer_points = 25;
	UINT Lagrange_degree = 1;
	bool cheat = false, add_suffix = false, compute_scattered_field = false;

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
			// not possible - handled above
		case CHEAT:
			cheat = true;
			break;
		case SCATTERED:
			compute_scattered_field = true;
			break;
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
		case DEGREE:
			Lagrange_degree = strtol(opt.arg, NULL, 10);
			break;
		case SUFFIX:
			temp_str = opt.arg;
			strcpy(suffix, temp_str.c_str());
			add_suffix = true;
			break;
		case UNKNOWN:
			// not possible because Arg::Unknown returns ARG_ILLEGAL
			// which aborts the parse with an error
			break;
		}
	}

	// Check for tests
	for (option::Option* o = options[TESTS]; o; o = o->next())
	{
		string this_arg = o->arg;

		if (this_arg == "computeConvolutions")
		{
			if (!computeConvolutions_example()){ return -1; }
		}
		else {
			printf("Not enough or invalid arguments, please try again.\n");
			return 1;
		}
		// tests completed
		return 0;
	}

	// Load the geometry
	printf("Opening geometry file: %s.mat", filename);
	GEOMETRY geometry; double dt, c, num_shapes_, dual_;
	string str1 = "./input/"; str1 = str1 + filename + ".mat";
	// read:  'boundary', 'dt', 'c', 'num_shapes', 'dual'
	if (!ReadGeometry(geometry, dt, c, num_shapes_, dual_, str1.c_str()))
	{
		fprintf(stderr, "\n\nError opening file. \n"); return false;
	}
	// if computing scattered field, read 'rho'
	GRID rho;
	if (compute_scattered_field)
	{
		if (!ReadScatteredFieldGeometry(rho, str1.c_str()))
		{
			fprintf(stderr, "\n\nError opening file for scattered field. \n");
			return false;
		}
	}
	UINT num_shapes = (UINT)num_shapes_;
	UINT dual = (UINT)dual_;
	printf("\t- done. \nNumber of distinct shapes = %i\n", num_shapes);
	printf("Using dual basis functions = %s\n", dual ? "true" : "false");

	// Output filename
	string str2 = "./results/";
	string str3 = str2 + filename;
	if (add_suffix)	str3 = str3 + suffix;
	str3 += ".mat";
	char result_file[124];
	strcpy(result_file, str3.c_str());

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
	if (dual_ == 1)
	{
		Z_matrices.basis_function_Z = BasisFunction.createDualSquare(geometry, false);
		Z_matrices.basis_function_S = BasisFunction.createDualHat(geometry, false, num_shapes);
		Z_matrices.test_function_Z = BasisFunction.createDualSquare(geometry, true);
		Z_matrices.test_function_S = BasisFunction.createDualHat(geometry, true, num_shapes);
	}
	else
	{
		Z_matrices.basis_function_Z = BasisFunction.createSquare(geometry, true);
		Z_matrices.basis_function_S = BasisFunction.createHat(geometry, false, num_shapes);
		Z_matrices.test_function_Z = BasisFunction.createSquare(geometry, true);
		Z_matrices.test_function_S = BasisFunction.createHat(geometry, true, num_shapes);
	}


	// Find out how many computing nodes we can use, then use them all (if not using cheat)
	setup_omp();


	// Perform computation
	if (compute_scattered_field)
		run_Zmatrices_calculation_scattered_field(Z_matrices, Lagrange_degree, rho, result_file);
	else
		run_Zmatrices_calculation(Z_matrices, Lagrange_degree, outer_points, result_file);


	return 0;
}


void run_Zmatrices_calculation(Zmatrices& Z_matrices, UINT Lagrange_degree,
	UINT outer_points, const char* result_file)
{
	// Simulation parameters
	double dt = Z_matrices.z_dt;
	double c = Z_matrices.z_c;
	UINT inner_points = outer_points + 1;
	UINT outer_points_sp = 5 * outer_points;
	UINT inner_points_sp = 5 * outer_points + 1;
	printf("\nSimulation parameters:"
		"\n\tc = %e"
		"\n\tdt = %e"
		"\n\tN_T = %i"
		"\n\touter_points = %i"
		"\n\tinner_points = %i"
		"\n\touter_points_sp = %i"
		"\n\tinner_points_sp = %i"
		"\n\tLagrange_degree = %i\n\n",
		c, dt, Z_matrices.z_N_T, outer_points, inner_points,
		outer_points_sp, inner_points_sp, Lagrange_degree);

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

	// Combine singular and hypersingular contributions of N operator
	cube N(Nh + Ns / (c*c));

	// Output MAT file that stores the operators
	printf("Compressing matrices into Matlab form...\n");
	mat_t *matfpZ = NULL;
	matvar_t *matvar = NULL;

	// Save matrices as separate matlab variables
	double NT = S.n_slices;
	if (CreateMatFile(&matfpZ, result_file) == -1)
	{
		return;
	}
	else
	{
		InsertVar(&matfpZ, "N_T", &NT);
		InsertCube(&matfpZ, "S", S);
		InsertCube(&matfpZ, "D", D);
		InsertCube(&matfpZ, "Dp", Dp);
		InsertCube(&matfpZ, "N", N);
		FinishMatFile(&matfpZ);

		// free memory
		S.clear(), D.clear(), Dp.clear(), Ns.clear(), Nh.clear();

		printf("\t done. \nOutput file location: %s\n\n", result_file);
	}

}

void run_Zmatrices_calculation_scattered_field(Zmatrices& Z_matrices, UINT Lagrange_degree,
	GRID& rho, const char* result_file)
{
	// Simulation parameters
	double dt = Z_matrices.z_dt;
	double c = Z_matrices.z_c;
	UINT inner_points = 1;
	printf("\nSimulation parameters for scattered field:"
		"\n\tNumber of grid vertices = %i"
		"\n\tc = %e"
		"\n\tdt = %e"
		"\n\tN_T = %i"
		"\n\tinner_points = %i"
		"\n\tLagrange_degree = %i\n\n",
		rho.size(), c, dt, Z_matrices.z_N_T, inner_points, Lagrange_degree);

	// Lagrange interpolators (temporal basis functions)
	CLagrange_interp timeBasis = CLagrange_interp(dt, Lagrange_degree);
	Z_matrices.timeBasis_D = timeBasis;
	CLagrange_interp timeBasis_Nh = timeBasis;
	timeBasis_Nh.integrate();
	Z_matrices.timeBasis_Nh = timeBasis_Nh;
	CLagrange_interp timeBasis_Ns = timeBasis;
	timeBasis_Ns.diff();
	Z_matrices.timeBasis_Ns = timeBasis_Ns;

	Z_matrices.z_inner_points = inner_points;

	// do main computation and time
	cube S, D;
	printf("\n%s\n\n", "Computing operators...");
#ifdef OS_WIN
	clock_t t;
#else
	struct timeval * t = new struct timeval;
#endif
	start_timing(t);
	Z_matrices.compute_fields(S, D, rho);
	finish_timing(t);

	// TODO: Factor of c?
	D /= c;

	// Output MAT file that stores the operators
	printf("Compressing matrices into Matlab form...\n");
	mat_t *matfpZ = NULL;
	matvar_t *matvar = NULL;

	// Save matrices as separate matlab variables
	double NT = S.n_slices;
	if (CreateMatFile(&matfpZ, result_file) == -1)
	{
		return;
	}
	else
	{
		InsertVar(&matfpZ, "N_T", &NT);
		InsertCube(&matfpZ, "S", S);
		InsertCube(&matfpZ, "D", D);
		FinishMatFile(&matfpZ);

		// free memory
		S.clear(), D.clear();

		printf("\t done. \nOutput file location: %s\n\n", result_file);
	}

}
