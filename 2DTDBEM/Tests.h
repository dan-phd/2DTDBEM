#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"
#include "Matlab_IO.h"
#include "Zmatrices.h"
#include "Common_functions.h"


bool computeConvolutions_example()
{
	printf("\nRunning computeConvolutions example...");
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
		tempconvs.compute(P1, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, vFh, vFs, vdF);

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

	CreateMatFile(&matfpT, "./results/computeConvolutions.mat");
	CreateStruct(&matvar, "TempConvs", fieldnames, nfields);

	InsertMatrixIntoStruct(&matvar, "Fh", Fh);
	InsertMatrixIntoStruct(&matvar, "Fs", Fs);
	InsertMatrixIntoStruct(&matvar, "dF", dF);

	FinishStruct(&matfpT, &matvar);
	FinishMatFile(&matfpT);

	return true;
}
