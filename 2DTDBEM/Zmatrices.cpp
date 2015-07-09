#include "Zmatrices.h"
#include "Lagrange_interp.h"
#include "BasisFunction.h"
#include "TempConvs.h"
#include <omp.h>

Zmatrices::Zmatrices(void)
{
}

Zmatrices::~Zmatrices(void)
{
}

//constructor
Zmatrices::Zmatrices(UINT N_T, double dt, GEOMETRY& ZGegeometry, double c)
{
	z_N_T = N_T;
	z_dt = dt;
	z_geom_obj = ZGegeometry;
	z_c = c;

	z_outer_points_sp = 50;
	z_inner_points_sp = 51;
	z_outer_points = 3;
	z_inner_points = 4;
	N_E = (UINT)ZGegeometry.size();

	cheat = false;
	status_int = 100;
}

void Zmatrices::use_cheat()
{
	cheat = true;
}

int Zmatrices::getMax(POLYMAT& obj)
{
	int res = obj[0].m_degree;
	for (int i = 1; i < obj.size(); i++)
	{
		if (res < obj[i].m_degree)
			res = obj[i].m_degree;
	}
	return res;
}

MATRIX Zmatrices::vertcat(POLYMAT& obj, UINT num_segments)
{
	MATRIX res((const arma::uword)obj.size(), num_segments, fill::zeros);

	for (int i = 0; i < obj.size(); i++)
	{
		//obj[i].m_coeffs.print("coeffs");
		res.row(i) = obj[i].m_coeffs;
	}
	return res;
}

MATRIX Zmatrices::bsxPower(VECTOR& bsxfun_A, VECTOR& bsxfun_B)
{
	MATRIX res;
	res.set_size(bsxfun_A.n_rows, bsxfun_B.size());

	for (UINT i = 0; i < bsxfun_B.size(); i++){
		res.col(i) = pow(bsxfun_A, bsxfun_B(i));
	}

	return res;
}

MATRIX Zmatrices::bsxPlus(POINT2D a, POINT2D b, VECTOR& s)
{
	MATRIX res;
	res.set_size(s.size(), 2);

	res.col(0) = a.x + s * (b.x - a.x);
	res.col(1) = a.y + s * (b.y - a.y);

	return res;
}

MATRIX Zmatrices::Gcoeffs_create(VECTOR& s0, VECTOR& w0, VECTOR&si, VECTOR& wi, int max_deg_test, int max_deg_basis)
{
	//Multiply quadrature points up to max polynomial degree
	VECTOR tmp = max_deg_test - linspace(0, max_deg_test, max_deg_test + 1);
	MATRIX Gcoeffs_test = bsxPower(s0, tmp);

	tmp = max_deg_basis - linspace(0, max_deg_basis, max_deg_basis + 1);
	MATRIX Gcoeffs_basis = bsxPower(si, tmp);

	//Multiply with quadrature weights
	Gcoeffs_test.each_col() %= w0;
	Gcoeffs_basis.each_col() %= wi;

	int test_points = s0.size();
	int basis_points = si.size();

	MATRIX G(test_points*(max_test_deg + 1), basis_points*(max_basis_deg + 1), fill::zeros);

	MATRIX tmp2(test_points, basis_points);
	for (int pm = 0; pm < max_deg_basis + 1; pm++)
	{
		int offset_m = pm*basis_points;

		for (int pn = 0; pn < max_deg_test + 1; pn++)
		{
			int offset_n = pn*test_points;

			for (int i = 0; i < basis_points; i++)
			{
				tmp2.col(i) = Gcoeffs_test.col(pn) * Gcoeffs_basis(i, pm);
			}

			G(span(offset_n, offset_n + test_points-1),
						 span(offset_m, offset_m + basis_points-1)) = tmp2;	//( span(rows),span(cols) )
		}
	}
	
	return G;


}

void Zmatrices::make_distances_lookup_table()
{
	//parameterise each segment for Gaussian quadrature integration
	VECTOR s0, w0;
	lgquad1(z_outer_points, s0, w0);
	VECTOR si, wi;
	lgquad1(z_inner_points, si, wi);

	//Gaussian quadrature integration for singularity points
	VECTOR s0_sp, w0_sp;
	lgquad1(z_outer_points_sp, s0_sp, w0_sp);
	VECTOR si_sp, wi_sp;
	lgquad1(z_inner_points_sp, si_sp, wi_sp);

	//container for lookup tables
	MATRIX **distances_ = CreateMatrix(N_E, N_E);
	MATRIX **over_dnp_ = CreateMatrix(N_E, N_E);
	MATRIX **over_dn_ = CreateMatrix(N_E, N_E);

	UINT k(0), m(0), n(0), max_n(N_E);
	for (m = 0; m < N_E; m++){

		cheat ? max_n = (m == 0) ? N_E : 1 : false;
		for (n = 0; n < max_n; n++){

			//When dealing with singularities at self patch and neighbouring edges, increase number of quadrature points
			if (n == (m + 1) % N_E || n == m || n == (m + N_E - 1) % N_E)
			{
				find_distances_between_2_edges(m,n,distances_[m][n], over_dnp_[m][n], over_dn_[m][n], si_sp, s0_sp);
			}
			else
			{
				//normal calculation when there are no singularities
				find_distances_between_2_edges(m,n,distances_[m][n], over_dnp_[m][n], over_dn_[m][n], si, s0);
			}

		} // end n
	} // end m

	distances = distances_;
	over_dnp = over_dnp_;
	over_dn = over_dn_;

}

void Zmatrices::find_distances_between_2_edges(UINT m, UINT n, MATRIX& P, MATRIX& dn_p, MATRIX& dn_,
	VECTOR& s_i, VECTOR& s_o)
{

	//Find geometry around observation point
	EDGE tmp_edge = z_geom_obj[m];
	POINT2D a_m = tmp_edge.a;
	POINT2D b_m = tmp_edge.b;
	POINT2D n_m = tmp_edge.n;

	//Find geometry around source point
	tmp_edge = z_geom_obj[n];
	POINT2D a_n = tmp_edge.a;
	POINT2D b_n = tmp_edge.b;
	POINT2D n_n = tmp_edge.n;

	// Find inner and outer Gaussian quadrature points
	int inner_quad_points(numel(s_i));
	int outer_quad_points(numel(s_o));

	// Make rho_m a vector of outer Gaussian quadrature points along observation edge
	MATRIX rho_m = bsxPlus(a_m, b_m, s_o);

	// Make rho_n a vector of inner Gaussian quadrature points along source edge
	MATRIX rho_n = bsxPlus(a_n, b_n, s_i);

	// P is an array of distances from rho_n to observation point, rho_m
	cube rho_mn;
	rho_mn.set_size(outer_quad_points, inner_quad_points, 2);
	for (int i = 0; i < 2; i++)
	{
		MATRIX tmp(outer_quad_points, inner_quad_points);
		for (int j = 0; j < inner_quad_points; j++)
		{
			tmp.col(j) = rho_m.col(i) - rho_n(j, i);
		}
		rho_mn.slice(i) = tmp;
	}
	cube t_rho_mn = square(rho_mn);
	MATRIX P_tmp( sqrt(t_rho_mn.slice(0) + t_rho_mn.slice(1)) );
	P = P_tmp / z_c;

	// n dot P/|P| and n' dot P/|P|
	cube unit_P(rho_mn);
	for (UINT i = 0; i < unit_P.n_slices; i++)
		unit_P.slice(i) /= P_tmp;
	
	dn_p = -(unit_P.slice(0) * n_n.x + unit_P.slice(1) * n_n.y);
	dn_ = (unit_P.slice(0) * n_m.x + unit_P.slice(1) * n_m.y);
}

MATRIX Zmatrices::Perform_quadrature(MATRIX& G_coeffs, MATRIX& F)
{
	MATRIX res(max_test_deg+1, max_basis_deg+1, fill::zeros);
	int basis_points(F.n_cols), test_points(F.n_rows);

	for (UINT pm = 0; pm < max_test_deg+1; pm++)
	{
		int offset_m(pm*basis_points);

		for (UINT pn = 0; pn < max_basis_deg+1; pn++)
		{
			int offset_n(pn*test_points);
			res(pn, pm) = accu(G_coeffs(span(offset_n, offset_n + test_points - 1),
										span(offset_m, offset_m + basis_points - 1)) % F);
		}
	}

	return res;
}

//computation of single elements of all operator matrices
void Zmatrices::Z_calc(const UINT& m, const UINT& n,
	const POINT2D& t_m, const POINT2D& t_n, const double& lmn,
	CLagrange_interp shiftedTB_D, CLagrange_interp shiftedTB_Nh, CLagrange_interp shiftedTB_Ns,
	const UINT& inner_quad_points, const UINT& outer_quad_points,
	MATRIX& G_dd, MATRIX& G_SZ, MATRIX& G_ZS, MATRIX& G_SS, MATRIX& G_ZZ,
	MATRIX& coeffs_nh, MATRIX& coeffs_ns, MATRIX& coeffs_d, MATRIX& coeffs_s, MATRIX& coeffs_dp)
{
	// Compute the temporal convolutions with an array of distances
	MATRIX P(distances[m][n]);
	MATRIX Fh(P.n_rows, P.n_cols, fill::zeros), F(Fh), dF(Fh);

	tempconvs.compute2(P, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, Fh, F, dF);
	dF /= z_c;

	MATRIX dF_dnp( over_dnp[m][n] % dF );
	MATRIX dF_dn( over_dn[m][n] % dF );

	//Perform quadrature
	coeffs_nh = Perform_quadrature(G_dd, Fh)	* lmn;
	coeffs_ns = Perform_quadrature(G_SS, F)		* lmn *(t_m.x*t_n.x + t_m.y*t_n.y);
	coeffs_d  = Perform_quadrature(G_ZS, dF_dnp)* lmn;
	coeffs_s  = Perform_quadrature(G_ZZ, F)		* lmn;
	coeffs_dp = Perform_quadrature(G_SZ, dF_dn) * lmn;
}

void Zmatrices::compute(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns)
{
	//parameterise each segment for Gaussian quadrature integration
	VECTOR s0, w0;
	lgquad1(z_outer_points, s0, w0);
	VECTOR si, wi;
	lgquad1(z_inner_points, si, wi);

	//Gaussian quadrature integration for singularity points
	VECTOR s0_sp, w0_sp;
	lgquad1(z_outer_points_sp, s0_sp, w0_sp);
	VECTOR si_sp, wi_sp;
	lgquad1(z_inner_points_sp, si_sp, wi_sp);

	//Make distances lookup table
	make_distances_lookup_table();

	//Find divergence of transverse basis and test functions
	CBasisFunction basis_function_d = basis_function_S.divergence(z_geom_obj);
	CBasisFunction test_function_d = test_function_S.divergence(z_geom_obj);

	//Number of basis and test segments for Z-directed, transverse and divergence fields
	UINT num_basis_segments_Z = (UINT)basis_function_Z.pol.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_Z = (UINT)test_function_Z.pol.size() / (UINT)z_geom_obj.size();
	UINT num_basis_segments_S = (UINT)basis_function_S.pol.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_S = (UINT)test_function_S.pol.size() / (UINT)z_geom_obj.size();
	UINT num_basis_segments_d = (UINT)basis_function_d.pol.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_d = (UINT)basis_function_d.pol.size() / (UINT)z_geom_obj.size();

	//Concatenate the basis and test coefficients and find the maximum degree
	int max_deg_basis_Z = getMax(basis_function_Z.pol);
	MATRIX basis_coeffs_Z = vertcat(basis_function_Z.pol, max_deg_basis_Z + 1);
	int max_deg_test_Z = getMax(test_function_Z.pol);
	MATRIX test_coeffs_Z = vertcat(test_function_Z.pol, max_deg_test_Z + 1);
	int max_deg_basis_S = getMax(basis_function_S.pol);
	MATRIX basis_coeffs_S = vertcat(basis_function_S.pol, max_deg_basis_S + 1);
	int max_deg_test_S = getMax(test_function_S.pol);
	MATRIX test_coeffs_S = vertcat(test_function_S.pol, max_deg_test_S + 1);
	int max_deg_basis_d = getMax(basis_function_d.pol);
	MATRIX basis_coeffs_d = vertcat(basis_function_d.pol, max_deg_basis_d + 1);
	int max_deg_test_d = getMax(test_function_d.pol);
	MATRIX test_coeffs_d = vertcat(test_function_d.pol, max_deg_test_d + 1);

	//Create Gaussian quadrature table to make ready to multiply with spatial basis functions
	max_basis_deg = max(max_deg_basis_Z, max_deg_basis_S);
	max_test_deg = max(max_deg_test_Z, max_deg_test_S);
	MATRIX Gcoeffs_ZZ = Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_Z);
	MATRIX Gcoeffs_ZS = Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_S);
	MATRIX Gcoeffs_SZ = Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_Z);
	MATRIX Gcoeffs_SS = Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_S);
	MATRIX Gcoeffs_dd = Gcoeffs_create(s0, w0, si, wi, max_deg_test_d, max_deg_basis_d);
	MATRIX Gcoeffs_ZZ_sp = Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_Z);
	MATRIX Gcoeffs_ZS_sp = Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_S);
	MATRIX Gcoeffs_SZ_sp = Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_Z);
	MATRIX Gcoeffs_SS_sp = Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_S);
	MATRIX Gcoeffs_dd_sp = Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_d, max_deg_basis_d);

	//Get index table to determine which edges interact
	umat idx_basis_Z = basis_function_Z.idx_table;
	umat idx_basis_S = basis_function_S.idx_table;
	umat idx_basis_d = basis_function_d.idx_table;
	umat idx_test_Z = test_function_Z.idx_table;
	umat idx_test_S = test_function_S.idx_table;
	umat idx_test_d = test_function_d.idx_table;
	int N_F = idx_basis_Z.n_rows;

	//Containers for operators (2D regions with 3rd dimension varying with time)
	S.zeros(N_F, N_F, z_N_T);
	D.zeros(N_F, N_F, z_N_T);
	Dp.zeros(N_F, N_F, z_N_T);
	Ns.zeros(N_F, N_F, z_N_T);
	Nh.zeros(N_F, N_F, z_N_T);

	//Pad the interpolators so the sizes match Nh
	CLagrange_interp tmp_Lag = timeBasis_Nh;
	timeBasis_Ns = tmp_Lag.pad_coeffs(timeBasis_Ns);
	timeBasis_D = tmp_Lag.pad_coeffs(timeBasis_D);
	tempconvs = CTempConvs();

	//Containers for operator coefficients for each edge
	MATRIX** coeffs_nh = CreateMatrix(N_E, N_E);
	MATRIX** coeffs_ns = CreateMatrix(N_E, N_E);
	MATRIX** coeffs_d = CreateMatrix(N_E, N_E);
	MATRIX** coeffs_s = CreateMatrix(N_E, N_E);
	MATRIX** coeffs_dp = CreateMatrix(N_E, N_E);

	//Variables to be thread-private
	int k(0), m(0), n(0), max_n(N_E), max_n2(N_F);
	POINT2D t_m, t_n;
	double l_m, l_n;
	CLagrange_interp shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns;

#pragma omp parallel default(shared) private(k,m,n,shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns,t_m,t_n,l_m,l_n)
	{
		//Loop over all segments so they all act as observation and source ( m and n) for all time steps
		for (k = 0; k < (int)z_N_T; k++){

			//Shifted time basis - shift and flip the time basis functions to get T(k*dt - t)
			shiftedTB_D = timeBasis_D;   shiftedTB_D.translate((k)* z_dt, -1);
			shiftedTB_Nh = timeBasis_Nh; shiftedTB_Nh.translate((k)* z_dt, -1);
			shiftedTB_Ns = timeBasis_Ns; shiftedTB_Ns.translate((k)* z_dt, -1);

#pragma omp for schedule(static)
			for (m = 0; m < (int)N_E; m++){

				//Find tangential vector at observation point
				t_m = z_geom_obj[m].t;
				l_m = z_geom_obj[m].l;

				cheat ? max_n = (m == 0) ? N_E : 1 : false;
				for (n = 0; n < max_n; n++){

					//Find tangential vector at source point
					t_n = z_geom_obj[n].t;
					l_n = z_geom_obj[n].l;

					double lmn(l_n*l_m);

					//When dealing with singularities at self patch and neighbouring edges, increase number of quadrature points
					if (n == (m + 1) % N_E || n == m || n == (m + N_E - 1) % N_E){

						//change si, s0, inner points, outer points, and Gcoeffs
						Z_calc(m, n, t_m, t_n, lmn, shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns,
							z_inner_points_sp, z_outer_points_sp,
							Gcoeffs_dd_sp, Gcoeffs_SZ_sp, Gcoeffs_ZS_sp, Gcoeffs_SS_sp, Gcoeffs_ZZ_sp,
							coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
					}
					else {

						//normal calculation when there are no singularities
						Z_calc(m, n, t_m, t_n, lmn, shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns,
							z_inner_points, z_outer_points,
							Gcoeffs_dd, Gcoeffs_SZ, Gcoeffs_ZS, Gcoeffs_SS, Gcoeffs_ZZ,
							coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
					}

				} // end n
			} // end m

#pragma omp single
		{
			if (cheat)
			{
				SPD_cheat_coeffs(coeffs_s, N_E, N_E);
				SPD_cheat_coeffs(coeffs_nh, N_E, N_E);
				SPD_cheat_coeffs(coeffs_ns, N_E, N_E);
				SPD_cheat_coeffs(coeffs_d, N_E, N_E);
				SPD_cheat_coeffs(coeffs_dp, N_E, N_E);
			}
		}

		// use integrated convolution values along with basis/test function coefficients
		// and index table to loop over basis functions and compute final matrix entries
#pragma omp for schedule(static)
		for (m = 0; m < (int)N_F; m++){

			//current test edges
			uvec index_test_Z = idx_test_Z.row(m).t();
			uvec index_test_S = idx_test_S.row(m).t();
			uvec index_test_d = idx_test_d.row(m).t();

			//initial function coefficients for current test edges
			uvec TFZ_start(index_test_Z * num_test_segments_Z);
			uvec TFS_start(index_test_S * num_test_segments_S);
			uvec TFd_start(index_test_d * num_test_segments_d);

			//create array of function coefficients
			MATRIX TCZ = create_function_coeffs_tbl(test_coeffs_Z, TFZ_start);
			MATRIX TCS = create_function_coeffs_tbl(test_coeffs_S, TFS_start);
			MATRIX TCd = create_function_coeffs_tbl(test_coeffs_d, TFd_start);

			cheat ? max_n2 = (m == 0) ? N_F : 1 : false;
			for (n = 0; n < max_n2; n++){

				//current test edges
				uvec index_basis_Z = idx_basis_Z.row(n).t();
				uvec index_basis_S = idx_basis_S.row(n).t();
				uvec index_basis_d = idx_basis_d.row(n).t();

				//initial function coefficients for current basis edges
				uvec BFZ_start(index_basis_Z * num_basis_segments_Z);
				uvec BFS_start(index_basis_S * num_basis_segments_S);
				uvec BFd_start(index_basis_d * num_basis_segments_d);

				//create array of function coefficients
				MATRIX BCZ = create_function_coeffs_tbl(basis_coeffs_Z, BFZ_start);
				MATRIX BCS = create_function_coeffs_tbl(basis_coeffs_S, BFS_start);
				MATRIX BCd = create_function_coeffs_tbl(basis_coeffs_d, BFd_start);

				//combine contributions
				Ns(m, n, k) = combine_contributions(coeffs_ns, index_test_S, index_basis_S, TCS, BCS);
				Nh(m, n, k) = combine_contributions(coeffs_nh, index_test_d, index_basis_d, TCd, BCd);
				S(m, n, k) = combine_contributions(coeffs_s, index_test_Z, index_basis_Z, TCZ, BCZ);
				Dp(m, n, k) = combine_contributions(coeffs_dp, index_test_S, index_basis_Z, TCS, BCZ);
				D(m, n, k) = combine_contributions(coeffs_d, index_test_Z, index_basis_S, TCZ, BCS);

			} // end n

		} // end m

		//#pragma omp flush

#pragma omp single
		{
			if (cheat)
			{
				SPD_cheat(S.slice(k));
				SPD_cheat(Nh.slice(k));
				SPD_cheat(Ns.slice(k));
				SPD_cheat(D.slice(k));
				SPD_cheat(Dp.slice(k));
			}

			//Status
			(k % status_int == 0) ? printf("%i ", k) : false;
			fflush(stdout);
		}

		} //end k
	}

	//Free memory
	deletepMat(coeffs_d, N_E);
	deletepMat(coeffs_s, N_E);
	deletepMat(coeffs_dp, N_E);
	deletepMat(coeffs_ns, N_E);
	deletepMat(coeffs_nh, N_E);

}

inline MATRIX Zmatrices::create_function_coeffs_tbl(MATRIX& coeff_tbl, uvec& start_elements)
{
	//start_elements.print("start");
	//coeff_tbl.print("coeff table");

	MATRIX coeffs = zeros(start_elements.n_elem, coeff_tbl.n_cols);
	for (int i = 0; i < (int)start_elements.n_elem; i++)
	{
		coeffs.row(i) = coeff_tbl.row(start_elements(i)+i);
	}

	//coeffs.print("coeffs");

	return coeffs;
}

void Zmatrices::lgquad1(UINT N, VECTOR& s, VECTOR& w)
{
	N -= 1;
	UINT N1 = N + 1;
	UINT N2 = N + 2;
	VECTOR xu = linspace(-1, 1, N1);

	//Initial guess
	VECTOR tmp = linspace(-0, N, N + 1);
	VECTOR y = cos((2 * tmp + 1) * PI / (2 * N + 2)) + (0.27 / N1) * sin(PI * xu * N / N2);

	// Legendre-Gauss Vandermonde Matrix
	MATRIX L(N1, N2);	 L.fill(0);

	// Derivative of LGVM
	MATRIX Lp(N1, N2);	 Lp.fill(0);

	// Compute the zeros of the N+1 Legendre Polynomial
	// using the recursion relation and the Newton-Raphson method
	VECTOR y0(y.size());
	y0.fill(2);

	//Iterate until new points are uniformly within epsilon of old points
	double eps(2.2204e-016);
	VECTOR _y = y - y0;
	VECTOR Lpp(N1);
	int c_cycle = 0;
	while (abs(max(_y)) > eps){
		if (c_cycle > 50)break;
		L.col(0).fill(1);
		Lp.col(0).fill(0);
		L.col(1) = y;
		Lp.col(1).fill(1);

		for (UINT k = 1; k < N1; k++){
			L.col(k + 1) = ((2 * k + 1) * (y % L.col(k)) - k * L.col(k - 1)) / (k + 1);
		}
		Lpp = N2 * (L.col(N1 - 1) - y % L.col(N2 - 1)) / (1 - square(y));
		y0 = y;
		y = y0 - L.col(N2 - 1) / Lpp;
		c_cycle++;
	}
	s = (1 + y) / 2;
	w = 1 / ((1 - square(y)) % square(Lpp)) * ((double)N2 / (double)N1) * ((double)N2 / (double)N1);
}

MATRIX ** Zmatrices::CreateMatrix(int i, int j)
{
	MATRIX **ppMat = (MATRIX**)new MATRIX*[j];
	for (int ii = 0; ii<i; ii++){
		ppMat[ii] = (MATRIX *)(calloc(sizeof(MATRIX) * j, 1));
	}
	return ppMat;
}

void Zmatrices::deletepMat(MATRIX** pMat, int count)
{
	for (int ii = 0; ii < count; ii++)
		delete pMat[ii];
	delete pMat;
}

double Zmatrices::combine_contributions(MATRIX** operator_coeffs, uvec& test_index, uvec& basis_index,
	MATRIX& test_coeffs, MATRIX& basis_coeffs)
{
	//loop through all interacting edges and multiply the basis and test polynomial coefficients with 
	//the integrated convolution ( operator ) coefficients
	double z(0);
	MATRIX tmp;

	/*printf("--------------------------------\n");
	test_index.print("test_index:");
	basis_index.print("basis_index:");
	test_coeffs.print("test_coeffs:");
	basis_coeffs.print("basis_coeffs:");//*/

	for (UINT alpha = 0; alpha < numel(test_index); alpha++)
	{
		for (UINT beta = 0; beta < numel(basis_index); beta++)
		{
			MATRIX operator_coeff = operator_coeffs[(int)test_index(alpha)][(int)basis_index(beta)];
			//operator_coeff.print("op coeff:");

			MATRIX polynomial_coeffs = test_coeffs.row(alpha).t() * basis_coeffs.row(beta);
			polynomial_coeffs.resize(max_test_deg+1, max_basis_deg+1);
			//polynomial_coeffs.print("polynomial_coeffs:");

			tmp = polynomial_coeffs % operator_coeff;

			z += accu(tmp);
		}
	}
	return z;
}

//CHEATS ONLY APPLICABLE FOR PEC CYLINDER (since Z matrix is SPD)
//Cheat to act on operator coeffs
void Zmatrices::SPD_cheat_coeffs(MATRIX** Z, UINT m, UINT n)
{
	//Using just the 1st row and column, copy elements diagonally down and right
	for (UINT i = 1; i < m; i++)
		for (UINT j = 1; j < n; j++)
			(Z[i][j]) = (Z[i - 1][j - 1]);
}
//Cheat to act on operators
void Zmatrices::SPD_cheat(MATRIX& Z)
{
	for (UINT i = 1; i < Z.n_cols; i++)
		for (UINT j = 1; j < Z.n_rows; j++)
			Z(i,j) = Z(i - 1, j - 1);
}
