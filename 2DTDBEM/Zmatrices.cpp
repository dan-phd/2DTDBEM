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

MATRIX Zmatrices::vertcat(POLYMAT& obj)
{
	MATRIX tmp = obj[0].m_coeffs;
	MATRIX res((const arma::uword)obj.size(), tmp.n_cols, fill::zeros);

	for (int i = 0; i < obj.size(); i++)
	{
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

void Zmatrices::Gcoeffs_create(VECTOR& s0, VECTOR& w0, VECTOR&si, VECTOR& wi, int max_deg_test, int max_deg_basis)
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
	m_data.type = 0;	m_data.m_cube.clear();	m_data.m_mat.clear();	m_data.m_poly.clear();

	if (max_deg_test == 0 && max_deg_basis == 0)
	{
		m_data.m_mat.reshape(test_points, basis_points);
		m_data.type = MAT;
		for (int i = 0; i < basis_points; i++)
		{
			m_data.m_mat.col(i) = Gcoeffs_test.col(0) * Gcoeffs_basis(i, 0);
		}
		return;
	}

	if (max_deg_test == 1 && max_deg_basis == 0)
	{
		m_data.m_cube.reshape(test_points, basis_points, max_deg_test + 1);
		m_data.type = CUBE;
		for (int i = 0; i < max_deg_test + 1; i++)
		{
			MATRIX tmp;		tmp.set_size(test_points, basis_points);
			for (int j = 0; j < basis_points; j++)
			{
				tmp.col(j) = Gcoeffs_test.col(i) * Gcoeffs_basis(j, 0);
			}
			m_data.m_cube.slice(i) = tmp;
		}
		return;
	}
	if (max_deg_test == -1 || max_deg_basis == -1)
	{
		m_data.m_mat.set_size(1, 1);
		m_data.m_mat.fill(0);
		m_data.type = EMPTY;
	}
	else{
		m_data.m_poly.clear();
		m_data.type = POLY;
		cube obj;
		obj.reshape(test_points, basis_points, max_deg_test + 1);

		for (int p = 0; p < max_deg_basis + 1; p++)
		{
			MATRIX tmp;		tmp.set_size(test_points, basis_points);
			for (int i = 0; i < max_deg_test + 1; i++)
			{
				for (int j = 0; j < basis_points; j++)
				{
					tmp.col(j) = Gcoeffs_test.col(i) * Gcoeffs_basis(j, p);
				}
				obj.slice(i) = tmp;
			}
			m_data.m_poly.push_back(obj);
		}
		return;
	}
}

// Create index table to determine which edges interact
MATRIX Zmatrices::create_index_table(int num_segments, int N_E)
{
	VECTOR tmp = linspace(1, num_segments, num_segments);
	MATRIX res = repmat(tmp.t(), N_E, 1);
	tmp = linspace(1, N_E, N_E);

	res.each_col() += tmp - num_segments - 1;
	res = res - floor(res/N_E)*N_E + 1;

	//res.print("idx:");

	return res;
}

MATRIX Zmatrices::Perform_quadrature(CustomData& m_data, MATRIX& m_mat)
{
	
	if (m_data.type == EMPTY)
	{
		MATRIX res(1, 1, fill::zeros);
		return res;
	}

	if (m_data.type == MAT)
	{
		MATRIX res(1, 1);
		res(0, 0) = accu(m_data.m_mat % m_mat);
		return res;
	}

	if (m_data.type == CUBE)
	{
		MATRIX res(m_data.m_cube.n_slices, 1);
		for (UINT i = 0; i < m_data.m_cube.n_slices; i++)
		{
			res(i) = sum(sum(m_data.m_cube.slice(i) % m_mat));
		}
		return res;
	}

	if (m_data.type == POLY)
	{
		MATRIX res(m_data.m_poly[0].n_slices, (const arma::uword)m_data.m_poly.size());
		for (UINT i = 0; i < m_data.m_poly[0].n_slices; i++){
			for (int j = 0; j < m_data.m_poly.size(); j++){
				res(i, j) = sum(sum(m_data.m_poly[j].slice(i) % m_mat));
			}
		}
		return res;
	}
}

//computation of single elements of all operator matrices
void Zmatrices::Z_calc(VECTOR& s_i, VECTOR& s_o, CustomData& G_dd, CustomData& G_SZ, CustomData& G_ZS, CustomData& G_SS, CustomData& G_ZZ,
	MATRIX& coeffs_nh, MATRIX& coeffs_ns, MATRIX& coeffs_d, MATRIX& coeffs_s, MATRIX& coeffs_dp)
{
	// Find inner and outer Gaussian quadrature points
	int inner_quad_points = numel(s_i);
	int outer_quad_points = numel(s_o);

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
	MATRIX tmp = t_rho_mn.slice(0) + t_rho_mn.slice(1);
	MATRIX P = sqrt(tmp);

	// Compute the temporal convolutions
	CTempConvs	tempconvs;
	VECTOR t_P = reshape(P, numel(P), 1) / z_c;
	VECTOR t_Fh(t_P.size(), fill::zeros), t_F(t_P.size(), fill::zeros), t_dF(t_P.size(), fill::zeros);

//#pragma omp parallel sections private(tempconvs,t_P)
	{
//#pragma omp section
		tempconvs.compute2(t_P, shiftedTB_Nh, shiftedTB_Ns, shiftedTB_D, t_Fh, t_F, t_dF);
	}

	MATRIX Fh = reshape(t_Fh, outer_quad_points, inner_quad_points);
	MATRIX F = reshape(t_F, outer_quad_points, inner_quad_points);
	MATRIX dF = reshape(t_dF, outer_quad_points, inner_quad_points) / z_c;

	// dF/dn = n dot P/|P| * dF
	cube unit_P = rho_mn;
	for (UINT i = 0; i < unit_P.n_slices; i++)
		unit_P.slice(i) = rho_mn.slice(i) / P;

	MATRIX dF_dnp = -(unit_P.slice(0) * n_n.x + unit_P.slice(1) * n_n.y) % dF;
	MATRIX dF_dn = (unit_P.slice(0) * n_m.x + unit_P.slice(1) * n_m.y) % dF;

	// Perform quadrature
	MATRIX nh_intF = Perform_quadrature(G_dd, Fh);
	MATRIX ns_intF = Perform_quadrature(G_SS, F);
	MATRIX dp_intF = Perform_quadrature(G_SZ, dF_dn);
	MATRIX d_intF = Perform_quadrature(G_ZS, dF_dnp);
	MATRIX s_intF = Perform_quadrature(G_ZZ, F);

	//scale using edge lengths
	coeffs_nh = nh_intF * l_n * l_m;
	coeffs_ns = ns_intF * l_n * l_m *(t_m.x*t_n.x + t_m.y*t_n.y);
	coeffs_d = d_intF * l_n * l_m;
	coeffs_s = s_intF* l_n * l_m;
	coeffs_dp = dp_intF * l_n * l_m;
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

	//Find divergence of transverse basis and test functions
	CBasisFunction	BasisFunction;
	POLYMAT	basis_function_d = basis_function_S;
	BasisFunction.divergence(basis_function_d, z_geom_obj);

	POLYMAT	test_function_d = test_function_S;
	BasisFunction.divergence(test_function_d, z_geom_obj);

	//Concatenate the basis and test coefficients and find the maximum degree
	POLYMAT tmp = basis_function_Z;
	int max_deg_basis_Z = getMax(tmp);
	MATRIX basis_coeffs_Z = vertcat(tmp);
	tmp = test_function_Z;
	int max_deg_test_Z = getMax(tmp);
	MATRIX test_coeffs_Z = vertcat(tmp);
	tmp = basis_function_S;
	int max_deg_basis_S = getMax(tmp);
	MATRIX basis_coeffs_S = vertcat(tmp);
	tmp = test_function_S;
	int max_deg_test_S = getMax(tmp);
	MATRIX test_coeffs_S = vertcat(tmp);
	int max_deg_basis_d = getMax(basis_function_d);
	MATRIX basis_coeffs_d = vertcat(basis_function_d);
	int max_deg_test_d = getMax(test_function_d);
	MATRIX test_coeffs_d = vertcat(test_function_d);

	//Number of basis and test segments for Z-directed, transverse and divergence fields
	UINT num_basis_segments_Z = (UINT)basis_function_Z.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_Z = (UINT)test_function_Z.size() / (UINT)z_geom_obj.size();
	UINT num_basis_segments_S = (UINT)basis_function_S.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_S = (UINT)test_function_S.size() / (UINT)z_geom_obj.size();
	UINT num_basis_segments_d = (UINT)basis_function_d.size() / (UINT)z_geom_obj.size();
	UINT num_test_segments_d = (UINT)basis_function_d.size() / (UINT)z_geom_obj.size();

	//make ready to multiply with spatial basis function
	CustomData Gcoeffs_ZZ, Gcoeffs_dd, Gcoeffs_ZZ_sp, Gcoeffs_dd_sp;		//	MATRIX ->0,0
	CustomData Gcoeffs_ZS, Gcoeffs_SS, Gcoeffs_ZS_sp, Gcoeffs_SS_sp;		//  poly->(0,1), (1,1)
	CustomData Gcoeffs_SZ, Gcoeffs_SZ_sp;									//	cube->1,0

	Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_Z);	Gcoeffs_ZZ = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_S);	Gcoeffs_ZS = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_Z);	Gcoeffs_SZ = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_S);	Gcoeffs_SS = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_d, max_deg_basis_d);	Gcoeffs_dd = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_Z);	Gcoeffs_ZZ_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_S);	Gcoeffs_ZS_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_Z);	Gcoeffs_SZ_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_S);	Gcoeffs_SS_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_d, max_deg_basis_d);	Gcoeffs_dd_sp = m_data;

	//Create index table to determine which edges interact
	MATRIX idx_basis_Z = create_index_table(num_basis_segments_Z, N_E);
	MATRIX idx_basis_S = create_index_table(num_basis_segments_S, N_E);
	MATRIX idx_basis_d = create_index_table(num_basis_segments_d, N_E);
	MATRIX idx_test_Z = create_index_table(num_test_segments_Z, N_E);
	MATRIX idx_test_S = create_index_table(num_test_segments_S, N_E);
	MATRIX idx_test_d = create_index_table(num_test_segments_d, N_E);

	//Containers for operators (2D regions with 3rd dimension varying with time
	S.resize(N_E, N_E, z_N_T);    S.fill(0);
	D.resize(N_E, N_E, z_N_T);    D.fill(0);
	Dp.resize(N_E, N_E, z_N_T);   Dp.fill(0);
	Ns.resize(N_E, N_E, z_N_T);   Ns.fill(0);
	Nh.resize(N_E, N_E, z_N_T);   Nh.fill(0);

	//Pad the interpolators so the sizes match Nh
	CLagrange_interp tmp_Lag = timeBasis_Nh;
	timeBasis_Ns = tmp_Lag.pad_coeffs(timeBasis_Ns);
	timeBasis_D = tmp_Lag.pad_coeffs(timeBasis_D);

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
	omp_set_num_threads(4);

	//Loop over all segments so they all act as observation and source ( m and n) for all time steps
	int k(0), m(0), n(0), max_n(N_E);
	for (k = 0; k < z_N_T; k++){

		//Shifted time basis - shift and flip the time basis functions to get T(k*dt - t)
		shiftedTB_D = timeBasis_D;
		shiftedTB_D.translate((k)* z_dt, -1);
		shiftedTB_Nh = timeBasis_Nh;
		shiftedTB_Nh.translate((k)* z_dt, -1);
		shiftedTB_Ns = timeBasis_Ns;
		shiftedTB_Ns.translate((k)* z_dt, -1);

		//Containers for operator coefficients
		MATRIX** coeffs_nh = CreateMatrix(N_E, N_E);
		MATRIX** coeffs_ns = CreateMatrix(N_E, N_E);
		MATRIX** coeffs_d = CreateMatrix(N_E, N_E);
		MATRIX** coeffs_s = CreateMatrix(N_E, N_E);
		MATRIX** coeffs_dp = CreateMatrix(N_E, N_E);

//#pragma omp parallel for default(shared) private(m,n) schedule(dynamic)
//#pragma omp parallel for private(n)
		for (m = 0; m < N_E; m++){

			//Find geometry around observation point
			EDGE tmp_edge = z_geom_obj[m];
			a_m = tmp_edge.a;
			b_m = tmp_edge.b;
			l_m = tmp_edge.l;
			t_m = tmp_edge.t;
			n_m = tmp_edge.n;

			cheat ? max_n = (m == 0) ? N_E : 1 : false;
			//Zmatrices zcopy = *this;
//#pragma omp parallel for default(shared) private(zcopy,tmp_edge)
			for (n = 0; n < max_n; n++){

				//Find geometry around source point
				tmp_edge = z_geom_obj[n];
				a_n = tmp_edge.a;
				b_n = tmp_edge.b;
				l_n = tmp_edge.l;
				t_n = tmp_edge.t;
				n_n = tmp_edge.n;

				//When dealing with singularities at self patch and neighbouring edges, increase number of quadrature points
				if (n == (m + 1) % N_E || n == m || n == (m + N_E - 1) % N_E){

					//change si, s0, inner points, outer points, and Gcoeffs
					Z_calc(si_sp, s0_sp, Gcoeffs_dd_sp, Gcoeffs_SZ_sp, Gcoeffs_ZS_sp, Gcoeffs_SS_sp, Gcoeffs_ZZ_sp,
						coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
				}
				else {

					//normal calculation when there are no singularities
					Z_calc(si, s0, Gcoeffs_dd, Gcoeffs_SZ, Gcoeffs_ZS, Gcoeffs_SS, Gcoeffs_ZZ,
						coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
				}

			} // end n
		} // end m

//#pragma omp flush

//#pragma omp single
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

		//coeffs_s[0][0].print("0:"); coeffs_s[1][1].print("1:");

		// use integrated convolution values along with basis/test function coefficients
		// and index table to compute final matrix entries
//#pragma omp parallel for
		for (m = 0; m < N_E; m++){

			//current test edges
			VECTOR index_test_Z = idx_test_Z.row(m).t();
			VECTOR test_v = idx_test_S.row(0).t();
			VECTOR index_test_S = idx_test_S.row(m).t();
			VECTOR index_test_d = idx_test_d.row(m).t();

			//current test function coefficients
			MATRIX TCZ, TCS, TCd;
			TCZ.set_size(num_test_segments_Z, test_coeffs_Z.n_cols);
			TCZ = test_coeffs_Z.rows(m * num_test_segments_Z, (m + 1) * num_test_segments_Z - 1);

			TCS.set_size(num_test_segments_S, test_coeffs_S.n_cols);
			TCS = test_coeffs_S.rows(m * num_test_segments_S, (m + 1) * num_test_segments_S - 1);

			TCd.set_size(num_test_segments_d, test_coeffs_d.n_cols);
			TCd = test_coeffs_d.rows(m * num_test_segments_d, (m + 1) * num_test_segments_d - 1);

			cheat ? max_n = (m == 0) ? N_E : 1 : false;
			for (n = 0; n < max_n; n++){
				//current basis edges
				VECTOR index_basis_Z = idx_basis_Z.row(n).t();
				VECTOR index_basis_S = idx_basis_S.row(n).t();
				VECTOR index_basis_d = idx_basis_d.row(n).t();
				//current basis function coefficients
				MATRIX BCZ, BCS, BCd;
				BCZ.set_size(num_basis_segments_Z, basis_coeffs_Z.n_cols);
				BCZ = basis_coeffs_Z.rows(n * num_basis_segments_Z, (n + 1) * num_basis_segments_Z - 1);

				BCS.set_size(num_basis_segments_S, basis_coeffs_S.n_cols);
				BCS = basis_coeffs_S.rows(n * num_basis_segments_S, (n + 1) * num_basis_segments_S - 1);

				BCd.set_size(num_basis_segments_d, basis_coeffs_d.n_cols);
				BCd = basis_coeffs_d.rows(n * num_basis_segments_d, (n + 1) * num_basis_segments_d - 1);

				//combine contributions
				Ns(m, n, k) = combine_contributions(coeffs_ns, index_test_S, index_basis_S, TCS, BCS);
				Nh(m, n, k) = combine_contributions(coeffs_nh, index_test_d, index_basis_d, TCd, BCd);
				S(m, n, k) = combine_contributions(coeffs_s, index_test_Z, index_basis_Z, TCZ, BCZ);
				Dp(m, n, k) = combine_contributions(coeffs_dp, index_test_S, index_basis_Z, TCS, BCZ);
				D(m, n, k) = combine_contributions(coeffs_d, index_test_Z, index_basis_S, TCZ, BCS);

			} // end n

		} // end m

//#pragma omp flush

//#pragma omp single
		{
			if (cheat)
			{
				SPD_cheat(S.slice(k));
				SPD_cheat(Nh.slice(k));
				SPD_cheat(Ns.slice(k));
				SPD_cheat(D.slice(k));
				SPD_cheat(Dp.slice(k));
			}
		}

//#pragma omp single
		{
			//Free memory
			deletepMat(coeffs_d, N_E);
			deletepMat(coeffs_s, N_E);
			deletepMat(coeffs_dp, N_E);
			deletepMat(coeffs_ns, N_E);
			deletepMat(coeffs_nh, N_E);

			//Status
			(k % status_int == 0) ? printf("%i ", k) : false;
			fflush(stdout);
		}
		
	} //end k

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

double Zmatrices::combine_contributions(MATRIX** operator_coeffs, VECTOR& test_index, VECTOR& basis_index,
	MATRIX& test_coeffs, MATRIX& basis_coeffs)
{
	//loop through all interacting edges and multiply the basis and test polynomial coefficients with 
	//the integrated convolution ( operator ) coefficients
	double z = 0;
	for (UINT alpha = 0; alpha < numel(test_index); alpha++){
		for (UINT beta = 0; beta < numel(basis_index); beta++){
			MATRIX operator_coeff = operator_coeffs[(int)test_index(alpha) - 1][(int)basis_index(beta) - 1];
			MATRIX polynomial_coeffs; polynomial_coeffs.set_size(test_coeffs.n_cols, basis_coeffs.n_cols);
			if (polynomial_coeffs.n_cols * polynomial_coeffs.n_rows == 0){ return 0; }
			for (UINT k = 0; k < test_coeffs.n_cols; k++){
				polynomial_coeffs.row(k) = test_coeffs(alpha, k) * basis_coeffs.row(beta);
			}
			z = z + accu(polynomial_coeffs % operator_coeff);
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
