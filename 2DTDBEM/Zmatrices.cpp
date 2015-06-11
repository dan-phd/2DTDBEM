#include "Zmatrices.h"
#include "Lagrange_interp.h"
#include "BasisFunction.h"
#include "TempConvs.h"

Zmatrices::Zmatrices(void)
{
}

Zmatrices::~Zmatrices(void)
{

}

//constructor
Zmatrices::Zmatrices(UINT N_T, double dt, ZGEOMETRY& ZGegeometry, double c)
{
	z_N_T = N_T;
	z_dt = dt;
	z_geom_obj = ZGegeometry;
	z_c = c;

	z_outer_points_sp = 50;
	z_inner_points_sp = 51;
	z_outer_points = 3;
	z_inner_points = 4;
	N_E = ZGegeometry.size();
}

double Zmatrices::getMax(vec& v)
{
	double res = v(0);
	for (int i = 0; i < v.size(); i++){
		if (res < v(i))
			res = v(i);
	}
	return res;
}

mat	 Zmatrices::vertcat(POLYMAT& obj){
	mat tmp = obj[0].m_coeffs;
	mat res = tmp;
	res.set_size(obj.size(), tmp.n_cols);
	res.fill(0);
	for (int i = 0; i < obj.size(); i++){
		res.row(i) = obj[i].m_coeffs;
	}
	return res;
}

double Zmatrices::getMax(POLYMAT& obj){
	double res = obj[0].m_degree;
	for (int i = 1; i < obj.size(); i++){
		if (res < obj[i].m_degree)
			res = obj[i].m_degree;
	}
	return res;
}

mat	Zmatrices::bsxfun(int n_MODE, vec& bsxfun_A, vec& bsxfun_B)
{
	mat res;
	res.set_size(bsxfun_A.n_rows, bsxfun_B.size());
	if (n_MODE == POWER)
	{
		for (int i = 0; i < bsxfun_B.size(); i++){
			res.col(i) = pow(bsxfun_A, bsxfun_B(i));
		}
	}

	return res;
}

mat	Zmatrices::bsxfun(int n_MODE, vec& bsxfun_A, mat& bsxfun_B)
{
	mat res;
	if (n_MODE == TIMES)
	{
		res.set_size(bsxfun_B.n_rows, bsxfun_B.n_cols);
		for (int i = 0; i < bsxfun_B.n_cols; i++){
			res.col(i) = bsxfun_A % bsxfun_B.col(i);
		}
	}

	return res;
}

mat	Zmatrices::bsxfun(int n_MODE, mat& bsxfun_A, vec& bsxfun_B)
{
	if (n_MODE == PLUS)
	{
		for (int i = 0; i < bsxfun_A.n_cols; i++){
			bsxfun_A.col(i) = bsxfun_A.col(i) + bsxfun_B;
		}
	}
	return bsxfun_A;
}

void Zmatrices::Gcoeffs_create(vec& s0, vec& w0, vec&si, vec& wi, int max_deg_test, int max_deg_basis)
{
	//Multiply quadrature points up to max polynomial degree				////	mat ->0,0
	vec tmp = max_deg_test - linspace(0, max_deg_test, max_deg_test + 1);	////    poly->(0,1), (1,1)
	mat Gcoeffs_test = bsxfun(POWER, s0, tmp);								////	cube->1,0

	tmp = max_deg_basis - linspace(0, max_deg_basis, max_deg_basis + 1);
	mat Gcoeffs_basis = bsxfun(POWER, si, tmp);

	//Multiply with quadrature weights
	Gcoeffs_test = bsxfun(TIMES, w0, Gcoeffs_test);
	Gcoeffs_basis = bsxfun(TIMES, wi, Gcoeffs_basis);
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
			mat tmp;		tmp.set_size(test_points, basis_points);
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
			mat tmp;		tmp.set_size(test_points, basis_points);
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

mat Zmatrices::create_index_table(int num_segments, int num_copies_per_row, int num_copies_per_col){
	mat res, res1;
	vec tmp = linspace(1, num_segments, num_segments);
	res1 = repmat(tmp, num_copies_per_col, num_copies_per_row);
	vec tmp1 = linspace(1, num_copies_per_row, num_copies_per_row);
	res = res1.t();
	res = bsxfun(PLUS, res, tmp1) - num_segments - 1;
	for (int row = 0; row < res.n_rows; row++)
		for (int col = 0; col < res.n_cols; col++)
			if (res(row, col) < 0){ res(row, col) = res(row, col) + num_copies_per_row; }

	res = res + 1;
	return res;
}

mat Zmatrices::Perform_quadrature(MyData& m_data, mat& m_mat)
{
	mat res;
	if (m_data.type == EMPTY)
	{
		res.set_size(1, 1);
		res.fill(0);
		return res;
	}

	if (m_data.type == MAT)
	{
		res.set_size(1, 1);
		res(0, 0) = accu(m_data.m_mat % m_mat);
		return res;
	}

	if (m_data.type == CUBE)
	{
		res.set_size(m_data.m_cube.n_slices, 1);
		for (int i = 0; i < m_data.m_cube.n_slices; i++)
		{
			res(i) = sum(sum(m_data.m_cube.slice(i) % m_mat));
		}
		return res;
	}

	if (m_data.type == POLY)
	{
		res.set_size(m_data.m_poly[0].n_slices, m_data.m_poly.size());
		for (int i = 0; i < m_data.m_poly[0].n_slices; i++){
			for (int j = 0; j < m_data.m_poly.size(); j++){
				res(i, j) = sum(sum(m_data.m_poly[j].slice(i) % m_mat));
			}
		}
		return res;
	}
}

//computation of single elements of all operator matrices
void Zmatrices::Z_calc(Zmatrices& obj, vec& s_i, vec& s_o, MyData& G_dd, MyData& G_SZ, MyData& G_ZS, MyData& G_SS, MyData& G_ZZ,
	mat& coeffs_nh, mat& coeffs_ns, mat& coeffs_d, mat& coeffs_s, mat& coeffs_dp)
{
	// Find inner and outer Gaussian quadrature points
	int inner_quad_points = numel(s_i);
	int outer_quad_points = numel(s_o);

	// Make rho_m a vector of outer Gaussian quadrature points along observation edge
	mat rho_m = bsxfun(PLUS, obj.a_m, s_o, minus(obj.b_m, obj.a_m)); //rho_m.print();

	// Make rho_n a vector of inner Gaussian quadrature points along source edge
	mat rho_n = bsxfun(PLUS, obj.a_n, s_i, minus(obj.b_n, obj.a_n)); //rho_n.print();

	// P is an array of distances from rho_n to observation point, rho_m
	cube rho_mn;
	rho_mn.set_size(outer_quad_points, inner_quad_points, 2);
	for (int i = 0; i < 2; i++)///loop slices
	{
		mat tmp;		tmp.set_size(outer_quad_points, inner_quad_points);
		for (int j = 0; j < inner_quad_points; j++)
		{
			tmp.col(j) = rho_m.col(i) - rho_n(j, i);
		}
		rho_mn.slice(i) = tmp;
	}
	cube t_rho_mn = square(rho_mn);
	mat tmp = t_rho_mn.slice(0) + t_rho_mn.slice(1);
	mat P = sqrt(tmp);

	// Compute the temporal convolutions
	vec t_Fh, t_F, t_dF;
	CTempConvs	tempconvs;
	vec t_P = reshape(P, numel(P), 1) / obj.z_c;
	tempconvs.TempConvs(t_P, obj.shiftedTB_Nh, obj.shiftedTB_Ns, obj.shiftedTB_D, t_Fh, t_F, t_dF);

	mat Fh = reshape(t_Fh, outer_quad_points, inner_quad_points);
	mat F = reshape(t_F, outer_quad_points, inner_quad_points);
	mat dF = reshape(t_dF, outer_quad_points, inner_quad_points) / obj.z_c;

	// dF/dn = n dot P/|P| * dF
	cube unit_P = rho_mn;
	for (int i = 0; i < unit_P.n_slices; i++)
		unit_P.slice(i) = rho_mn.slice(i) / P;		// P / |P| : unit_P = bsxfun(@rdivide, rho_mn, P);    

	mat dF_dnp = -(unit_P.slice(0) * obj.n_n.x + unit_P.slice(1) * obj.n_n.y) % dF;
	mat dF_dn = (unit_P.slice(0) * obj.n_m.x + unit_P.slice(1) * obj.n_m.y) % dF;

	//	Alternatively:
	//	rdgt = bsxfun(@times, rho_mn, dF ./ P);                      % dF * P/|P|
	//	dF_dnp = sum(bsxfun(@times, reshape(-n_n,[1,1,2]), rdgt),3);
	//	dF_dn = sum(bsxfun(@times, reshape(n_m,[1,1,2]), rdgt),3);

	// Perform quadrature
	mat nh_intF = Perform_quadrature(G_dd, Fh);
	mat ns_intF = Perform_quadrature(G_SS, F);
	mat dp_intF = Perform_quadrature(G_SZ, dF_dn);
	mat d_intF = Perform_quadrature(G_ZS, dF_dnp);
	mat s_intF = Perform_quadrature(G_ZZ, F);

	//scale using edge lengths
	coeffs_nh = nh_intF * obj.l_n * obj.l_m;
	coeffs_ns = ns_intF * obj.l_n * obj.l_m *(obj.t_m.x*obj.t_n.x + obj.t_m.y*obj.t_n.y);
	coeffs_d = d_intF * obj.l_n * obj.l_m;
	coeffs_s = s_intF* obj.l_n * obj.l_m;
	coeffs_dp = dp_intF * obj.l_n * obj.l_m;
}

bool Zmatrices::compute(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns, Zmatrices& z_Zmatrices)
{
	//parameterise each segment for Gaussian quadrature integration
	vec s0, w0;
	lgquad1(z_Zmatrices.z_outer_points, s0, w0);
	vec si, wi;
	lgquad1(z_Zmatrices.z_inner_points, si, wi);

	//Gaussian quadrature integration for singularity points
	vec s0_sp, w0_sp;
	lgquad1(z_Zmatrices.z_outer_points_sp, s0_sp, w0_sp);
	vec si_sp, wi_sp;
	lgquad1(z_Zmatrices.z_inner_points_sp, si_sp, wi_sp);

	//Find divergence of transverse basis and test functions
	CBasisFunction	BasisFunction;
	POLYMAT	basis_function_d = z_Zmatrices.basis_function_S;
	BasisFunction.divergence(basis_function_d, z_Zmatrices.z_geom_obj);

	POLYMAT	test_function_d = z_Zmatrices.test_function_S;
	BasisFunction.divergence(test_function_d, z_Zmatrices.z_geom_obj);
	//Concatenate the basis and test coefficients and find the maximum degree
	POLYMAT tmp = z_Zmatrices.basis_function_Z;
	double max_deg_basis_Z = getMax(tmp);
	mat basis_coeffs_Z = vertcat(tmp);

	tmp = z_Zmatrices.test_function_Z;
	double max_deg_test_Z = getMax(tmp);
	mat test_coeffs_Z = vertcat(tmp);

	tmp = z_Zmatrices.basis_function_S;
	double max_deg_basis_S = getMax(tmp);
	mat basis_coeffs_S = vertcat(tmp);

	tmp = z_Zmatrices.test_function_S;
	double max_deg_test_S = getMax(tmp);
	mat test_coeffs_S = vertcat(tmp);

	double max_deg_basis_d = getMax(basis_function_d);
	mat basis_coeffs_d = vertcat(basis_function_d);

	double max_deg_test_d = getMax(test_function_d);
	mat test_coeffs_d = vertcat(test_function_d);
	//Number of basis and test segments for Z-directed, transverse and divergence fields
	UINT num_basis_segments_Z = z_Zmatrices.basis_function_Z.size() / z_geom_obj.size();
	UINT num_test_segments_Z = z_Zmatrices.test_function_Z.size() / z_geom_obj.size();

	UINT num_basis_segments_S = z_Zmatrices.basis_function_S.size() / z_geom_obj.size();
	UINT num_test_segments_S = z_Zmatrices.test_function_S.size() / z_geom_obj.size();

	UINT num_basis_segments_d = basis_function_d.size() / z_geom_obj.size();
	UINT num_test_segments_d = basis_function_d.size() / z_geom_obj.size();
	//make ready to multiply with spatial basis function
	MyData Gcoeffs_ZZ, Gcoeffs_dd, Gcoeffs_ZZ_sp, Gcoeffs_dd_sp;			////	mat ->0,0
	MyData Gcoeffs_ZS, Gcoeffs_SS, Gcoeffs_ZS_sp, Gcoeffs_SS_sp;		////    poly->(0,1), (1,1)
	MyData Gcoeffs_SZ, Gcoeffs_SZ_sp;										////	cube->1,0

	Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_Z);				 	Gcoeffs_ZZ = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_Z, max_deg_basis_S);				 	Gcoeffs_ZS = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_Z);				 	Gcoeffs_SZ = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_S, max_deg_basis_S);				 	Gcoeffs_SS = m_data;
	Gcoeffs_create(s0, w0, si, wi, max_deg_test_d, max_deg_basis_d);				 	Gcoeffs_dd = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_Z); 	Gcoeffs_ZZ_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_Z, max_deg_basis_S); 	Gcoeffs_ZS_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_Z); 	Gcoeffs_SZ_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_S, max_deg_basis_S); 	Gcoeffs_SS_sp = m_data;
	Gcoeffs_create(s0_sp, w0_sp, si_sp, wi_sp, max_deg_test_d, max_deg_basis_d); 	Gcoeffs_dd_sp = m_data;
	//Create index table to determine which edges interact
	mat idx_basis_Z = create_index_table(num_basis_segments_Z, z_Zmatrices.N_E, 1);
	mat idx_basis_S = create_index_table(num_basis_segments_S, z_Zmatrices.N_E, 1);
	mat idx_basis_d = create_index_table(num_basis_segments_d, z_Zmatrices.N_E, 1);
	mat idx_test_Z = create_index_table(num_test_segments_Z, z_Zmatrices.N_E, 1);
	mat idx_test_S = create_index_table(num_test_segments_S, z_Zmatrices.N_E, 1);
	mat idx_test_d = create_index_table(num_test_segments_d, z_Zmatrices.N_E, 1);
	//Containers for operators (2D regions with 3rd dimension varying with time
	S.resize(z_Zmatrices.N_E, z_Zmatrices.N_E, z_Zmatrices.z_N_T);    S.fill(0);
	D.resize(z_Zmatrices.N_E, z_Zmatrices.N_E, z_Zmatrices.z_N_T);    D.fill(0);
	Dp.resize(z_Zmatrices.N_E, z_Zmatrices.N_E, z_Zmatrices.z_N_T);   Dp.fill(0);
	Ns.resize(z_Zmatrices.N_E, z_Zmatrices.N_E, z_Zmatrices.z_N_T);   Ns.fill(0);
	Nh.resize(z_Zmatrices.N_E, z_Zmatrices.N_E, z_Zmatrices.z_N_T);   Nh.fill(0);
	//Pad the interpolators so the sizes mathc Nh
	CLagrange_interp tmp_Lag = z_Zmatrices.timeBasis_Nh;
	z_Zmatrices.timeBasis_Ns = tmp_Lag.pad_coeffs(z_Zmatrices.timeBasis_Ns);
	z_Zmatrices.timeBasis_D = tmp_Lag.pad_coeffs(z_Zmatrices.timeBasis_D);

	//Loop over all segments so they all act as observation and source ( m and n) for all time steps
	///set status bar
	for (int k = 0; k < z_Zmatrices.z_N_T; k++){
		//Shifted time basis - shift and flip the time basis functions to get T(k*dt - t)
		z_Zmatrices.shiftedTB_D = z_Zmatrices.timeBasis_D;
		z_Zmatrices.shiftedTB_D.translate((k)* z_Zmatrices.z_dt, -1);

		z_Zmatrices.shiftedTB_Nh = z_Zmatrices.timeBasis_Nh;
		z_Zmatrices.shiftedTB_Nh.translate((k)* z_Zmatrices.z_dt, -1);

		z_Zmatrices.shiftedTB_Ns = z_Zmatrices.timeBasis_Ns;
		z_Zmatrices.shiftedTB_Ns.translate((k)* z_Zmatrices.z_dt, -1);

		//Containers for operator coefficients
		mat** coeffs_nh = CreateMatrix(z_Zmatrices.N_E, z_Zmatrices.N_E);
		mat** coeffs_ns = CreateMatrix(z_Zmatrices.N_E, z_Zmatrices.N_E);
		mat** coeffs_d = CreateMatrix(z_Zmatrices.N_E, z_Zmatrices.N_E);
		mat** coeffs_s = CreateMatrix(z_Zmatrices.N_E, z_Zmatrices.N_E);
		mat** coeffs_dp = CreateMatrix(z_Zmatrices.N_E, z_Zmatrices.N_E);
		//		deletepMat(coeffs_d, z_Zmatrices.N_E);

		for (int m = 0; m < z_Zmatrices.N_E; m++){
			//Find geometry around observation point
			ZEDGE tmp_edge = z_Zmatrices.z_geom_obj[m];
			z_Zmatrices.a_m = tmp_edge.a;
			z_Zmatrices.b_m = tmp_edge.b;
			z_Zmatrices.l_m = tmp_edge.l;
			z_Zmatrices.t_m = tmp_edge.t;
			z_Zmatrices.n_m = tmp_edge.n;

			for (int n = 0; n < z_Zmatrices.N_E; n++){
				//Find geometry around source point
				ZEDGE tmp_edge = z_Zmatrices.z_geom_obj[n];
				z_Zmatrices.a_n = tmp_edge.a;
				z_Zmatrices.b_n = tmp_edge.b;
				z_Zmatrices.l_n = tmp_edge.l;
				z_Zmatrices.t_n = tmp_edge.t;
				z_Zmatrices.n_n = tmp_edge.n;

				//When dealing with singularities at self patch and neighbouring
				//edges, increase number of quadrature points or use split
				//integral routine
				if (n == (m + 1) % z_Zmatrices.N_E || n == m || n == (m + z_Zmatrices.N_E - 1) % z_Zmatrices.N_E){
					//change si, s0, inner points, outer points, and Gcoeffs
					Z_calc(z_Zmatrices, si_sp, s0_sp, Gcoeffs_dd_sp, Gcoeffs_SZ_sp, Gcoeffs_ZS_sp, Gcoeffs_SS_sp, Gcoeffs_ZZ_sp,
						coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
				}
				else{
					//normal calculation when there are no singularities
					Z_calc(z_Zmatrices, si, s0, Gcoeffs_dd, Gcoeffs_SZ, Gcoeffs_ZS, Gcoeffs_SS, Gcoeffs_ZZ,
						coeffs_nh[m][n], coeffs_ns[m][n], coeffs_d[m][n], coeffs_s[m][n], coeffs_dp[m][n]);
				}
			}//  end n

		}// end m

		//		use integrated convolution values along with basis/test function coefficients
		//		and index table to compute final matrix entries
		for (int m = 0; m < z_Zmatrices.N_E; m++){
			//current test edges
			vec index_test_Z = idx_test_Z.row(m).t();
			vec test_v = idx_test_S.row(0).t();
			vec index_test_S = idx_test_S.row(m).t();
			vec index_test_d = idx_test_d.row(m).t();

			//current test function coefficients
			mat TCZ, TCS, TCd;
			TCZ.set_size(num_test_segments_Z, test_coeffs_Z.n_cols);
			TCZ = test_coeffs_Z.rows(m * num_test_segments_Z, (m + 1) * num_test_segments_Z - 1);

			TCS.set_size(num_test_segments_S, test_coeffs_S.n_cols);
			TCS = test_coeffs_S.rows(m * num_test_segments_S, (m + 1) * num_test_segments_S - 1);

			TCd.set_size(num_test_segments_d, test_coeffs_d.n_cols);
			TCd = test_coeffs_d.rows(m * num_test_segments_d, (m + 1) * num_test_segments_d - 1);

			for (int n = 0; n < z_Zmatrices.N_E; n++){
				//current basis edges
				vec index_basis_Z = idx_basis_Z.row(n).t();
				vec index_basis_S = idx_basis_S.row(n).t();
				vec index_basis_d = idx_basis_d.row(n).t();
				//current basis function coefficients
				mat BCZ, BCS, BCd;
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

			}// end n

		}// end m
		deletepMat(coeffs_d, z_Zmatrices.N_E);
		deletepMat(coeffs_s, z_Zmatrices.N_E);
		deletepMat(coeffs_dp, z_Zmatrices.N_E);
		deletepMat(coeffs_ns, z_Zmatrices.N_E);
		deletepMat(coeffs_nh, z_Zmatrices.N_E);

		// indicator
		(k % 100 == 0) ? printf("%i ", k) : false;
		//cout << "Cycle K = " << k << ";" << endl;
	}//end k

	return true;
}

void Zmatrices::lgquad1(UINT N, vec& s, vec& w)
{
	N = N - 1;
	UINT N1 = N + 1;
	UINT N2 = N + 2;
	vec xu = linspace(-1, 1, N1);

	//Initial guess
	vec tmp = linspace(-0, N, N + 1);
	vec y = cos((2 * tmp + 1) * PI / (2 * N + 2)) + (0.27 / N1) * sin(PI * xu * N / N2);

	// Legendre-Gauss Vandermonde Matrix
	mat L(N1, N2);	 L.fill(0);

	// Derivative of LGVM
	mat Lp(N1, N2);	 Lp.fill(0);

	// Compute the zeros of the N+1 Legendre Polynomial
	// using the recursion relation and the Newton-Raphson method
	vec y0(y.size());
	y0.fill(2);

	//Iterate until new points are uniformly within epsilon of old points
	vec _y = y - y0;
	vec Lpp(N1);
	int c_cycle = 0;
	while (abs(getMax(_y)) > eps){
		if (c_cycle > 50)break;
		L.col(0).fill(1);
		Lp.col(0).fill(0);
		L.col(1) = y;
		Lp.col(1).fill(1);

		for (int k = 1; k < N1; k++){
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
mat ** Zmatrices::CreateMatrix(int i, int j)
{
	//mat **ppMat = (mat**)new mat*[j];
	//for( int ii = 0; ii<i; ii++ )
	//	ppMat[ii] = new mat[j];
	//return ppMat;

	mat **ppMat = (mat**)new mat*[j];
	for (int ii = 0; ii<i; ii++){
		ppMat[ii] = (mat *)(calloc(sizeof(mat) * j, 1));
	}
	return ppMat;
}

void Zmatrices::deletepMat(mat** pMat, int count)
{
	for (int ii = 0; ii < count; ii++)
		delete pMat[ii];
	delete pMat;
}

POINT2D  Zmatrices::minus(POINT2D a, POINT2D b)
{
	POINT2D res;
	res.x = a.x - b.x;
	res.y = a.y - b.y;
	return res;
}

mat	Zmatrices::bsxfun(int MODE, POINT2D a, vec& s, POINT2D b)
{
	if (MODE == PLUS){
		mat res; res.set_size(s.size(), 2);
		res.col(0) = s * b.x + a.x;
		res.col(1) = s * b.y + a.y;
		return res;
	}
}

double Zmatrices::combine_contributions(mat** operator_coeffs, vec& test_index, vec& basis_index,
	mat& test_coeffs, mat& basis_coeffs)
{
	//loop through all interacting edges and multiply the basis and test polynomial coefficients with 
	//the integrated convolution ( operator ) coefficients
	double z = 0;
	for (int alpha = 0; alpha < numel(test_index); alpha++){
		for (int beta = 0; beta < numel(basis_index); beta++){
			mat operator_coeff = operator_coeffs[(int)test_index(alpha) - 1][(int)basis_index(beta) - 1];
			mat polynomial_coeffs; polynomial_coeffs.set_size(test_coeffs.n_cols, basis_coeffs.n_cols);
			if (polynomial_coeffs.n_cols * polynomial_coeffs.n_rows == 0){ return 0; }
			for (int k = 0; k < test_coeffs.n_cols; k++){
				polynomial_coeffs.row(k) = test_coeffs(alpha, k) * basis_coeffs.row(beta);
			}
			z = z + accu(polynomial_coeffs % operator_coeff);
		}
	}
	return z;
}
