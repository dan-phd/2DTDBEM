#pragma once
#include "Define.h"
#include "Lagrange_interp.h"

enum MYDATA_TYPE{
	MAT = 1,
	CUBE = 2,
	POLY = 3,
	EMPTY = 0
};

class MyData
{
public:
	MyData(){
		type = 0;
		m_mat.clear();
		m_cube.clear();
		m_poly.clear();
	};
	~MyData(){};

	mat m_mat;
	POLYCUBE m_poly;
	cube m_cube;
	int type;
	void getType(int n_type){
		if (n_type = MAT){ return; }

	}
};

class Zmatrices
{
public:
	//constructor
	Zmatrices(void);
	Zmatrices(UINT N_T, double dt, ZGEOMETRY& ZGegeometry, double c);
	~Zmatrices(void);

	//functions
	bool	 compute(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns, Zmatrices& z_Zmatrices);
	void	 lgquad1(UINT N, vec& s, vec& w);
	void	 Gcoeffs_create(vec& s0, vec& w0, vec&si, vec& wi, int max_deg_test, int max_deg_basis);
	void	 deletepMat(mat** pMat, int count);
	double	 getMax(vec& v);
	double	 getMax(POLYMAT& obj);
	double	 combine_contributions(mat** operator_coeffs, vec& test_index, vec& basis_index,
		mat& test_coeffs, mat& basis_coeffs);
	mat		 vertcat(POLYMAT& obj);
	mat		 bsxfun(int n_MODE, vec& bsxfun_A, vec& bsxfun_B);
	mat		 bsxfun(int n_MODE, vec& bsxfun_A, mat& bsxfun_B);
	mat		 bsxfun(int n_MODE, mat& bsxfun_A, vec& bsxfun_B);
	mat		 bsxfun(int n_MODE, POINT2D a, vec& s, POINT2D b);
	mat		 create_index_table(int num_segments, int num_copies_per_row, int num_copies_per_col);
	mat		 Perform_quadrature(MyData& m_data, mat& m_mat);
	mat **  CreateMatrix(int i, int j);

	//define operator- between both POINT2D structure
	POINT2D  minus(POINT2D a, POINT2D b);

	//properties:
	UINT z_N_T;
	double z_dt;
	double z_c;
	CLagrange_interp timeBasis_D;
	CLagrange_interp timeBasis_Nh;
	CLagrange_interp timeBasis_Ns;
	POLYMAT basis_function_Z;
	POLYMAT basis_function_S;
	POLYMAT test_function_Z;
	POLYMAT test_function_S;
	int	z_outer_points_sp;
	int z_inner_points_sp;
	int z_outer_points;
	int z_inner_points;
	MyData m_data;

private:

	//functions
	void	 Z_calc(Zmatrices& obj, vec& s_i, vec& s_o, MyData& G_dd, MyData& G_SZ, MyData& G_ZS, MyData& G_SS, MyData& G_ZZ, mat& coeffs_nh, mat& coeffs_ns, mat& coeffs_d, mat& coeffs_s, mat& coeffs_dp);

	//properties:
	POINT2D a_m, b_m, t_m, n_m, a_n, b_n, t_n, n_n;
	double l_m, l_n;
	ZGEOMETRY z_geom_obj;
	CLagrange_interp shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns;
	UINT N_E;
};
