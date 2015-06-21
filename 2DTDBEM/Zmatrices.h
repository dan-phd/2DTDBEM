#pragma once
#include "Define.h"
#include "Lagrange_interp.h"

enum CUSTOMDATA_TYPE{
	MAT = 1,
	CUBE = 2,
	POLY = 3,
	EMPTY = 0
};

class CustomData
{
public:
	CustomData(){
		type = 0;
		m_mat.clear();
		m_cube.clear();
		m_poly.clear();
	};
	~CustomData(){};

	MATRIX m_mat;
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
	Zmatrices(UINT N_T, double dt, GEOMETRY& ZGegeometry, double c);
	~Zmatrices(void);

	//functions
	void	compute(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns);
	void	use_cheat();

	//properties:
	UINT z_N_T, status_int;
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
	CustomData m_data;

private:

	//functions
	void	Z_calc(VECTOR& s_i, VECTOR& s_o,
				CustomData& G_dd, CustomData& G_SZ, CustomData& G_ZS, CustomData& G_SS, CustomData& G_ZZ,
				MATRIX& coeffs_nh, MATRIX& coeffs_ns, MATRIX& coeffs_d, MATRIX& coeffs_s, MATRIX& coeffs_dp);
	void	SPD_cheat_coeffs(MATRIX** Z, UINT m, UINT n);
	void	SPD_cheat(MATRIX& Z);
	void	lgquad1(UINT N, VECTOR& s, VECTOR& w);
	void	Gcoeffs_create(VECTOR& s0, VECTOR& w0, VECTOR&si, VECTOR& wi, int max_deg_test, int max_deg_basis);
	void	deletepMat(MATRIX** pMat, int count);
	int		getMax(POLYMAT& obj);
	double	combine_contributions(MATRIX** operator_coeffs, VECTOR& test_index, VECTOR& basis_index,
				MATRIX& test_coeffs, MATRIX& basis_coeffs);
	MATRIX	vertcat(POLYMAT& obj);
	MATRIX	bsxPower(VECTOR& bsxfun_A, VECTOR& bsxfun_B);
	MATRIX	bsxPlus(POINT2D a, POINT2D b, VECTOR& s);
	MATRIX	create_index_table(int num_segments, int N_E);
	MATRIX	Perform_quadrature(CustomData& m_data, MATRIX& m_mat);
	MATRIX** CreateMatrix(int i, int j);

	//properties:
	bool cheat;
	POINT2D a_m, b_m, t_m, n_m, a_n, b_n, t_n, n_n;
	double l_m, l_n;
	GEOMETRY z_geom_obj;
	CLagrange_interp shiftedTB_D, shiftedTB_Nh, shiftedTB_Ns;
	UINT N_E;
};
