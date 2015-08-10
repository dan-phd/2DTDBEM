#pragma once
#include "Define.h"
#include "BasisFunction.h"
#include "Lagrange_interp.h"
#include "TempConvs.h"

class Zmatrices
{
public:
	//constructors
	Zmatrices(void);
	Zmatrices(UINT N_T, double dt, GEOMETRY& ZGegeometry, double c);
	~Zmatrices(void);

	//functions
	void	compute(cube& S, cube& D, cube& Dp, cube& Nh, cube& Ns);
	void	use_cheat();

	//properties:
	UINT z_N_T;
	double z_dt;
	double z_c;
	CLagrange_interp timeBasis_D;
	CLagrange_interp timeBasis_Nh;
	CLagrange_interp timeBasis_Ns;
	CBasisFunction basis_function_Z;
	CBasisFunction basis_function_S;
	CBasisFunction test_function_Z;
	CBasisFunction test_function_S;
	int	z_outer_points_sp;
	int z_inner_points_sp;
	int z_outer_points;
	int z_inner_points;

private:

	//functions
	void	Z_calc(const UINT& m, const UINT& n, const POINT2D& t_m, const POINT2D& t_n, const double& lmn,
					CLagrange_interp shiftedTB_D, CLagrange_interp shiftedTB_Nh, CLagrange_interp shiftedTB_Ns,
					const UINT& inner_quad_points, const UINT& outer_quad_points,
					MATRIX& G_dd, MATRIX& G_SZ, MATRIX& G_ZS, MATRIX& G_SS, MATRIX& G_ZZ,
					MATRIX& coeffs_nh, MATRIX& coeffs_ns, MATRIX& coeffs_d, MATRIX& coeffs_s, MATRIX& coeffs_dp);
	void	SPD_cheat_coeffs(MATRIX** Z, UINT m, UINT n);
	void	SPD_cheat(MATRIX& Z);
	void	lgquad1(UINT N, VECTOR& s, VECTOR& w);
	void	deletepMat(MATRIX** pMat, int count);
	void	find_distances_between_2_edges(UINT m, UINT n, MATRIX& P, MATRIX& dn_p, MATRIX& dn_,
					VECTOR& s_i, VECTOR& s_o);
	void	make_distances_lookup_table();
	int		getMax(POLYMAT& obj);
	double	combine_contributions(MATRIX** operator_coeffs, uvec& test_index, uvec& basis_index,
				MATRIX& test_coeffs, MATRIX& basis_coeffs);
	MATRIX	vertcat(POLYMAT& obj, UINT num_segments);
	MATRIX	bsxPower(VECTOR& bsxfun_A, VECTOR& bsxfun_B);
	MATRIX	bsxPlus(POINT2D a, POINT2D b, VECTOR& s);
	MATRIX	Perform_quadrature(MATRIX& G_coeffs, MATRIX& F);
	MATRIX** CreateMatrix(int i, int j);
	MATRIX	Gcoeffs_create(VECTOR& s0, VECTOR& w0, VECTOR& si, VECTOR& wi, int max_deg_test, int max_deg_basis);
	inline MATRIX create_function_coeffs_tbl(MATRIX& coeff_tbl, uvec& start_elements);

	//properties:
	bool cheat;
	GEOMETRY z_geom_obj;
	CTempConvs tempconvs;
	UINT N_E, max_basis_deg, max_test_deg;
	MATRIX **distances, **over_dnp, **over_dn;
};
