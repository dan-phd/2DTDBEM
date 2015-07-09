#include "BasisFunction.h"

CBasisFunction::CBasisFunction(void)
{
}

CBasisFunction::~CBasisFunction(void)
{
}

//Create index table to determine which edges interact
umat create_index_table(UINT N_E, UINT N_S, UINT start, UINT iterate)
{
	VECTOR a = linspace(1, N_E, N_E)-2 + start +1;
	umat b(N_E, N_S, fill::zeros);

	for (int i = 0; i < (int)N_S; i++)
	{
		for (int j = 0; j < (int)N_E; j++)
		{
			b(j, i) = (int)(a(j)+i) % N_E;
		}
	}

	int num_functions = (int)(N_E / iterate);
	umat res(num_functions, N_S, fill::zeros);
	for (int i = 0; i < (int)N_E; i += iterate)
	{
		res.row(i/iterate) = b.row(i);
	}

	//res.print("res:");

	return res;
}

//Create hat basis, "scale" set to true scales by 1/length
CBasisFunction CBasisFunction::createHat(GEOMETRY& geometry, bool scale)
{
	UINT N_E = (UINT)geometry.size();

	//2 polynomial contributions per segment
	UINT N_S = 2;
	POLYMAT pol_tmp(N_S * N_E);

	m_degree = 1;

	VECTOR partition(2);
	partition(0) = 0; partition(1) = 1;

	MATRIX coeffs(N_S, m_degree + 1);
	coeffs(0, 0) = 1; coeffs(0, 1) = 0;			// upslope = s
	coeffs(1, 0) = -1; coeffs(1, 1) = 1;		// downslope = -s+1

	for (int edge = 0; edge<(int)N_E; edge++)
	{
		double	l = 1;
		if (scale == true) l = geometry.at(edge).l;

		// upslope
		pol_tmp.at(N_S * edge) = CPiecewisePol(partition, coeffs.row(0) / l, m_degree);

		// downslope
		pol_tmp.at(N_S * edge + 1) = CPiecewisePol(partition, coeffs.row(1) / l, m_degree);
	}

	CBasisFunction obj = CBasisFunction();
	obj.pol = pol_tmp;
	obj.idx_table = create_index_table(N_E, N_S, N_E-1, 1);

	return obj;
}

// //Create square basis, "scale" set to true scales by length
CBasisFunction CBasisFunction::createSquare(GEOMETRY& geometry, bool scale)
{
	UINT N_E = (UINT)geometry.size();

	//1 polynomial contributions per segment
	UINT N_S = 1;
	POLYMAT pol_tmp(N_S * N_E);

	m_degree = 0;

	VECTOR partition(2);
	partition(0) = 0; partition(1) = 1;
	
	MATRIX coeffs(N_S, m_degree + 1);
	coeffs(0, 0) = 1;				// constant

	for (int edge = 0; edge<(int)N_E; edge++)
	{
		double	l = 1;
		if (scale == true) l = geometry.at(edge).l;

		pol_tmp.at(edge) = CPiecewisePol(partition, coeffs / l, m_degree);
	}

	CBasisFunction obj = CBasisFunction();
	obj.pol = pol_tmp;
	obj.idx_table = create_index_table(N_E, 1, 0, 1);

	return obj;
}

//Create dual hat basis, "scale" set to true scales by 1/length
CBasisFunction CBasisFunction::createDualHat(GEOMETRY& dual_geometry, bool scale)
{
	UINT N_E = (UINT)dual_geometry.size();

	//4 polynomial contributions per segment
	UINT N_S = 4;
	POLYMAT pol_tmp(N_S * N_E);

	m_degree = 1;

	VECTOR partition(2);
	partition(0) = 0; partition(1) = 1;

	MATRIX coeffs(N_S, m_degree + 1);
	coeffs(0, 0) = 0.5; coeffs(0, 1) = 0;		// upslope 1 = s/2
	coeffs(1, 0) = 0.5; coeffs(1, 1) = 0.5;		// upslope 2 = s/2 + 1/2
	coeffs(2, 0) = -0.5; coeffs(2, 1) = 1;		// downslope 3 = -s/2 + 1
	coeffs(3, 0) = -0.5; coeffs(3, 1) = 0.5;	// downslope 4 = -s/2 + 1/2

	// zero polynomial
	CPiecewisePol zero_pol = CPiecewisePol(zeros(2, 1), zeros(1, 2), 0);

	for (int edge = 0; edge<(int)N_E; edge+=2)
	{
		double	l = 1;
		if (scale == true) l = dual_geometry.at(edge).l;

		// first edge has upslope 2 and downslope 4 polynomials
		pol_tmp.at(N_S * edge + 0) = zero_pol;
		pol_tmp.at(N_S * edge + 1) = CPiecewisePol(partition, coeffs.row(1) / l, m_degree);
		pol_tmp.at(N_S * edge + 2) = zero_pol;
		pol_tmp.at(N_S * edge + 3) = CPiecewisePol(partition, coeffs.row(3) / l, m_degree);

		// second edge has upslope 1 and downslope 3 polynomials
		pol_tmp.at(N_S * edge + 4) = CPiecewisePol(partition, coeffs.row(0) / l, m_degree);
		pol_tmp.at(N_S * edge + 5) = zero_pol;
		pol_tmp.at(N_S * edge + 6) = CPiecewisePol(partition, coeffs.row(2) / l, m_degree);
		pol_tmp.at(N_S * edge + 7) = zero_pol;
	}

	CBasisFunction obj = CBasisFunction();
	obj.pol = pol_tmp;
	obj.idx_table = create_index_table(N_E, N_S, N_E-1, 2);

	return obj;
}

//Create dual hat basis, "scale" set to true scales by 1/length
CBasisFunction CBasisFunction::createDualSquare(GEOMETRY& dual_geometry, bool scale)
{
	UINT N_E = (UINT)dual_geometry.size();

	//4 polynomial contributions per segment
	UINT N_S = 2;
	POLYMAT pol_tmp(N_S * N_E);

	m_degree = 0;

	VECTOR partition(2);
	partition(0) = 0; partition(1) = 1;
	
	MATRIX coeffs(1, m_degree + 1);
	coeffs(0, 0) = 1;				// constant

	// zero polynomial
	CPiecewisePol zero_pol = CPiecewisePol(zeros(2, 1), zeros(1, 1), 0);

	for (int edge = 0; edge<(int)N_E; edge += 2)
	{
		double	l = 1;
		//Each edge is split into 2, thus the pulse will need to be doubled
		//to get back to original length
		if (scale == true) l = dual_geometry.at(edge).l * 2;

		// first edge has first half of square, second edge has the second half
		pol_tmp.at(N_S * edge + 0) = CPiecewisePol(partition, coeffs / l, m_degree);
		pol_tmp.at(N_S * edge + 1) = zero_pol;
		pol_tmp.at(N_S * edge + 2) = zero_pol;
		pol_tmp.at(N_S * edge + 3) = CPiecewisePol(partition, coeffs / l, m_degree);

	}

	CBasisFunction obj = CBasisFunction();
	obj.pol = pol_tmp;
	obj.idx_table = create_index_table(N_E, N_S, 0, 2);

	return obj;
}

CBasisFunction CBasisFunction::divergence(GEOMETRY& geometry)
{
	POLYMAT polys = pol;

	UINT numBasisFunctions = (UINT)geometry.size();

	UINT numElements = (UINT)polys.size() / (UINT)geometry.size();
	for (int segment = 0; segment<(int)numBasisFunctions; segment++)
	{
		double	l = geometry.at(segment).l;

		for (int e = 0; e<(int)numElements; e++)
		{
			polys.at(e + numElements * segment).diff();
			polys.at(e + numElements * segment).m_coeffs = polys.at(e + numElements * segment).m_coeffs / l;
		}
	}

	CBasisFunction obj = CBasisFunction();
	obj.pol = polys;
	obj.idx_table = idx_table;

	return	obj;
}
