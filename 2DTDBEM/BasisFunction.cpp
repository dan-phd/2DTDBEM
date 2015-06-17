#include "BasisFunction.h"

CBasisFunction::CBasisFunction(void)
{
}

CBasisFunction::~CBasisFunction(void)
{
}

//Create hat basis, "scale" set to true scales by 1/length
POLYMAT	CBasisFunction::createHat(GEOMETRY& geometry, bool scale)
{
	UINT numBasisFunctions = (UINT)geometry.size();

	//2 polynomial contributions per segment
	POLYMAT obj(2 * numBasisFunctions);
	if (numBasisFunctions<2 && scale == false)
		return	obj;
	m_degree = 1;
	for (int segment = 0; segment<(int)numBasisFunctions; segment++)
	{
		double	l = 1;
		int next = segment + 1;
		if (next == numBasisFunctions)
			next = 0;
		if (scale == true)
			l = geometry.at(segment).l;
		// upslope = s
		VECTOR	partition(2);
		partition(0) = 0; partition(1) = 1;
		MATRIX	coeffs(1, 2);
		coeffs(0, 0) = 1; coeffs(0, 1) = 0;
		obj.at(2 * segment) = CPiecewisePol(partition, coeffs / l, m_degree);
		// downslope = -s+1
		coeffs(0, 0) = -1; coeffs(0, 1) = 1;
		obj.at(1 + 2 * segment) = CPiecewisePol(partition, coeffs / l, m_degree);
	}
	return	obj;
}

// //Create square basis, "scale" set to true scales by length
POLYMAT	CBasisFunction::createSquare(GEOMETRY& geometry, bool scale)
{
	UINT numBasisFunctions = (UINT)geometry.size();
	POLYMAT obj(numBasisFunctions);
	if (numBasisFunctions<2 && scale == false)
		return	obj;
	m_degree = 0;

	// 1 polynomial contribution per segment
	for (int segment = 0; segment<(int)numBasisFunctions; segment++)
	{
		double	l = 1;
		int next = segment + 1;
		if (next == numBasisFunctions)
			next = 0;
		if (scale == true)
			l = geometry.at(segment).l;
		// constant
		VECTOR	partition(2);
		partition(0) = 0; partition(1) = 1;
		MATRIX	coeffs(1, 1);
		coeffs(0, 0) = 1;
		obj.at(segment) = CPiecewisePol(partition, coeffs / l, m_degree);
	}
	return	obj;
}

void CBasisFunction::divergence(POLYMAT& obj, GEOMETRY& geometry)
{
	UINT numBasisFunctions = (UINT)geometry.size();

	UINT numElements = (UINT)obj.size() / (UINT)geometry.size();
	for (int segment = 0; segment<(int)numBasisFunctions; segment++)
	{
		double	l = geometry.at(segment).l;

		for (int e = 0; e<(int)numElements; e++)
		{
			obj.at(e + numElements * segment).diff();
			obj.at(e + numElements * segment).m_coeffs = obj.at(e + numElements * segment).m_coeffs / l;

		}
	}
}
