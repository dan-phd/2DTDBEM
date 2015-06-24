#pragma once
#include "PiecewisePol.h"

class CLagrange_interp :
	// inheritance
	public CPiecewisePol
{
public:
	CLagrange_interp(void);
	CLagrange_interp(const double dt, int degree);
	~CLagrange_interp(void);

	//Properties
	double	m_dt;

	//Functions
	CLagrange_interp	pad_coeffs(CLagrange_interp& varargin);
	MATRIX				padarray(MATRIX& mt, const int padsize, const char* direction);
	VECTOR				padarray(VECTOR& vt, const int padsize);

	CLagrange_interp*	operator=(const CLagrange_interp& other);
};
