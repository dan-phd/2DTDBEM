#pragma once
#include "piecewisepol.h"

class CLagrange_interp :
	public CPiecewisePol
{
public:
	CLagrange_interp(void);
	CLagrange_interp(const double dt, unsigned int degree);
	~CLagrange_interp(void);

	//Properties
	double	m_dt;

	//Function
	CLagrange_interp	pad_coeffs(CLagrange_interp& varargin);
	MATRIX				padarray(MATRIX& mt, const int padsize, const char* direction);
	VECTOR				padarray(VECTOR& vt, const int padsize);

	void	operator=(CLagrange_interp& other);
};
