#include "Lagrange_interp.h"

CLagrange_interp::CLagrange_interp(void)
{
}

CLagrange_interp::~CLagrange_interp(void)
{
}

CLagrange_interp::CLagrange_interp(const double dt, unsigned int degree)
{
	//Set the partitions for which each polynomial will act upon
	m_degree = degree;
	m_dt = dt;
	m_partition = makevec(-dt, degree*dt, dt);

	//Get the coefficients for each polynomial
	MATRIX	coeffs = zeromat(degree + 1, degree + 1);
	for (int i = 0; i<(int)degree + 1; i++)
	{
		VECTOR f(1);
		f(0) = 1;
		for (int phi = 1; phi <= i; phi++)
		{
			VECTOR vec1(2);
			vec1(0) = (double)-1 / phi / dt;
			vec1(1) = 1;
			f = conv(f, vec1);
		}

		VECTOR g(1);
		g(0) = 1;
		for (int phi = 1; phi <= (int)degree - i; phi++)
		{
			VECTOR vec1(2);
			vec1(0) = (double)1 / phi / dt;
			vec1(1) = 1;
			g = conv(g, vec1);
		}
		setrow(coeffs, i, conv(f, g));
	}
	m_coeffs = coeffs;
}

CLagrange_interp CLagrange_interp::pad_coeffs(CLagrange_interp& varargin)
{
	CLagrange_interp	varargout = varargin;
	// pad columns
	varargout.m_coeffs = padarray(varargin.m_coeffs,
		m_coeffs.n_cols - varargin.m_coeffs.n_cols, "pre");

	// pad rows
	varargout.m_coeffs = padarray(varargout.m_coeffs,
		m_coeffs.n_rows - varargin.m_coeffs.n_rows, "post");

	// pad partition
	varargout.m_partition = padarray(varargout.m_partition,
		m_coeffs.n_rows - varargin.m_coeffs.n_rows);
	//MATRIX	ret;
	return	varargout;
}

MATRIX CLagrange_interp::padarray(MATRIX& mt, int padsize, const char* direction)
{
	int i = 0, j = 0;
	MATRIX ret;
	if (strcmp(direction, "pre") == 0)
	{
		ret = zeromat(mt.n_rows, mt.n_cols + padsize);
		for (i = 0; i<(int)ret.n_rows; i++)
			for (j = padsize; j<(int)ret.n_cols; j++)
				ret(i, j) = mt(i, j - padsize);
	}
	else if (strcmp(direction, "post") == 0)
	{
		ret = zeromat(mt.n_rows + padsize, mt.n_cols);
		for (i = 0; i<(int)mt.n_rows; i++)
			for (j = 0; j<(int)ret.n_cols; j++)
				ret(i, j) = mt(i, j);
	}
	else
		return mt;
	return	ret;
}

VECTOR CLagrange_interp::padarray(VECTOR& vt, const int padsize)
{
	VECTOR	ret = zerovec(vt.size() + padsize);
	for (int i = 0; i<(int)vt.size(); i++)
		ret(i) = vt(i);
	return	ret;
}

void CLagrange_interp::operator=(CLagrange_interp& other)
{
	CPiecewisePol::operator=(other);

	this->m_dt = other.m_dt;
}
