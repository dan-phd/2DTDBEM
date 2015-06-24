#pragma once
#include "Define.h"

class CPiecewisePol
{
public:
	CPiecewisePol(void);
	CPiecewisePol(VECTOR partition, MATRIX coeffs, UINT degree);
	~CPiecewisePol(void);

	//Properties
	VECTOR		m_partition;
	MATRIX		m_coeffs;
	int			m_degree;	//0:Constant 1:Linear 2:Quadric 3:Cubic 4:Quartic

	//Function
	VECTOR		eval(VECTOR& t);
	void		translate(const double k, const int p);
	void		diff();
	void		integrate();

	VECTOR				conv(const VECTOR& u, const VECTOR& v);
	void				flipdim(VECTOR& vec1);
	void				flipdim(MATRIX& mt, UINT dim);
	VECTOR				polyder(VECTOR& poly);
	double				polyval(VECTOR& poly, const double x);
	static	VECTOR		makevec(const double begin, const double end, const double step);

	CPiecewisePol*		operator=(const CPiecewisePol& other);
};
