#pragma once
#include "Define.h"
#include "Lagrange_interp.h"

class CTempConvs
{
public:
	CTempConvs(void);
	~CTempConvs(void);

	void	compute(MATRIX& distances,
				CLagrange_interp& intTB,
				CLagrange_interp& TB,
				CLagrange_interp& dTB,
				MATRIX& Fh,
				MATRIX& Fs,
				MATRIX& dF);

private:

	VECTOR	dotdiv(VECTOR& a, VECTOR& b);
	VECTOR	dotpow(VECTOR& a, UINT pow);
	
	double	max(const double a, const double b);
	double	heaviside(const double a);
	double	dotdiv(double& a, double& b);
	double	dotpow(double& a, UINT power);
};
