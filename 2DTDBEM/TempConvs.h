#pragma once
#include "Define.h"
#include "Lagrange_interp.h"

class CTempConvs
{
public:
	CTempConvs(void);
	~CTempConvs(void);

	VECTOR	dotdiv(VECTOR& a, VECTOR& b);
	void	compute(VECTOR& distances,
				CLagrange_interp& intTB,
				CLagrange_interp& TB,
				CLagrange_interp& dTB,
				VECTOR& Fh,
				VECTOR& Fs,
				VECTOR& dF);
};
