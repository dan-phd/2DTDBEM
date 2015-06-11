#pragma once
#include "piecewisepol.h"

class CBasisFunction :
	public CPiecewisePol
{
public:
	CBasisFunction(void);
	~CBasisFunction(void);

	POLYMAT			createHat(GEOMETRY& geometry, BOOL scale);
	POLYMAT			createSquare(GEOMETRY& geometry, BOOL scale);
	POLYMAT			createHat(ZGEOMETRY& zgeometry, BOOL scale);
	POLYMAT			createSquare(ZGEOMETRY& zgeometry, BOOL scale);
	void			divergence(POLYMAT& obj, GEOMETRY& geometry);
	void			divergence(POLYMAT& obj, ZGEOMETRY& zgeometry);

	double			distance(POINT2D&	pt1, POINT2D&	pt2)
	{
		double	ret = sqrt((pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y));
		return	ret;
	}
};

