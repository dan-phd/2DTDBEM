#pragma once
#include "piecewisepol.h"

class CBasisFunction :
	public CPiecewisePol
{
public:
	CBasisFunction(void);
	~CBasisFunction(void);

	POLYMAT		createHat(GEOMETRY& geometry, bool scale);
	POLYMAT		createSquare(GEOMETRY& geometry, bool scale);
	void		divergence(POLYMAT& obj, GEOMETRY& geometry);

	double		distance(POINT2D& pt1, POINT2D& pt2)
	{
		double	ret = sqrt((pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y));
		return	ret;
	}
};

