#pragma once
#include "PiecewisePol.h"

class CBasisFunction :
	// inheritance
	public CPiecewisePol
{
public:
	CBasisFunction(void);
	~CBasisFunction(void);

	CBasisFunction	createHat(GEOMETRY& geometry, bool scale, UINT num_shapes);
	CBasisFunction	createSquare(GEOMETRY& geometry, bool scale);
	CBasisFunction	createDualHat(GEOMETRY& dual_geometry, bool scale, UINT num_shapes);
	CBasisFunction	createDualSquare(GEOMETRY& dual_geometry, bool scale);
	CBasisFunction	divergence(GEOMETRY& geometry);

	inline double distance(POINT2D& pt1, POINT2D& pt2) {
		return sqrt((pt1.x - pt2.x)*(pt1.x - pt2.x) + (pt1.y - pt2.y)*(pt1.y - pt2.y));
	}

	//properties
	umat idx_table;		//index table to define which polynomial segments to apply to each edge 
	POLYMAT pol;		//polynomial basis function(s) for each edge
};

