// Definition of types

#pragma once		// this should work for most compilers

#ifndef DEFINE_H
#define DEFINE_H

// Check if we're running on Windows
#if defined(_WIN32) || defined(WIN32)
#define OS_WIN
#endif

#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <vector>

using namespace std;
using namespace arma;

#define	VECTOR			vec
#define	MATRIX			mat
#define	GEOMETRY		vector<EDGE>
#define	POLYMAT			vector<CPiecewisePol>
#define POLYCUBE		vector<cube>

#define INF				datum::inf
#define NaN				1<<30
#define PI				datum::pi

typedef	unsigned int	UINT;


struct POINT2D
{
	POINT2D()
	{
		x = 0; y = 0;
	}
	POINT2D(double ax, double ay)
	{
		x = ax;
		y = ay;
	}
	double	x;
	double	y;
};

struct EDGE
{
	UINT he_idx;
	POINT2D	a;
	POINT2D	b;
	double	l;
	POINT2D t;
	POINT2D n;
	double shape;
};


#endif