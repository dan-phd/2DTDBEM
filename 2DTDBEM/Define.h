// Definition of types

#pragma once

#include <stdio.h>
#include <math.h>
#include <armadillo>
#include <vector>
using namespace std;
using namespace arma;

struct	POINT2D
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
	POINT2D	a;
	POINT2D	b;
	double	l;
};

struct ZEDGE
{
	double he_idx;
	POINT2D	a;
	POINT2D	b;
	double	l;
	POINT2D t;
	POINT2D n;
};
#define	VECTOR			vec
#define	PTVECTOR		vector<POINT2D>
#define	GEOMETRY		vector<EDGE>
#define	ZGEOMETRY		vector<ZEDGE>
#define	POLYMAT			vector<CPiecewisePol>
#define POLYCUBE		vector<cube>
#define	MATRIX			mat
#define INF				1<<30
#define NaN				1<<30
#define	MSG( msg )		printf( "%s",msg );
#define END( ary )		ary( ary.size()-1 )
#define PI				3.1415926539
#define eps				2.2204e-016
typedef	unsigned int	UINT;
typedef	int				BOOL;
#define TRUE			1
#define FALSE			0
#define POWER			10
#define TIMES			100 
#define PLUS			1000
enum
{
	NO_ERR = 0,
	MAT_ERR = 10,
};
typedef	int	ErrorCode;
