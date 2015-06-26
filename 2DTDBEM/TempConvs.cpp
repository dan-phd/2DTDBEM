#include "TempConvs.h"
#include "PiecewisePol.h"
#include <math.h>
#include <omp.h>

CTempConvs::CTempConvs(void)
{}
CTempConvs::~CTempConvs(void)
{}

VECTOR CTempConvs::dotdiv(VECTOR& a, VECTOR& b)
{
	
	///*
	VECTOR div(a.size());
	for (int i = 0; i < (int)a.size(); i++)
	{
		// Armadillo elementwise division gives NaN for 0/0
		if ((a(i) == 0))
			div(i) = 0;
		else if (b(i) == 0)
			div(i) = INF;
		else div(i) = a(i) / b(i);
	}
	//*/

	// Alternatively, the following gives identicle results to Matlab
	/*
	VECTOR div(a.size());
	for (int i = 0; i<(int)a.size(); i++){
		div(i) = 0;
	}
	if (a.size() != b.size())
		return	div;

	for (int i = 0; i<(int)a.size(); i++)
	{
		if (a(i) == NaN)
			div(i) = NaN;
		if (a(i) == INF)
		{
			if (b(i) == INF || b(i) == NaN)
				div(i) = NaN;
			else
				div(i) = INF;
		}
		if (a(i) == 0)
		{
			if (b(i) != 0 || b(i) == INF)
				div(i) = 0;
			else
				div(i) = NaN;
		}
		else
		{
			if (b(i) == 0 || b(i) == NaN)
				div(i) = NaN;
			if (b(i) == INF)
				div(i) = 0;
			else
				div(i) = a(i) / b(i);
		}
	}
	//*/

	//a.print("a:"); b.print("b:"); div.print("div:");

	return div;

}
VECTOR CTempConvs::dotpow(VECTOR& a, UINT power)
{
	VECTOR tmp = a;
	for (UINT p = 1; p < power; p++)
		tmp = tmp % a;

	return tmp;
}

void CTempConvs::compute(VECTOR& P, CLagrange_interp& intTB, CLagrange_interp& TB, CLagrange_interp& dTB,
	VECTOR& Fh, VECTOR& Fs, VECTOR& dF)
{
	UINT num_partitions = intTB.m_partition.size() - 1;
	UINT num_degree = intTB.m_coeffs.n_cols - 1;

	MATRIX Ia(P.size(), num_degree + 1, fill::zeros);
	MATRIX Ib(P.size(), num_degree + 1, fill::zeros);
	MATRIX dIa(Ia), dIb(Ib);

	for (int i = 0; i<(int)num_partitions; i++)
	{

		//Current partition properties
		double partition_start = intTB.m_partition(i);		//(k-i)dt  
		double partition_end = intTB.m_partition(i + 1);    //(k-i+1)dt

		//Integration limits
		VECTOR a = clamp(P, partition_start, INF);
		VECTOR b = clamp(P, partition_end, INF);
		
		//Differentiated limits wrt P (using Heaviside)
		VECTOR da = P;
		da.elem(find(da > partition_start)).ones();
		da.elem(find(da == partition_start)).fill(0.5);
		da.elem(find(da < partition_start)).zeros();
		VECTOR db = P;
		db.elem(find(db > partition_end)).ones();
		db.elem(find(db == partition_end)).fill(0.5);
		db.elem(find(db < partition_end)).zeros();

		//Compute dI_1 term and set this term to 0 after threshold point
		VECTOR tmp1 = a % da - P;
		VECTOR tmp2 = sqrt(square(a) - square(P));
		VECTOR dI1a = dotdiv(tmp1, tmp2);
		dI1a.elem(find(dI1a > partition_start)).zeros();

		VECTOR tmp3 = b % db - P;
		VECTOR tmp4 = sqrt(square(b) - square(P));
		VECTOR	dI1b = dotdiv(tmp3, tmp4);
		dI1b.elem(find(dI1b > partition_end)).zeros();

		//First integral equation, I_0 and dI_0
		tmp2 += a;
		tmp4 += b;
		Ia.col(0) = log(tmp2);
		Ib.col(0) = log(tmp4);
		tmp1 = da + dI1a;
		tmp3 = db + dI1b;
		dIa.col(0) = tmp1 / tmp2;
		dIb.col(0) = tmp3 / tmp4;

		//Convolve polynomial coefficients with solved integral
		tmp1 = Ib.col(0) - Ia.col(0);
		Fh += tmp1 * intTB.m_coeffs(i, num_degree);
		Fs += tmp1 * TB.m_coeffs(i, num_degree);
		tmp1 = dIb.col(0) - dIa.col(0);
		dF += tmp1* dTB.m_coeffs(i, num_degree);

		//Second integral equation, I_1
		if (num_degree>0)
		{
			//I_1
			Ia.col(1) = sqrt(square(a) - square(P));
			Ib.col(1) = sqrt(square(b) - square(P));

			//dI_1
			dIa.col(1) = dI1a;
			dIb.col(1) = dI1b;

			//Convolve
			tmp1 = Ib.col(1) - Ia.col(1);
			Fh += tmp1 * intTB.m_coeffs(i, num_degree - 1);
			Fs += tmp1 * TB.m_coeffs(i, num_degree - 1);
			tmp1 = dIb.col(1) - dIa.col(1);
			dF += tmp1* dTB.m_coeffs(i, num_degree - 1);

			//All other integral equations, I_2...I_p
			for (int n = 2; n <= (int)num_degree; n++)
			{
				// I_2...I_p
				tmp1 = pow(a, n - 1) % Ia.col(1);
				tmp2 = square(P);
				tmp2 = (tmp1 + (n - 1) * (tmp2 % Ia.col(n - 2))) / n;
				Ia.col(n) = tmp2;
				tmp1 = dotpow(b, n - 1);
				tmp1 = tmp1 % Ib.col(1);
				tmp2 = square(P);
				tmp2 = (tmp1 + (n - 1)*(tmp2 % Ib.col(n - 2))) / n;
				Ib.col(n) = tmp2;

				//dI_2 ... dI_p
				tmp1 = dotpow(da, n - 1);
				tmp1 = tmp1 % Ia.col(1);
				tmp2 = dotpow(a, n - 1);
				tmp2 = tmp2 % dIa.col(1);
				tmp3 = 2 * P;
				tmp3 = tmp3 % Ia.col(n - 2);
				tmp4 = square(P) % dIa.col(n - 2);
				tmp4 = (n - 1)*(tmp3 + tmp4);
				tmp4 = (tmp1 + tmp2 + tmp4) / (double)n;
				dIa.col(n) = tmp4;

				tmp1 = dotpow(db, n - 1);
				tmp1 = tmp1 % Ib.col(1);
				tmp2 = dotpow(b, n - 1);
				tmp2 = tmp2 % dIb.col(1);
				tmp3 = 2 * P;
				tmp3 = tmp3 % Ib.col(n - 2);
				tmp4 = square(P) % dIb.col(n - 2);
				tmp4 = (n - 1)*(tmp3 + tmp4);
				tmp4 = (tmp1 + tmp2 + tmp4) / (double)n;
				dIb.col(n) = tmp4;

				// Convolve
				tmp1 = Ib.col(n) - Ia.col(n);
				Fh += tmp1 * intTB.m_coeffs(i, num_degree - n);
				Fs += tmp1 * TB.m_coeffs(i, num_degree - n);
				tmp1 = dIb.col(n) - dIa.col(n);
				dF += tmp1* dTB.m_coeffs(i, num_degree - n);
			}
		}
	}

	Fh = Fh / (2 * PI);
	Fs = Fs / (2 * PI);
	dF = dF / (2 * PI);
}


double CTempConvs::max(const double a, const double b)
{
	if (a >= b)
		return a;
	else
		return b;
}
double CTempConvs::heaviside(const double a)
{
	if (a > 0.0)
		return 1.0;
	else if (a==0.0)
		return 0.5;
	else
		return 0.0;
}
double CTempConvs::dotdiv(double& a, double& b)
{
	if (b == 0) return INF;
	else if (a == 0) return 0;
	else return ( a / b );
}
double CTempConvs::dotpow(double& a, UINT power)
{
	double tmp(a);
	for (UINT p = 1; p < power; p++)
		tmp*=a;

	return tmp;
}

void CTempConvs::compute2(VECTOR& P,
	CLagrange_interp& intTB, CLagrange_interp& TB, CLagrange_interp& dTB,
	VECTOR& Fh, VECTOR& Fs, VECTOR& dF )
{
	const UINT num_partitions( intTB.m_partition.size() - 1 );
	const UINT num_degree( intTB.m_coeffs.n_cols - 1 );

	double *Ia, *Ib, *dIa, *dIb;

	int i, n, ii;
#pragma omp parallel default(shared) private(i, ii, n, Ia, Ib, dIa, dIb)
	{
		Ia = (double*)(calloc(sizeof(double), num_degree+1));
		Ib = (double*)(calloc(sizeof(double), num_degree + 1));
		dIa = (double*)(calloc(sizeof(double), num_degree + 1));
		dIb = (double*)(calloc(sizeof(double), num_degree + 1));

#pragma omp for schedule(static)
		for (ii = 0; ii < (int)P.size(); ii++)
		{

			//Current distance
			double PP(P(ii));

			for (i = 0; i<(int)num_partitions; i++)
			{
				//Current partition properties
				double partition_start(intTB.m_partition(i));		 //(k-i)dt  
				double partition_end(intTB.m_partition(i + 1));    //(k-i+1)dt

				//Integration limits
				double a(max(PP, partition_start));
				double b(max(PP, partition_end));

				//Differentiated limits wrt P (using Heaviside)
				double da(heaviside(PP - partition_start));
				double db(heaviside(PP - partition_end));

				//Compute dI_1 term and set this term to 0 after threshold point
				double tmp1(a*da - PP);
				double tmp2(sqrt(a*a - PP*PP));
				double dI1a(dotdiv(tmp1, tmp2));
				if (dI1a>partition_start) dI1a = 0;

				double tmp3(b*db - PP);
				double tmp4(sqrt(b*b - PP*PP));
				double dI1b(dotdiv(tmp3, tmp4));
				if (dI1b > partition_end) dI1b = 0;

				//First integral equation, I_0 and dI_0
				tmp2 += a;
				tmp4 += b;
				Ia[0] = log(tmp2);
				Ib[0] = log(tmp4);
				dIa[0] = (da + dI1a) / tmp2;
				dIb[0] = (db + dI1b) / tmp4;

				//Convolve polynomial coefficients with solved integral
				tmp1 = Ib[0] - Ia[0];
				Fh.at(ii) += tmp1 * intTB.m_coeffs(i, num_degree);
				Fs.at(ii) += tmp1 * TB.m_coeffs(i, num_degree);
				dF.at(ii) += (dIb[0] - dIa[0]) * dTB.m_coeffs(i, num_degree);

				//Second integral equation, I_1
				if (num_degree > 0)
				{
					//I_1
					Ia[1] = sqrt(a*a - PP*PP);
					Ib[1] = sqrt(b*b - PP*PP);

					//dI_1
					dIa[1] = dI1a;
					dIb[1] = dI1b;

					//Convolve
					tmp1 = Ib[1] - Ia[1];
					Fh.at(ii) += tmp1 * intTB.m_coeffs(i, num_degree - 1);
					Fs.at(ii) += tmp1 * TB.m_coeffs(i, num_degree - 1);
					dF.at(ii) += (dIb[1] - dIa[1]) * dTB.m_coeffs(i, num_degree - 1);

					//All other integral equations, I_2...I_p
					for (n = 2; n <= (int)num_degree; n++)
					{
						double over_n(1 / (double)n);

						// I_2...I_p
						Ia[n] = (dotpow(a, n - 1) * Ia[1] + (n - 1) * PP*PP * Ia[n - 2]) * over_n;
						Ib[n] = (dotpow(b, n - 1) * Ib[1] + (n - 1) * PP*PP * Ib[n - 2]) * over_n;

						//dI_2 ... dI_p
						tmp1 = dotpow(da, n - 1) * Ia[1];
						tmp2 = dotpow(a, n - 1) * dIa[1];
						tmp3 = (n - 1)*(2 * PP * Ia[n - 2] + PP*PP * dIa[n - 2]);
						dIa[n] = (tmp1 + tmp2 + tmp3) * over_n;

						tmp1 = dotpow(db, n - 1) * Ib[1];
						tmp2 = dotpow(b, n - 1) * dIb[1];
						tmp3 = (n - 1)*(2 * PP * Ib[n - 2] + PP*PP * dIb[n - 2]);
						dIb[n] = (tmp1 + tmp2 + tmp3) * over_n;

						// Convolve
						tmp1 = Ib[n] - Ia[n];
						Fh.at(ii) += tmp1 * intTB.m_coeffs(i, num_degree - n);
						Fs.at(ii) += tmp1 * TB.m_coeffs(i, num_degree - n);
						dF.at(ii) += (dIb[n] - dIa[n]) * dTB.m_coeffs(i, num_degree - n);

					} // n -> num_degree
				} // if num_degree>0
			} // i -> num_partitions
		} // ii -> size(P)

		delete[] Ia;
		delete[] Ib;
		delete[] dIa;
		delete[] dIb;

	} // omp parallel

	const double over_2pi(1 / (2 * PI));
	Fh *= over_2pi;
	Fs *= over_2pi;
	dF *= over_2pi;
}
