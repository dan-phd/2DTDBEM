#include "PiecewisePol.h"
#include "Lagrange_interp.h"

CPiecewisePol::CPiecewisePol(void)
{
}

CPiecewisePol::~CPiecewisePol(void)
{
}

//Constructor
CPiecewisePol::CPiecewisePol(VECTOR partition, MATRIX coeffs, UINT degree)
{
	m_partition = partition;
	m_coeffs = coeffs;
	m_degree = degree;
}

//Evaluates the piecewise polynomial at times t
VECTOR CPiecewisePol::eval(VECTOR& t)
{
	VECTOR	y = zerovec(t.size());

	for (int i = 0; i<(int)t.size(); i++)
	{
		for (int interval = 0; interval<(int)m_partition.size() - 1; interval++)
		{
			if (t(i) >= m_partition(interval) && t(i) <= m_partition(interval + 1))
				y(i) = polyval(getrow(m_coeffs, interval), t(i));
		}
	}
	return	y;
}

//Shift function using k (includes dt) and scale using p, to obtain T(k+pt)
void CPiecewisePol::translate(const double k, const int p)
{
	if (p == 0)
		return;
	int i, j;
	UINT	eff_degree = m_coeffs.n_cols;       //effective number of degrees
	//Calculate new coefficients from transform and shift procedure
	MATRIX	new_coeffs = zeromat(m_coeffs.n_rows, eff_degree);
	for (i = 0; i<(int)m_coeffs.n_rows; i++)
	{
		VECTOR	q(1);
		q(0) = m_coeffs(i, 0);
		for (j = 1; j<(int)eff_degree; j++)
		{
			VECTOR	vec1(2);
			vec1(0) = p; vec1(1) = k;
			q = conv(vec1, q);
			q(q.size() - 1) += m_coeffs(i, j);
		}
		setrow(new_coeffs, i, q);
	}
	//Sort out partitioning for new object
	VECTOR new_partition(m_partition.size());
	for (i = 0; i<(int)m_partition.size(); i++)
		new_partition(i) = (m_partition(i) - k) / p;
	//Mirror the function if needed
	if (p < 0)
	{
		flipdim(new_partition);
		flipdim(new_coeffs, 1);
	}
	m_partition = new_partition;
	m_coeffs = new_coeffs;
}

void CPiecewisePol::diff()
{
	MATRIX	coeffs = zeromat(m_coeffs.n_rows, m_coeffs.n_cols - 1);
	// 	//Resize the coefficient matrix so there is 1 less coefficient in each partition
	for (int i = 0; i<(int)coeffs.n_rows; i++)
		for (int j = 0; j<(int)coeffs.n_cols; j++)
			coeffs(i, j) = m_coeffs(i, j);

	//Derivative for each polynomial row of the coefficient matrix
	for (int i = 0; i<(int)coeffs.n_rows; i++)
		setrow(coeffs, i, polyder(getrow(m_coeffs, i)));

	m_degree = m_degree - 1;
	m_coeffs = coeffs;
}

void CPiecewisePol::integrate()
{
	//Divide each coefficient by its degree+1
	MATRIX	coeffs = zeromat(m_coeffs.n_rows + 1, m_coeffs.n_cols + 1);
	VECTOR	partition = zerovec(m_partition.size() + 1);
	for (int i = 0; i<(int)m_partition.size(); i++)
		partition(i) = m_partition(i);
	int i, j;
	for (i = 0; i<(int)m_coeffs.n_rows; i++)
		for (j = 0; j<(int)m_coeffs.n_cols; j++)
			coeffs(i, j) = m_coeffs(i, j) / (m_coeffs.n_cols - j);

	//Make sure the end partition goes to infinity
	coeffs(m_degree + 1, m_degree + 1) = 0;

	//Since this function will be smoother, compute the extra
	//(constant) coefficient for each polynomial using initial condition, f_x
	double	f_x = 0;
	for (int i = 0; i<(int)m_partition.size(); i++)
	{
		coeffs(i, coeffs.n_cols - 1) = f_x - polyval(getrow(coeffs, i), partition(i));
		f_x = polyval(getrow(coeffs, i), partition(i + 1));
	}
	m_degree++;
	m_coeffs = coeffs;
	m_partition = partition;
}

//Convolution and polynomial multiplication
VECTOR	CPiecewisePol::conv(const VECTOR& u, const VECTOR& v)
{
	int j, k;
	UINT	m = u.size();
	UINT	n = v.size();
	VECTOR	ret = zerovec(m + n - 1);
	for (k = 0; k<(int)(m + n - 1); k++)
		for (j = 0; j <= k; j++)
		{
			if (j<(int)m && (k - j)<(int)n)
				ret(k) += u(j)*v(k - j);
		}
	return	ret;
}

//Create zero vector
VECTOR	CPiecewisePol::zerovec(const UINT size)
{
	VECTOR	ret(size);
	for (int i = 0; i<(int)size; i++)
		ret(i) = 0;
	return	ret;
}

//Get a row vector from matrix 
VECTOR CPiecewisePol::getrow(const MATRIX& mt, const UINT	row)
{
	VECTOR	ret(mt.n_cols);
	for (int i = 0; i<(int)ret.size(); i++)
		ret(i) = mt(row, i);
	return	ret;
}

//Set a row vector to matrix
bool CPiecewisePol::setrow(MATRIX& mt, const UINT row, const VECTOR& rowvec)
{
	if (mt.n_cols != rowvec.size())
	{
		printf("Error setting row!");
		return false;
	}
	for (int i = 0; i<(int)rowvec.size(); i++)
		mt(row, i) = rowvec(i);
	return true;
}


//Create zero matrix
MATRIX CPiecewisePol::zeromat(const UINT size1, const UINT size2)
{
	MATRIX	ret(size1, size2);
	for (int i = 0; i<(int)size1; i++)
		for (int j = 0; j<(int)size2; j++)
			ret(i, j) = 0;
	return	ret;
}

//Flip vector ,for example [1 2 3]-->[3 2 1]
void CPiecewisePol::flipdim(VECTOR& vec1)
{
	for (int i = 0; i<(int)vec1.size() / 2; i++)
	{
		double tmp = vec1(i);
		vec1(i) = vec1(vec1.size() - i - 1);
		vec1(vec1.size() - i - 1) = tmp;
	}
}

void CPiecewisePol::flipdim(MATRIX& mt, UINT dim)
{
	if (dim == 1)//row
	{
		for (int i = 0; i<(int)mt.n_rows / 2; i++)
			for (int j = 0; j<(int)mt.n_cols; j++)
			{
				double tmp = mt(i, j);
				mt(i, j) = mt(mt.n_rows - i - 1, j);
				mt(mt.n_rows - i - 1, j) = tmp;
			}
	}
	else if (dim == 2)//column
	{
		for (int i = 0; i<(int)mt.n_rows; i++)
			for (int j = 0; j<(int)mt.n_cols / 2; j++)
			{
				double tmp = mt(i, j);
				mt(i, j) = mt(i, mt.n_cols - j - 1);
				mt(i, mt.n_cols - j - 1) = tmp;
			}
	}
	else
		return;
}

//Get polynomial derivation
VECTOR CPiecewisePol::polyder(VECTOR& poly)
{
	VECTOR	ret = zerovec(poly.size() - 1);
	for (int i = 0; i<(int)poly.size() - 1; i++)
		ret(i) = poly(i)*(poly.size() - i - 1);
	return	ret;
}

//Get polynomial value
double CPiecewisePol::polyval(VECTOR& poly, const double x)
{
	double ret = 0;
	for (int i = 0; i<(int)poly.size(); i++)
		ret += poly(poly.size() - 1 - i)*pow(x, i);
	return	ret;
}

//Create vector from begin,end and step values
VECTOR CPiecewisePol::makevec(const double begin, const double end, const double step)
{
	int size = int((end - begin) / step + 1);
	VECTOR	ret(size);
	for (int i = 0; i<size; i++)
		ret(i) = begin + step*i;
	return	ret;
}

void CPiecewisePol::operator=(CPiecewisePol& other)
{
	this->m_partition = other.m_partition;
	this->m_coeffs = other.m_coeffs;
	this->m_degree = other.m_degree;
}
