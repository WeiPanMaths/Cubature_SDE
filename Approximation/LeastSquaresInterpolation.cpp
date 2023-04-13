#include "LeastSquaresInterpolation.h"
#include "mkl.h"



typedef std::vector<double> DVector;

double getIntervalLength_(size_t bit)
{
	//assert(bit_max >= bit);

	int half_bit_rep[dlength] = { 0 };
	double v = 0.5;		// reference value

	// Boilerplate to circumvent the fact bitwise operators can't be applied to double
	union {
		double value;
		char   array[dlength];
	};

	value = v;
	size_t bit_max = 52;

	for (int i = 0; i < dlength; ++i)
	{
		int relativeToByte = i % CHAR_BIT;
		bool isBitSet = (array[sizeof(double) - 1 - i / CHAR_BIT] &
			(1 << (CHAR_BIT - relativeToByte - 1))) == (1 << (CHAR_BIT - relativeToByte - 1));
		half_bit_rep[i] = (isBitSet ? 1 : 0);
	}

	for (int i = 0; i < bit_max - bit; ++i)
		half_bit_rep[12 + i] = 1;

	//return 1. - bitstring_to_double(std::begin(half_bit_rep), std::end(half_bit_rep));
	return 1. - bitstring_to_double(std::begin(half_bit_rep), std::end(half_bit_rep));
}


double asup_norm_col_major(const DVector & mymat, const int m, const int n)
{
	double norm(0);
	for (int i = 0; i < m; ++i)		// row
	{
		double _temp(0);

		for (int j = 0; j < n; ++j) // column
		{
			_temp += mymat[i + j*m];
		}
		if (_temp > norm)
			norm = _temp;
	}
	return norm;
}

double compute_condition_number(const DVector & a, int no_of_rows, int no_of_cols)
{
	double anorm = asup_norm_col_major(a, no_of_rows, no_of_cols);
	DVector _mymatrix(a);

	// call dgetrf to compute the LU factorization of A
	lapack_int m = no_of_rows, n = no_of_cols, lda = no_of_rows;
	std::vector<lapack_int> ipiv;
	ipiv.resize(no_of_rows);
	auto info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, m, n, &_mymatrix[0], lda, &ipiv[0]);		// column major!

	double rcond(-1);

	if (info != 0)
	{
		std::cout << "The diagonal element " << info << " of the triangular factor " << std::endl;
		std::cout << "of A is zero, so that A does not have full rank;\n";
		std::cout << "the least squares solution could not be computed.\n";
		exit(1);
	}
	else
	{
		// success - copy result into local variable
		//std::cout << "\n success " << std::endl;
		char charnorm('I');
		LAPACKE_dgecon(LAPACK_COL_MAJOR, charnorm, n, &_mymatrix[0], lda, anorm, &rcond);
		//std::cout << "reciprocal condition number: " << rcond << std::endl;
		//std::cout << "condition number: " << 1. / rcond << std::endl;
	}
	return 1. / rcond;
}

#if dim_ ==2 

bool InterpolationMorton::interpolate_svd(const std::vector<point_2d> & pts, const std::vector<double> & values)
{
	int no_of_rows = pts.size(), no_of_cols = myNoOfMonomials;
	//int no_of_cols = (2+myDegree)*(1+myDegree)/2;	// number of monomials
	char charN('N');
	int m = no_of_rows, n = no_of_cols, nrhs = 1, lda = no_of_rows, ldb = std::max(1, std::max(no_of_rows, no_of_cols));

	std::vector<double> work;
	std::vector<double> a(no_of_rows*no_of_cols);
	std::vector<double> b(values);
	double rcond(-1);  // if rcond < 0 then machine precision is used. used to determine the effective rank of A.
	int rank(0);
	std::vector<double> s;
	s.resize(std::min(m, n));

	// generate matrix A in column-major
	this->generate_monomial_matrix(a, pts);
	{
		//std::cout << "calculating sup norm condition number: " << compute_condition_number(a, no_of_rows, no_of_cols) << "\n";
	}

	auto info = LAPACKE_dgelsd(LAPACK_COL_MAJOR, m, n, nrhs, &a[0], lda, &b[0], ldb, &s[0], rcond, &rank);

	if (info != 0)
	{
		std::cout << "The diagonal element " << info << " of the triangular factor " << std::endl;
		std::cout << "of A is zero, so that A does not have full rank;\n";
		std::cout << "the least squares solution could not be computed.\n";
		exit(1);
	}
	else
	{
		//// success - copy result into local variable
		//std::cout << "l2 condition number " << s[0] / s[std::min(m, n) - 1] << std::endl;
		///// could do data search here ///////

		////std::cout << *(s.begin()) / *(--s.end()) << std::endl;
		////std::cout << s[0] / s[std::min(m, n)-1] << std::endl;
		//std::cout << "\n singular values of A: " << std::endl;
		//for (auto itr = s.begin(); itr != s.end(); ++itr)
		//	std::cout << *itr << ", ";
		//std::cout << "\n effective rank: " << rank << std::endl << std::endl;
		for (int i = 0; i<myNoOfMonomials; ++i)
			myApproximationCoefficients[i] = b[i];
		return true;
	}
}

#endif
