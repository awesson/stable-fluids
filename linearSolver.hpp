// Karen's CGD

#define MAX_STEPS 100000


/*****************************************************************************
 * Matrix class the solver will accept
 *****************************************************************************/
class ImplicitMatrix
{
public:
	virtual ~ImplicitMatrix()
	{
	}

	virtual void MatVecMult(const double x[], double r[]) const = 0;
};


/*****************************************************************************
 * This implicit matrix represents a discrete Laplacian operator.
 *****************************************************************************/
class ImplicitMatrixLap : public ImplicitMatrix
{
public:
	explicit ImplicitMatrixLap(int gridSize) : m_GridSize(gridSize)
	{
	}

	/*************************************************************************
	 * Calculates the result of an implicit matrix multiply with the vector x
	 * and stores it in r.
	 * This implicit matrix multiply calculates the Laplacian of x
	 * assuming x is a grid of m_NumCells by m_NumCells. In general,
	 * the Laplacian at x[i] will be -4*x[i] plus each of its neighboring
	 * cells in the grid. This would then be divided by the size of a cell, but
	 * I assume that the cell size is unit length.
	 *
	 * @param x The grid to take the Laplacian of.
	 * @param r Stores the calculated Laplacian.
	 *************************************************************************/
	virtual void MatVecMult(const double x[], double r[]) const override;

private:
	// The width and length of grid of cells.
	int m_GridSize;
};

/*****************************************************************************
 * Solves Ax = b for a symmetric, positive definite matrix A.
 * A is represented implicitly by the function MatVecMult()
 * which performs a matrix vector multiply, Av, and places the result in r.
 * 
 * @param n             The length of the vectors x and b.
 * @param A             Symmetric, positive definite matrix.
 * @param epsilon       The error tolerance.
 * @param steps[in,out] When passed in, steps is the maximum number of
 *                      iterations to take to solve. If not set, or set to
 *                      less than 1, MAX_STEPS is used instead.
 *                      Upon completion, steps contains the actual number of
 *                      iterations taken.
 *****************************************************************************/
double ConjGrad(int n,
				const ImplicitMatrix *A,
				double x[],
				const double b[],
				double epsilon,
				int *steps);

/** Some vector helper functions. */
void VecAddEqual(int n, double r[], const double v[]);
void VecDiffEqual(int n, double r[], const double v[]);
void VecAssign(int n, double v1[], const double v2[]);
void VecTimesScalar(int n, double v[], double s);
double VecDot(int n, const double v1[], const double v2[]);
double VecSqrLen(int n, const double v[]);
