#include <stdio.h>
#include <string.h>

#include "linearSolver.hpp"


#define IX_DIM(i, j) ((i) + (m_GridSize * (j)))
#define ITER_DIM    for(int i = 1; i < m_GridSize - 1; i++){      \
						for(int j = 1; j < m_GridSize - 1; j++){
#define ENDITER_DIM }}


void ImplicitMatrixLap::MatVecMult(const double x[], double r[]) const
{
	// Set the boundaries explicitly to zero
	for (int i = 0; i < m_GridSize; ++i)
	{
		// Bottom
		r[IX_DIM(i, 0)] = 0;
		// Left
		r[IX_DIM(0, i)] = 0;
		// Top
		r[IX_DIM(i, m_GridSize - 1)] = 0;
		// Right
		r[IX_DIM(m_GridSize - 1, i)] = 0;
	}

	/** Interior corners */
	// Bottom left
	r[IX_DIM(1, 1)] = x[IX_DIM(1, 2)] + x[IX_DIM(2, 1)] - 2.0f * x[IX_DIM(1, 1)];
	// Top right
	r[IX_DIM(m_GridSize - 2, m_GridSize - 2)] = x[IX_DIM(m_GridSize - 2, m_GridSize - 3)]
											  + x[IX_DIM(m_GridSize - 3, m_GridSize - 2)]
											  - (2.0 * x[IX_DIM(m_GridSize - 2, m_GridSize - 2)]);
	// Bottom right
	r[IX_DIM(m_GridSize - 2, 1)] = x[IX_DIM(m_GridSize - 2, 2)]
								 + x[IX_DIM(m_GridSize - 3, 1)]
								 - (2.0 * x[IX_DIM(m_GridSize - 2, 1)]);
	// Top left
	r[IX_DIM(1, m_GridSize - 2)] = x[IX_DIM(2, m_GridSize - 2)]
								 + x[IX_DIM(1, m_GridSize - 3)]
								 - 2.0 * x[IX_DIM(1, m_GridSize - 2)];

	/** Interior edges */
	for (int i = 2; i < m_GridSize - 2; ++i)
	{
		// Left column
		r[IX_DIM(1, i)] = x[IX_DIM(1, i-1)] + x[IX_DIM(1, i+1)]
						+ x[IX_DIM(2, i)] - 3.0 * x[IX_DIM(1, i)];
		// Right column
		r[IX_DIM(m_GridSize - 2, i)] = x[IX_DIM(m_GridSize - 2, i-1)]
									 + x[IX_DIM(m_GridSize - 2, i+1)]
									 + x[IX_DIM(m_GridSize - 3, i)]
									 - (3.0 * x[IX_DIM(m_GridSize - 2, i)]);
		// Bottom row
		r[IX_DIM(i, 1)] = x[IX_DIM(i-1, 1)] + x[IX_DIM(i+1, 1)]
						+ x[IX_DIM(i, 2)] - 3.0 * x[IX_DIM(i, 1)];
		// Top row
		r[IX_DIM(i, m_GridSize - 2)] = x[IX_DIM(i-1, m_GridSize - 2)]
									 + x[IX_DIM(i+1, m_GridSize - 2)]
									 + x[IX_DIM(i, m_GridSize - 3)]
									 - 3.0 * x[IX_DIM(i, m_GridSize - 2)];
	}

	/** Middle grid cells */
	for (int i = 2; i < m_GridSize - 2; i++)
	{
		for (int j = 2; j < m_GridSize - 2; j++)
		{
			r[IX_DIM(i, j)] = x[IX_DIM(i+1, j)]
							+ x[IX_DIM(i-1, j)]
							+ x[IX_DIM(i, j+1)]
							+ x[IX_DIM(i, j-1)]
							- (4.0 * x[IX_DIM(i, j)]);
		}
	}
}

double ConjGrad(int n,
				const ImplicitMatrix *A,
				double x[],
				const double b[],
				double epsilon, // how low should we go?
				int *steps)
{
	int i = 0;
	int iMax;
	double alpha, beta, rSqrLen, rSqrLenOld, u;

	double *r = new double[n];
	double *d = new double[n];
	double *t = new double[n];
	double *temp = new double[n];

	VecAssign(n, x, b);

	VecAssign(n, r, b);
	A->MatVecMult(x, temp);
	VecDiffEqual(n, r, temp);

	rSqrLen = VecSqrLen(n, r);

	VecAssign(n, d, r);
	
	if (steps && *steps > 0)
	{
		iMax = *steps;
	}
	else
	{
		iMax = MAX_STEPS;
	}

	if (rSqrLen > epsilon)
	{
		while (i < iMax)
		{
			i++;
			A->MatVecMult(d, t);
			u = VecDot(n, d, t);

			if (u == 0)
			{
				printf("(SolveConjGrad) d'Ad = 0\n");
				break;
			}

			// How far should we go?
			alpha = rSqrLen / u;

			// Take a step along direction d
			VecAssign(n, temp, d);
			VecTimesScalar(n, temp, alpha);
			VecAddEqual(n, x, temp);

			if (i & 0x3F)
			{
				VecAssign(n, temp, t);
				VecTimesScalar(n, temp, alpha);
				VecDiffEqual(n, r, temp);
			}
			else
			{
				// For stability, correct r every 64th iteration
				VecAssign(n, r, b);
				A->MatVecMult(x, temp);
				VecDiffEqual(n, r, temp);
			}

			rSqrLenOld = rSqrLen;
			rSqrLen = VecSqrLen(n, r);

			// Converged! Let's get out of here
			if (rSqrLen <= epsilon)
			{
				break;
			}

			// Change direction: d = r + beta * d
			beta = rSqrLen / rSqrLenOld;
			VecTimesScalar(n, d, beta);
			VecAddEqual(n, d, r);
		}
	}

	delete[] r;
	delete[] d;
	delete[] t;
	delete[] temp;

	if (steps)
	{
		*steps = i;
	}

	return rSqrLen;
}

/** Vector helper functions */
void VecAddEqual(int n, double r[], const double v[])
{
	for (int i = 0; i < n; i++)
	{
		r[i] = r[i] + v[i];
	}
}

void VecDiffEqual(int n, double r[], const double v[])
{
	for (int i = 0; i < n; i++)
	{
		r[i] = r[i] - v[i];
	}
}

void VecAssign(int n, double v1[], const double v2[])
{
	memcpy(v1, v2, n * sizeof(double));
}

void VecTimesScalar(int n, double v[], double s)
{
	for (int i = 0; i < n; i++)
	{
		v[i] *= s;
	}
}

double VecDot(int n, const double v1[], const double v2[])
{
	double dot = 0;
	for (int i = 0; i < n; i++)
	{
		dot += v1[i] * v2[i];
	}
	return dot;
}

double VecSqrLen(int n, const double v[])
{
	return VecDot(n, v, v);
}
