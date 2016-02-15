#include "linearSolver.hpp"

#define IX_DIM(i,j) ((i)+(m_NumCells+2)*(j))
#define ITER_DIM    for(int i = 1; i <= m_NumCells; i++){      \
                        for(int j = 1; j <= m_NumCells; j++){
#define ENDITER_DIM }}

implicitMatrixLap::implicitMatrixLap(int i_NumCells): m_NumCells(i_NumCells){ }

implicitMatrixLap::~implicitMatrixLap(){ }

/*******************************************************************************
 * Calculates the result of an implicit matrix multiply with the vector x and
 * stores it in r.
 * This implicit matrix multiply calculates the laplacian of x
 * assuming x is a grid of m_NumCells by m_NumCells. In general,
 * the laplacian at x[i] will be -4*x[i] plus each of its neighboring
 * cells in the grid. This would then be divided by the size of a cell, but
 * I assume that the cell size is unit length.
 *
 * @param x The grid to take the laplacian of.
 * @param r Stores the calculated laplacian.
 ******************************************************************************/
void implicitMatrixLap::matVecMult(double x[], double r[])
{
    /**************
     * boundaries *
     **************/
    
    // set the boundaries explicitly to zero
    for(int i = 0; i <= m_NumCells + 1; ++i){
        // bottom
        r[IX_DIM(i, 0)] = 0;
        // left
        r[IX_DIM(0, i)] = 0;
        // top
        r[IX_DIM(i, m_NumCells+1)] = 0;
        // right
        r[IX_DIM(m_NumCells+1, i)] = 0;
    }
    
    /********************
     * interior corners *
     ********************/

    // bottom left
    r[IX_DIM(1, 1)] = x[IX_DIM(1, 2)] + x[IX_DIM(2, 1)] - 2.0f*x[IX_DIM(1, 1)];
    // top right
    r[IX_DIM(m_NumCells, m_NumCells)] = x[IX_DIM(m_NumCells, m_NumCells - 1)]
                                      + x[IX_DIM(m_NumCells - 1, m_NumCells)]
                                      - 2.0*x[IX_DIM(m_NumCells, m_NumCells)];
    // bottom right
    r[IX_DIM(m_NumCells, 1)] = x[IX_DIM(m_NumCells, 2)]
                             + x[IX_DIM(m_NumCells - 1, 1)]
                             - 2.0*x[IX_DIM(m_NumCells, 1)];
    // top left
    r[IX_DIM(1, m_NumCells)] = x[IX_DIM(2, m_NumCells)]
                             + x[IX_DIM(1, m_NumCells - 1)]
                             - 2.0*x[IX_DIM(1, m_NumCells)];

    /******************
     * interior edges *
     ******************/
    for(int i = 2; i < m_NumCells; ++i){
        // left column
        r[IX_DIM(1, i)] = x[IX_DIM(1, i-1)] + x[IX_DIM(1, i+1)]
                        + x[IX_DIM(2, i)] - 3.0*x[IX_DIM(1, i)];
        // right column
        r[IX_DIM(m_NumCells, i)] = x[IX_DIM(m_NumCells, i-1)]
                                 + x[IX_DIM(m_NumCells, i+1)]
                                 + x[IX_DIM(m_NumCells-1, i)] 
                                 - 3.0*x[IX_DIM(m_NumCells, i)];
        // bottom row
        r[IX_DIM(i, 1)] = x[IX_DIM(i-1, 1)] + x[IX_DIM(i+1, 1)]
                        + x[IX_DIM(i, 2)] - 3.0*x[IX_DIM(i, 1)];
        // top row
        r[IX_DIM(i, m_NumCells)] = x[IX_DIM(i-1, m_NumCells)]
                                 + x[IX_DIM(i+1, m_NumCells)]
                                 + x[IX_DIM(i, m_NumCells-1)]
                                 - 3.0*x[IX_DIM(i, m_NumCells)];
    }

    /*********************
     * middle grid cells *
     *********************/
    for(int i = 2; i < m_NumCells; i++){
        for(int j = 2; j < m_NumCells; j++){
            r[IX_DIM(i, j)] = x[IX_DIM(i+1,j)]
                            + x[IX_DIM(i-1,j)]
                            + x[IX_DIM(i,j+1)]
                            + x[IX_DIM(i,j-1)]
                            - 4.0*x[IX_DIM(i,j)];
        }
    }
}

// vector helper functions

void vecAddEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] + v[i];
}

void vecDiffEqual(int n, double r[], double v[])
{
  for (int i = 0; i < n; i++)
    r[i] = r[i] - v[i];
}

void vecAssign(int n, double v1[], double v2[])
{
  for (int i = 0; i < n; i++)
    v1[i] = v2[i];
}

void vecTimesScalar(int n, double v[], double s)
{
  for (int i = 0; i < n; i++)
    v[i] *= s;
}

double vecDot(int n, double v1[], double v2[])
{
  double dot = 0;
  for (int i = 0; i < n; i++)
    dot += v1[i] * v2[i];
  return dot;
}

double vecSqrLen(int n, double v[])
{
  return vecDot(n, v, v);
}

double ConjGrad(int n, implicitMatrix *A, double x[], double b[],
                double epsilon, // how low should we go?
                int    *steps)
{
  int       i, iMax;
  double    alpha, beta, rSqrLen, rSqrLenOld, u;

  double *r = (double *) malloc(sizeof(double) * n);
  double *d = (double *) malloc(sizeof(double) * n);
  double *t = (double *) malloc(sizeof(double) * n);
  double *temp = (double *) malloc(sizeof(double) * n);

  vecAssign(n, x, b);

  vecAssign(n, r, b);
  A->matVecMult(x, temp);
  vecDiffEqual(n, r, temp);

  rSqrLen = vecSqrLen(n, r);

  vecAssign(n, d, r);

  i = 0;
  if (*steps)
    iMax = *steps;
  else
    iMax = MAX_STEPS;

  if (rSqrLen > epsilon)
    while (i < iMax) {
      i++;
      A->matVecMult(d, t);
      u = vecDot(n, d, t);

      if (u == 0) {
        printf("(SolveConjGrad) d'Ad = 0\n");
        break;
      }

      // How far should we go?
      alpha = rSqrLen / u;

      // Take a step along direction d
      vecAssign(n, temp, d);
      vecTimesScalar(n, temp, alpha);
      vecAddEqual(n, x, temp);

      if (i & 0x3F) {
        vecAssign(n, temp, t);
        vecTimesScalar(n, temp, alpha);
        vecDiffEqual(n, r, temp);
      } else {
        // For stability, correct r every 64th iteration
        vecAssign(n, r, b);
        A->matVecMult(x, temp);
        vecDiffEqual(n, r, temp);
      }

      rSqrLenOld = rSqrLen;
      rSqrLen = vecSqrLen(n, r);

      // Converged! Let's get out of here
      if (rSqrLen <= epsilon)
        break;

      // Change direction: d = r + beta * d
      beta = rSqrLen/rSqrLenOld;
      vecTimesScalar(n, d, beta);
      vecAddEqual(n, d, r);
    }

  // free memory

  free(r);
  free(d);
  free(t);
  free(temp);

  *steps = i;
  return(rSqrLen);
}


