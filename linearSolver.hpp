#include <math.h> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// Karen's CGD

#define MAX_STEPS 100000

// Matrix class the solver will accept
class implicitMatrix
{
 public:
    virtual ~implicitMatrix() { }
    virtual void matVecMult(double x[], double r[]) = 0;

};

/*******************************************************************************
 * This implicit matrix represents a discrete laplacian operator.
 ******************************************************************************/
class implicitMatrixLap : public implicitMatrix
{
public:
    implicitMatrixLap(int m_NumCells);
    virtual ~implicitMatrixLap();
    virtual void matVecMult(double x[], double r[]);

private:
    // the width and length of cells of the grid
    int m_NumCells;
};



// Solve Ax = b for a symmetric, positive definite matrix A
// A is represented implicitly by the function "matVecMult"
// which performs a matrix vector multiple Av and places result in r
// "n" is the length of the vectors x and b
// "epsilon" is the error tolerance
// "steps", as passed, is the maximum number of steps, or 0 (implying MAX_STEPS)
// Upon completion, "steps" contains the number of iterations taken
double ConjGrad(int n, implicitMatrix *A, double x[], double b[], 
        double epsilon, // how low should we go?
        int    *steps);

// Some vector helper functions
void vecAddEqual(int n, double r[], double v[]);
void vecDiffEqual(int n, double r[], double v[]);
void vecAssign(int n, double v1[], double v2[]);
void vecTimesScalar(int n, double v[], double s);
double vecDot(int n, double v1[], double v2[]);
double vecSqrLen(int n, double v[]);
