/*
-------------------------------------------------------------------------------
VectorField.cpp : Implementation of the vector field class.
-------------------------------------------------------------------------------
*/

#include "BaseFluidField.hpp"
#include "LinearSolver.hpp"
#include "VectorField.hpp"


// Macros for looping and indexing into the grid
#define ITER_DIM    for(int i = 1; i < m_GridSize - 1; i++){     \
						for(int j = 1; j < m_GridSize - 1; j++){
#define ENDITER_DIM }}

#define VORT_SCALE 0.00001

// Flags for bilinear interpolation
#define X_COMP 0
#define Y_COMP 1


void VectorField::TimeStep(const VectorField &accelerationField, const VectorField &velocityField)
{
	AddField(accelerationField);
	Advection(velocityField);
	VorticityConfinement();
	Projection();
}

void VectorField::AddVectorField(const Vec2 sourceField[])
{
	for (int i = 0; i < m_TotalNumCells; i++)
	{
		m_Field[i] += m_Dt * sourceField[i];
	}
}

void VectorField::Advection(const VectorField &oldField)
{
	Vec2 prevPos, curPos;
	Vec2 vel;
	// Temporary field
	Vec2 *newField = new Vec2[m_TotalNumCells];

	ITER_DIM
		// Determine the velocity at the grid centers
		curPos = Vec2(i - 0.5, j - 0.5);
		vel = oldField.Interpolate(curPos);

		// Backtrack with the interpolated velocity
		prevPos = curPos - (vel * m_Dt);

		// Apply clipping at boundary if necessary
		ClipPos(prevPos, curPos, vel, m_GridSize);

		// Find the new grid centered velocities based on the previous position
		newField[Idx2DTo1D(i, j)] = oldField.Interpolate(prevPos);
	ENDITER_DIM

	ITER_DIM
		// Interpolate the new grid centered velocities to the grid edges
		if (i == 1)
		{
			m_Field[Idx2DTo1D(i, j)][0] = 0;
		}
		else
		{
			m_Field[Idx2DTo1D(i, j)][0] = 0.5 * (newField[Idx2DTo1D(i, j)][0]
										+ newField[Idx2DTo1D(i - 1, j)][0]);
		}

		if (j == 1)
		{
			m_Field[Idx2DTo1D(i, j)][1] = 0;
		}
		else
		{
			m_Field[Idx2DTo1D(i, j)][1] = 0.5 * (newField[Idx2DTo1D(i, j)][1]
										+ newField[Idx2DTo1D(i, j - 1)][1]);
		}
	ENDITER_DIM

	delete[] newField;
}

void VectorField::VorticityConfinement()
{
	// The curl
	double *curl = new double[m_TotalNumCells]();
	// The normalized gradient of the curl
	Vec2 *gradCurl = new Vec2[m_TotalNumCells]();
	// The forces which will counteract numerical diffusion
	Vec2 *forces = new Vec2[m_TotalNumCells]();

	// Calculate the magnitude of the curl of v at grid centers
	ITER_DIM
		double bottom, right, top, left;

		// Check for boundaries
		if (i == 1)
		{
			left = Interpolate(Vec2(i - 0.5, j - 0.5))[1];
		}
		else
		{
			left = Interpolate(Vec2(i - 1.5, j - 0.5))[1];
		}

		if (i == m_GridSize - 2)
		{
			right = Interpolate(Vec2(i - 0.5, j - 0.5))[1];
		}
		else
		{
			right = Interpolate(Vec2(i + 0.5, j - 0.5))[1];
		}

		if (j == 1)
		{
			bottom = Interpolate(Vec2(i - 0.5, j - 0.5))[0];
		}
		else
		{
			bottom = Interpolate(Vec2(i - 0.5, j - 1.5))[0];
		}

		if (j == m_GridSize - 2)
		{
			top = Interpolate(Vec2(i - 0.5, j - 0.5))[0];
		}
		else
		{
			top = Interpolate(Vec2(i - 0.5, j + 0.5))[0];
		}

		curl[Idx2DTo1D(i, j)] = (right - left - top + bottom) / 2;
	ENDITER_DIM

	Vec2 grad_w;
	// Calculate the normalized gradient of the curls on the faces
	ITER_DIM
		// Check for boundaries
		if (i == 1)
		{
			if (j == 1)
			{
				grad_w = Vec2(0, 0);
			}
			else
			{
				grad_w = Vec2(0, curl[Idx2DTo1D(i, j)] - curl[Idx2DTo1D(i, j - 1)]);
			}
		}
		else if (j == 1)
		{
			grad_w = Vec2(curl[Idx2DTo1D(i, j)] - curl[Idx2DTo1D(i - 1, j)], 0);
		}
		else
		{
			grad_w = Vec2(curl[Idx2DTo1D(i, j)] - curl[Idx2DTo1D(i - 1, j)],
						  curl[Idx2DTo1D(i, j)] - curl[Idx2DTo1D(i, j - 1)]);
		}

		if (norm(grad_w) > EPSILON)
		{
			gradCurl[Idx2DTo1D(i, j)] = grad_w / norm(grad_w);
		}
		else
		{
			gradCurl[Idx2DTo1D(i, j)] = Vec2(0, 0);
		}
	ENDITER_DIM

	// Average the normals to the grid centers
	ITER_DIM
		gradCurl[Idx2DTo1D(i, j)] = Vec2(0.5 * (gradCurl[Idx2DTo1D(i, j)][0] + gradCurl[Idx2DTo1D(i + 1, j)][0]),
										 0.5 * (gradCurl[Idx2DTo1D(i, j)][1] + gradCurl[Idx2DTo1D(i, j + 1)][1]));
	ENDITER_DIM

	// Compute the needed forces at the centers of each cell
	ITER_DIM
		forces[Idx2DTo1D(i, j)] = VORT_SCALE * Vec2(gradCurl[Idx2DTo1D(i, j)][1] * curl[Idx2DTo1D(i, j)],
													-gradCurl[Idx2DTo1D(i, j)][0] * curl[Idx2DTo1D(i, j)]);
	ENDITER_DIM

	// Average the forces to the faces of each cell
	ITER_DIM
		if (i == 1)
		{
			if (j == 1)
			{
				forces[Idx2DTo1D(i, j)] = Vec2(0, 0);
			}
			else
			{
				double forceY = 0.5 * (forces[Idx2DTo1D(i, j)][1] + forces[Idx2DTo1D(i, j - 1)][1]);
				forces[Idx2DTo1D(i, j)] = Vec2(0, forceY);
			}
		}
		else if (j == 1)
		{
			double forceX = 0.5 * (forces[Idx2DTo1D(i, j)][0] + forces[Idx2DTo1D(i - 1, j)][0]);
			forces[Idx2DTo1D(i, j)] = Vec2(forceX, 0);
		}
		else
		{
			double forceX = 0.5 * (forces[Idx2DTo1D(i, j)][0] + forces[Idx2DTo1D(i - 1, j)][0]);
			double forceY = 0.5 * (forces[Idx2DTo1D(i, j)][1] + forces[Idx2DTo1D(i, j - 1)][1]);
			forces[Idx2DTo1D(i, j)] = Vec2(forceX, forceY);
		}
	ENDITER_DIM

	// Apply forces to the field
	AddVectorField(forces);

	delete[] curl;
	delete[] gradCurl;
	delete[] forces;
}

void VectorField::Projection()
{
	double *pressure = new double[m_TotalNumCells]();
	double *divergence = new double[m_TotalNumCells]();

	// Calculate the divergence of the current velocity field
	ITER_DIM
		divergence[Idx2DTo1D(i, j)] = m_Field[Idx2DTo1D(i + 1, j)][0]
									- m_Field[Idx2DTo1D(i, j)][0]
									+ m_Field[Idx2DTo1D(i, j + 1)][1]
									- m_Field[Idx2DTo1D(i, j)][1];
	ENDITER_DIM

	// Initialize Laplacian matrix
	ImplicitMatrixLap laplacian(m_GridSize);

	// Solve for pressure using conjugate gradient
	int steps = MAX_STEPS;
	double err = ConjGrad(m_TotalNumCells,
						  &laplacian,
						  pressure,
						  divergence,
						  EPSILON,
						  &steps);
	if (err > EPSILON)
	{
		fprintf(stderr, "solution for projection did not converge!\n");
		exit(1);
	}

	// Calculate the new velocity field as v - grad(p)
	ITER_DIM
		// Checks for boundaries
		if (i == 1)
		{
			if (j == 1)
			{
				m_Field[Idx2DTo1D(i, j)] = Vec2(0, 0);
			}
			else
			{
				double yPressureDelta = pressure[Idx2DTo1D(i, j)]
									  - pressure[Idx2DTo1D(i, j - 1)];
				m_Field[Idx2DTo1D(i, j)] = m_Field[Idx2DTo1D(i, j)]
										 - Vec2(0, yPressureDelta);
			}
		}
		else
		{
			if (j == 1)
			{
				m_Field[Idx2DTo1D(i, j)] = m_Field[Idx2DTo1D(i, j)]
										 - Vec2(pressure[Idx2DTo1D(i, j)]
										 - pressure[Idx2DTo1D(i - 1, j)], 0);
			}
			else
			{
				double xPressureDelta = pressure[Idx2DTo1D(i, j)]
									  - pressure[Idx2DTo1D(i - 1, j)];
				double yPressureDelta = pressure[Idx2DTo1D(i, j)]
									  - pressure[Idx2DTo1D(i, j - 1)];
				m_Field[Idx2DTo1D(i, j)] = m_Field[Idx2DTo1D(i, j)]
										 - Vec2(xPressureDelta, yPressureDelta);
			}
		}
	ENDITER_DIM

	delete[] pressure;
	delete[] divergence;
}

Vec2 VectorField::Interpolate(const Vec2 &pos) const
{
	double Vx, Vy;

	// Interpolate the x component
	int low_x = static_cast<int>(floor(pos[0] + 1) + EPSILON);
	int low_y = static_cast<int>(pos[1] + 0.5 + EPSILON);
	Vx = BilinearInterpolation(pos, low_x, low_y, X_COMP);

	// Interpolate the y component
	low_x = static_cast<int>(pos[0] + 0.5);
	low_y = static_cast<int>(floor(pos[1] + 1) + EPSILON);
	Vy = BilinearInterpolation(pos, low_x, low_y, Y_COMP);

	return Vec2(Vx, Vy);
}

double VectorField::BilinearInterpolation(const Vec2 &pos,
										  int low_x,
										  int low_y,
										  int comp) const
{
	double x_low_inter, x_high_inter;

	// Interpolate between y values at the lower x positions
	x_low_inter = (((low_y + 0.5 * (1 - comp)) - pos[1])
					* m_Field[Idx2DTo1D(low_x, low_y)][comp])
				+ ((pos[1] - (low_y - 0.5 - 0.5 * comp))
					* m_Field[Idx2DTo1D(low_x, low_y + 1)][comp]);
	// Interpolate between y values at the higher x positions
	x_high_inter = (((low_y + 0.5 * (1 - comp)) - pos[1])
					* m_Field[Idx2DTo1D(low_x + 1, low_y)][comp])
				 + ((pos[1] - (low_y - 0.5 - 0.5 * comp))
					* m_Field[Idx2DTo1D(low_x + 1, low_y + 1)][comp]);

	// Interpolate between the averaged values at the low and high x positions
	return ((low_x + 0.5 * comp - pos[0]) * x_low_inter)
		+ ((pos[0] - (low_x - 0.5 - 0.5 * (1 - comp))) * x_high_inter);
}
