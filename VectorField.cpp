/*
-------------------------------------------------------------------------------
VectorField.cpp : Implementation of the vector field class.
-------------------------------------------------------------------------------
*/

#include "FieldHelper.hpp"
#include "VectorField.hpp"

// macros for allocating memory for the field
#define CREATE_DIM1 (new Vec2[(numCells+2)*(numCells+2)])
#define CREATE_DIM2 (new Vec2[(copyField->m_NumCells+2)       \
							   *(copyField->m_NumCells+2)])
#define CREATE_DIM  (new Vec2[(m_NumCells+2)*(m_NumCells+2)])

// macros for looping and indexing into the grid
#define IX_DIM(i,j) ((i)+(m_NumCells+2)*(j))
#define ITER_DIM    for(int i = 1; i <= m_NumCells; i++){     \
						for(int j = 1; j <= m_NumCells; j++){
#define ENDITER_DIM }}

#define VORT_SCALE 0.00001

// flags for bilinear interpolation
#define X_COMP 0
#define Y_COMP 1
#include "linearSolver.hpp"


	/***********************************************************************
 * Standard constructor.
 *
 * @param numCells The width and height of the grid in cells.
 * @param viscosity The viscosity of the fluid. (currently not used)
 * @param dt The time step of the simulation.
 ***********************************************************************/
VectorField::VectorField(int numCells, double viscosity, double dt) : m_NumCells(numCells)
																	, m_Viscosity(viscosity)
																	, m_Dt(dt)
{
	int gridDimention = numCells + 2;
	int gridSize = gridDimention * gridDimention;
	m_Field = new Vec2[gridSize];

	for (int i = 0; i < (m_NumCells + 2); i++)
	{
		for (int j = 0; j < (m_NumCells + 2); j++)
		{
			m_Field[IX_DIM(i,j)] = 0.0;
		}
	}
}


/*****************************************************************
 * Copy constructor.
 *
 * @param CopyField The field to copy.
 *****************************************************************/
VectorField::VectorField(VectorField *copyField) : m_NumCells(copyField->m_NumCells)
												 , m_Field(CREATE_DIM2)
												 , m_Viscosity(copyField->m_Viscosity)
												 , m_Dt(copyField->m_Dt)
{
	for (int i = 0; i < (m_NumCells + 2) * (m_NumCells + 2); i++)
	{
		m_Field[i][0] = (*copyField)[i][0];
		m_Field[i][1] = (*copyField)[i][1];
	}
}


/*******************************************************************
 * The deconstructor
 *******************************************************************/
VectorField::~VectorField(void)
{
	if (m_Field) delete m_Field;
}


/******************************************************************************
 * Advances the field by one timestep.
 * Adds forces from the user input and applies the advection,
 * vorticity confinement and projection steps.
 *
 * @param accelerationField A vector field representing forces from user input to be applied.
 * @param VelocityField		The velocity field to advect through (usually itself).
 ******************************************************************************/
void VectorField::TimeStep(VectorField *accelerationField, VectorField *velocityField)
{
	AddField(accelerationField);
	Advection(velocityField);
	VorticityConfinement();
	Projection();
}


/*******************************************************************
 * Adds another vector field to this one.
 *
 * @param a_SrcField The vector field to add to this one.
 *******************************************************************/
void VectorField::AddField(VectorField *srcField)
{
	for (int i = 0; i < ((m_NumCells + 2) * (m_NumCells + 2)); i++)
	{
		m_Field[i][0] += m_Dt * ((*srcField)[i][0]);
		m_Field[i][1] += m_Dt * ((*srcField)[i][1]);
	}
}


/*******************************************************************
 * Adds another vector field, stored as an array, to this one.
 *
 * @param a_SrcField The vector field to add to this one.
 *******************************************************************/
void VectorField::AddField(Vec2 *srcField)
{
	for (int i = 0; i < ((m_NumCells + 2) * (m_NumCells + 2)); i++)
	{
		m_Field[i][0] += m_Dt * srcField[i][0];
		m_Field[i][1] += m_Dt * srcField[i][1];
	}
}


/******************************************************************************
 * Applies Lagrangian advection to the field. This backtracks from the current
 * position using the velocity at the current point and then uses the value of
 * the field at the previous point as the new value.
 *
 * @param oldField The vector field before the current time step.
 *****************************************************************************/
void VectorField::Advection(VectorField *oldField)
{
	Vec2 prevPos, curPos;
	Vec2 vel;
	Vec2 *new_u = CREATE_DIM; // temporary field

	ITER_DIM
		// determine the velocity at the grid centers
		curPos = Vec2(i - 0.5, j - 0.5);
		vel = oldField->Interpolate(curPos);

		// backtrace with the interpolated velocity
		prevPos = curPos - (vel * m_Dt);

		// apply clipping at boundary if necessary
		FieldHelper::ClipPos(prevPos, curPos, vel, m_NumCells);

		// find the new grid centered velocities based on the previous position
		new_u[IX_DIM(i,j)] = oldField->Interpolate(prevPos);
	ENDITER_DIM

	ITER_DIM
		// interpolate the new grid centered velocities to the grid edges
		if (i == 1)
		{
			m_Field[IX_DIM(i,j)][0] = 0;
		}
		else
		{
			m_Field[IX_DIM(i,j)][0] = 0.5 * (new_u[IX_DIM(i,j)][0]
									+ new_u[IX_DIM(i-1,j)][0]);
		}

		if (j == 1)
		{
			m_Field[IX_DIM(i,j)][1] = 0;
		}
		else
		{
			m_Field[IX_DIM(i,j)][1] = 0.5 * (new_u[IX_DIM(i,j)][1]
									+ new_u[IX_DIM(i,j-1)][1]);
		}
	ENDITER_DIM

	delete new_u;
}


/*****************************************************************************
 * Applies vorticity confinement to the field. This finds a measure of the
 * vorticity at each point and then reinforces this to counter act numerical
 * diffusion.
 *****************************************************************************/
void VectorField::VorticityConfinement()
{
	// the curl
	double *curl = static_cast<double*>(malloc((m_NumCells + 2) * (m_NumCells + 2) * sizeof(double)));
	// the normalized gradient of the curl
	Vec2 *gradCurl = static_cast<Vec2*>(malloc((m_NumCells + 2) * (m_NumCells + 2) * sizeof(Vec2)));
	// the forces which will counteract numerical diffusion
	Vec2 *forces = static_cast<Vec2*>(malloc((m_NumCells + 2) * (m_NumCells + 2) * sizeof(Vec2)));

	if (!curl || !gradCurl || !forces)
	{
		fprintf(stderr, "malloc failed\n");
		exit(1);
	}

	// initialize to zero
	for (int i = 0; i < m_NumCells + 2; ++i)
	{
		curl[IX_DIM(i, 0)] = curl[IX_DIM(0, i)] = 0;
		curl[IX_DIM(i, m_NumCells+1)] = 0;
		curl[IX_DIM(m_NumCells+1, i)] = 0;
		gradCurl[IX_DIM(i, 0)] = gradCurl[IX_DIM(0, i)] = Vec2(0, 0);
		gradCurl[IX_DIM(i, m_NumCells+1)] = gradCurl[IX_DIM(m_NumCells+1, i)] = Vec2(0, 0);
		forces[IX_DIM(i, 0)] = forces[IX_DIM(0, i)] = Vec2(0, 0);
		forces[IX_DIM(i, m_NumCells+1)] = forces[IX_DIM(m_NumCells+1, i)] = Vec2(0, 0);
	}

	// calculate the magnitude of the curl of v at grid centers
	ITER_DIM
		double bottom, right, top, left;

		// check for boundaries
		if (i == 1)
		{
			left = Interpolate(Vec2(i - 0.5, j - 0.5))[1];
		}
		else
		{
			left = Interpolate(Vec2(i - 1.5, j - 0.5))[1];
		}

		if (i == m_NumCells)
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

		if (j == m_NumCells)
		{
			top = Interpolate(Vec2(i - 0.5, j - 0.5))[0];
		}
		else
		{
			top = Interpolate(Vec2(i - 0.5, j + 0.5))[0];
		}

		curl[IX_DIM(i, j)] = (right - left - top + bottom) / 2;
	ENDITER_DIM

	Vec2 grad_w;
	// calculate the normalized gradient of the curls on the faces
	ITER_DIM
		// check for boundaries
		if (i == 1)
		{
			if (j == 1)
			{
				grad_w = Vec2(0, 0);
			}
			else
			{
				grad_w = Vec2(0, curl[IX_DIM(i,j)] - curl[IX_DIM(i,j-1)]);
			}
		}
		else if (j == 1)
		{
			grad_w = Vec2(curl[IX_DIM(i,j)] - curl[IX_DIM(i-1,j)], 0);
		}
		else
		{
			grad_w = Vec2(curl[IX_DIM(i,j)] - curl[IX_DIM(i-1,j)],
						  curl[IX_DIM(i,j)] - curl[IX_DIM(i,j-1)]);
		}

		if (norm(grad_w) > EPSILON)
		{
			gradCurl[IX_DIM(i,j)] = grad_w / norm(grad_w);
		}
		else
		{
			gradCurl[IX_DIM(i,j)] = Vec2(0, 0);
		}
	ENDITER_DIM

	// average the normals to the grid centers
	ITER_DIM
		gradCurl[IX_DIM(i, j)] = Vec2(0.5 * (gradCurl[IX_DIM(i, j)][0] + gradCurl[IX_DIM(i+1, j)][0]),
									  0.5 * (gradCurl[IX_DIM(i, j)][1] + gradCurl[IX_DIM(i, j+1)][1]));
	ENDITER_DIM

	// compute the needed forces at the centers of each cell
	ITER_DIM
		forces[IX_DIM(i, j)] = VORT_SCALE * Vec2(gradCurl[IX_DIM(i, j)][1] * curl[IX_DIM(i, j)],
												-gradCurl[IX_DIM(i, j)][0] * curl[IX_DIM(i, j)]);
	ENDITER_DIM

	// average the forces to the faces of each cell
	ITER_DIM
		if (i == 1)
		{
			if (j == 1)
			{
				forces[IX_DIM(i,j)] = Vec2(0, 0);
			}
			else
			{
				double forceY = 0.5 * (forces[IX_DIM(i, j)][1] + forces[IX_DIM(i, j-1)][1]);
				forces[IX_DIM(i, j)] = Vec2(0, forceY);
			}
		}
		else if (j == 1)
		{
			double forceX = 0.5 * (forces[IX_DIM(i, j)][0] + forces[IX_DIM(i-1, j)][0]);
			forces[IX_DIM(i, j)] = Vec2(forceX, 0);
		}
		else
		{
			double forceX = 0.5 * (forces[IX_DIM(i, j)][0] + forces[IX_DIM(i-1, j)][0]);
			double forceY = 0.5 * (forces[IX_DIM(i, j)][1] + forces[IX_DIM(i, j-1)][1]);
			forces[IX_DIM(i, j)] = Vec2(forceX, forceY);
		}
	ENDITER_DIM

	// apply forces to the field
	AddField(forces);

	free(curl);
	free(gradCurl);
	free(forces);
}


/******************************************************************************
 * Applies the projection step to the field. The projection step solves for a
 * field which has zero divergence. This is done by solving a linear system of
 * laplacian(pressure) = divergence(field) for the pressure. Then the zero
 * divergence field is the current field - gradient(pressure).
 ******************************************************************************/
void VectorField::Projection()
{
	// pressure
	double *pressure = static_cast<double*>(malloc((m_NumCells + 2) * (m_NumCells + 2) * sizeof(double)));
	// divergence of the field
	double *div_v = static_cast<double*>(malloc((m_NumCells + 2) * (m_NumCells + 2) * sizeof(double)));

	// calculate divergence of the current velocity field
	// and initializes the pressure
	ITER_DIM
		pressure[IX_DIM(i, j)] = 0;
		div_v[IX_DIM(i, j)] = m_Field[IX_DIM(i+1, j)][0]
							- m_Field[IX_DIM(i, j)][0]
							+ m_Field[IX_DIM(i, j+1)][1]
							- m_Field[IX_DIM(i, j)][1];
	ENDITER_DIM

	// set boundaries explicitly to zero
	for (int i = 0; i <= m_NumCells + 1; ++i)
	{
		// bottom
		pressure[IX_DIM(i, 0)] = div_v[IX_DIM(i, 0)] = 0;
		// left
		pressure[IX_DIM(0, i)] = div_v[IX_DIM(0, i)] = 0;
		// top
		pressure[IX_DIM(i, m_NumCells+1)] = 0;
		div_v[IX_DIM(i, m_NumCells+1)] = 0;
		// right
		pressure[IX_DIM(m_NumCells+1, i)] = 0;
		div_v[IX_DIM(m_NumCells+1, i)] = 0;
	}

	// initialize laplacian matrix
	ImplicitMatrix *Lap_p = new ImplicitMatrixLap(m_NumCells);

	// solve for pressure using conjugate gradient
	int steps = MAX_STEPS;
	double err = ConjGrad((m_NumCells + 2) * (m_NumCells + 2),
						  Lap_p,
						  pressure,
						  div_v,
						  EPSILON,
						  &steps);
	if (err > EPSILON)
	{
		fprintf(stderr, "solution for projection did not converge\n");
		exit(1);
	}

	// calculate the new velocity field as v-grad(p)
	ITER_DIM
		// checks for boundaries
		if (i == 1)
		{
			if (j == 1)
			{
				m_Field[IX_DIM(i, j)] = Vec2(0, 0);
			}
			else
			{
				m_Field[IX_DIM(i, j)] = m_Field[IX_DIM(i, j)]
									  - Vec2(0, pressure[IX_DIM(i, j)]
									  - pressure[IX_DIM(i, j-1)]);
			}
		}
		else
		{
			if (j == 1)
			{
				m_Field[IX_DIM(i, j)] = m_Field[IX_DIM(i, j)]
									  - Vec2(pressure[IX_DIM(i, j)]
									  - pressure[IX_DIM(i-1, j)], 0);
			}
			else
			{
				m_Field[IX_DIM(i, j)] = m_Field[IX_DIM(i, j)]
									  - Vec2(pressure[IX_DIM(i, j)] - pressure[IX_DIM(i-1, j)],
											 pressure[IX_DIM(i, j)] - pressure[IX_DIM(i, j-1)]);
			}
		}
	ENDITER_DIM

	free(pressure);
	free(div_v);
	delete(Lap_p);
}


/******************************************************************************
 * Interpolate the velocity at the given position using bilinear interpolation.
 * Each component is interpolated independently as a weighted average of the
 * velocities stored at the 4 nearest edge centers.
 *
 * @param pos The position to get the velocity at.
 * @return The velocity at the given position.
 ******************************************************************************/
Vec2 VectorField::Interpolate(Vec2 pos)
{
	double Vx, Vy;

	// interpolate the x component
	int low_x = static_cast<int>(floor(pos[0] + 1) + EPSILON);
	int low_y = static_cast<int>(pos[1] + 0.5 + EPSILON);
	Vx = BilinearInterpolation(pos, low_x, low_y, X_COMP);

	// interpolate the y component
	low_x = static_cast<int>(pos[0] + 0.5);
	low_y = static_cast<int>(floor(pos[1] + 1) + EPSILON);
	Vy = BilinearInterpolation(pos, low_x, low_y, Y_COMP);

	return Vec2(Vx, Vy);
}


/******************************************************************************
 * Interpolate the given component of velocity using bilinear interpolation.
 *
 * @param pos The position to interpolate the value at.
 * @param low_x The x index of the cell whose x component of velocity is the
 *              closest on the left to the given position.
 * @param low_y The y index of the cell whose y component of velocity is the
 *              closest on the bottom to the given position.
 * @param comp Flag for which component is being interpolated.
 *             0 for x and 1 for y.
 * @return The interpolated value.
 *****************************************************************************/
double VectorField::BilinearInterpolation(Vec2 pos, int low_x,
										  int low_y, int comp)
{
	double x_low_inter, x_high_inter;

	// interpolate between y values at the lower x positions
	x_low_inter = ((low_y + 0.5 * (1 - comp)) - pos[1])
		* m_Field[IX_DIM(low_x,low_y)][comp]
		+ (pos[1] - (low_y - 0.5 - 0.5 * comp))
		* m_Field[IX_DIM(low_x, low_y+1)][comp];
	// interpolate between y values at the higher x positions
	x_high_inter = ((low_y + 0.5 * (1 - comp)) - pos[1])
		* m_Field[IX_DIM(low_x+1,low_y)][comp]
		+ (pos[1] - (low_y - 0.5 - 0.5 * comp))
		* m_Field[IX_DIM(low_x+1, low_y+1)][comp];

	// interpolate between the averaged values at the low and high x positions
	return (low_x + 0.5 * comp - pos[0]) * x_low_inter
		+ (pos[0] - (low_x - 0.5 - 0.5 * (1 - comp))) * x_high_inter;
}
