/*
----------------------------------------------------------------------
ScalarField.cpp : Implementation of the scalar field class.
----------------------------------------------------------------------
*/

#include "FieldHelper.hpp"
#include "ScalarField.hpp"
#include "VectorField.hpp"

#include <math.h>

// macros for allocating memory for the field
#define CREATE_DIM1  (new double[(numCells + 2)*(numCells + 2)])
#define CREATE_DIM2  (new double[(copyField->m_NumCells + 2)        \
								 *(copyField->m_NumCells + 2)])
#define CREATE_DIM   (new double[(m_NumCells + 2)*(m_NumCells + 2)])

// macros for looping and indexing into the grid
#define IX_DIM(i, j) ((i)+(m_NumCells+2)*(j))
#define ITER_DIM     for(int i = 1; i <= m_NumCells; i++){         \
						 for(int j = 1; j <= m_NumCells; j++){
#define ENDITER_DIM  }}

/**********************************************************************
 * Standard constructor.
 * @param a_NumCells The width and height of the grid in cells.
 * @param a_Viscosity The viscosity of the fluid. (currently not used)
 * @param a_Dt The time step of the simulation.
 **********************************************************************/
ScalarField::ScalarField(int numCells, double viscosity, double dt) : m_NumCells(numCells)
																	, m_Field(CREATE_DIM1)
																	, m_Viscosity(viscosity)
																	, m_Dt(dt)
{
	for (int i = 0; i < (m_NumCells + 2) * (m_NumCells + 2); i++)
	{
		m_Field[i] = 0.f;
	}
}

/************************************************************
 * Constructor used to make a copy of another scalar field.
 * @param The field to copy.
 *************************************************************/
ScalarField::ScalarField(ScalarField *copyField) : m_NumCells(copyField->m_NumCells)
												 , m_Field(CREATE_DIM2)
												 , m_Viscosity(copyField->m_Viscosity)
												 , m_Dt(copyField->m_Dt)
{
	for (int i = 0; i < (m_NumCells + 2) * (m_NumCells + 2); i++)
	{
		m_Field[i] = (*copyField)[i];
	}
}

/***********************
 * The deconstructor
 ***********************/
ScalarField::~ScalarField(void)
{
	delete m_Field;
}

/******************************************************************************
 * Advances the field by one timestep. Adds any new scalar values from the user
 * input and applies the advection step.
 *
 * @param srcField		A scalar field representing the change in scalar values from user
 *						input to be added.
 * @param VelocityField The velocity field to advect the scalar through.
 ******************************************************************************/
void ScalarField::TimeStep(ScalarField *srcField, VectorField *velocityField)
{
	AddField(srcField);
	Advection(velocityField);
}

/******************************************************************************
 * Adds another scalar field to this one over one timestep.
 *
 * @param a_SrcField A scalar field representing the change in scalar values.
 ******************************************************************************/
void ScalarField::AddField(ScalarField *srcField)
{
	for (int i = 0; i < ((m_NumCells + 2) * (m_NumCells + 2)); i++)
	{
		m_Field[i] += m_Dt * ((*srcField)[i]);
	}
}


/******************************************************************************
 * Applies Lagrangian advection to the field. This backtracks from the current
 * position using the velocity at the current point and then uses the value of
 * the field at the previous point as the new value.
 *
 * @param u The vector field to advect the scalar field through.
 ******************************************************************************/
void ScalarField::Advection(VectorField *u)
{
	double *newField = CREATE_DIM; // temporary field

	ITER_DIM
		// interpolate the velocity at the center of the current grid cell
		Vec2 curPos = Vec2(i - 0.5, j - 0.5);
		Vec2 vel = u->Interpolate(curPos);

		// take a step backward
		Vec2 prevPos = curPos - vel * m_Dt;

		// clip the trajectory at the boundary if necessary.
		FieldHelper::ClipPos(prevPos, curPos, vel, m_NumCells);

		// set the scalar to the interpolated scalar at the earlier position
		newField[IX_DIM(i,j)] = Interpolate(prevPos);
	ENDITER_DIM

	/* set boundaries explicitly to the interior edges of the grid so that
	   interpolation at the boundary returns the correct scalar value and
	   conserves the scalar quantity. */
	for (int i = 1; i <= m_NumCells; ++i)
	{
		// bottom
		newField[IX_DIM(i, 0)] = newField[IX_DIM(i, 1)];
		// left
		newField[IX_DIM(0, i)] = newField[IX_DIM(1, i)];
		// top
		newField[IX_DIM(i, m_NumCells+1)] = newField[IX_DIM(i, m_NumCells)];
		// right
		newField[IX_DIM(m_NumCells+1, i)] = newField[IX_DIM(m_NumCells, i)];
	}
	// bottom-left
	newField[IX_DIM(0, 0)] = newField[IX_DIM(1, 1)];
	// top-left
	newField[IX_DIM(0, m_NumCells+1)] = newField[IX_DIM(1, m_NumCells)];
	// bottom-right
	newField[IX_DIM(m_NumCells+1, 0)] = newField[IX_DIM(m_NumCells, 1)];
	// top-right
	newField[IX_DIM(m_NumCells+1, m_NumCells+1)] = newField[IX_DIM(m_NumCells, m_NumCells)];

	// update scalar field
	ITER_DIM
		this->m_Field[IX_DIM(i,j)] = newField[IX_DIM(i,j)];
	ENDITER_DIM

	delete newField;
}

/******************************************************************************
 * Interpolates the scalar value at the given position on the grid using
 * bilinear interpolation.
 *
 * @param pos The position to get the value at.
 * @return The interpolated scalar value at the given position
 ******************************************************************************/
double ScalarField::Interpolate(const Vec2 &pos)
{
	// calculate the indices of the cells corresponding to the position
	int low_x = static_cast<int>(pos[0] + 0.5);
	int low_y = static_cast<int>(pos[1] + 0.5);

	double scalar_low_x, scalar_high_x;
	// interpolate between y values at the lower x positions
	scalar_low_x = ((low_y + 0.5f) - pos[1]) * m_Field[IX_DIM(low_x,low_y)]
				 + (pos[1] - (low_y - 0.5f)) * m_Field[IX_DIM(low_x, low_y+1)];
	// interpolate between y values at the higher x positions
	scalar_high_x = ((low_y + 0.5f) - pos[1]) * m_Field[IX_DIM(low_x+1,low_y)]
				  + (pos[1] - (low_y - 0.5f)) * m_Field[IX_DIM(low_x+1, low_y+1)];

	// interpolate between the averaged values at the low and high x positions
	return ((low_x + 0.5f) - pos[0]) * scalar_low_x
		+ (pos[0] - (low_x - 0.5f)) * scalar_high_x;
}
