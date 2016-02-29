/*
----------------------------------------------------------------------
ScalarField.cpp : Implementation of the scalar field class.
----------------------------------------------------------------------
*/

#include "ScalarField.hpp"
#include "VectorField.hpp"

#include <math.h>


// macros for looping and indexing into the grid
#define IX_DIM(i, j) ((i) + (m_GridSize * (j)))
#define ITER_DIM     for(int i = 1; i < m_GridSize - 1; i++){         \
						 for(int j = 1; j < m_GridSize - 1; j++){
#define ENDITER_DIM  }}


void ScalarField::TimeStep(const ScalarField& srcField, const VectorField& velocityField)
{
	AddField(srcField);
	Advection(velocityField);
}

void ScalarField::Advection(const VectorField& vectorField)
{
	// Temporary field
	double *newField = new double[m_TotalNumCells];

	ITER_DIM
		// Interpolate the velocity at the center of the current grid cell
		Vec2 curPos = Vec2(i - 0.5, j - 0.5);
		Vec2 vel = vectorField.Interpolate(curPos);

		// Take a step backward
		Vec2 prevPos = curPos - vel * m_Dt;

		// Clip the trajectory at the boundary if necessary.
		ClipPos(prevPos, curPos, vel, m_GridSize);

		// Set the scalar to the interpolated scalar at the earlier position
		newField[IX_DIM(i, j)] = Interpolate(prevPos);
	ENDITER_DIM

	// Set boundaries explicitly to the interior edges of the grid so that
	// interpolation at the boundary returns the correct scalar value and
	// conserves the scalar quantity.
	for (int i = 1; i < m_GridSize - 1; ++i)
	{
		// Bottom
		newField[IX_DIM(i, 0)] = newField[IX_DIM(i, 1)];
		// Left
		newField[IX_DIM(0, i)] = newField[IX_DIM(1, i)];
		// Top
		newField[IX_DIM(i, m_GridSize - 1)] = newField[IX_DIM(i, m_GridSize - 2)];
		// Right
		(m_GridSize - 1, i) = (m_GridSize - 2, i);
	}
	// Bottom-left
	newField[IX_DIM(0, 0)] = newField[IX_DIM(1, 1)];
	// Top-left
	newField[IX_DIM(0, m_GridSize - 1)] = newField[IX_DIM(1, m_GridSize - 2)];
	// Bottom-right
	newField[IX_DIM(m_GridSize - 1, 0)] = newField[IX_DIM(m_GridSize - 2, 1)];
	// Top-right
	newField[IX_DIM(m_GridSize - 1, m_GridSize - 1)] = newField[IX_DIM(m_GridSize - 2, m_GridSize - 2)];

	// Update scalar field
	ITER_DIM
		m_Field[IX_DIM(i,j)] = newField[IX_DIM(i,j)];
	ENDITER_DIM

	delete[] newField;
}

double ScalarField::Interpolate(const Vec2 &pos) const
{
	// Calculate the indices of the cells corresponding to the position
	int low_x = static_cast<int>(pos[0] + 0.5);
	int low_y = static_cast<int>(pos[1] + 0.5);

	double scalar_low_x, scalar_high_x;
	// Interpolate between y values at the lower x positions
	scalar_low_x = ((low_y + 0.5f) - pos[1]) * m_Field[IX_DIM(low_x, low_y)]
				 + (pos[1] - (low_y - 0.5f)) * m_Field[IX_DIM(low_x, low_y+1)];
	// Interpolate between y values at the higher x positions
	scalar_high_x = ((low_y + 0.5f) - pos[1]) * m_Field[IX_DIM(low_x+1, low_y)]
				  + (pos[1] - (low_y - 0.5f)) * m_Field[IX_DIM(low_x+1, low_y+1)];

	// Interpolate between the averaged values at the low and high x positions
	return ((low_x + 0.5f) - pos[0]) * scalar_low_x
		+ (pos[0] - (low_x - 0.5f)) * scalar_high_x;
}
