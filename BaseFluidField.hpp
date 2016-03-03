#pragma once

#include <cstring>
#include <gfx/vec2.h>

#define EPSILON 1.0e-10


/*****************************************************************************
 * Helper function for 2D accesses into a field.
 *****************************************************************************/
inline int Idx2DTo1D(int size, int x, int y)
{
	return x + (size * y);
}

/*****************************************************************************
 * Abstract class representing a fluid field.
 * It uses a MAC grid representation with values defined
 * at the center of grid points. Has a boundary of one cell width in order to
 * enforce the boundary conditions when interpolating.
 *****************************************************************************/
template<class T>
class FluidField
{
public:
	/*************************************************************************
	 * Standard constructor.
	 *
	 * @param size      The width and height of the grid in cells.
	 * @param viscosity The viscosity of the fluid. (currently not used)
	 * @param dt        The time step of the simulation.
	 *************************************************************************/
	FluidField(int size, double viscosity, double dt): m_GridSize(size)
													 , m_TotalNumCells(m_GridSize * m_GridSize)
													 , m_Viscosity(viscosity)
													 , m_Dt(dt)
	{
		m_Field = new T[m_TotalNumCells]();
	}

	virtual ~FluidField()
	{
		delete[] m_Field;
	}

	/*************************************************************************
	 * Adds another field to this one.
	 * If they are not the same size, nothing happens.
	 *
	 * @param sourceField The vector field to add to this one.
	 *************************************************************************/
	virtual void AddField(const FluidField& sourceField)
	{
		if(sourceField.m_GridSize != m_GridSize)
		{
			printf("Tried to add two fields with different sizes, %d and %d",
				   m_GridSize,
				   sourceField.m_GridSize);
			return;
		}

		for (int i = 0; i < m_TotalNumCells; i++)
		{
			m_Field[i] += m_Dt * sourceField[i];
		}
	}

	/*************************************************************************
	 * Index accessor to the ith grid point.
	 *************************************************************************/
	bool IsPositionWithinBounds(int x, int y)
	{
		return ((x > 0) && (x < m_GridSize - 1)
			&& (y > 0) && (y < m_GridSize - 1));
	}

	/*************************************************************************
	 * Index accessor to the ith grid point.
	 *************************************************************************/
	T& operator[](int i) const
	{
		return m_Field[i];
	}

	int Idx2DTo1D(int x, int y) const
	{
		return ::Idx2DTo1D(m_GridSize, x, y);
	}

	/*************************************************************************
	 * Coordinate accessor to the grid point at (x, y).
	 *************************************************************************/
	T& operator()(int x, int y) const
	{
		return m_Field[Idx2DTo1D(x, y)];
	}

	/*************************************************************************
	 * Interpolates the field at the given position on the grid using
	 * bilinear interpolation.
	 *
	 * @param pos The position to get the value at.
	 * @return The interpolated field value at the given position
	 *************************************************************************/
	virtual T Interpolate(const Vec2 &pos) const = 0;
	
protected:
	/*************************************************************************
	 * Helper function for the field's advection step that calculates
	 * the position for which the new value will be interpolated at by
	 * clipping the trajectory at the boundary if necessary. Clipping at the
	 * boundary means only values strictly within and not on the boundary
	 * are valid.
	 * 
	 * @param prevPosGuess[in,out]  The position to be clipped. It is expected
	 *								to be -vel*dt from the current position.
	 * @param curPos                The current position being considered.
	 * @param vel                   The velocity at the current position.
	 * @param gridSize              The width and length of the grid in cells.
	 *
	 * @return The clipped position.
	 *************************************************************************/
	void ClipPos(Vec2 &prevPosGuess, const Vec2 &curPos, const Vec2 &vel, int gridSize) const
	{
		if (prevPosGuess[0] > 0 && prevPosGuess[0] <= gridSize
			&& prevPosGuess[1] > 0 && prevPosGuess[1] <= gridSize)
		{
			// The guess position is inbounds so no need to clip
			return;
		}

		// Check for components of velocity being zero since we will divide by them
		bool xVelIsZero = fabs(vel[0]) < EPSILON;
		bool yVelIsZero = fabs(vel[1]) < EPSILON;

		if (xVelIsZero && yVelIsZero)
		{
			// If the velocity has zero magnitude, we can't clip it
			return;
		}

		// Initialize to max to handle cases where the velocity in that axis is 0
		double timeToBoundaryInX = std::numeric_limits<double>::max();
		double timeToBoundaryInY = timeToBoundaryInX;

		if (!xVelIsZero)
		{
			timeToBoundaryInX = CalcTimeToBoundaryAlongAxis(0, vel, curPos, gridSize);
		}

		if (!yVelIsZero)
		{
			timeToBoundaryInY = CalcTimeToBoundaryAlongAxis(1, vel, curPos, gridSize);
		}

		// Correct for the direction that goes out of bounds first
		double timeToBoundary = timeToBoundaryInX;
		if (timeToBoundaryInY < timeToBoundaryInX)
		{
			timeToBoundary = timeToBoundaryInY;
		}

		prevPosGuess = curPos - (vel * timeToBoundary);
	}

	double CalcTimeToBoundaryAlongAxis(int axis,
									   const Vec2 &vel,
									   const Vec2 &curPos,
									   int gridSize) const
	{
		// Out of high bounds (the position we want is at -vel)
		if (vel[axis] < 0)
		{
			return (curPos[axis] - gridSize) / vel[axis];
		}
		else // Out of low bounds
		{
			return curPos[axis] / vel[axis];
		}
	}

	int m_GridSize;
	int m_TotalNumCells;
	T *m_Field;
	double m_Viscosity;
	double m_Dt;
};
