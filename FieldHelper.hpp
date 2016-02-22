#pragma once

#include <gfx/vec2.h>

#define EPSILON 1.0e-10

/*************************************************************************************************
 * Class wrapper for a helper function used in both the vector and scalar field classes.
 * This makes sure multiple instances of this function are not created each time it is included,
 * such as in StableFluids.cpp as well as ScalarField.cpp.
 *************************************************************************************************/
class FieldHelper
{
public:
	/*********************************************************************************************
	 * Helper function for the field's advection step that calculates the position for which the
	 * new value will be interpolated at by clipping the trajectory at the boundary if necessary.
	 * The positions are relative to the bottom left corner of the cell (1,1)
	 * 
	 * @param prevPosGuess[in,out]  The position to be clipped.
	 *                              It is expected to be -vel*dt from the current position.
	 * @param curPos                The current position being considered.
	 * @param vel                   The velocity at the current position.
	 * @param numCells              The width and length of the grid in cells.
	 *
	 * @return The clipped position.
	 *********************************************************************************************/
	static void ClipPos(Vec2 &prevPosGuess, const Vec2 &curPos, const Vec2 &vel, int numCells)
	{
		if (prevPosGuess[0] >= 0 && prevPosGuess[0] <= numCells
			&& prevPosGuess[1] >= 0 && prevPosGuess[1] <= numCells)
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
			timeToBoundaryInX = CalcTimeToBoundaryAlongAxis(0, vel, curPos, numCells);
		}

		if (!yVelIsZero)
		{
			timeToBoundaryInY = CalcTimeToBoundaryAlongAxis(1, vel, curPos, numCells);
		}

		// Correct for the direction that goes out of bounds first
		double timeToBoundary = timeToBoundaryInX;
		if (timeToBoundaryInY < timeToBoundaryInX)
		{
			timeToBoundary = timeToBoundaryInY;
		}

		prevPosGuess = curPos - (vel * timeToBoundary);
	}

private:
	static double CalcTimeToBoundaryAlongAxis(int axis,
											  const Vec2 &vel,
											  const Vec2 &curPos,
											  int numCells)
	{
		if (vel[axis] < 0) // Out of high bounds
		{
			return (curPos[axis] - numCells) / vel[axis];
		}
		else // Out of low bounds
		{
			return curPos[axis] / vel[axis];
		}
	}
};
