#pragma once

#include "BaseFluidField.hpp"

/******************************************************************************
 * MAC grid representation of a vector field. The vector components are defined
 * on the bottom and left edges of each cell. Has a boundary of one cell width
 * in order to have the velocity defined at the top and right edges and to
 * enforce the boundary conditions when interpolating.
 ******************************************************************************/
class VectorField : public FluidField<Vec2>
{
public:
	/***********************************************************************
	 * Standard constructor.
	 *
	 * @param numCells The width and height of the grid in cells.
	 * @param viscosity The viscosity of the fluid. (currently not used)
	 * @param dt The time step of the simulation.
	 ***********************************************************************/
	VectorField(int dimentionSize, double viscosity, double dt)
		: FluidField(dimentionSize, viscosity, dt)
	{
	}

	/******************************************************************************
	 * Advances the field by one timestep.
	 * Adds forces from the user input and applies the advection,
	 * vorticity confinement and projection steps.
	 *
	 * @param accelerationField A vector field representing forces from user input to be applied.
	 * @param VelocityField		The velocity field to advect through (usually itself).
	 ******************************************************************************/
	void TimeStep(const VectorField& accelerationField, const VectorField& velocityField);

	/******************************************************************************
	 * Interpolate the velocity at the given position using bilinear interpolation.
	 * Each component is interpolated independently as a weighted average of the
	 * velocities stored at the 4 nearest edge centers.
	 *
	 * @param pos The position to get the velocity at.
	 * 
	 * @return The velocity at the given position.
	 ******************************************************************************/
	virtual Vec2 Interpolate(const Vec2& pos) const override;

private:
	/*************************************************************************
	 * Adds another vector field, stored as an array, to this one.
	 * It assumes that the size of the array is m_TotalNumCells.
	 *
	 * @param sourceField The vector field to add to this one.
	 *************************************************************************/
	void AddVectorField(const Vec2 sourceField[]);
	
	/******************************************************************************
	 * Applies Lagrangian advection to the field. This backtracks from the current
	 * position using the velocity at the current point and then uses the value of
	 * the field at the previous point as the new value.
	 *
	 * @param oldField The vector field before the current time step.
	 *****************************************************************************/
	void Advection(const VectorField& oldField);

	/*****************************************************************************
	 * Applies vorticity confinement to the field. This finds a measure of the
	 * vorticity at each point and then reinforces this to counter act numerical
	 * diffusion.
	 *****************************************************************************/
	void VorticityConfinement();

	/*************************************************************************
	 * Applies the projection step to the field.
	 * The projection step solves for a field which has zero divergence.
	 * This is done by solving a linear system of
	 * Laplacian(pressure) = divergence(field) for the pressure. Then the zero
	 * divergence field is the current field - gradient(pressure).
	 *************************************************************************/
	void Projection();
	
	/*************************************************************************
	 * Interpolate the given component of velocity using bilinear interpolation.
	 *
	 * @param pos   The position to interpolate the value at.
	 * @param low_x The x index of the cell whose x component of velocity is the
	 *              closest on the left to the given position.
	 * @param low_y The y index of the cell whose y component of velocity is the
	 *              closest on the bottom to the given position.
	 * @param comp  Flag for which component is being interpolated.
	 *              0 for x and 1 for y.
	 *
	 * @return The interpolated value.
	 *************************************************************************/
	double BilinearInterpolation(const Vec2& pos, int lowX, int lowY, int comp) const;
};
