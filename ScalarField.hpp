#pragma once

#include "BaseFluidField.hpp"

class VectorField;

/*****************************************************************************
 * MAC grid representation of a scalar field. Values are defined at the center
 * of grid points. Has a boundary of one cell width in order to enforce the
 * boundary conditions when interpolating.
 *****************************************************************************/
class ScalarField : public FluidField<double>
{
public:
	/*************************************************************************
	 * Standard constructor.
	 *
	 * @param size      The width and height of the grid in cells.
	 * @param viscosity The viscosity of the fluid. (currently not used)
	 * @param dt        The time step of the simulation.
	 *************************************************************************/
	ScalarField(int size, double viscosity, double dt)
		: FluidField(size, viscosity, dt)
	{
	}

	/*************************************************************************
	 * Advances the field by one timestep. Adds any new scalar values from
	 * the user input and applies the advection step.
	 *
	 * @param srcField		A scalar field representing the change in
	 *						scalar values from user input to be added.
	 * @param velocityField The velocity field to advect the scalar through.
	 *************************************************************************/
	void TimeStep(const ScalarField& srcField, const VectorField& velocityField);

	/*************************************************************************
	 * Interpolates the scalar value at the given position on the grid using
	 * bilinear interpolation.
	 *
	 * @param pos The position to get the value at.
	 * @return The interpolated scalar value at the given position
	 *************************************************************************/
	virtual double Interpolate(const Vec2 &pos) const override;

private:
	/*************************************************************************
	 * Applies Lagrangian advection to the field.
	 * This backtracks from the current position using the velocity at
	 * the current point and then uses the value of the field at
	 * the previous point as the new value.
	 *
	 * @param vectorField The vector field to advect the scalar field through.
	 *************************************************************************/
	void Advection(const VectorField& vectorField);
};
