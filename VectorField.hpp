#pragma once

#include <gfx/vec2.h>

/******************************************************************************
 * MAC grid representation of a vector field. The vector components are defined
 * on the bottom and left edges of each cell. Has a boundary of one cell width
 * in order to have the velocity defined at the top and right edges and to
 * enforce the boundary conditions when interpolating.
 ******************************************************************************/
class VectorField
{
public:
	VectorField(int dimentionSize, double viscosity, double dt);
	VectorField(VectorField *copyField);
	virtual ~VectorField(void);

	void TimeStep(VectorField *sourceField, VectorField *velocityField);
	void AddField(VectorField *a_SrcField);
	void AddField(Vec2 *srouceField);
	void Advection(VectorField *oldField);
	void VorticityConfinement();
	void Projection();
	Vec2 Interpolate(Vec2 pos);
	double BilinearInterpolation(Vec2 pos, int lowX, int lowY, int comp);

	Vec2& operator[](int i)
	{
		return m_Field[i];
	}

	int m_NumCells;
	int m_DimentionSize;
	Vec2 *m_Field;
	double m_Viscosity;
	double m_Dt;
};
