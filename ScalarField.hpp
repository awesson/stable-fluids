#pragma once
#include "VectorField.hpp"

/******************************************************************************
 * MAC grid representation of a scalar field. Values are defined at the center
 * of grid points. Has a boundary of one cell width in order to enforce the
 * boundary conditions when interpolating.
 ******************************************************************************/
class ScalarField
{
public:
    ScalarField(int a_NumCells, double a_Viscosity, double a_Dt);
    ScalarField(ScalarField *CopyField);
    virtual ~ScalarField(void);

    void    TimeStep(ScalarField *a_SrcField, VectorField *VelocityField);
    void    AddField(ScalarField *a_SrcField);
    void    Advection(VectorField *u);
    double  Interpolate(Vec2 pos);

    double  &operator[](int i) { return m_Field[i]; }

    int     m_NumCells;
    double  *m_Field;
    double  m_Viscosity;
    double  m_Dt;
};
