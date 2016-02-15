#pragma once
#include <gfx/vec2.h>
#include "linearSolver.hpp"

#define EPSILON 1.0e-10

/******************************************************************************
 * MAC grid representation of a vector field. The vector components are defined
 * on the bottom and left edges of each cell. Has a boundary of one cell width
 * in order to have the velocity defined at the top and right edges and to
 * enforce the boundary conditions when interpolating.
 ******************************************************************************/
class VectorField
{
public:
    VectorField(int a_NumCells, double a_Viscosity, double a_Dt);
    VectorField(VectorField *CopyField);
    virtual ~VectorField(void);

    void    TimeStep(VectorField *a_SrcField, VectorField *VelocityField);
    void    AddField(VectorField *a_SrcField);
    void    AddField(Vec2 *a_SrcField);
    void    Advection(VectorField *oldField);
    void    VorticityConfinement();
    void    Projection();
    Vec2    Interpolate(Vec2 pos);
    double  BilinearInterpolation(Vec2 pos, int low_x, int low_y, int comp);

    Vec2    &operator[](int i)       { return m_Field[i]; }

    int     m_NumCells;
    Vec2    *m_Field;
    double  m_Viscosity;
    double  m_Dt;
};

/******************************************************************************
 * Class wrapper for a helper function used in both the vector and scalar field
 * classes. This makes sure multiple instances of this function are not created
 * each time it is included, such as in StableFluids.cpp as well as
 * ScalarField.cpp.
 ******************************************************************************/
class Field_Helper
{
public:
    /**************************************************************************
     * Helper function for the field's advection step that calculates
     * the position for which the new value will be interpolated at
     * by clipping the trajectory at the boundary if necessary.
     * The positions are relative to the bottom left corner of the cell (1,1)
     * 
     * @param prevPos_guess The position that is -vel*dt from the current
     *                      position.
     * @param curPos The current position being considered.
     * @param vel The velocity at the current position.
     * @param NumCells The width and length of the grid in cells.
     * @return The clipped position
     **************************************************************************/
    static Vec2 clipped_pos(Vec2 prevPos_guess, Vec2 curPos,
        Vec2 vel, int NumCells)
    {
        double time_to_boundary_in_x, time_to_boundary_in_y;
        
        if(prevPos_guess[0] < 0 || prevPos_guess[0] > NumCells
            || prevPos_guess[1] < 0 || prevPos_guess[1] > NumCells)
        {
            // check for components of velocity being zero since we divide by it
            if(fabs(vel[0]) < EPSILON){
                // velocity in x direction is 0 so only account for y velocity
                
                if(vel[1] < 0) // out of high bounds
                    time_to_boundary_in_y = (curPos[1] - NumCells) / vel[1];
                else // out of low bounds
                    time_to_boundary_in_y = curPos[1] / vel[1];
                    
                return curPos - vel * time_to_boundary_in_y;
            } else{
                if(fabs(vel[1]) < EPSILON){
                    // vel in the y direction is 0 so only account for x vel
                    
                    if(vel[0] < 0) // out of high bounds
                        time_to_boundary_in_x = (curPos[0] - NumCells) / vel[0];
                    else // out of low bounds
                        time_to_boundary_in_x = curPos[0] / vel[0];
                        
                    return curPos - vel * time_to_boundary_in_x;
                } else{
                    // neither x nor y velocity is zero
                    
                    if(vel[0] < 0) // out of high bounds
                        time_to_boundary_in_x = (curPos[0] - NumCells) / vel[0];
                    else // out of low bounds
                        time_to_boundary_in_x = curPos[0] / vel[0];
                        
                    if(vel[1] < 0) // out of high bounds
                        time_to_boundary_in_y = (curPos[1] - NumCells) / vel[1];
                    else // out of low bounds
                        time_to_boundary_in_y = curPos[1] / vel[1];
                    
                    // correct for the direction that goes out of bounds first
                    if(time_to_boundary_in_x < time_to_boundary_in_y)
                        return curPos - vel * time_to_boundary_in_x;
                    else
                        return curPos - vel * time_to_boundary_in_y;
                }
            }
        }
        
        // the guess position is inbounds so no need to clip
        return prevPos_guess;
    }
};
