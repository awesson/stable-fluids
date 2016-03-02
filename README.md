#Stable Fluids
Resurrected old project from college. 2D fluid simulation of incompressible fluids based on Jos Stam's Stable Fluids paper.

## TO COMPILE:
Useful reference for installing openGL/GLUT:
http://web.eecs.umich.edu/~sugih/courses/eecs487/glut-howto/

####linux/unix:
You will need to have installed...
- g++
- openGL
- GLUT
- libpng

Make sure that it is in the include path (CPATH and LIBRARY_PATH).
Then just run make in the top level directory of the repo.

####Windows:
There is a Visual Studio 2010 solution in StableFluidsVSSln.
You will need to install openGL and GLUT and put them in the correct folder
depending on your version of Visual Studio (see above for help).
Once that is done, you should be able to open the .sln file and
run in either Debug or Release.

## HOW TO USE PROGRAM:
- Add density and temperature with the right mouse button.
- Change color of density by toggling the r, g, and b keys.
- Add velocities with the left mouse button and dragging the mouse.
- Toggle density/velocity display with the 'v' key.
- Clear the simulation by pressing the 'c' key. Quit by pressing the 'q' key.

## FILES:
ScalarField.h - Holds definition of scalar field class. Has functionality to
                advance the field by one time step, add one field to another,
                and interpolate values anywhere on the grid.
    
ScalarField.cpp - Implementation of the scalar field class.

VectorField.h - Holds definition of vector field class and a class with helper
                functions used by both the vector and scalar field classes.
                Has functionality to advance the field by one time step,
                add one field to another, and interpolate values anywhere
                on the grid.
    
VectorField.cpp - Implementation of the vector field class.

StableFluids.cpp - Defines the entry point for the console application.
                   Initializes and frees memory for the fields,
                   calls the functions which advance the simulation,
                   handles input, and renders the density and velocity.
    
LinearSolver.cpp/LinearSolver.h - Has an implicit matrix class which defines
                                  the matrix by a function which applies a
                                  matrix-vector multiply. Also holds the
                                  implementation of conjugate gradient solver.
    
imageio.h/imageio.cpp - Has functionality for loading images
                        and saving frames as images. (Not used on Windows)
