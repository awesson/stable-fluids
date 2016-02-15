TO COMPILE:
You will need to install g++ and libpng and make sure that it is in the include path (CPATH and LIBRARY_PATH)

HOW TO USE PROGRAM:
Add density and temperature with the right mouse button.
Change color of density by toggling the r, g, and b keys.
Add velocities with the left mouse button and dragging the mouse.
Toggle density/velocity display with the 'v' key.
Clear the simulation by pressing the 'c' key. Quit by pressing the 'q' key.

FILES:
ScalarField.h - Holds definition of scalar field class. Has functionality to
    advance the field by one time step, add one field to another, and
    interpolate values anywhere on the grid.
ScalarField.cpp - Implementation of the scalar field class.
VectorField.h - Holds definition of vector field class and a class with helper
    functions used by both the vector and scalar field classes. Has 
    functionality to advance the field by one time step, add one field to 
    another, and interpolate values anywhere on the grid.
VectorField.cpp - Implementation of the vector field class.
StableFluids.cpp - Defines the entry point for the console application.
    Initializes and frees memory for the fields, calls the functions
    which advance the simulation, and renders the density and velocity.
linearSolver.cpp/linearSolver.h - Has an implicit matrix class which defines the
    matrix by a function which applies a matrix-vector multiply. Also holds the
    implementation of a conjugate gradient solver.
imageio.h/imageio.cpp - Has functionality for loading images and saving frames
    as images.