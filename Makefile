# $Id: gfx-config.in 343 2008-09-13 18:34:59Z garland $

CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H
OBJS = imageio.o StableFluids.o ScalarField.o VectorField.o linearSolver.o

StableFluids: $(OBJS)
	$(CXX) -o $@ $^ -lpng -framework GLUT -framework OpenGL
clean:
	rm $(OBJS) StableFluids
