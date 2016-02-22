# $Id: gfx-config.in 343 2008-09-13 18:34:59Z garland $

CXX = g++
CXXFLAGS = -g -O2 -Wall -Wno-sign-compare -Iinclude -DHAVE_CONFIG_H
OBJS = imageio.o StableFluids.o ScalarField.o VectorField.o LinearSolver.o

# Linux
EXE = StableFluids
LDFLAGS = -lpng -lGL -lGLU -lglut

# OS X
ifeq ($(shell uname), Darwin)
LDFLAGS = -lpng -framework OpenGL -framework GLUT
endif

$(EXE) : $(OBJS)
	$(CXX) -o $@ $^ $(CFLAGS) $(LDFLAGS)
clean:
	rm $(OBJS) $(EXE)
