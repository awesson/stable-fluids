/*
-------------------------------------------------------------------------------
StableFluids.cpp : Defines the entry point for the console application.
-------------------------------------------------------------------------------
*/

#include "ScalarField.hpp"
#include "VectorField.hpp"
#ifndef _MSC_BUILD
#include "imageio.hpp"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <gfx/vec3.h>
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

/* Macros */

#define IX(i,j) ((i) + (GridSize * (j)))

// Weights for temperature and density
#define alpha .0005
#define beta .003

/* Global variables */

// Grid size in each dimension
static int GridSize;

// Total number of cells in the grid
static int TotalGridCells;

// Simulation time step
static double Dt;

// Diffusion rate of the fluid
static double Diff;

// Viscosity of the fluid
static double Visc;

// Amount that is added to velocity, density and temperature from user input
static double Force, Source, Temp;

// Flag for whether to draw velocity or density field
static bool DisplayVel;

// Flag for whether to dump the frames
static bool DumpFrames;
static int FrameNumber;

static VectorField *VelocityField, *AccelerationField;
static ScalarField *DensityFieldRed, *DeltaDensityFieldRed;
static ScalarField *DensityFieldGreen, *DeltaDensityFieldGreen;
static ScalarField *DensityFieldBlue, *DeltaDensityFieldBlue;
static ScalarField *TemperatureField, *DeltaTemperatureField;

static int WinId;
static int WinX, WinY;

static int PrevFpsTakenTime, FpsFrameNumber;
static char FpsString[11];

/* Variables for mouse interaction */
static int MouseDown[3];
static int ClickOriginMouseX, ClickOriginMouseY, MouseX, MouseY;

static bool PlaceRed, PlaceGreen, PlaceBlue;


static void InitValues()
{
	GridSize = 100;
	TotalGridCells = GridSize * GridSize;
	Dt = 2;
	Diff = 0.0;
	Visc = 0.0;
	Force = 1;
	Source = 10;
	Temp = 5;

	DisplayVel = false;
	DumpFrames = false;

	PrevFpsTakenTime = 0;
	FpsString[0] = '\0';

	PlaceRed = false;
	PlaceGreen = true;
	PlaceBlue = false;
}

/*
-------------------------------------------------------------------------------
free/clear/allocate simulation data
-------------------------------------------------------------------------------
*/
static void FreeData()
{
	if (VelocityField) delete (VelocityField);
	if (AccelerationField) delete (AccelerationField);
	if (DensityFieldRed) delete (DensityFieldRed);
	if (DeltaDensityFieldRed) delete (DeltaDensityFieldRed);
	if (DensityFieldGreen) delete (DensityFieldGreen);
	if (DeltaDensityFieldGreen) delete (DeltaDensityFieldGreen);
	if (DensityFieldBlue) delete (DensityFieldBlue);
	if (DeltaDensityFieldBlue) delete (DeltaDensityFieldBlue);
	if (TemperatureField) delete (TemperatureField);
	if (DeltaTemperatureField) delete (DeltaTemperatureField);
}

static void ClearData()
{
	for (int i = 0; i < GridSize; i++)
	{
		for (int j = 0; j < GridSize; ++j)
		{
			(*VelocityField)[IX(i,j)] = 0.0;
			(*AccelerationField)[IX(i,j)] = 0.0;
			(*DensityFieldRed)[IX(i,j)] = 0.0;
			(*DeltaDensityFieldRed)[IX(i,j)] = 0.0;
			(*DensityFieldGreen)[IX(i,j)] = 0.0;
			(*DeltaDensityFieldGreen)[IX(i,j)] = 0.0;
			(*DensityFieldBlue)[IX(i,j)] = 0.0;
			(*DeltaDensityFieldBlue)[IX(i,j)] = 0.0;
			(*DeltaTemperatureField)[IX(i,j)] = 0.0;
			(*TemperatureField)[IX(i,j)] = 0.0;
		}
	}
}

static int AllocateData()
{
	VelocityField = new VectorField(GridSize, Visc, Dt);
	AccelerationField = new VectorField(GridSize, Visc, Dt);
	DensityFieldRed = new ScalarField(GridSize, Diff, Dt);
	DeltaDensityFieldRed = new ScalarField(GridSize, Diff, Dt);
	DensityFieldGreen = new ScalarField(GridSize, Diff, Dt);
	DeltaDensityFieldGreen = new ScalarField(GridSize, Diff, Dt);
	DensityFieldBlue = new ScalarField(GridSize, Diff, Dt);
	DeltaDensityFieldBlue = new ScalarField(GridSize, Diff, Dt);
	TemperatureField = new ScalarField(GridSize, Diff, Dt);
	DeltaTemperatureField = new ScalarField(GridSize, Diff, Dt);

	if (!VelocityField || !AccelerationField
		|| !DensityFieldRed || !DeltaDensityFieldRed
		|| !DensityFieldGreen || !DeltaDensityFieldGreen
		|| !DensityFieldBlue || !DeltaDensityFieldBlue
		|| !TemperatureField || !DeltaTemperatureField)
	{
		fprintf(stderr, "cannot allocate data\n");
		return 0;
	}

	return 1;
}


/*
-------------------------------------------------------------------------------
OpenGL specific drawing routines
-------------------------------------------------------------------------------
*/

void SetOrthographicProjection()
{
	// Switch to projection mode
	glMatrixMode(GL_PROJECTION);

	// Save previous matrix which contains the
	// settings for the perspective projection
	glPushMatrix();

	// Reset matrix
	glLoadIdentity();

	// Set a 2D orthographic projection
	gluOrtho2D(0, WinX, WinY, 0);

	// Switch back to modelview mode
	glMatrixMode(GL_MODELVIEW);
}

void RestorePerspectiveProjection()
{
	glMatrixMode(GL_PROJECTION);
	// Restore previous projection matrix
	glPopMatrix();

	// Get back to modelview mode
	glMatrixMode(GL_MODELVIEW);
}

void RenderBitmapString(float x, float y, void *font, char *string)
{
	glRasterPos3f(x, y, 0.0f);
	for (char *c = string; *c != '\0'; c++)
	{
		glutBitmapCharacter(font, *c);
	}
}

static void RenderQuad(float x, float y, float width, float height)
{
	glPushMatrix();
	glTranslatef(x, y, 0.f);

	glBegin(GL_QUADS);

	glVertex3f(-0.5f * width, -0.5f * height, 0.0f);
	glVertex3f(-0.5f * width, 0.5f * height, 0.0f);
	glVertex3f(0.5f * width, 0.5f * height, 0.0f);
	glVertex3f(0.5f * width, -0.5f * height, 0.0f);

	glEnd();

	glPopMatrix();
}

static void DrawFps()
{
	glPushMatrix();
	glLoadIdentity();
	SetOrthographicProjection();

	glColor3f(0.0f, 0.0f, 0.0f);
	RenderQuad(25, 6, 50, 12);

	glColor3f(1.0f, 1.0f, 1.0f);
	RenderBitmapString(5, 10, GLUT_BITMAP_TIMES_ROMAN_10, FpsString);

	glPopMatrix();
	RestorePerspectiveProjection();
}

static void PreDisplay()
{
	glViewport(0, 0, WinX, WinY);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 1.0, 0.0, 1.0);
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
}

static void PostDisplay()
{
#ifndef _MSC_BUILD
	// Write frames if necessary.
	if (DumpFrames)
	{
		const int FRAME_INTERVAL = 1;
		if ((FrameNumber % FRAME_INTERVAL) == 0)
		{
			const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
			const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
			
			unsigned char * buffer = new char[w * h * 4];
			glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);

			char filename[13];
			sprintf(filename, "img%.5i.png", FrameNumber / FRAME_INTERVAL);
			printf("Dumped %s.\n", filename);

			saveImageRGBA(filename, buffer, w, h);

			delete[] buffer;
		}
		
		FrameNumber++;
	}
#endif

	DrawFps();

	glutSwapBuffers();
}

/**
  * Draws vectors scaled by the magnitude of the velocity,
  * representing the velocity field.
  **/
static void DrawVelocity()
{
	int i, j;
	float x, y, h;

	h = 1.0f / (GridSize - 2);

	glColor3f(1.0f, 1.0f, 1.0f);
	glLineWidth(1.0f);

	glBegin(GL_LINES);

	for (i = 1; i < GridSize - 1; i++)
	{
		x = (i - 0.5f) * h;

		for (j = 1; j < GridSize - 1; j++)
		{
			y = (j - 0.5f) * h;

			glVertex2f(x, y);
			glVertex2f(x + (*VelocityField)[IX(i, j)][0],
					   y + (*VelocityField)[IX(i, j)][1]);
		}
	}

	glEnd();
}

/**
  * Draws the color of each density field scaled
  * by the magnitude of the density at that point.
  **/
static void DrawDensity()
{
	int i, j;
	float x, y, h;

	h = 1.0f / (GridSize - 2);

	glBegin(GL_QUADS);

	for (i = 1; i < GridSize - 1; i++)
	{
		x = (i - 0.5f) * h;

		for (j = 1; j < GridSize - 1; j++)
		{
			y = (j - 0.5f) * h;

			Vec3f bl = Vec3(0, 0, 0); // bottom left color
			Vec3f br = Vec3(0, 0, 0); // bottom right color
			Vec3f tr = Vec3(0, 0, 0); // top right color
			Vec3f tl = Vec3(0, 0, 0); // top left color

			bl[0] += (*DensityFieldRed)[IX(i, j)];
			br[0] += (*DensityFieldRed)[IX(i+1, j)];
			tl[0] += (*DensityFieldRed)[IX(i, j+1)];
			tr[0] += (*DensityFieldRed)[IX(i+1, j+1)];

			bl[1] += (*DensityFieldGreen)[IX(i, j)];
			br[1] += (*DensityFieldGreen)[IX(i+1, j)];
			tl[1] += (*DensityFieldGreen)[IX(i, j+1)];
			tr[1] += (*DensityFieldGreen)[IX(i+1, j+1)];

			bl[2] += (*DensityFieldBlue)[IX(i, j)];
			br[2] += (*DensityFieldBlue)[IX(i+1, j)];
			tl[2] += (*DensityFieldBlue)[IX(i, j+1)];
			tr[2] += (*DensityFieldBlue)[IX(i+1, j+1)];

			glColor3f(bl[0], bl[1], bl[2]);
			glVertex2f(x, y);
			glColor3f(br[0], br[1], br[2]);
			glVertex2f(x + h, y);
			glColor3f(tr[0], tr[1], tr[2]);
			glVertex2f(x + h, y + h);
			glColor3f(tl[0], tl[1], tl[2]);
			glVertex2f(x, y + h);
		}
	}

	glEnd();
}

/**
 * Relates mouse movements to forces/sources.
 * Adds a force to the velocity field for single click dragging.
 * Adds density and temperature when double clicking.
 * The color will be a combination of the bool values set by the r, g, and b keys.
 *
 * @param[out] d    The red fluid density
 * @param[out] d2   The green fluid density
 * @param[out] d3   The blue fluid density
 * @param[out] t    The temperate field
 * @param[out] vel  The velocity field representing the effects of the user forces
 **/
static void GetFromUI(ScalarField *d,
					  ScalarField *d2,
					  ScalarField *d3,
					  ScalarField *temp,
					  VectorField *vel)
{
	int i, j, size = TotalGridCells;

	// Initialize fields
	for (i = 0; i < size; i++)
	{
		(*vel)[i] = (*d)[i] = (*d2)[i] = (*d3)[i] = (*temp)[i] = 0.0;
	}

	// If there is no mouse input, then everything stays empty
	if (!MouseDown[0] && !MouseDown[2]) return;

	// Convert the mouse position to a grid position
	i = static_cast<int>((MouseX / static_cast<float>(WinX)) * GridSize);
	j = static_cast<int>(((WinY - MouseY) / static_cast<float>(WinY)) * GridSize);

	// Make sure the position is within the boundaries
	if ((i < 1) || (i > GridSize - 2) || (j < 1) || (j > GridSize - 2)) return;

	// Force on velocity field
	if (MouseDown[0])
	{
		(*vel)[IX(i, j)][0] = Force * (MouseX - ClickOriginMouseX);
		(*vel)[IX(i, j)][1] = Force * (ClickOriginMouseY - MouseY);
	}

	// Add fluid
	if (MouseDown[2])
	{
		if (PlaceRed)
		{
			(*d)[IX(i,j)] = Source;
		}

		if (PlaceGreen)
		{
			(*d2)[IX(i,j)] = Source;
		}

		if (PlaceBlue)
		{
			(*d3)[IX(i,j)] = Source;
		}

		(*temp)[IX(i,j)] = Temp;
	}

	ClickOriginMouseX = MouseX;
	ClickOriginMouseY = MouseY;

	return;
}

/*
-------------------------------------------------------------------------------
GLUT callback routines
-------------------------------------------------------------------------------
*/

static void KeyFunc(unsigned char key, int x, int y)
{
	switch (key)
	{
		case 'c':
		case 'C':
			ClearData();
			break;

		case 'd':
		case 'D':
			FrameNumber = 0;
			DumpFrames = !DumpFrames;
			break;

		case 'q':
		case 'Q':
			FreeData();
			exit(0);
			break;

		case 'v':
		case 'V':
			DisplayVel = !DisplayVel;
			break;

		case 'r':
		case 'R':
			PlaceRed = !PlaceRed;
			break;
		case 'g':
		case 'G':
			PlaceGreen = !PlaceGreen;
			break;
		case 'b':
		case 'B':
			PlaceBlue = !PlaceBlue;
			break;
	}
}

static void MouseFunc(int button, int state, int x, int y)
{
	ClickOriginMouseX = MouseX = x;
	ClickOriginMouseY = MouseY = y;

	if (glutGetModifiers() & GLUT_ACTIVE_SHIFT)
	{
		button = 2;
	}

	MouseDown[button] = state == GLUT_DOWN;
}

static void MotionFunc(int x, int y)
{
	MouseX = x;
	MouseY = y;
}

static void ReshapeFunc(int width, int height)
{
	glutSetWindow(WinId);
	glutReshapeWindow(width, height);

	WinX = width;
	WinY = height;
}

static void IdleFunc()
{
	// Get user input
	GetFromUI(DeltaDensityFieldRed,
			  DeltaDensityFieldGreen,
			  DeltaDensityFieldBlue,
			  DeltaTemperatureField,
			  AccelerationField);

	// Add gravity and temperature forces to velocity field.
	// Assumes the ambient temperature is zero.
	for (int i = 0; i < TotalGridCells; ++i)
	{
		double tot_den = (*DensityFieldRed)[i]
					   + (*DensityFieldGreen)[i]
					   + (*DensityFieldBlue)[i];
		double T = (*TemperatureField)[i];
		(*AccelerationField)[i] += Vec2(0, -(alpha * tot_den) + (beta * T));
	}

	// Advance sim
	VelocityField->TimeStep(*AccelerationField, *VelocityField);
	DensityFieldRed->TimeStep(*DeltaDensityFieldRed, *VelocityField);
	DensityFieldGreen->TimeStep(*DeltaDensityFieldGreen, *VelocityField);
	DensityFieldBlue->TimeStep(*DeltaDensityFieldBlue, *VelocityField);
	TemperatureField->TimeStep(*DeltaTemperatureField, *VelocityField);

	glutSetWindow(WinId);
	glutPostRedisplay();

	// Calculate fps
	int cur_time = glutGet(GLUT_ELAPSED_TIME);
	if (cur_time - PrevFpsTakenTime > 1000)
	{
		float fps = (1000.0f * FpsFrameNumber) / static_cast<float>(cur_time - PrevFpsTakenTime);
		sprintf(FpsString, "FPS: %4.2f", fps);
		PrevFpsTakenTime = cur_time;
		FpsFrameNumber = 0;
	}
	FpsFrameNumber++;
}

static void DisplayFunc()
{
	PreDisplay();

	if (DisplayVel)
	{
		DrawVelocity();
	}
	else
	{
		DrawDensity();
	}

	PostDisplay();
}


/**
 * Open a glut compatible window and set callbacks
 **/
static void OpenGlutWindow()
{
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);

	glutInitWindowPosition(0, 0);
	glutInitWindowSize(WinX, WinY);
	WinId = glutCreateWindow("Stable Fluids!");

	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();
	glClear(GL_COLOR_BUFFER_BIT);
	glutSwapBuffers();

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);

	PreDisplay();

	glutKeyboardFunc(KeyFunc);
	glutMouseFunc(MouseFunc);
	glutMotionFunc(MotionFunc);
	glutReshapeFunc(ReshapeFunc);
	glutIdleFunc(IdleFunc);
	glutDisplayFunc(DisplayFunc);
}


int main(int argc, char **argv)
{
	glutInit(&argc, argv);

	if (argc != 1 && argc != 8)
	{
		fprintf(stderr, "usage : %s N Dt Diff Visc Force Source\n", argv[0]);
		fprintf(stderr, "where:\n");
		fprintf(stderr, "\t GridSize : Grid resolution\n");
		fprintf(stderr, "\t Dt       : Time step\n");
		fprintf(stderr, "\t Diff     : Diffusion rate of the density\n");
		fprintf(stderr, "\t Visc     : Viscosity of the fluid\n");
		fprintf(stderr, "\t Force    : Scales the mouse movement that generate a force\n");
		fprintf(stderr, "\t Source   : Amount of density that will be deposited\n");
		fprintf(stderr, "\t Temp     : The temperature of the fluid\n");
		exit(1);
	}

	InitValues();

	if (argc == 1)
	{
		fprintf(stderr, "Using defaults : GridSize=%d Dt=%g Diff=%g Visc=%g Force=%g Source=%g Temp=%g\n",
					  GridSize, Dt, Diff, Visc, Force, Source, Temp);
	}
	else
	{
		GridSize = atoi(argv[1]);
		TotalGridCells = GridSize * GridSize;
		Dt = atof(argv[2]);
		Diff = atof(argv[3]);
		Visc = atof(argv[4]);
		Force = atof(argv[5]);
		Source = atof(argv[6]);
		Temp = atof(argv[7]);
	}

	printf("\n\nHow to use this demo:\n\n");
	printf("\t Add density and temperature with the right mouse button\n");
	printf("\t Change color of density being by toggling the r, g, and b keys\n");
	printf("\t Add velocities with the left mouse button and dragging the mouse\n");
	printf("\t Toggle density/velocity display with the 'v' key\n");
	printf("\t Clear the simulation by pressing the 'c' key\n");
	printf("\t Quit by pressing the 'q' key\n");
	fflush(stdout);

	if (!AllocateData())
	{
		exit(1);
	}

	ClearData();

	WinX = 512;
	WinY = 512;
	OpenGlutWindow();

	glutMainLoop();

	exit(0);
}
