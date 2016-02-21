/*
-------------------------------------------------------------------------------
StableFluids.cpp : Defines the entry point for the console application.
-------------------------------------------------------------------------------
*/

#include "ScalarField.hpp"
#include "VectorField.hpp"
#include "imageio.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <gfx/vec3.h>
#include <GLUT/glut.h>

/* macros */

#define IX(i,j) ((i)+(N+2)*(j))

// weights for temperature and density
#define alpha .0005
#define beta .003

/* global variables */

static int N; // grid size
static double dt, diff, visc;

// amount that is added to velocity, density and temperature from user input
static float force, source, temp;

static int dvel; // bool for whether to draw velocity or density field
static int dump_frames; // flag for whether to dump the frames
static int frame_number;

static VectorField *VelocityField, *PrevVelocityField;
static ScalarField *DensityFieldRed, *PrevDensityFieldRed;
static ScalarField *DensityFieldGreen, *PrevDensityFieldGreen;
static ScalarField *DensityFieldBlue, *PrevDensityFieldBlue;
static ScalarField *TemperatureField, *PrevTemperatureField;

static int win_id;
static int win_x, win_y;

static int prev_fps_taken_time, fpsFrameNumber;
static char fpsString[10];

/* variables for mouse interaction */
static int mouse_down[3];
static int omx, omy, mx, my;

static bool red, green, blue;


/*
-------------------------------------------------------------------------------
free/clear/allocate simulation data
-------------------------------------------------------------------------------
*/


static void free_data ( void )
{
    if( VelocityField ) delete ( VelocityField );
    if( PrevVelocityField ) delete ( PrevVelocityField );
    if( DensityFieldRed ) delete ( DensityFieldRed );
    if( PrevDensityFieldRed ) delete ( PrevDensityFieldRed );
    if( DensityFieldGreen ) delete ( DensityFieldGreen );
    if( PrevDensityFieldGreen ) delete ( PrevDensityFieldGreen );
    if( DensityFieldBlue ) delete ( DensityFieldBlue );
    if( PrevDensityFieldBlue ) delete ( PrevDensityFieldBlue );
    if( TemperatureField ) delete ( TemperatureField );
    if( PrevTemperatureField ) delete ( PrevTemperatureField );
}

/**
  * Reset all values to zero.
  **/
static void clear_data ( void )
{
    for (int i=0 ; i< N+2 ; i++ ) {
        for(int j=0; j < N+2; ++j) {
            (*VelocityField)[IX(i,j)][0] = (*VelocityField)[IX(i,j)][1] = 0.0;
            (*PrevVelocityField)[IX(i,j)][0] = 0.0;
            (*PrevVelocityField)[IX(i,j)][1] = 0.0;
            (*DensityFieldRed)[IX(i,j)] = 0.0;
            (*PrevDensityFieldRed)[IX(i,j)] = 0.0;
            (*DensityFieldGreen)[IX(i,j)] = 0.0;
            (*PrevDensityFieldGreen)[IX(i,j)] = 0.0;
            (*DensityFieldBlue)[IX(i,j)] = 0.0;
            (*PrevDensityFieldBlue)[IX(i,j)] = 0.0;
            (*PrevTemperatureField)[IX(i,j)] = 0.0;
            (*TemperatureField)[IX(i,j)] = 0.0;
        }
    }
}

static int allocate_data ( void )
{
    VelocityField         = new VectorField(N, visc, dt);
    PrevVelocityField     = new VectorField(N, visc, dt);
    DensityFieldRed       = new ScalarField(N, diff, dt);
    PrevDensityFieldRed   = new ScalarField(N, diff, dt);
    DensityFieldGreen     = new ScalarField(N, diff, dt);
    PrevDensityFieldGreen = new ScalarField(N, diff, dt);
    DensityFieldBlue      = new ScalarField(N, diff, dt);
    PrevDensityFieldBlue  = new ScalarField(N, diff, dt);
    TemperatureField      = new ScalarField(N, diff, dt);
    PrevTemperatureField  = new ScalarField(N, diff, dt);

    if ( !VelocityField || !PrevVelocityField ||
         !DensityFieldRed || !PrevDensityFieldRed ||
         !DensityFieldGreen || !PrevDensityFieldGreen ||
         !DensityFieldBlue || !PrevDensityFieldBlue ||
         !TemperatureField || !PrevTemperatureField) {
        fprintf ( stderr, "cannot allocate data\n" );
        return ( 0 );
    }
    
    prev_fps_taken_time = 0;
    fpsString[0] = '\0';

    // bool for which color fluid will be added
    red = false;
    green = true;
    blue = false;

    return ( 1 );
}


/*
-------------------------------------------------------------------------------
OpenGL specific drawing routines
-------------------------------------------------------------------------------
*/

void setOrthographicProjection()
{
	// switch to projection mode
	glMatrixMode(GL_PROJECTION);

	// save previous matrix which contains the
	//settings for the perspective projection
	glPushMatrix();

	// reset matrix
	glLoadIdentity();

	// set a 2D orthographic projection
	gluOrtho2D(0, win_x, win_y, 0);

	// switch back to modelview mode
	glMatrixMode(GL_MODELVIEW);
}

void restorePerspectiveProjection()
{
	glMatrixMode(GL_PROJECTION);
	// restore previous projection matrix
	glPopMatrix();

	// get back to modelview mode
	glMatrixMode(GL_MODELVIEW);
}

void renderBitmapString(float x,
		                float y,
		                void *font,
		                char *string)
{
	glRasterPos3f(x, y, 0.0f);
	for (char *c = string; *c != '\0'; c++)
	{
		glutBitmapCharacter(font, *c);
	}
}

static void renderQuad(float x, float y, float width, float height)
{
    glPushMatrix();
    glTranslatef(x, y, 0.f);
    
    glBegin(GL_QUADS); // Start drawing a quad primitive  
  
    glVertex3f(-0.5f * width, -0.5f * height, 0.0f);
    glVertex3f(-0.5f * width, 0.5f * height, 0.0f);
    glVertex3f(0.5f * width, 0.5f * height, 0.0f);
    glVertex3f(0.5f * width, -0.5f * height, 0.0f);
  
    glEnd();  
    
    glPopMatrix();
}

static void draw_fps()
{
    glPushMatrix();
    glLoadIdentity();
    setOrthographicProjection();
    glColor3f(0.0f, 0.0f, 0.0f);
    renderQuad(25, 6, 50, 12);
    glColor3f(1.0f, 1.0f, 1.0f);
    renderBitmapString(5, 10,  GLUT_BITMAP_TIMES_ROMAN_10, fpsString);
    glPopMatrix();
    restorePerspectiveProjection();
}

static void pre_display ( void )
{
    glViewport ( 0, 0, win_x, win_y );
    glMatrixMode ( GL_PROJECTION );
    glLoadIdentity ();
    gluOrtho2D ( 0.0, 1.0, 0.0, 1.0 );
    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
}

static void post_display ( void )
{
    // Write frames if necessary.
    if (dump_frames) {
        const int FRAME_INTERVAL = 1;
        if ((frame_number % FRAME_INTERVAL) == 0) {
            const unsigned int w = glutGet(GLUT_WINDOW_WIDTH);
            const unsigned int h = glutGet(GLUT_WINDOW_HEIGHT);
            unsigned char * buffer = (unsigned char *) malloc(w * h * 4 * sizeof(unsigned char));
            if (!buffer)
                exit(-1);
            glReadPixels(0, 0, w, h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
            char filename[13];
            sprintf(filename, "img%.5i.png", frame_number / FRAME_INTERVAL);
            printf("Dumped %s.\n", filename);
            saveImageRGBA(filename, buffer, w, h);
            free(buffer);
        }
        frame_number++;
    }
    
    draw_fps ();

    glutSwapBuffers ();
}

/**
  * Draws vectors scaled by the magnitude of the velocity,
  * representing the velocity field.
  **/
static void draw_velocity ( void )
{
    int i, j;
    float x, y, h;

    h = 1.0f/N;

    glColor3f ( 1.0f, 1.0f, 1.0f );
    glLineWidth ( 1.0f );

    glBegin ( GL_LINES );

    for ( i=1 ; i<=N ; i++ ) {
        x = (i-0.5f)*h;
        for ( j=1 ; j<=N ; j++ ) {
            y = (j-0.5f)*h;

            glVertex2f ( x, y );
            glVertex2f ( x+(*VelocityField)[IX(i,j)][0],
                         y+(*VelocityField)[IX(i,j)][1] );
        }
    }

    glEnd ();
}

/**
  * Draws the color of each density field scaled
  * by the magnitude of the density at that point.
  **/
static void draw_density ( void )
{
    int i, j;
    float x, y, h;

    h = 1.0f/N;

    glBegin ( GL_QUADS );

    for ( i=1; i<=N ; i++ ) {
        x = (i-0.5f)*h;
        for ( j=1; j<=N ; j++ ) {
            y = (j-0.5f)*h;
            Vec3f bl = Vec3(0,0,0); // bottom left color
            Vec3f br = Vec3(0,0,0); // bottom right color
            Vec3f tr = Vec3(0,0,0); // top right color
            Vec3f tl = Vec3(0,0,0); // top left color

            bl[0] += (*DensityFieldRed)[IX(i,j)];
            br[0] += (*DensityFieldRed)[IX(i+1,j)];
            tl[0] += (*DensityFieldRed)[IX(i,j+1)];
            tr[0] += (*DensityFieldRed)[IX(i+1,j+1)];

            bl[1] += (*DensityFieldGreen)[IX(i,j)];
            br[1] += (*DensityFieldGreen)[IX(i+1,j)];
            tl[1] += (*DensityFieldGreen)[IX(i,j+1)];
            tr[1] += (*DensityFieldGreen)[IX(i+1,j+1)];

            bl[2] += (*DensityFieldBlue)[IX(i,j)];
            br[2] += (*DensityFieldBlue)[IX(i+1,j)];
            tl[2] += (*DensityFieldBlue)[IX(i,j+1)];
            tr[2] += (*DensityFieldBlue)[IX(i+1,j+1)];

            glColor3f ( bl[0], bl[1], bl[2] ); glVertex2f ( x, y );
            glColor3f ( br[0], br[1], br[2] ); glVertex2f ( x+h, y );
            glColor3f ( tr[0], tr[1], tr[2] ); glVertex2f ( x+h, y+h );
            glColor3f ( tl[0], tl[1], tl[2] ); glVertex2f ( x, y+h );
        }
    }

    glEnd ();
}

/*
-------------------------------------------------------------------------------
relates mouse movements to forces/sources
-------------------------------------------------------------------------------
*/

/**
 * Adds a force to the velocity field for single click dragging. Adds density
 * and temperature when double clicking. The color will be a combination of the
 * bool values set by the r, g, and b keys.
 *
 * @param d The red fluid density
 * @param d2 The green fluid density
 * @param d3 The blue fluid density
 * @param T The temperate field
 * @param u_v The velocity field representing the effects of the user forces
 **/
static void get_from_UI( ScalarField * d,
                         ScalarField * d2,
                         ScalarField * d3,
                         ScalarField * T,
                         VectorField * u_v )
{
    int i, j, size = (N+2)*(N+2);

    // initialize fields
    for ( i=0 ; i<size ; i++ ) {
        (*u_v)[i][0] = (*u_v)[i][1] = (*d)[i] = 0.0;
        (*d2)[i] = (*d3)[i] = (*T)[i] = 0.0;
    }

    if ( !mouse_down[0] && !mouse_down[2] ) return;

    i = (int)((       mx /(float)win_x)*N+1);
    j = (int)(((win_y-my)/(float)win_y)*N+1);

    if ( i<1 || i>N || j<1 || j>N ) return;

    if ( mouse_down[0] ) { // force on velocity field
            (*u_v)[IX(i,j)][0] = force * (mx-omx);
            (*u_v)[IX(i,j)][1] = force * (omy-my);
    }

    if ( mouse_down[2] ) { // add fluid
        if(red)
            (*d)[IX(i,j)] = source;
        if(green)
            (*d2)[IX(i,j)] = source;
        if(blue)
            (*d3)[IX(i,j)] = source;
        (*T)[IX(i,j)] = temp;
    }

    omx = mx;
    omy = my;

    return;
}

/*
-------------------------------------------------------------------------------
GLUT callback routines
-------------------------------------------------------------------------------
*/

static void key_func ( unsigned char key, int x, int y )
{
    switch ( key )
    {
    case 'c':
    case 'C':
        clear_data ();
        break;

    case 'd':
    case 'D':
        frame_number = 0;
        dump_frames = !dump_frames;
        break;

    case 'q':
    case 'Q':
        free_data ();
        exit ( 0 );
        break;

    case 'v':
    case 'V':
        dvel = !dvel;
        break;

    case 'r':
    case 'R':
        red = !red;
        break;
    case 'g':
    case 'G':
        green = !green;
        break;
    case 'b':
    case 'B':
        blue = !blue;
        break;
    }
}

static void mouse_func ( int button, int state, int x, int y )
{
    omx = mx = x;
    omx = my = y;

    if (glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
        button = 2;
    }

    mouse_down[button] = state == GLUT_DOWN;
}

static void motion_func ( int x, int y )
{
    mx = x;
    my = y;
}

static void reshape_func ( int width, int height )
{
    glutSetWindow ( win_id );
    glutReshapeWindow ( width, height );

    win_x = width;
    win_y = height;
}

static void idle_func ( void )
{   
    // get user input
    get_from_UI( PrevDensityFieldRed, PrevDensityFieldGreen,
                 PrevDensityFieldBlue, PrevTemperatureField,
                 PrevVelocityField );

    /* add gravity and temperature forces to velocity field.
       Assumes the ambient temperature is just zero */
    for(int i = 0; i < (N+2)*(N+2); ++i){
        double tot_den = DensityFieldRed->m_Field[i]
                       + DensityFieldGreen->m_Field[i]
                       + DensityFieldBlue->m_Field[i];
        double T = TemperatureField->m_Field[i];
        PrevVelocityField->m_Field[i] += Vec2(0, -alpha*tot_den + beta*T);
    }

    // advance sim
    VelocityField->TimeStep( PrevVelocityField, VelocityField );
    DensityFieldRed->TimeStep( PrevDensityFieldRed, VelocityField );
    DensityFieldGreen->TimeStep( PrevDensityFieldGreen, VelocityField );
    DensityFieldBlue->TimeStep( PrevDensityFieldBlue, VelocityField );
    TemperatureField->TimeStep(PrevTemperatureField, VelocityField );

    glutSetWindow ( win_id );
    glutPostRedisplay ();
    
    // calculate fps
	int cur_time = glutGet(GLUT_ELAPSED_TIME);
	if(cur_time - prev_fps_taken_time > 1000)
	{
		float fps = (1000.0f * fpsFrameNumber) / ((float)(cur_time - prev_fps_taken_time));
		sprintf(fpsString, "FPS: %4.2f", fps);
		prev_fps_taken_time = cur_time;
		fpsFrameNumber = 0;
	}
	fpsFrameNumber++;
}

static void display_func()
{
    pre_display ();

    if ( dvel ) draw_velocity ();
    else    draw_density ();

    post_display ();
}


/**
 * Open a glut compatible window and set callbacks
 **/
static void open_glut_window ( void )
{
    glutInitDisplayMode ( GLUT_RGBA | GLUT_DOUBLE );

    glutInitWindowPosition ( 0, 0 );
    glutInitWindowSize ( win_x, win_y );
    win_id = glutCreateWindow ( "Stable Fluids!" );

    glClearColor ( 0.0f, 0.0f, 0.0f, 1.0f );
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();
    glClear ( GL_COLOR_BUFFER_BIT );
    glutSwapBuffers ();

    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);

    pre_display ();

    glutKeyboardFunc ( key_func );
    glutMouseFunc ( mouse_func );
    glutMotionFunc ( motion_func );
    glutReshapeFunc ( reshape_func );
    glutIdleFunc ( idle_func );
    glutDisplayFunc ( display_func );
}


int main ( int argc, char ** argv )
{
    glutInit ( &argc, argv );

    if ( argc != 1 && argc != 8 ) {
        fprintf ( stderr, "usage : %s N dt diff visc force source\n", argv[0] );
        fprintf ( stderr, "where:\n" );\
        fprintf ( stderr, "\t N      : grid resolution\n" );
        fprintf ( stderr, "\t dt     : time step\n" );
        fprintf ( stderr, "\t diff   : diffusion rate of the density\n" );
        fprintf ( stderr, "\t visc   : viscosity of the fluid\n" );
        fprintf ( stderr, "\t force  : scales the mouse movement that generate a force\n" );
        fprintf ( stderr, "\t source : amount of density that will be deposited\n" );
        fprintf ( stderr, "\t temp : the temperature of the fluid\n" );
        exit ( 1 );
    }

    if ( argc == 1 ) {
        N = 89;
        dt = 2;
        diff = 0.0;
        visc = 0.0;
        force = 1;
        source = 10;
        temp = 5;
        fprintf (stderr, "Using defaults : N=%d dt=%g diff=%g visc=%g force=%g source=%g temp=%g\n",
                 N, dt, diff, visc, force, source, temp);
    } else {
        N = atoi(argv[1]);
        dt = atof(argv[2]);
        diff = atof(argv[3]);
        visc = atof(argv[4]);
        force = atof(argv[5]);
        source = atof(argv[6]);
        temp = atof(argv[7]);
    }

    printf ( "\n\nHow to use this demo:\n\n" );
    printf ( "\t Add density and temperature with the right mouse button\n" );
    printf ( "\t Change color of density being by toggling the r, g, and b keys\n");
    printf ( "\t Add velocities with the left mouse button and dragging the mouse\n" );
    printf ( "\t Toggle density/velocity display with the 'v' key\n" );
    printf ( "\t Clear the simulation by pressing the 'c' key\n" );
    printf ( "\t Quit by pressing the 'q' key\n" );
    fflush(stdout);

    dvel = 0;

    if ( !allocate_data () ) exit ( 1 );
    clear_data ();

    win_x = 512;
    win_y = 512;
    open_glut_window ();

    glutMainLoop ();

    exit ( 0 );
}


