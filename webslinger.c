//-------------------------------------------------------------------------
// This program checks the rigidity of various structures by building them
// as mass-spring systems and seeing if they stay rigid under various
// forces and gravity.

#ifdef _WIN32
#include <windows.h>
#endif

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/gl.h>
#include <glut.h>
#include "massmesh.h"
#include "graphics.h"

// File name to store CSV output points into.
char      *g_csv_outfile_name = NULL;
unsigned  g_csv_frame_number = 0;

const double	GACCEL	=      0;	//< Acceleration due to gravity
const double	TIME	=	(1.0/5000.0);
double	medium_dx = 0.0, medium_dy = 0, medium_dz = 0;  //< How fast the medium around the masses is moving

/* The maximum velocity of the wind, how fast the wind changes, and the
   fraction of the screen in which the string can move without changing
   the wind. */

#define	MAXWIND 10.0
#define	WINDCHANGE 0.2
#define	SCREENPART (1.0/3.0)

// Parameters that affect the graphics
const double	SCALE = 0.004;
const int	SIMS_PER_DRAW = 800;
const int	DRAWS_PER_FRAME = 1;

/* =====================================================================
 * Section of code dealing with making and simulating a 2D fibrin-like mesh. */

const int	BCOUNT	=	10;
const double	BLENGTH =	40.0;
const double	BOFFSETX =	-200.0;
const double	BOFFSETY =	-200.0;
const double	BMINSFORCE	=	0.3;
const double	BMAXSFORCE	=	3.0;
const double	BMASS	=	5.0;
const double	BDAMP	=	5.0;
const double	BBREAK	=	50.0;

static  MASS_node	*mlist = NULL;	/* Points to first node in mass list */
static  SPRING_node	*slist = NULL;	/* Points to first node in spring list */
static  GENERAL_SPRING_node	*gslist = NULL;	/* Points to first node in general spring list */
static  HINGE_node      *hlist = NULL;  /* Points to the first hings in the hinge list */

// Pointers dealing with mouse-based interaction.
static MASS_node *mouse_mlist = NULL;
static double mouse_x, mouse_y;

double	uniform_rand(void)
{
  return rand() / ( (double)RAND_MAX );
}

double	random_spring(double min, double max)
{
  return min + uniform_rand() / (max-min);
}

void	simulate_and_draw_meshes(void)
{
  // Clear color depends on whether we just saved a frame.  If so, then
  // flash as feedback.
  static unsigned last_frame = 0;
  if (last_frame != g_csv_frame_number) {
    glClearColor(0.7, 0.7, 0.7, 1.0);
    glClear(GL_COLOR);
    last_frame = g_csv_frame_number;
  } else {
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR);
  }

  // Draw the masses as spheres and the springs as lines.
  // General springs in yellow, normal ones in green.
  int j,k;
  for (j = 0; j < DRAWS_PER_FRAME; j++) {
      glColor3f(0.7,0.7,0);
      if (draw_general_springs(gslist)) { exit(-2); }
      glColor3f(0,1,0);
      if (draw_springs(slist)) { exit(-1); }
      if (draw_masses(mlist)) { exit(-3); }
      // XXX Draw hinges?

      for (k = 0; k < SIMS_PER_DRAW; k++) {
	apply_springs(&slist);
	apply_general_springs(&gslist);
        apply_hinges(hlist);

        // If we've grabbed a node from either list, drag it around with us.
        if (mouse_mlist) {
          mouse_mlist->x = mouse_x;
          mouse_mlist->y = mouse_y;
        }

	step_masses(mlist, 0.0,GACCEL,0.0, TIME, medium_dx,medium_dy,medium_dz);
      }
  }

} /* end of simulation and drawing 2D fibrin-like mesh. */

void	simulate_and_draw(void)
{
    setup_graphics_frame();

    /* This section is the code being developed with an eye towards
       simulating biological fibers. */
    simulate_and_draw_meshes();

} /* end of simulation and drawing loop */

void display_func(void)
{
  simulate_and_draw();
  glutPostRedisplay();
  glutSwapBuffers();
}

void keyboardCallbackForGLUT(unsigned char key, int x, int y)
{
  switch (key) {
    case 'q':
    case 'Q':
    case 27: // Escape key
      exit(0);
      break;

    case 's':
    case 'S':
      // Save the coordinates of all of the points into the CSV output file.
      // If this works, set a flag to cause the background to flash the next time a frame
      // is drawn, as feedback.
      FILE *csv_file = fopen(g_csv_outfile_name, "a+");
      if (csv_file != NULL) {
        // If this is the first frame, put a header.
        if (g_csv_frame_number == 0) {
          fprintf(csv_file, "FrameNumber,Spot ID,X,Y,Z,Radius,Orientation (if meaningful),Length (if meaningful), Fit Background (for FIONA), Gaussian Summed Value (for FIONA), Mean Background (FIONA), Summed Value (for FIONA)\n");
        }

        // Store the this frame for all existing masses.
        MASS_node *m = mlist;
        unsigned beadno = 0;
        while (m != NULL) {
          fprintf(csv_file, "%d, %d, %lg,%lg,%lg, %lg, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0\n",
            g_csv_frame_number, beadno, m->x, m->y, m->z, m->radius);
          m = m->next;
          beadno++;
        }

        fclose(csv_file);
        g_csv_frame_number++;

      } else {
        perror("Could not write frame to output file");
      }

      break;
  }
}

// Return the squared distance between the points, useful for
// finding the closest points but not requiring square roots.
inline double dist2(double x1, double y1, double x2, double y2)
{
  return (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1);
}

void mouseCallbackForGLUT(int button, int state, int raw_x, int raw_y)
{
  double x = mouse_x_in_world(raw_x);
  double y = mouse_y_in_world(raw_y);

  // Store this state for others to use.
  mouse_x = x;
  mouse_y = y;

  switch(button) {
    // The left mouse button selects a single mass from each of the
    // meshes (the closest to the mouse position) and records it
    // for use in the motion callback.  It also moves the mass to
    // mouse location.
    case GLUT_LEFT_BUTTON:
      if (state == GLUT_DOWN) {

        // Find the closest mass in the first mass list.
        MASS_node *p = mlist;
        double min_d = dist2(x, y, p->x, p->y);
        MASS_node *min_p = p;
        while (p != NULL) {
          if (dist2(x, y, p->x, p->y) < min_d) {
            min_d = dist2(x, y, p->x, p->y);
            min_p = p;
          }
          p = p->next;
        }

        // Say that we've grabbed this point, then move it to our location.
        mouse_mlist = min_p;
        min_p->x = x;
        min_p->y = y;

      } else {
        // Let go of the nodes.
        mouse_mlist = NULL;
      }
      break;

    case GLUT_MIDDLE_BUTTON:
      if (state == GLUT_DOWN) {
      }
      break;

    case GLUT_RIGHT_BUTTON:
      if (state == GLUT_DOWN) {
      }
      break;
  }
}

void motionCallbackForGLUT(int raw_x, int raw_y)
{
  double x = mouse_x_in_world(raw_x);
  double y = mouse_y_in_world(raw_y);

  // Store this state for others to use.
  mouse_x = x;
  mouse_y = y;

  return;
}

void Usage(const char *s)
{
  fprintf(stderr,"Usage: %s [infile_name]\n", s);
  fprintf(stderr,"       infile_name: Configuration file (default webslinger.cfg)\n");
  exit(-1);
}

int main(unsigned argc, const char *argv[])
{
  // File name to read, defaults to webslinger.cfg
  const char  *infile_name = "webslinger.cfg";

  // Command-line parsing
  unsigned  i;
  unsigned  real_params = 0;  //< Number of non-flag parameters
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-help")) {
      Usage(argv[0]);
    } else if (argv[i][0] == '-') {
      Usage(argv[0]);
    } else {
      switch (real_params) {
      case 0:
        infile_name = argv[i];
        break;
      default:
        Usage(argv[0]);
      }
    }
  }

  // Read the structure description from the configuration file.
  FILE *f = fopen(infile_name, "r");
  if (f == NULL) {
    char s[1024];
    sprintf(s, "Could not open configuration file '%s' for reading", infile_name);
    perror(s);
    exit(-3);
  }
  if (!parse_structure_from_file(f, &mlist, &slist, &gslist, &hlist)) {
    fprintf(stderr,"Cannot parse %s\n", infile_name);
    exit(-4);
  }
  fclose(f);

  // Make the name of the output file by adding ".csv" to the input file name
  // Then delete the output file, so it will be created empty if written.
  g_csv_outfile_name = new char[strlen(infile_name) + 10];
  if (g_csv_outfile_name == NULL) {
    fprintf(stderr,"Out of memory creating CSV outfile name\n");
    exit(-5);
  }
  sprintf(g_csv_outfile_name, "%s.csv", infile_name);
  unlink(g_csv_outfile_name);


  init_graphics("Webslinger v2.3: Pull with mouse, Q quit, S save coords", display_func);
  glutMotionFunc(motionCallbackForGLUT);
  glutMouseFunc(mouseCallbackForGLUT);
  glutKeyboardFunc(keyboardCallbackForGLUT);
  glutMainLoop();

  // Never gets here, but to keep the compiler happy...
  return 0;
}
