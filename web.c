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

//const double	GACCEL	=      -9.8;	//< Acceleration due to gravity
const double	GACCEL	=      0;	//< Acceleration due to gravity
const double	TIME	=	(1.0/5000.0);
//double	dx = 1.0, dy = 0, dz = 0;
double	dx = 0.0, dy = 0, dz = 0;

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

static  MASS_node *masses1[BCOUNT][BCOUNT];
static  MASS_node *masses2[BCOUNT][BCOUNT];
static  MASS_node	*mlist1 = NULL;	/* Points to first node in mass list for the first mesh */
static  MASS_node	*mlist2 = NULL;	/* Points to first node in mass list for the second mesh */
static  SPRING_node	*slist1 = NULL;	/* Points to first node in spring list for first mesh */
static  SPRING_node	*slist2 = NULL;	/* Points to first node in spring list for second mesh */


// Pointers dealing with mouse-based interaction.
static MASS_node *mouse_mlist1 = NULL;
static MASS_node *mouse_mlist2 = NULL;
static double mouse_x, mouse_y;

static bool UNTIE = false;

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
    static bool	first_time = true;
    static int	loop = 0;
    int i,j,k;

    if (first_time) {
	first_time = false;

	// Make two 2D meshes of masses, overlaid on top of each other.
	for (i = 0; i < BCOUNT; i++) {
	  for (j = 0; j < BCOUNT; j++) {
	    add_mass(&mlist1, BMASS, i*BLENGTH + BOFFSETX,j*BLENGTH + BOFFSETY,0, BDAMP);
	    masses1[i][j] = mlist1;
	    add_mass(&mlist2, BMASS, i*BLENGTH + BOFFSETX,j*BLENGTH + BOFFSETY,0, BDAMP);
	    masses2[i][j] = mlist2;
	  }
	}

	// Add the horizontal springs, picking a new spring constant for each line of them.
	// Make the rest length slightly less than the distance between the points.
	for (j = 0; j < BCOUNT; j++) {
	  double kspring1 = random_spring(BMINSFORCE,BMAXSFORCE);
	  double kspring2 = random_spring(BMINSFORCE,BMAXSFORCE);
	  for (i = 0; i < BCOUNT-1; i++) {
	    add_spring(&slist1, 0.8*BLENGTH, kspring1, masses1[i][j], masses1[i+1][j], BBREAK);
	    add_spring(&slist2, 0.8*BLENGTH, kspring2, masses2[i][j], masses2[i+1][j], BBREAK);
	  }
	}

	// Add the vertical springs, picking a new spring constant for each line of them.
	// Make the rest length slightly less than the distance between the points.
	for (i = 0; i < BCOUNT; i++) {
	  double kspring1 = random_spring(BMINSFORCE,BMAXSFORCE);
	  double kspring2 = random_spring(BMINSFORCE,BMAXSFORCE);
	  for (j = 0; j < BCOUNT-1; j++) {
	    add_spring(&slist1, 0.8*BLENGTH, kspring1, masses1[i][j], masses1[i][j+1], BBREAK);
	    add_spring(&slist2, 0.8*BLENGTH, kspring2, masses2[i][j], masses2[i][j+1], BBREAK);
	  }
	}
    }

    //------------------------------------------------------------
    // The following happens each time through the loop.

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	glColor3f(0,1,0);
	if (draw_springs(slist1)) { exit(-1); }

	glColor3f(1,1,0);
	if (draw_springs(slist2)) { exit(-1); }

	for (k = 0; k < SIMS_PER_DRAW; k++) {
	  apply_springs(&slist1);
	  apply_springs(&slist2);

	  if (!UNTIE) {
/*
	    // Tie down individual masses.
	    i = 0; j = 0;
	    masses1[i][j]->x = i*BLENGTH + BOFFSETX; masses1[i][j]->y = j*BLENGTH + BOFFSETY; masses1[i][j]->z = 0;
	    masses2[i][j]->x = i*BLENGTH + BOFFSETX; masses2[i][j]->y = j*BLENGTH + BOFFSETY; masses2[i][j]->z = 0;
	    i = 7; j = 5;
	    masses1[i][j]->x = i*BLENGTH + BOFFSETX; masses1[i][j]->y = j*BLENGTH + BOFFSETY; masses1[i][j]->z = 0;
	    masses2[i][j]->x = i*BLENGTH + BOFFSETX; masses2[i][j]->y = j*BLENGTH + BOFFSETY; masses2[i][j]->z = 0;
*/
            // Tie down lines of masses.
	    for (i = 0; i < BCOUNT; i++) {
	      // Hold the masses along the top edge still for both meshes
	      masses1[i][BCOUNT-1]->x = i*BLENGTH + BOFFSETX; masses1[i][BCOUNT-1]->y = (BCOUNT-1)*BLENGTH + BOFFSETY; masses1[i][BCOUNT-1]->z = 0;
	      masses2[i][BCOUNT-1]->x = i*BLENGTH + BOFFSETX; masses2[i][BCOUNT-1]->y = (BCOUNT-1)*BLENGTH + BOFFSETY; masses2[i][BCOUNT-1]->z = 0;
	    }
	  }

          // If we've grabbed a node from either list, drag it around with us.
          if (mouse_mlist1) {
            mouse_mlist1->x = mouse_x;
            mouse_mlist1->y = mouse_y;
          }
          if (mouse_mlist2) {
            mouse_mlist2->x = mouse_x;
            mouse_mlist2->y = mouse_y;
          }

	  step_masses(mlist1, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
	  step_masses(mlist2, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
	  loop++;
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

  case 'u':
  case 'U':
    // Untie/tie the points.
    UNTIE = !UNTIE;
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
        MASS_node *p = mlist1;
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
        mouse_mlist1 = min_p;
        min_p->x = x;
        min_p->y = y;

        // Find the closest mass in the first mass list.
        p = mlist2;
        min_d = dist2(x, y, p->x, p->y);
        min_p = p;
        while (p != NULL) {
          if (dist2(x, y, p->x, p->y) < min_d) {
            min_d = dist2(x, y, p->x, p->y);
            min_p = p;
          }
          p = p->next;
        }

        // Say that we've grabbed this point, then move it to our location.
        mouse_mlist2 = min_p;
        min_p->x = x;
        min_p->y = y;

      } else {
        // Let go of the nodes.
        mouse_mlist1 = NULL;
        mouse_mlist2 = NULL;
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


void main()
{
  init_graphics("Pull with mouse, U to untie", display_func);
  glutMotionFunc(motionCallbackForGLUT);
  glutMouseFunc(mouseCallbackForGLUT);
  glutKeyboardFunc(keyboardCallbackForGLUT);
  glutMainLoop();
}
