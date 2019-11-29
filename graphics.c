#ifdef _WIN32
#include <windows.h>
#include <float.h>
#define finite(x) _finite(x)
#endif

#include <stdio.h>
#include <math.h>

#ifdef _WIN32
  #include <GL/gl.h>
  #include <GL/glu.h>
  #include <glut.h>
#elif __APPLE__
  //#include <OpenGL/gl.h>
  //#include <OpenGL/glu.h>
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif

#include "massmesh.h"
#include "graphics.h"

const int SCREEN_X=512;
const int SCREEN_Y=512;

int getmaxx(void) { return SCREEN_X; }
int getmaxy(void) { return SCREEN_Y; }

double	SCREEN_SCALE = 0.004;
extern SPRING_node* remove_nd;

void resize_window(int w, int h)
{
}

static int vertex(const MASS_node *mn)
{
	// Returns error if the points are not finite.
	if ( !finite(mn->x) || !finite(mn->y) || !finite(mn->z) ) {
		fprintf(stderr,"Infinite points can't be drawn...\n");
		return -1;
	}
	glVertex3d(mn->x, mn->y, mn->z);

	return 0;
}

int	draw_masses(const MASS_node *mn)
{
  static bool first = true;
  static GLuint sphereDlist;

  if (first) {
    GLUquadric	*myQuadric = gluNewQuadric();
    sphereDlist=glGenLists(1);
    glNewList(sphereDlist, GL_COMPILE);
    gluSphere(myQuadric, 1.0, 10,10);
    glEndList();

    first = false;
  }
    while (mn != NULL) {
      glPushMatrix();
	glTranslated(mn->x, mn->y, mn->z);
	glScaled(mn->radius, mn->radius, mn->radius);
	glCallList(sphereDlist);
	mn = mn->next;
      glPopMatrix();
    }
    return 0;
}

int	draw_springs(const SPRING_node *sn)
{
    if (sn != NULL) {
	float currentColor[4];
	glGetFloatv(GL_CURRENT_COLOR, currentColor);
	glBegin(GL_LINES);
	while (sn != NULL) {
		if (sn->flag == INVALID) {
			sn = sn->next;
			continue;
		}
		else if (sn->flag == SELECTED) {
			glColor3f(1, 0.5, 0);
		}
		vertex(sn->m1);
		vertex(sn->m2);
		if (sn->flag == SELECTED) {
			glColor3f(currentColor[0], currentColor[1], currentColor[2]);
		}
		sn = sn->next;
	}
	glEnd();
    }
    return 0;
}

int	draw_general_springs(const GENERAL_SPRING_node *sn)
{
    if (sn != NULL) {
	glBegin(GL_LINES);
	while (sn != NULL) {
		vertex(sn->m1);
		vertex(sn->m2);
		sn = sn->next;
	}
	glEnd();
    }
    return 0;
}

int	setup_graphics_frame(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glScaled(SCREEN_SCALE, SCREEN_SCALE, SCREEN_SCALE);

	return 0;
}

int	init_graphics(const char *name, Display_Func df)
{
	static	int	no_args = 1;
	static	char	*fargv[1] = {"Glut_program"};
	int	glut_window;

	//---------------------------------------------------------
	// Ask for a window, using Glut.

	static int glutAttributeList = GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH;
	glutInit(&no_args, fargv);
	glutInitDisplayMode(glutAttributeList);
	glutInitWindowSize(SCREEN_X, SCREEN_Y);
	glut_window = glutCreateWindow(name);
	glutSetWindow(glut_window);
	glutDisplayFunc(df);
	glutIdleFunc(df);
	glutReshapeFunc(resize_window);

	//---------------------------------------------------------
	// Set up the graphics defaults.

	glClearColor(0.0,0.0,0.0,0.0);
	glClearDepth(1.0);
//	glEnable(GL_LINE_SMOOTH);

	return 0;
}

double mouse_x_in_world(int mouse_raw_x)
{
  // Convert from window coordinates to the range 0..1
  double norm = mouse_raw_x / static_cast<double>(SCREEN_X-1);

  // Convert from that to -1..1
  double ndc = (norm-0.5) * 2;

  // Convert from that to world coordinates
  return ndc / SCREEN_SCALE;
}

double mouse_y_in_world(int mouse_raw_y)
{
  // Convert from window coordinates to the range 0..1
  double norm = mouse_raw_y / static_cast<double>(SCREEN_Y-1);

  // Convert from that to -1..1
  double ndc = (norm-0.5) * 2;

  // Swap the sign; Y axis is inverted.
  ndc *= -1;

  // Convert from that to world coordinates
  return ndc / SCREEN_SCALE;
}
