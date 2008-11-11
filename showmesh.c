/*	This program will show the effects of springs on a connected
 * string of mass points.  The springs may or may not be damped.  It uses
 * the massmesh.c routines to access the mass mesh and to perform the
 * position and force updates for it.
 *	This program currently shows a string being blown about on the
 * screen by winds that change.  The winds tend to keep the string near
 * the center of the screen.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "massmesh.h"
#include "graphics.h"

#define	NUM	70
#define	REST	0.01
#define	SFORCE	4.0
#define	HFORCE	0.0
#define	GACCEL	-9.8
#define	MASS	0.1
#define	DAMP	0.1
#define	TIME	(1.0/2000.0)
#define	ENDMULTIPLE 3.0
#define	HEADMULTIPLE 1.0

/* Either AMP or MOVE can be nonzero.  Both should not be nonzero at the
 * same time.  If UNTIE is 1, then both of the others are ignored.  */

#define	UNTIE	0
#define	AMP	50.0
#define	ANGLE	(3.14159 * 1)
#define	MOVE	0.0

/* The maximum velocity of the wind, how fast the wind changes, and the
   fraction of the screen in which the string can move without changing
   the wind. */

#define	MAXWIND 40.0
#define	WINDCHANGE 2.0
#define	SCREENPART (1.0/3.0)

int	SIMS_PER_DRAW = 100;
int	DRAWS_PER_FRAME = 1;

MASS_node	*mh;	/* Points to first node in mass list */
SPRING_node	*sh;	/* Points to first node in spring list */
HINGE_node	*hh;	/* Points to first node in hinge list */
int		loop;
MASS_node	*cm;
SPRING_node	*cs;
double	dx,dy,dz;

void	simulate_and_draw(void)
{
    static int	first_time = 1;
    int j,k;

    if (first_time) {
	first_time = 0;

	if (make_string(NUM, &mh, &sh, &hh, MASS,DAMP, REST,SFORCE, HFORCE)) {
		fprintf(stderr,"Error creating the string.\n");
		exit(-1);
	}

	/* Make the last one heavier */

	cm = mh;
	cm->mass *=HEADMULTIPLE;
	while (cm->next != NULL)
		cm = cm->next;
	cm->mass *= ENDMULTIPLE;


	dx = 5*WINDCHANGE;
	dy = dz = 0.0;
	loop = 0;
    }

	//------------------------------------------------------------
	// The following happens each time through the loop

	setup_graphics_frame();

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	cm = mh;
	if (draw_masses(mh)) {
		exit(-1);
	}
      for (k = 0; k < SIMS_PER_DRAW; k++) {
	apply_springs(sh);
	apply_hinges(hh);
	if (UNTIE == 0) {
		if (AMP != 0.0) {
			mh->x = AMP * sin(ANGLE*loop*TIME);
			mh->vx= AMP * cos(ANGLE*loop*TIME);
			mh->y = mh->z = mh->vy = mh->vz = 0.0;
		}
		else {
			mh->x += MOVE*TIME;
			mh->vx = MOVE*TIME;
			mh->y = mh->z = mh->vy = mh->vz = 0.0;
		}
	}
	if (mh->x > (getmaxx()/2*SCREENPART)) dx -= WINDCHANGE;
	if (dx > MAXWIND) dx = MAXWIND;
	if (mh->x < -(getmaxx()/2*SCREENPART)) dx += WINDCHANGE;
	if (dx < -MAXWIND) dx = -MAXWIND;
	if (mh->y > (getmaxy()/2*SCREENPART)) dy -= WINDCHANGE;
	if (dy > MAXWIND) dy = MAXWIND;
	if (mh->y < -(getmaxy()/2*SCREENPART)) dy += WINDCHANGE;
	if (dy < 0) dy = 0;
	step_masses(mh, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
        clear_forces(mh);
	loop++;
      }
    }

} /* end of simulation and drawing loop */

void display_func(void)
{
	simulate_and_draw();
	glutPostRedisplay();
	glutSwapBuffers();
}

// Initialize graphics window

main()
{
	init_graphics("ShowMesh", display_func);
	glutMainLoop();
}

