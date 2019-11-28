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
#ifdef _WIN32
#include <glut.h>
#elif __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#include "massmesh.h"
#include "graphics.h"

const double	GACCEL	=	-9.8;
const double	TIME	=	(1.0/5000.0);
double	dx = 1.0, dy = 0, dz = 0;

/* The maximum velocity of the wind, how fast the wind changes, and the
   fraction of the screen in which the string can move without changing
   the wind. */

#define	MAXWIND 10.0
#define	WINDCHANGE 0.2
#define	SCREENPART (1.0/3.0)

// Parameters that affect the graphics
const double	SCALE = 0.004;
const int	SIMS_PER_DRAW = 400;
const int	DRAWS_PER_FRAME = 1;

/* =====================================================================
 * Section of code dealing with making and simulating hinged and unhinged string. */

const int	COUNT	=	10;
const double	REST	=	20.0;
const double	SFORCE	=	1.0;
const double	HFORCE	=	500.0;
const double	MASS	=	1.0;
const double	DAMP	=	1.0;

MASS_node	*smh = NULL;	/* Points to first node in mass list */
MASS_node	*smt = NULL;	/* Points to the last node in the mass list */
SPRING_node	*ssh = NULL;	/* Points to first node in spring list */
HINGE_node	*shh = NULL;	/* Points to first node in hinge list */
MASS_node	*smh2 = NULL;	/* Points to first node in mass list */
MASS_node	*smt2 = NULL;	/* Points to the last node in the mass list */
SPRING_node	*ssh2 = NULL;	/* Points to first node in spring list */
HINGE_node	*shh2 = NULL;	/* Points to first node in hinge list */

#define	UNTIE	0

void	simulate_and_draw_strings(void)
{
    static int	first_time = 1;
    static int	loop = 0;
    int j,k;

    if (first_time) {
	first_time = 0;

	make_string(COUNT, &smh, &ssh, &shh, MASS, DAMP, REST, SFORCE, HFORCE);
	smt = smh; while (smt->next != NULL) { smt = smt->next; }
	make_string(COUNT, &smh2, &ssh2, &shh2, MASS, DAMP, REST, SFORCE, 0.0);
	smt2 = smh2; while (smt2->next != NULL) { smt2 = smt2->next; }
    }

    //------------------------------------------------------------
    // The following happens each time through the loop

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	if (draw_springs(ssh)) { exit(-1); }
	if (draw_springs(ssh2)) { exit(-1); }

	for (k = 0; k < SIMS_PER_DRAW; k++) {
	  apply_springs(&ssh);
	  apply_hinges(shh);
	  apply_springs(&ssh2);
	  apply_hinges(shh2);

	  if (UNTIE == 0) {
	    smh->x = smh->y = smh->z = 0;
	    smt->x = REST * (COUNT-1); smt->y = smt->z = 0;
	    smh2->x = smh2->y = smh2->z = 0;
	    smt2->x = REST * (COUNT-1); smt2->y = smt2->z = 0;
	  }

	  step_masses(smh, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(smh);
	  step_masses(smh2, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(smh2);
	  loop++;
	}
    }

} /* end of simulation and drawing strings loop */

/* =====================================================================
 * Section of code dealing with making and simulating breaking string. */

const int	BCOUNT	=	10;
const double	BREST	=	20.0;
const double	BSFORCE	=	2.0;
const double	BHFORCE	=	0.0;
const double	BMASS	=	5.0;
const double	BDAMP	=	1.0;
const double	BBREAK	=	24.0;

MASS_node	*bmh = NULL;	/* Points to first node in mass list */
MASS_node	*bmt = NULL;	/* Points to the last node in the mass list */
SPRING_node	*bsh = NULL;	/* Points to first node in spring list */
HINGE_node	*bhh = NULL;	/* Points to first node in hinge list */

#define	BUNTIE	0

void	simulate_and_draw_breaking(void)
{
    static int	first_time = 1;
    static int	loop = 0;
    int j,k;

    if (first_time) {
	first_time = 0;

	make_string(BCOUNT, &bmh, &bsh, &bhh, BMASS, BDAMP, BREST, BSFORCE, BHFORCE, BBREAK);
	bmt = bmh; while (bmt->next != NULL) { bmt = bmt->next; }
    }

    //------------------------------------------------------------
    // The following happens each time through the loop.

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	if (draw_springs(bsh)) { exit(-1); }

	for (k = 0; k < SIMS_PER_DRAW; k++) {
	  apply_springs(&bsh);
	  apply_hinges(bhh);

	  // Stick a plane at the bottom to keep the string from dropping
	  // out of the viewing volume.
	  apply_plane(bmh, 0, 1, 0, 200.0, 50.0);

	  if (BUNTIE == 0) {
	    bmh->x = bmh->y = bmh->z = 0;
	    bmt->x = BREST * (BCOUNT-1); bmt->y = bmt->z = 0;
	  }

	  step_masses(bmh, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(bmh);
	  loop++;
	}
    }

} /* end of simulation and drawing breaking string loop */

/* =====================================================================
 * Section of code dealing with making and simulating a fat string. */

const int	FCOUNT	=	20;
const double	FREST	=	10.0;
const double	FSFORCE	=	2.0;
const double	FHFORCE	=	0.0;
const double	FMASS	=	5.0;
const double	FDAMP	=	1.0;
const double	FRADIUS	=	20.0;
const double	FBREAK	=	1e100;

MASS_node	*fmh = NULL;	/* Points to first node in mass list */
MASS_node	*fmt = NULL;	/* Points to the last node in the mass list */
SPRING_node	*fsh = NULL;	/* Points to first node in spring list */
HINGE_node	*fhh = NULL;	/* Points to first node in hinge list */
MASS_node	*fmh2 = NULL;	/* Second mass list to be collided with */

void	simulate_and_draw_fat(void)
{
    static int	first_time = 1;
    static int	loop = 0;
    int j,k;

    if (first_time) {
	first_time = 0;

	// Make a string for the first object
	make_string(FCOUNT, &fmh, &fsh, &fhh, FMASS, FDAMP, FREST, FSFORCE, FHFORCE, FBREAK, FRADIUS);
	fmt = fmh; while (fmt->next != NULL) { fmt = fmt->next; }

	// Make a single mass for the other object
	add_mass(&fmh2, FMASS, FCOUNT*FREST/2, -200 + FRADIUS*2, 0.000001, FDAMP, FRADIUS*2);

	// Move both of these over to provide variety
	translate_masses(fmh, -200, 0, 0);
	translate_masses(fmh2, -200, 0, 0);
    }

    //------------------------------------------------------------
    // The following happens each time through the loop.

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	if (draw_masses(fmh)) { exit(-1); }
	if (draw_masses(fmh2)) { exit(-1); }

	for (k = 0; k < SIMS_PER_DRAW; k++) {
	  apply_springs(&fsh);
	  apply_hinges(fhh);

	  // Stick a plane at the bottom to keep the string from dropping
	  // out of the viewing volume.
	  apply_plane(fmh, 0, 1, 0, 200.0, 50.0);
	  apply_plane(fmh2, 0, 1, 0, 200.0, 50.0);

	  // Make the two strings bounce off each other.
	  apply_collisions(fmh, fmh2, 50.0);

	  // Make the first string not self-intersect, except for
	  // masses within four springs of each other.
	  apply_self_collisions(fmh, 50.0, 4);

	  step_masses(fmh, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(fmh);
	  step_masses(fmh2, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(fmh2);
	  loop++;
	}
    }

} /* end of simulation and drawing breaking string loop */

/* =====================================================================
 * Section of code dealing with making and simulating capped cubes. */

#define	NUM_X 3
#define NUM_Y 2
#define NUM_Z 2

const double	CREST	=	40.0;
const double	CSFORCE	=	3.5;
const double	CHFORCE	=	0.0;
const double	CMASS	=	0.1;
const double	CDAMP	=	0.1;

/* Either AMP or MOVE can be nonzero.  Both should not be nonzero at the
 * same time.  If UNTIE is 1, then both of the others are ignored.  */

#define	CUNTIE	0
#define	AMP	30.0
#define	ANGLE	(3.14159 * 0.07)
#define	MOVE	0.0

CAPPED_CUBE	cap_cubes[NUM_X][NUM_Y][NUM_Z];
CUBE		cubes[NUM_X][NUM_Y][NUM_Z];

MASS_node	*cmh = NULL;	/* Points to first node in mass list */
SPRING_node	*csh = NULL;	/* Points to first node in spring list */
HINGE_node	*chh = NULL;	/* Points to first node in hinge list */

void	simulate_and_draw_cubes(void)
{
    static int	first_time = 1;
    static int	loop = 0;
    int j,k;
    const double height = CREST / sqrt(2.0);
    const double FULLSTEP = CREST + 2*height;
    int x,y,z;
    MASS_node *xm, *ym, *zm;

    if (first_time) {
	first_time = 0;

	// Make the structure that is to be simulated.
	for (x = 0; x < NUM_X; x++) {
	 for (y = 0; y < NUM_Y; y++) {
	  for (z = 0; z < NUM_Z; z++) {
//XXX Only the x one seems to work, for some reason...
		if (x > 0) {
			xm = (cap_cubes[x-1][y][z])[1];
		} else {
			xm = NULL;
		}
		if (y > 0) {
			ym = (cap_cubes[x][y-1][z])[3];
		} else {
			ym = NULL;
		}
		if (z > 0) {
			zm = (cap_cubes[x][y][z-1])[5];
		} else {
			zm = NULL;
		}
		make_capped_cube(&cmh, &csh,
			x*FULLSTEP, y*FULLSTEP, z*FULLSTEP,
			CMASS, CDAMP, CREST, CSFORCE,
			xm, ym, zm,
			cap_cubes[x][y][z], cubes[x][y][z]);
	  }
	 }
	}
    }

    //------------------------------------------------------------
    // The following happens each time through the loop

    for (j = 0; j < DRAWS_PER_FRAME; j++) {
	if (draw_springs(csh)) { exit(-1); }

	for (k = 0; k < SIMS_PER_DRAW; k++) {
	  apply_springs(&csh);
	  apply_hinges(chh);

	  if (CUNTIE == 0) {
		  MASS_node *mn = cubes[0][0][0][0];
		  if (AMP != 0.0) {
			  mn->x = AMP * sin(ANGLE*loop*TIME);
			  mn->vx= AMP * cos(ANGLE*loop*TIME);
			  mn->y = mn->z = mn->vy = mn->vz = 0.0;
		  }
		  else {
			  mn->x += MOVE*TIME;
			  mn->vx = MOVE*TIME;
			  mn->y = mn->z = mn->vy = mn->vz = 0.0;
		  }
	  }
	  step_masses(cmh, 0.0,GACCEL,0.0, TIME, dx,dy,dz);
          clear_forces(cmh);
	  loop++;
	}
    }

} /* end of simulation and drawing loop for capped cubes */

void	simulate_and_draw_snotmatrix(void)
{
    static int	first_time = 1;
    static int	loop = 0;
    const int MATRIX_COUNT = 10;
    const unsigned MATRIX_STRAND_LENGTH = 100;
    const double MMASS = 0.1;
    const double MDAMP = 0.1;
    const double MREST = 0.0;
    const double MSFORCE = 0.1;
    const double MHFORCE = 0.0;
    const double MRADIUS = 15.0;
    const double MSPACING = 45.0;
    const double MY = 250.0;
    static MASS_node	*mmh[MATRIX_COUNT];	/* First node in each the matrix mass list. */
    static SPRING_node	*msh[MATRIX_COUNT];	/* First node in each of the matrix spring lists. */
    static HINGE_node	*mhh[MATRIX_COUNT];
    int i;
    int j,k;

    if (first_time) {
	first_time = 0;

	for (i = 0; i < MATRIX_COUNT; i++) {
	  make_string(MATRIX_STRAND_LENGTH, &mmh[i], &msh[i], &mhh[i], MMASS + ((i * 1047) % 10)*0.3*MMASS, MDAMP, MREST, MSFORCE, MHFORCE, 1e100, MRADIUS);
	  translate_masses(mmh[i], MSPACING*(i-MATRIX_COUNT/2), MY,0);
	}
    }

    //------------------------------------------------------------
    // The following happens each time through the loop

    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (j = 0; j < DRAWS_PER_FRAME; j++) {
      for (i = 0; i < MATRIX_COUNT; i++) {
	if (draw_masses(mmh[i])) { exit(-1); }
	if (draw_springs(msh[i])) { exit(-1); }
      }

      for (k = 0; k < SIMS_PER_DRAW; k++) {
	for (i = 0; i < MATRIX_COUNT; i++) {
	  apply_springs(&msh[i]);

	  if (UNTIE == 0) {
	    mmh[i]->x = MSPACING*(i - MATRIX_COUNT/2);
	    mmh[i]->y = MY;
	    mmh[i]->z = 0;
	  }

	  step_masses(mmh[i], 0.0,GACCEL,0.0, TIME, 0,0,0);
          clear_forces(mmh[i]);
	}
	loop++;
      }
    }

} /* end of simulation and drawing snotmatrix loop */


void	simulate_and_draw(void)
{
    setup_graphics_frame();

    /* This section is the code being developed with an eye towards
       simulating biological fibers.
    glColor3f(1,0,0);
    simulate_and_draw_cubes();
    glColor3f(0,1,0);
    simulate_and_draw_strings();
    glColor3f(0,0,1);
    simulate_and_draw_breaking();
    glColor3f(1,1,0);
    simulate_and_draw_fat();
    */

    /* This section is being developed with an eye towards a
       stupid Matrix knock-off.
       */
    glColor3f(0,1,0);
    simulate_and_draw_snotmatrix();

} /* end of simulation and drawing loop */

void display_func(void)
{
	simulate_and_draw();
	glutPostRedisplay();
	glutSwapBuffers();
}

int main()
{
	init_graphics("Various simulations", display_func);
	glutMainLoop();
	return 0;
}
