/*	This program will show the effects of springs on a connected
 * string of mass points.  The springs may or may not be damped.  It uses
 * the massmesh.c routines to access the mass mesh and to perform the
 * position and force updates for it.
 *	This routine uses the Turbo C version 1.5 graphics routines to
 * display the structure.  It is therefore non-portable.
 */

#include <stdio.h>
#include <math.h>
#include <time.h>
#include "massmesh.h"

#define	NUM	12
#define	REST	5.0
#define	FORCE	1.0
#define	MASS	1.0
#define	DAMP	0.0

#define	AMP	0.0
#define	ANGLE	(3.14159/15.0)

main()
{
	MASS_node	*mh;	/* Points to first node in mass list */
	SPRING_node	*sh;	/* Points to first node in spring list */
	HINGE_node	*hh;	/* Points to first node in hinge list */
	int	loop;
	MASS_node	*cm;
	SPRING_node	*cs;
	long		t1,t2;

	printf("Making the string\n");
	if (make_string(NUM, &mh, &sh, &hh, MASS,DAMP, REST,FORCE*2, 0.0)) {
		fprintf(stderr,"Error creating the string.\n");
		exit(-1);
	}

	printf("Starting the timing loop\n");
	loop = 0;
	time(&t1);
	while (loop < 50000) {
		apply_springs(sh);
		mh->x = AMP * sin(ANGLE*loop);
		mh->vx= AMP * cos(ANGLE*loop);
		mh->y = mh->z = mh->vy = mh->vz = 0.0;
		step_masses(mh, 0.0,-FORCE,0.0, 1.0/120.0, 0.0,0.0,0.0);
		loop++;
	}
	time(&t2);
	if (t2 != t1) {
		printf("Rate: %d per second\n", loop/(t2-t1));
	} else {
		printf("No time passed...\n");
	}

} /* end of main */

