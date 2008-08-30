/************************************************************************
 *	This file contains the routines needed to simulate 		*
 * objects as masses tied together by springs and hinges.  It includes 	*
 * the routines that set up the mass spring, and hinge lists, as well as*
 * the routines that control the motion simulation.  The file		*
 * massmesh.h needs to be included by other files that use these	*
 * routines.								*
 *	The routines in this file can be used to model 1, 2, or 	*
 * 3-dimensional objects in this manner.  The difference in		*
 * dimension is only a function of how the springs and hinges are 	*
 * attatched to the masses in the list.					*
 *	The normal usage of these routines is to build a mesh of masses	*
 * connected by springs and/or hinges.  Then a loop is entered in which	*
 * the forces acting on the masses are computed and then the forces are	*
 * applied to the masses.  The application of the forces clears the old	*
 * force values so that new ones can be computed for the next iteration.*
 ************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "massmesh.h"
#include <vector>

using namespace std;

#ifndef	M_PI
const double M_PI = 2*acos(0.0);
#endif

#pragma warning( disable : 4996 )

inline	double	sqr(double a) { return a*a; }

//#define	DEBUG

/****************************************************************
 *	These are utility routines that allow applications to	*
 * build the mass list, the hinge list and the spring list.	*
 ****************************************************************/

/*	This routine will add a mass to the head of a list of masses.  It gets
 * a pointer to the pointer to the head of the list.  If the head of the
 * list is null, then it is initialized.  In all cases, the new mass
 * becomes the head of the mass list. The initial force and velocity on the
 * new mass are set to zero.
 *	The function returns -1 on failure, 0 on success.
 */

int	add_mass(MASS_node **headptr, double mass, double x, double y, double z,
	double damping, double radius)
{
	MASS_node *newm; 	/* Pointer to the new mass */
	if ( (newm = (MASS_node *)malloc(sizeof(MASS_node))) == NULL ) {
		perror("add_mass (cannot allocate more memory)");
		return(-1);	/* FAIL */
	}

	newm->next = *headptr;	/* Insert before head of the list */
	*headptr = newm;	/* Point the head to this entry */

	newm->mass = mass;	/* Fill in the values for this mass */
	newm->x = x;
	newm->y = y;
	newm->z = z;
	newm->damping = damping;
	newm->radius = radius;
	newm->vx = newm->vy = newm->vz = 0.0;
	newm->fx = newm->fy = newm->fz = 0.0;

	return(0);		/* SUCCESS */

} /* end of add_mass() */


/*	This routine will add a spring to the head of a list of springs.
 * It gets a pointer to the pointer to the head of the list.  It will put
 * the new spring at the head of the list.  It works properly if the head
 * of the list is NULL to begin with.  m1 and m2 should point to the masses
 * at either end of the spring.
 *	The function returns -1 on failure and 0 on success.
 */

int	add_spring(SPRING_node **headptr, double rest, double k,
	MASS_node *m1, MASS_node *m2, double fbreak)
{
	SPRING_node	*news;	/* The new spring node location */

	if ( (m1 == NULL) || (m2 == NULL) ) {
		perror("add_spring (bad mass)");
		return(-1);	/* BAD MASS */
	}
	if ( (news = (SPRING_node *)malloc(sizeof(SPRING_node))) == NULL) {
		perror("add_spring (cannot allocate more memory)");
		return(-1);	/* FAIL */
	}

	news->next = *headptr;	/* Place new entry at head of list */
	*headptr = news;	/* Make the head point to the new one */

	news->rest = rest;	/* Fill in the parameter values */
	news->k = k;
	news->fbreak = fbreak;
	news->m1 = m1;
	news->m2 = m2;

	return(0);	/* SUCCESS */

} /* end of add_spring() */


/*	This routine will add a general spring to the head of a list of
 * general springs. It gets a pointer to the pointer to the head of the list.
 * It will put the new spring at the head of the list.  It works properly if
 * the head of the list is NULL to begin with.  m1 and m2 should point to
 * the masses at either end of the spring.
 *	The function returns -1 on failure and 0 on success.
 */

int	add_general_spring(GENERAL_SPRING_node **headptr, double rest, double k,
                           GENERAL_SPRING_TABLE_ENTRY  *table,
                           MASS_node *m1, MASS_node *m2,
                           double fbreak)
{
	GENERAL_SPRING_node	*news;	/* The new spring node location */

	if ( (m1 == NULL) || (m2 == NULL) ) {
		perror("add_general_spring (bad mass)");
		return(-1);	/* BAD MASS */
	}
	if ( (news = (GENERAL_SPRING_node *)malloc(sizeof(GENERAL_SPRING_node))) == NULL) {
		perror("add_general_spring (cannot allocate more memory)");
		return(-1);	/* FAIL */
	}

	news->next = *headptr;	/* Place new entry at head of list */
	*headptr = news;	/* Make the head point to the new one */

	news->rest = rest;	/* Fill in the parameter values */
	news->k = k;
        news->table = table,
	news->fbreak = fbreak;
	news->m1 = m1;
	news->m2 = m2;

	return(0);	/* SUCCESS */

} /* end of add_general_spring() */


/*	This routine will add a hinge to the head of a list of hinges.
 * It gets a pointer to the pointer to the head of the list.  It will put
 * the new hinge at the head of the list.  It works properly if the head
 * of the list is NULL to begin with.  m1, m2 and m3 should point to the
 * masses connected by the hinge.  m2 is in the middle.
 *	The function returns -1 on failure and 0 on success.
 */

int	add_hinge(HINGE_node **headptr, double k,
	MASS_node *m1, MASS_node *m2, MASS_node *m3)
{
	HINGE_node	*newh;	/* The new hinge node location */

	if ( (m1 == NULL) || (m2 == NULL) || (m3 == NULL) ) {
		perror("add_hinge (bad mass)");
		return(-1);	/* BAD MASS */
	}
	if ( (newh = (HINGE_node *)malloc(sizeof(HINGE_node))) == NULL) {
		perror("add_hinge (cannot allocate more memory)");
		return(-1);	/* FAIL */
	}

	newh->next = *headptr;	/* Place new entry at head of list */
	*headptr = newh;	/* Make the head point to the new one */

	newh->k = k;		/* Fill in the parameter values */
	newh->m1 = m1;
	newh->m2 = m2;
	newh->m3 = m3;

	return(0);	/* SUCCESS */

} /* end of add_hinge() */


/************************************************************************
 *	These routines work with the completed structure consisting of	*
 * the mass list, spring list, hinge list, and any external forces acting*
 * on the masses. The spring list is assumed to be linked into the mass	*
 * list via the fields in the spring nodes.  The routines work if passed*
 * empty lists, but this serves no purpose.				*
 ************************************************************************/

/*	This routine will zero the forces on all masses in the mass list
 * that is passed to it.  If the step_masses() function is used, it will
 * automatically clear the forces on all masses as it computes the new
 * position, so this routine is not normally needed.
 *	The forces must be cleared before the new forces acting on each mass
 * can be accumulated.
 *	This routine is passed a pointer to the head of the mass list.
 * Note that this is different from the pointer that is passed to the
 * add_mass() function.
 */

void clear_forces(MASS_node *mlist)
{
  while (mlist != NULL) {
    mlist->fx = mlist->fy = mlist->fz = 0.0;
    mlist = mlist ->next;
  }
}


/*	This routine will compute the forces applied by all of the springs
 * in the spring list.  It adds the force applied by each spring to the
 * force values of the masses on the ends of the spring.  Once the entire
 * spring list has been traversed, the effects off all springs will have been
 * added to the other forces acting on the masses in the mass list.
 *	This routine is called with a handle so that it can affect even the
 * first node in a list (to remove a broken spring).  It returns the number
 * of springs that broke, -1 on failure.
 *	If a spring breaks while it is being simulated, it is removed from
 * the list and its memory is freed.  The masses it points to will not have
 * any forces added to them by the spring.
 */

int apply_springs(SPRING_node **sh)
{
  double  dx,dy,dz;	/* The distance vector for the two ends */
  double  dist;		/* Distance between the two end masses */
  double  scale;	/* Scale on the distance vector making force */
  int	  num_broke = 0;/* How many springs broke during the run? */
  SPRING_node *last, *cur;

  last = cur = *sh;
  while (cur != NULL) {	/* Do for each spring in the list */

    /* Find the distance vector from one to the other */
    dx = (cur->m2->x) - (cur->m1->x);
    dy = (cur->m2->y) - (cur->m1->y);
    dz = (cur->m2->z) - (cur->m1->z);

    /* Find the current length of the spring */
    dist = sqrt( dx*dx + dy*dy + dz*dz );

    /* Find the scale factor for the distance vector that will
     * turn it into a force value.  Then divide the factor by
     * two so that half of this force is applied at either end
     * of the spring.  Apply the scale to the distance vector
     * to make it into a force vector.
     * If there is enough force to break the spring, then we
     * delete this one from the list and don't apply and force
     * to its masses.
     */

    scale = (dist - cur->rest) * cur->k / dist;
    if (scale >= cur->fbreak) {
      SPRING_node *newnext = cur->next;
      last->next = newnext;
      delete cur;
      num_broke++;
      if (last == cur) { /* Head of the list; make the change stick. */
	*sh = newnext;
      }
      cur = newnext;
      continue;
    }
    scale /= 2.0;
    dx *= scale;
    dy *= scale;
    dz *= scale;

    /* Apply the force to each end of the spring.  Note that
     * the force is opposite at each end.
     */

    cur->m1->fx += dx;
    cur->m1->fy += dy;
    cur->m1->fz += dz;
    cur->m2->fx -= dx;
    cur->m2->fy -= dy;
    cur->m2->fz -= dz;

    last = cur;
    cur = cur->next;	/* Next spring in the list */
  }

  return num_broke;
} /* end of apply_springs() */


/*	This routine will compute the forces applied by all of the general
 * springs in the spring list.  It adds the force applied by each spring to the
 * force values of the masses on the ends of the spring.  Once the entire
 * spring list has been traversed, the effects off all springs will have been
 * added to the other forces acting on the masses in the mass list.
 *	This routine is called with a handle so that it can affect even the
 * first node in a list (to remove a broken spring).  It returns the number
 * of springs that broke, -1 on failure.
 *	If a spring breaks while it is being simulated, it is removed from
 * the list and its memory is freed.  The masses it points to will not have
 * any forces added to them by the spring.
 */

int apply_general_springs(GENERAL_SPRING_node **sh)
{
  double  dx,dy,dz;	/* The distance vector for the two ends */
  double  dist;		/* Distance between the two end masses */
  double  force;        /* Force applied by the spring */
  int	  num_broke = 0;/* How many springs broke during the run? */
  GENERAL_SPRING_node *last, *cur;

  last = cur = *sh;
  while (cur != NULL) {	/* Do for each spring in the list */

    /* Make sure we have a valid table.  If not, no force. */
    if (cur->table == NULL) {
      fprintf(stderr,"Empty general spring list!\n");
      continue;
    }

    /* Find the distance vector from one to the other */
    dx = (cur->m2->x) - (cur->m1->x);
    dy = (cur->m2->y) - (cur->m1->y);
    dz = (cur->m2->z) - (cur->m1->z);

    /* Find the current length of the spring */
    dist = sqrt( dx*dx + dy*dy + dz*dz );

    /* Find the force based on the distance vector that will
     * turn it into a force value.  This is done by calculating
     * the strain (change in length over rest length) of the
     * spring and then looking up the corresponding force
     * (or linearly interpolating) in the force/strain table.
     * The result is multiplied by the spring constant, k.
     */

    /* Initialize the force so that if the distance is below
     * the first table entry we still interpolate correctly.
     */
    double strain = dist / cur->rest;
    GENERAL_SPRING_TABLE_ENTRY *prev = cur->table;
    GENERAL_SPRING_TABLE_ENTRY *next = cur->table->next;
    if (strain < prev->strain) {
      force = prev->force * cur->k;
    } else do {
      if (next == NULL) {
        force = prev->force * cur->k;
        break;
      }
      if (strain > next->strain) {
        prev = next;
        next = next->next;
      } else {
        double frac = (strain - prev->strain) / (next->strain - prev->strain);
        force = (frac*next->force + (1.0 - frac)*prev->force) * cur->k;
        break;
      }
    } while (true);

    /* If there is enough force to break the spring, then we
     * delete this one from the list and don't apply and force
     * to its masses.
     */
    if (force >= cur->fbreak) {
      GENERAL_SPRING_node *newnext = cur->next;
      last->next = newnext;
      delete cur;
      num_broke++;
      if (last == cur) { /* Head of the list; make the change stick. */
	*sh = newnext;
      }
      cur = newnext;
      continue;
    }

    /* Divide the factor by
     * two so that half of this force is applied at either end
     * of the spring.
     */
    force /= 2.0;

    /* Normalize the distance vector and multiply by force */
    double fx = dx / dist * force;
    double fy = dy / dist * force;
    double fz = dz / dist * force;

    /* Apply the force to each end of the spring.  Note that
     * the force is opposite at each end.
     */

    cur->m1->fx += fx;
    cur->m1->fy += fy;
    cur->m1->fz += fz;
    cur->m2->fx -= fx;
    cur->m2->fy -= fy;
    cur->m2->fz -= fz;

    last = cur;
    cur = cur->next;	/* Next spring in the list */
  }

  return num_broke;
} /* end of apply_general_springs() */


/*	This routine will add the values of the forces applied by the
 * hinges to the forces acting on the masses.  It acts in the plane
 * passing through the three masses, and perpendicular to the line
 * from the center mass to each end mass.  It applies a torque between
 * these two lines (the torque is applied independently on the two
 * lines).  The torque is converted to a force by the arm lengths, and
 * acts equally on the center mass and the far mass for each arm,
 * so it results in a force on each of the end masses acting perpendicular
 * to the line between them and the center mass, in the plane of the three
 * masses.  It also results in a pair of forces acting on the center mass,
 * pushing it towards the line between the two end masses.
 *	Note that this routine is called with a pointer to the first node
 * in the hinge list, not with a handle as in add_hinge().
 */

void apply_hinges(HINGE_node *hh)
{
  while (hh != NULL) {	/* Do for each hinge in the list */

    /* Compute the vectors from the center mass to each end mass. */
    double v1x = hh->m1->x - hh->m2->x;
    double v1y = hh->m1->y - hh->m2->y;
    double v1z = hh->m1->z - hh->m2->z;
    double v3x = hh->m3->x - hh->m2->x;
    double v3y = hh->m3->y - hh->m2->y;
    double v3z = hh->m3->z - hh->m2->z;

    /* Normalize these vectors, unless one has zero length. */
    double d1 = sqrt(v1x*v1x + v1y*v1y + v1z*v1z);
    double d3 = sqrt(v3x*v3x + v3y*v3y + v3z*v3z);
    if ( (d1 == 0) || (d3 == 0) ) {
      hh = hh->next;
      continue;
    }
    double n1x = v1x / d1;
    double n1y = v1y / d1;
    double n1z = v1z / d1;
    double n3x = v3x / d3;
    double n3y = v3y / d3;
    double n3z = v3z / d3;

    /* Compute the torque magnitude, based on the angle between the
     * two arms and the hinge constant.  The angle is found by taking
     * the arc-cosine of the dot product between the two normalized
     * vectors and then subtracting from 180; it is computed in degrees.
     * Half of this will be used on each arm. */
    double angle = 180.0 - 180.0 / M_PI * acos( n1x*n3x + n1y*n3y + n1z*n3z );
#ifdef	DEBUG
    fprintf(stderr,"Hinge dot = %lg\n", n1x*n3x + n1y*n3y + n1z*n3z);
#endif
    double halftorque = angle * hh->k / 2.0;
#ifdef	DEBUG
    fprintf(stderr,"Hinge Halftorque = %lg\n", halftorque);
#endif

    /* If there is no torque, then we're done with this hinge.  Quit now
     * to avoid singularities in the calculations below.  */
    if (halftorque != 0.0) {

      /* Find the vector perpendicular to the plane containing the
       * two arms.  This is used below to compute each of the torques 
       * on the arms.  */
      double planex = n1y*n3z - n1z*n3y;
      double planey = n1z*n3x - n1x*n3z;
      double planez = n1x*n3y - n1y*n3x;

      /*===================================================================*/
      /* Arm 1-2: Apply half of the force resulting from the torque to the
       * arm between the first and center mass.  The point of application
       * will be halfway down the arm, resulting in the forces acting
       * equally on the first and center masses; each will get half of
       * the force, and in opposite directions.  The direction is in the
       * plane of the three masses and perpendicular to the arm between
       * the first and center mass.  Find the direction by crossing the
       * two arms to find the out-of-plane normal (above) and then crossing that
       * with the line between the first and center masses to get the
       * line in the plane that is perpendicular to it.  Normalize the
       * vector and scale it by half its length, and then apply half of
       * the torque based on that.  */
      double vec1x = planey*n1z - planez*n1y;
      double vec1y = planez*n1x - planex*n1z;
      double vec1z = planex*n1y - planey*n1x;
      double vec1len = sqrt( vec1x*vec1x + vec1y*vec1y + vec1z*vec1z );
#ifdef	DEBUG
      fprintf(stderr,"Hinge Vec1len = %lg\n", vec1len); //XXX BUG HERE!
#endif
      double f1x = halftorque * vec1x / (2.0 * d1 * vec1len);
      double f1y = halftorque * vec1y / (2.0 * d1 * vec1len);
      double f1z = halftorque * vec1z / (2.0 * d1 * vec1len);

#ifdef	DEBUG
      fprintf(stderr,"Hinge Force1 = %lg, %lg, %lg\n", f1x, f1y, f1z);
#endif

      hh->m1->fx -= f1x;
      hh->m1->fy -= f1y;
      hh->m1->fz -= f1z;

      hh->m2->fx += f1x;
      hh->m2->fy += f1y;
      hh->m2->fz += f1z;

      /*===================================================================*/
      /* Arm 2-3: Apply half of the force resulting from the torque to the
       * arm between the third and center mass.  The point of application
       * will be halfway down the arm, resulting in the forces acting
       * equally on the third and center masses; each will get half of
       * the force, and in opposite directions.  The direction is in the
       * plane of the three masses and perpendicular to the arm between
       * the third and center mass.  Find the direction by crossing the
       * two arms to find the out-of-plane normal (above) and then crossing that
       * with the line between the first and center masses to get the
       * line in the plane that is perpendicular to it.  Normalize the
       * vector and scale it by half its length, and then apply half of
       * the torque based on that.  */
      double vec3x = planey*n3z - planez*n3y;
      double vec3y = planez*n3x - planex*n3z;
      double vec3z = planex*n3y - planey*n3x;
      double vec3len = sqrt( vec3x*vec3x + vec3y*vec3y + vec3z*vec3z );
      double f3x = vec3x / (2.0 * d1 * halftorque);
      double f3y = vec3y / (2.0 * d1 * halftorque);
      double f3z = vec3z / (2.0 * d1 * halftorque);

      hh->m3->fx += f3x;
      hh->m3->fy += f3y;
      hh->m3->fz += f3z;

      hh->m2->fx -= f3x;
      hh->m2->fy -= f3y;
      hh->m2->fz -= f3z;
    }

    hh = hh->next;	/* Go to the next hinge in the list */
  }

} /* end of apply_hinges() */


/*	This routine will compare each mass within the list against a
 * plane and apply forces to keep the masses from penetrating the plane.
 * Forces are proportional to the distance of penetration and to the spring
 * constant (k) passed in to the function.
 *	The plane is defined by Ax+By+Cz+D = 0.  Negative values are inside
 * the plane.  The vector (A,B,C) is normal to the plane and must be of unit
 * length for k to be in the correct units.  D offsets the plane from the
 * origin along the vector (A,B,C).
 *	This routine takes into account the radius of the masses.
 *	This routine is passed a pointer to the head of the mass list.
 * Note that this is different from the pointer that is passed to the
 * add_mass() function.
 */

void apply_plane(MASS_node *mlist, double A, double B, double C, double D, double k)
{
  while (mlist != NULL) {
    // Compute signed distance into half-space cut by plane
    double dist = A*mlist->x + B*mlist->y + C*mlist->z + D - mlist->radius;

    // Apply force normal to the plane if the mass is inside the plane.
    // Apply along the vector (A,B,C).
    if (dist < 0) {
      mlist->fx -= dist*A*k;
      mlist->fy -= dist*B*k;
      mlist->fy -= dist*C*k;
    }

    mlist = mlist ->next;
  }
}


/*	This routine will compare each mass within the first list against
 * each in the second and apply forces to keep the masses from penetrating
 * each other.  Forces are proportional to the distance of penetration and
 * to the spring constant (k) passed in to the function.  Half of the force
 * is applied to each colliding mass.  The force acts along the vector
 * between the centers of the masses.
 *	This routine takes into account the radius of the masses.
 *	This routine is passed a pointer to the head of the mass lists.
 * These CANNOT both point to the same list to check for self-intersections
 * (use apply_self_collisions for this).
 * Note that this is different from the pointer that is passed to the
 * add_mass() function.
 */

void apply_collisions(MASS_node *m1, MASS_node *m2, double k)
{
  MASS_node *mlist = m2;  //< Passes multiple times through this list

  // Double-indexed loop, first by m1 mass and then by m2 mass.  This
  // checks all masses in the first list against all masses in the second.
  while (m1 != NULL) {
    mlist = m2;
    while (mlist != NULL) {
      // Compute square of distance compared to radii
      double dx = m1->x - mlist->x; // Points from mlist towards m1
      double dy = m1->y - mlist->y;
      double dz = m1->z - mlist->z;
      double radsum = m1->radius + mlist->radius;
      double dist2 = sqr(dx) + sqr(dy) + sqr(dz);

      // Apply force along vector between them, based on depth of
      // penetration, half to each mass.
      if (dist2 < sqr(radsum)) {
	double dist = sqrt(dist2);    // Only do sqrt when we have to
	double pen = radsum - dist;   // Positive number
	double force = pen * k / 2;
	double nx = dx / dist;	      // Points from mlist towards m1
	double ny = dy / dist;
	double nz = dz / dist;

	m1->fx += nx * force;
	m1->fy += ny * force;
	m1->fz += nz * force;

	mlist->fx -= nx * force;
	mlist->fy -= ny * force;
	mlist->fz -= nz * force;
      }

      mlist = mlist ->next;
    }
    m1 = m1->next;
  }
}


/*	This routine will compare each mass within the list against
 * each other mass and apply forces to keep the masses from penetrating
 * each other.  Forces are proportional to the distance of penetration and
 * to the spring constant (k) passed in to the function.  Half of the force
 * is applied to each colliding mass.  The force acts along the vector
 * between the centers of the masses.
 *	This routine takes into account the radius of the masses.
 *	Because some masses are connected together into strings using
 * springs, and because the rest length of the springs may be less than
 * the diameter of the masses, the skip parameter specifies how many
 * adjacent masses to skip when doing this test; a nonzero value for
 * this parameter only makes sense for mass lists created using the
 * make_string function.
 *	This routine is passed a pointer to the head of the mass list.
 * Note that this is different from the pointer that is passed to the
 * add_mass() function.
 */

void apply_self_collisions(MASS_node *mlist, double k, int skip)
{
  // Double-indexed loop, first by m1 mass and then by all masses after
  // the current mass.
  while (mlist != NULL) {
    MASS_node *m1 = mlist->next;
    for (int i = 0; (i < skip) && (m1 != NULL); i++) {
      m1 = m1->next;
    }
    while (m1 != NULL) {
      // Compute square of distance compared to radii
      double dx = m1->x - mlist->x; // Points from mlist towards m1
      double dy = m1->y - mlist->y;
      double dz = m1->z - mlist->z;
      double radsum = m1->radius + mlist->radius;
      double dist2 = sqr(dx) + sqr(dy) + sqr(dz);

      // Apply force along vector between them, based on depth of
      // penetration, half to each mass.
      if (dist2 < sqr(radsum)) {
	double dist = sqrt(dist2);    // Only do sqrt when we have to
	double pen = radsum - dist;   // Positive number
	double force = pen * k / 2;
	double nx = dx / dist;	      // Points from mlist towards m1
	double ny = dy / dist;
	double nz = dz / dist;

	m1->fx += nx * force;
	m1->fy += ny * force;
	m1->fz += nz * force;

	mlist->fx -= nx * force;
	mlist->fy -= ny * force;
	mlist->fz -= nz * force;
      }

      m1 = m1->next;
    }
    mlist = mlist ->next;
  }
}


/*	This routine applies the forces that are acting on each mass to the
 * masses.  It then clears the forces for the next iteration.  The position
 * and velocity for each mass is updated by one time step.  The time step is
 * assumed to be a constant times the unit time step.
 *	The acceleration that is sent as a parameter is added to the
 * acceleration that is caused by the combined force that is
 * already acting on each mass point in the list before the movement for
 * the time step is computed.  This allows for gravity, etc, to be applied
 * uniformly to all masses on the list.
 *	The dx,dy,dz parameters specify the motion of the medium in which
 * the masses reside.  This allows viscosity to push the masses along with
 * any currents.
 *	The force is assumed to be constant over the time step.
 *	Note that the routine is passed a pointer to the head of the mass
 * list, not a handle as in add_mass().
 */

void	step_masses(MASS_node *mh,
	double iax, double iay, double iaz,
	double time,
	double dx, double dy, double dz)
{
	double	ax,ay,az;	/* Acceleration in x,y, and z */
	double	vx,vy,vz;	/* Average velocity over the interval */
	double	mass;		/* Mass of a given mass point */
	double	damping;	/* Damping factor of a given mass point */

	while (mh != NULL) {	/* Perform action on each mass */

		/* Forces acting on the object
		 *         = Force there - damping * velocity
		 */

		damping = mh->damping;
		mh->fx += damping*(dx - mh->vx);
		mh->fy += damping*(dy - mh->vy);
		mh->fz += damping*(dz - mh->vz);

		/* Acceleration = Force/mass + new acceleration */

		mass = mh->mass;
		ax = (mh->fx)/mass + iax;
		ay = (mh->fy)/mass + iay;
		az = (mh->fz)/mass + iaz;

		/* Average velocity = Initial + accel*time/2 */

		vx = mh->vx + ax*time/2;
		vy = mh->vy + ay*time/2;
		vz = mh->vz + az*time/2;

		/* Final velocity = Initial + accel*time */

		mh->vx += ax*time;
		mh->vy += ay*time;
		mh->vz += az*time;

		/* New position = old + Average velocity * time */

		mh->x += time*vx;
		mh->y += time*vy;
		mh->z += time*vz;

		/* Clear forces for the next go round */

		mh->fx = mh->fy = mh->fz = 0.0;

		mh = mh->next;	/* Go to the next mass */
	}

} /* end of step_masses() */

/*	This routine will make a string of masses tied together by
 * springs and hinges.  The head mass in the string will be located at
 * the origin, and the others will run along the positive x axis, with
 * each one separated from its predecessor by the rest length of the
 * spring connecting them.
 *	It is assumed that all masses, springs, and hinges are the same
 * and that their characteristics are passed in the function call.  If
 * the spring constant is zero, then no spring list is made.  If the
 * hinge constant is zero, then no hinge list is made.  The initial forces
 * and velocity of the masses is zero.
 *	Any masses, springs, or hinges in the lists passed to this routine
 * are lost when the new string is created.  The space allocated to them is
 * not freed when this is done.  It is assumed that this routine is acting
 * on a new list.
 *	The routine returns 0 on success, -1 on failure.
 */

int	make_string(int num, MASS_node **mh, SPRING_node **sh, HINGE_node **hh,
	double mass, double damping, double rest, double springk,
	double hingek, double sbreak, double radius)
{
	register int	loop;
	int	result;		/* Holds return values of called routines */
	MASS_node	*m1,*m2;/* The last two masses in the list */

	*mh = (MASS_node *)NULL;	/* Force creation of new lists */
	*sh = (SPRING_node *)NULL;
	*hh = (HINGE_node *)NULL;

	/*   Take care of the boundry cases, where there are fewer than three
	 * masses in the string.  The rest of the masses will have joints
	 * and springs attatched to them in the while loop.
	 *   Masses are added starting at the far end of the string, so that
	 * the head will end up at the origin.
	 */

	if (num <= 0) return(0);

	if (num >= 1)
		if (result = add_mass(mh, mass, rest*(num-1),0.0,0.0, damping, radius))
			return(result);
	m1 = *mh;

	if (num >= 2) {
		if ( (num==2) && (springk==0.0) ) {
			perror("make_string (two masses and no springs)");
			return(-1);
		}
		if (result = add_mass(mh,mass, rest*(num-2),0.0,0.0, damping, radius))
			return(result);
		if ( springk != 0.0 )	/* Put in a spring if they are there */
			if ( result = add_spring(sh, rest,springk, m1,*mh, sbreak) )
				return(result);
		m2 = m1;
		m1 = *mh;
	}

	/*   At this point, there are two masses possibly tied together by
	 * a spring.  The two masses are pointed to by m1 and m2.  m1 will
	 * be between m2 and any new mass.  There are also three masses
	 * possibly tied together by a hinge.
	 *   Note that either springs or hinges may be absent, but not both.
	 * This is checked for explicitly before the loop is entered.
	 */

	if (num >= 3) {
		if ( (springk == 0.0) && (hingek == 0.0) ) {
			perror("make_string (no springs or hinges)");
			return(-1);
		}
                loop = 3;
		while (loop <= num) {
                        if ( result = add_mass(mh, mass,
			     rest*(num-loop),0.0,0.0,damping, radius) )
				return(result);
			if ( springk != 0.0 )
			  if ( result = add_spring(sh, rest,springk, m1,*mh, sbreak) )
				return(result);
			if ( hingek != 0.0 )
			  if ( result = add_hinge(hh, hingek, m2,m1,*mh) )
				return(result);
			m2 = m1;
			m1 = *mh;
			loop++;
		} /* end of while loop adding masses */

	} /* More than two masses */

	return(0);	/* The creation has been a success */

} /* end of make_string() */


/*	This routine will make a cube of masses tied together by
 * springs and hinges.  One mass in the cube will be located at
 * (x,y,z), and the others will run along the positive x,y,z axes, with
 * each one separated from its predecessor by the rest length of the
 * spring connecting them.
 *	It is assumed that all masses, springs, and hinges are the same
 * and that their characteristics are passed in the function call.  If
 * the spring constant is zero, then no spring list is made.  If the
 * hinge constant is zero, then no hinge list is made.  The initial forces
 * and velocity of the masses is zero.
 *	The routine returns 0 on success, -1 on failure.
 */

int	make_cube(MASS_node **mh, SPRING_node **sh, HINGE_node **hh,
	double x, double y, double z,
	double mass, double damping, double rest, double springk,
	double hingek,
	CUBE c)
{
	int	result;		/* Holds return values of called routines */

	//---------------------------------------------------------------
	// Add the 8 masses to the cube, also storing them into the cube
	// structure.  The zeroeth is at the origin, the first at +x, the
	// second +y, the fourth +x,+y, then the next four are the same with
	// +z on each of them.

	if (result = add_mass(mh, mass, x     ,y     ,z     , damping)) {
		return(result);
	}
	c[0] = *mh;

	if (result = add_mass(mh, mass, x+rest,y     ,z     , damping)) {
		return(result);
	}
	c[1] = *mh;

	if (result = add_mass(mh, mass, x     ,y+rest,z     , damping)) {
		return(result);
	}
	c[2] = *mh;

	if (result = add_mass(mh, mass, x+rest,y+rest,z     , damping)) {
		return(result);
	}
	c[3] = *mh;

	if (result = add_mass(mh, mass, x     ,y     ,z+rest, damping)) {
		return(result);
	}
	c[4] = *mh;

	if (result = add_mass(mh, mass, x+rest,y     ,z+rest, damping)) {
		return(result);
	}
	c[5] = *mh;

	if (result = add_mass(mh, mass, x     ,y+rest,z+rest, damping)) {
		return(result);
	}
	c[6] = *mh;

	if (result = add_mass(mh, mass, x+rest,y+rest,z+rest, damping)) {
		return(result);
	}
	c[7] = *mh;

	//---------------------------------------------------------------
	// Add springs between the corners of the cube, unless the
	// spring constant is zero.  It adds the four on the bottom, then
	// the four on the top, then the four up the sides.

	if ( springk != 0.0 ) {
	  if ( add_spring(sh, rest,springk, c[0],c[1]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[0],c[2]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[1],c[3]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[2],c[3]) ) { return(-1); }

	  if ( add_spring(sh, rest,springk, c[4],c[5]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[4],c[6]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[5],c[7]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[6],c[7]) ) { return(-1); }

	  if ( add_spring(sh, rest,springk, c[0],c[4]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[1],c[5]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[2],c[6]) ) { return(-1); }
	  if ( add_spring(sh, rest,springk, c[3],c[7]) ) { return(-1); }
	}

	return(0);	/* The creation has been a success */

} /* end of make_cube() */

/*	This routine will make a pyramid of masses tied together by
 * springs.  The masses to used are passed in to the function, which
 * attaches the springs to them.  This is because the routine is
 * intended to be connected to an existing structure, such as a cube
 * face.  Actually, it only attaches the springs between the apex
 * and the four corners -- the ones between the corners are assumed
 * to be in place already.
 *	The routine returns 0 on success, -1 on failure.
 */

int	make_pyramid(SPRING_node **sh, MASS_node *apex,
		MASS_node *b1, MASS_node *b2, MASS_node *b3, MASS_node *b4,
		double rest, double springk)
{
	// Make sure the parameters are valid
	if ( (apex == NULL) || (b1 == NULL) || (b2 == NULL) ||
	     (b3 == NULL) || (b4 == NULL) ) {
		fprintf(stderr,"make_pyramid(): NULL mass pointer\n");
		return -1;
	}
	if ( (rest <= 0) || (springk < 0) ) {
		fprintf(stderr,"make_pyramid(): Invalid parameter\n");
		return -1;
	}

	// Construct the springs.
	if ( add_spring(sh, rest,springk, apex,b1) ) { return(-1); }
	if ( add_spring(sh, rest,springk, apex,b2) ) { return(-1); }
	if ( add_spring(sh, rest,springk, apex,b3) ) { return(-1); }
	if ( add_spring(sh, rest,springk, apex,b4) ) { return(-1); }

	return 0;
}

/*	This routine will make a pyramid-capped cube of masses tied
 * together by springs and hinges.  One mass in the cube will
 * be located at (x,y,z), and the others will run along the positive
 * x,y,z axes, with each one separated from its predecessor by the
 * rest length of the spring connecting them.  (This is done using
 * the make_cube() function listed above).
 *	The make_pyramid() function is used to place caps on each
 * of the faces of the cube.
 *	The masses at the apexes of the six cubes are placed into
 * the CAPPED_CUBE structure, in the order: -x, +x, -y, +y, -z, +z.
 *	You can cause the structure to be linked into an existing
 * structure by passing in pointers to the nodes that it is to be
 * connected to.  This enables another routine to form a matrix of
 * these structures.  Passing NULL for these pointers will cause this
 * routine to create them.
 *	The routine returns 0 on success, -1 on failure.
 */

int	make_capped_cube(MASS_node **mh, SPRING_node **sh, 
		double x, double y, double z,
		double mass, double damping,
		double rest, double springk,
		MASS_node *minusx, MASS_node *minusy, MASS_node *minusz,
		CAPPED_CUBE cc, CUBE cube)
{
	// Height of the pyramid above the cube face.
	double height = rest / sqrt(2.0);

	// Make a cube
	make_cube(mh, sh, NULL, x,y,z, mass, damping, rest, springk,
		0, cube);

	// Add a pyramid to the -X face; either link it to an existing
	// node that was passed in or else create a new node.
	if (minusx == NULL) {
		add_mass(mh, mass, x-height, y+rest/2, z+rest/2, damping);
		cc[0] = *mh;
	} else {
		cc[0] = minusx;
	}
	make_pyramid(sh, cc[0], cube[0], cube[2], cube[4], cube[6],
		rest, springk);

	// Add a pyramid to the +X face
	add_mass(mh, mass, x+rest+height, y+rest/2, z+rest/2, damping);
	cc[1] = *mh;
	make_pyramid(sh, cc[1], cube[1], cube[3], cube[5], cube[7],
		rest, springk);

	// Add a pyramid to the -Y face; either link it to an existing
	// node that was passed in or else create a new node.
	if (minusy == NULL) {
		add_mass(mh, mass, x+rest/2, y-height, z+rest/2, damping);
		cc[2] = *mh;
	} else {
		cc[2] = minusy;
	}
	make_pyramid(sh, cc[2], cube[0], cube[1], cube[4], cube[5],
		rest, springk);

	// Add a pyramid to the +Y face
	add_mass(mh, mass, x+rest/2, y+rest+height, z+rest/2, damping);
	cc[3] = *mh;
	make_pyramid(sh, cc[3], cube[2], cube[3], cube[6], cube[7],
		rest, springk);

	// Add a pyramid to the -Z face; either link it to an existing
	// node that was passed in or else create a new node.
	if (minusz == NULL) {
		add_mass(mh, mass, x+rest/2, y+rest/2, z-height, damping);
		cc[4] = *mh;
	} else {
		cc[4] = minusz;
	}
	make_pyramid(sh, cc[4], cube[0], cube[1], cube[2], cube[3],
		rest, springk);

	// Add a pyramid to the +Z face
	add_mass(mh, mass, x+rest/2, y+rest/2, z+rest+height, damping);
	cc[5] = *mh;
	make_pyramid(sh, cc[5], cube[4], cube[5], cube[6], cube[7],
		rest, springk);

	return 0;
}

void  translate_masses(MASS_node *mlist, double x, double y, double z)
{
  while (mlist != NULL) {
    mlist->x += x;
    mlist->y += y;
    mlist->y += z;

    mlist = mlist ->next;
  }
}

// Helper classes and functions for the next function.
class MASS_NAME {
public:
  MASS_NAME(const char *name, MASS_node *node) {
    strncpy(d_name, name, sizeof(d_name)-1);
    d_name[sizeof(d_name)-1] = '\0';
    d_node = node;
  }
  const char *name() const { return d_name; }
  MASS_node *node() const { return d_node; }

protected:
  char        d_name[1024];
  MASS_node   *d_node;
};

// Look up the mass node with the specified name in the vector of
// mass nodes.  Return a pointer to the mass node, or NULL if there
// was not one with this name.
MASS_node *lookup_mass_node(const vector<MASS_NAME> n, const char *name)
{
  vector<MASS_NAME>::const_iterator i;
  for (i = n.begin(); i != n.end(); i++) {
    if (strcmp(name, i->name()) == NULL) {
      // Found it!
      return i->node();
    }
  }

  // Didn't find it.
  return NULL;
}

class FORCE_CURVE_NAME {
public:
  FORCE_CURVE_NAME(const char *name, GENERAL_SPRING_TABLE_ENTRY *table) {
    strncpy(d_name, name, sizeof(d_name)-1);
    d_name[sizeof(d_name)-1] = '\0';
    d_table = table;
  }
  const char *name() const { return d_name; }
  GENERAL_SPRING_TABLE_ENTRY *table() const { return d_table; }

protected:
  char        d_name[1024];
  GENERAL_SPRING_TABLE_ENTRY   *d_table;
};

// Look up the force curve with the specified name in the vector of
// curves.  Return a pointer to the table in the curve, or NULL if there
// was not one with this name.
GENERAL_SPRING_TABLE_ENTRY *lookup_force_curve(const vector<FORCE_CURVE_NAME> n, const char *name)
{
  vector<FORCE_CURVE_NAME>::const_iterator i;
  for (i = n.begin(); i != n.end(); i++) {
    if (strcmp(name, i->name()) == NULL) {
      // Found it!
      return i->table();
    }
  }

  // Didn't find it.
  return NULL;
}

// Read in a description of a single structure from a file; this includes
// reading in named masses, springs to connect masses, and hinges to connect
// triples of masses.  The format of the data in the file is:
//    structure {
//      mass_damping DAMPING  {default 5}
//      mass_radius RADIUS {default 0}
//      mass NAME MASS X Y Z DAMPING          (or)
//      mass NAME MASS X Y Z DAMPING RADIUS
//      (however many masses, each with a unique name
//      spring_constant_over_length VALUE   {Defaults to 1.0}
//      rest_length_fraction VALUE   {Defaults to 1.0, 0.5 means extended to twice rest length}
//      spring NAME1 NAME2  {Computes based on spring_constant_over_length and rest_length_fraction and mass distance} (or)
//      spring NAME1 NAME2 REST K                 (or)
//      spring NAME1 NAME2 REST K BREAKING_FORCE
//      (however many springs)
//      generic_strain_force_curve NAME {
//        strain force
//        ...
//      }
//      generic_spring NAME1 NAME2 STRAIN_FORCE_CURVE_NAME {Computes based on spring_constant_over_length and rest_length_fraction and mass distance} (or)
//      generic_spring NAME1 NAME2 STRAIN_FORCE_CURVE_NAME REST K (or)
//      generic_spring NAME1 NAME2 STRAIN_FORCE_CURVE_NAME REST K BREAKING_FORCE (or)
//      (however many generic springs)
//      hinge NAME1 NAME2 NAME3 K
//      (however many hinges)
//    }
// Returns true if the parsing went okay, leaving the file pointer at
// the beginning of the line following the closing brace.  Returns pointers
// to the described structure's masses, springs, and hinges in the handles,
// NULL for each if there were no entries.
bool    parse_structure_from_file(FILE *f, MASS_node **mh, SPRING_node **sh, GENERAL_SPRING_node **gsh, HINGE_node **hh)
{
  // Variables used to help parse lines.  The line itself is NULL-terminated in case
  // we get a line in the file that is too long and fgets() doesn't read the NULL
  // character.
  char  line[1024];
  line[sizeof(line)-1] = '\0';
  char s1[1024], s2[1024];

  // The linked lists used to build up the structure.  Also a name lists
  // that hold names and pointers to the objects to which they belong.
  
  vector<MASS_NAME> mass_names;
  MASS_node   *masses = NULL;
  SPRING_node *springs = NULL;

  vector<FORCE_CURVE_NAME> force_curve_names;
  GENERAL_SPRING_node *general_springs = NULL;
  
  HINGE_node  *hinges = NULL;

  // Default values for constants that can be overridden in the file
  double mass_radius = 0.0; //< Invisible masses by default.
  double mass_damping = 5.0;  //< Damping applied to masses that don't specify it.
  double spring_constant_over_length = 1.0; //< If spring constant calculated based on distance between masses.
  double rest_length_fraction = 1.0; //< Fraction of initial length that is rest length, if sccbodbm

  // The first line in the file must have two strings in it, "structure" and "{".
  if (fgets(line, sizeof(line)-1, f) == NULL) {
    fprintf(stderr,"parse_structure_from_file: Could not read structure line\n");
    return false;
  }
  if (sscanf(line, "%s %s", s1, s2) != 2) {
    fprintf(stderr,"parse_structure_from_file: Too short structure line\n");
    return false;
  }
  if ( (strcmp(s1, "structure") != 0) || (strcmp(s2, "{") != 0) ) {
    fprintf(stderr,"parse_structure_from_file: Bad structure line:\n(%s)", line);
    return false;
  }

  // Read lines until we get to the "}" line, handling each according to its
  // initial word.
  if (fgets(line, sizeof(line)-1, f) == NULL) {
    fprintf(stderr,"parse_structure_from_file: End of structure not found\n");
    return false;
  }
  while (line[0] != '}') {
    
    // Check the first word in the line to see what kind of line it is,
    // then parse the rest of the line based on that word and take action
    // as appropriate.
    if (sscanf(line, "%s", s1) != 1) {
      // Empty line, which is okay.
      // Get the next line in the file, then head back around the loop.
      if (fgets(line, sizeof(line)-1, f) == NULL) {
        fprintf(stderr,"parse_structure_from_file: End of structure not found\n");
        return false;
      }
      continue;

    } else if (strcmp(s1, "mass_damping") == 0) {

      // Find the value for this constant.
      if (sscanf(line, "%s %lg", s1, &mass_damping) != 2) {
        fprintf(stderr,"parse_structure_from_file: Bad mass_damping line: %s\n", line);
        return false;
      }

    } else if (strcmp(s1, "mass_radius") == 0) {

      // Find the value for this constant.
      if (sscanf(line, "%s %lg", s1, &mass_radius) != 2) {
        fprintf(stderr,"parse_structure_from_file: Bad mass_radius line: %s\n", line);
        return false;
      }

    } else if (strcmp(s1, "mass") == 0) {

      // Got a mass line.  Read its name and other parameters, then create a
      // new mass and add it to the mass list and to the mass_name list.  Remember
      // that the first string in the line will be the keyword itself, which should
      // be skipped.  Remember that radius is an optional parameter.
      char name[1024];
      double mass, x,y,z, damping;
      double radius = 0.0;
      int ret;
      if ( (ret = sscanf(line, "%s %s %lg %lg %lg %lg %lg %lg",
                 s1, name, &mass, &x, &y, &z, &damping, &radius)) < 6) {
        fprintf(stderr,"parse_structure_from_file: Bad mass line: %s\n", line);
        return false;
      }

      // If the damping wasn't specified, use the default
      if (ret <= 6) {
        damping = mass_damping;
      }

      // If the radius wasn't specified, use the default
      if (ret <= 7) {
        radius = mass_radius;
      }

      // Make sure we don't already have a mass with this name.
      if (lookup_mass_node(mass_names, name) != NULL) {
        fprintf(stderr, "parse_structure_from_file: Mass name repeated: %s\n", name);
        return false;
      }

      // Good name, go ahead and add the new mass.
      if (add_mass(&masses, mass, x,y,z, damping, radius) != 0) {
        fprintf(stderr,"parse_structure_from_file: Could not add mass\n");
        return false;
      }
      mass_names.push_back(MASS_NAME(name, masses));

    } else if (strcmp(s1, "spring_constant_over_length") == 0) {

      // Find the value for this constant.
      if (sscanf(line, "%s %lg", s1, &spring_constant_over_length) != 2) {
        fprintf(stderr,"parse_structure_from_file: Bad spring_constant_over_length line: %s\n", line);
        return false;
      }

    } else if (strcmp(s1, "rest_length_fraction") == 0) {

      // Find the value for this constant.
      if (sscanf(line, "%s %lg", s1, &rest_length_fraction) != 2) {
        fprintf(stderr,"parse_structure_from_file: Bad rest_length_fraction line: %s\n", line);
        return false;
      }

    } else if (strcmp(s1, "spring") == 0) {

      // Got a spring line.  Read its parameters, then create a
      // new spring and add it to the spring list.  Remember
      // that the first string in the line will be the keyword itself, which should
      // be skipped.  Remember that breaking_force is an optional parameter.
      char mass1[1024], mass2[1024];
      double rest_len, k, breaking_force = 1e100;
      int ret;
      if ( (ret = sscanf(line, "%s %s %s %lg %lg %lg",
                 s1, mass1, mass2, &rest_len, &k, &breaking_force)) < 3) {
        fprintf(stderr,"parse_structure_from_file: Bad spring line: %s\n", line);
        return false;
      }
      if (ret == 4) {
        fprintf(stderr,"parse_structure_from_file: Bad spring line: %s\n", line);
        return false;
      }
      
      // Locate the named masses.
      MASS_node *m1, *m2;
      m1 = lookup_mass_node(mass_names, mass1);
      m2 = lookup_mass_node(mass_names, mass2);
      if ( (m1 == NULL) || (m2 == NULL) ) {
        fprintf(stderr,"parse_structure_from_file: Could not find masses: %s %s\n", mass1, mass2);
        return false;
      }

      // Compute the spring constant and rest length based on mass distance if we weren't
      // told the distance.
      if (ret == 3) {
        rest_len = mass_distance(m1, m2) * rest_length_fraction;
        k = spring_constant_over_length / mass_distance(m1, m2);
      }

      // Add the specified spring.
      if (add_spring(&springs, rest_len, k, m1, m2, breaking_force) != 0) {
        fprintf(stderr,"parse_structure_from_file: Could not add spring\n");
        return false;
      }

    } else if (strcmp(s1, "general_strain_force_curve") == 0) {

      // Got a general strain force curve line.  Create a new list and fill
      // it with the strain and force entries found in the following lines.  Remember
      // that the first string in the line will be the keyword itself, which should
      // be skipped.
      char name[1024];
      int ret;
      if ( (ret = sscanf(line, "%s %s",
                 s1, name)) < 2) {
        fprintf(stderr,"parse_structure_from_file: Bad general_strain_force_curve line: %s\n", line);
        return false;
      }

      // Read in lines that have "strain force" on them until we get to a line
      // that just has a closing brace.  Add each to the end of the current force
      // curve table.  Tack each new entry onto the end of the current list.
      GENERAL_SPRING_TABLE_ENTRY *curve = NULL;
      GENERAL_SPRING_TABLE_ENTRY **next = &curve;
      while (true) {
        if (fgets(line, sizeof(line)-1, f) == NULL) {
          fprintf(stderr,"parse_structure_from_file: End of general_strain_force_curve not found\n");
          return false;
        }

        if (strchr(line, '}') != NULL) {
          // End of the list!
          break;
        }

        double strain, force;
        if (sscanf(line, "%lg %lg", &strain, &force) != 2) {
          fprintf(stderr,"parse_structure_from_file: Bad general_strain_force_curve entry: %s\n", line);
          return false;
        }

	if ( (*next = (GENERAL_SPRING_TABLE_ENTRY *)malloc(sizeof(GENERAL_SPRING_TABLE_ENTRY))) == NULL ) {
		perror("parse_structure_from_file (cannot allocate more memory)");
		return(false);	/* FAIL */
	}
        (*next)->force = force;
        (*next)->strain = strain;
        (*next)->next = NULL;
        next = &((*next)->next);
      }
      
      // Make sure we don't already have a curve with this name.
      if (lookup_force_curve(force_curve_names, name) != NULL) {
        fprintf(stderr, "parse_structure_from_file: Force curve name repeated: %s\n", name);
        return false;
      }

      // Good name, go ahead and add the new curve.
      force_curve_names.push_back(FORCE_CURVE_NAME(name, curve));

    } else if (strcmp(s1, "general_spring") == 0) {

      // Got a general spring line.  Read its parameters, then create a
      // new spring and add it to the spring list.  Remember
      // that the first string in the line will be the keyword itself, which should
      // be skipped.  Remember that breaking_force is an optional parameter.
      char mass1[1024], mass2[1024], force_curve[1024];
      double rest_len, k, breaking_force = 1e100;
      int ret;
      if ( (ret = sscanf(line, "%s %s %s %s %lg %lg %lg",
                 s1, mass1, mass2, force_curve, &rest_len, &k, &breaking_force)) < 4) {
        fprintf(stderr,"parse_structure_from_file: Bad general spring line: %s\n", line);
        return false;
      }
      if (ret == 5) {
        fprintf(stderr,"parse_structure_from_file: Bad general spring line: %s\n", line);
        return false;
      }
      
      // Locate the named masses.
      MASS_node *m1, *m2;
      m1 = lookup_mass_node(mass_names, mass1);
      m2 = lookup_mass_node(mass_names, mass2);
      if ( (m1 == NULL) || (m2 == NULL) ) {
        fprintf(stderr,"parse_structure_from_file: Could not find masses: %s %s\n", mass1, mass2);
        return false;
      }

      // Locate the force curve
      GENERAL_SPRING_TABLE_ENTRY *fc = lookup_force_curve(force_curve_names, force_curve);
      if (fc == NULL) {
        fprintf(stderr,"parse_structure_from_file: Could not find force curve: %s\n", force_curve);
        return false;
      }

      // Compute the rest length based on distance if we weren't told the distance.
      // Set the spring constant to the spring constant per unit distance because the
      // strain/force curve already takes into account the distance division, it is
      // an intrinsic property of the material.
      if (ret == 4) {
        rest_len = mass_distance(m1, m2) * rest_length_fraction;
        k = spring_constant_over_length;
      }

      // Add the specified spring.
      if (add_general_spring(&general_springs, rest_len, k, fc, m1, m2, breaking_force) != 0) {
        fprintf(stderr,"parse_structure_from_file: Could not add general spring\n");
        return false;
      }

    } else {
      // Didn't recognize the keyword starting this line.
      fprintf(stderr,"parse_structure_from_file: Unrecognized keyword on line: %s\n", line);
      return false;
    }

    // Get the next line in the file.
    if (fgets(line, sizeof(line)-1, f) == NULL) {
      fprintf(stderr,"parse_structure_from_file: End of structure not found\n");
      return false;
    }
  }

  // Put the array values into the return parameters.
  *mh = masses;
  *sh = springs;
  *gsh = general_springs;
  *hh = hinges;
  return true;
}
