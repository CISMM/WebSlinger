/************************************************************************
 *	This file contains constants and type definitions required by	*
 * the massmesh.c file and all other files using its routines.		*
 ************************************************************************/

#ifndef	MASSMESH_H
#define MASSMESH_H

#ifndef	NULL
#define	NULL	0
#endif

#include  <stdio.h>

/*	The masses are stored in a linked list.  The last mass has a null
 * next pointer.  Each mass points stores its mass, its location, its
 * velocity when last computed, the damping factor (the percent of velocity
 * countered during one time step due to friction, etc), and a pointer to the
 * next mass in the list.  Each point also includes an area where forces on
 * the mass point can be accumulated by the objects (springs) acting on this
 * point.  This allows an arbitrary number of springs to be attatched to one
 * point and they will have their forces summed appropriately.
 *	The ordering of the masses in the list is not important, as any
 * spring attached to a mass point will have a pointer that will allow it
 * to reference the point.  Note that once the list is created and springs
 * have been attached, no mass point can be deleted or moved if it has a
 * spring pointing to it, or the spring will act incorrectly.
 */

typedef	struct	MASS {	/* Used to allow a pointer to this type in struct */
	  double mass;		/* The mass of this point */
	  double radius;	/* The radius of the mass for collision purposes */
	  double x,y,z;		/* The coordinates in space of this point */
	  double vx,vy,vz;	/* The velocity of this point */
	  double fx,fy,fz;	/* The forces acting on this point */
	  double damping;	/* The damping of movement of this point */
	  struct MASS *next;	/* The next mass in the list */
	} MASS_node;


/*	This struct defines a node in the list of springs.  A spring goes
 * between two masses and applies forces to them based on the distance and
 * the spring constant of the spring.  If the distance between the masses
 * is equal to the rest length of the spring, then there is no force generated
 * by the spring.  If the distance is too short, then a force that tends to
 * separate the masses is generated.  If the distance is too long, then a
 * force pulling them together will be generated.
 *	The order of springs in the list does not matter.  The last node has
 * a null next pointer.
 */

typedef struct	SPRING { /* Used to allow a pointer to this type in struct */
	  double rest;		/* The rest length of the spring */
	  double k;		/* The spring constant of the spring */
	  double fbreak;	/* The force which will break the spring */
	  MASS_node *m1;	/* One mass attatched to the spring */
	  MASS_node *m2;	/* The other mass attatched to the spring */
	  struct SPRING	*next;	/* The next spring in the list */
	} SPRING_node;


/*      This structure defines a generalized spring node.  This is like a
 * spring node, but rather than being linear, its force is based on a table
 * look-up based on the strain of the spring and its rest length.  The first
 * structure is one that is used to make the table in the second one.
 */

typedef struct  GENERAL_SPRING_ENTRY { /* Used to allow a pointer to this type in struct */
          double strain;    /* Strain (change in length over rest length) */
          double force;     /* Force applied at the above strain */
          struct GENERAL_SPRING_ENTRY *next;  /* The next entry in the table */
        } GENERAL_SPRING_TABLE_ENTRY;

typedef struct  GENERAL_SPRING {  /* Used to allow a pointer to this type in struct */
	  double rest;		/* The rest length of the spring */
	  double k;		/* The spring constant of the spring */
          GENERAL_SPRING_TABLE_ENTRY  *table; /* Pointer to the table of force/strain */
	  double fbreak;	/* The force which will break the spring */
	  MASS_node *m1;	/* One mass attatched to the spring */
	  MASS_node *m2;	/* The other mass attatched to the spring */
	  struct GENERAL_SPRING	*next;	/* The next general spring in the list */
	} GENERAL_SPRING_node;

/*	This structure defines a hinge node.  A hinge relates three masses
 * and tends to keep the second on a line between the first and third.  A
 * hinge will apply torque to the arms between the masses to make this so.
 * The hinge constant k, multiplied by the angle between the arms, determines
 * the strength of this torque.  The torque is applied both to the center
 * node and to the end nodes, being converted to force by the arm lengths.
 *	The order of hinges in the hinge list does not matter.  The last node
 * has a null next pointer.
 */

typedef	struct	HINGE {	/* Used to allow a pointer to this type in struct */
	 double	k;		/* The hinge constant of this hinge */
	 MASS_node *m1,*m2,*m3;	/* The masses this hinge is connected to */
	 struct HINGE *next;	/* The next hinge in the list */
        } HINGE_node;

/*	The return types of all functions in massmesh.c that are made to be
 * used by other routines can be found here.  The first four return -1 on
 * failure.
 */

extern	int	add_mass(MASS_node **headptr,
			double mass,
			double x, double y, double z,
			double damping,
			double radius = 0.0);

extern	int	add_spring(SPRING_node **headptr,
			double rest,
			double k,
			MASS_node *m1, MASS_node *m2,
			double fbreak = 1e100);

extern	int	add_general_spring(GENERAL_SPRING_node **headptr,
			double rest,
			double k,
                        GENERAL_SPRING_TABLE_ENTRY  *table,
			MASS_node *m1, MASS_node *m2,
			double fbreak = 1e100);

extern	int	add_hinge(HINGE_node **headptr,
			double k,
			MASS_node *m1, MASS_node *m2, MASS_node *m3);

extern	void	clear_forces(MASS_node *mlist);
extern	int	apply_springs(SPRING_node **sh);
extern	int	apply_general_springs(GENERAL_SPRING_node **sh);
extern	void	apply_hinges(HINGE_node *hh);
extern	void	apply_plane(MASS_node *mlist, double A, double B, double C, double D, double k);
extern	void	apply_collisions(MASS_node *m1, MASS_node *m2, double k);
extern	void	apply_self_collisions(MASS_node *mlist, double k, int skip = 0);

// Return the distance between the two mass nodes.
inline double mass_distance(const MASS_node *m1, const MASS_node *m2)
{
  if ( (m1 == NULL) || (m2 == NULL) ) { return 0.0; }
  double dx = m1->x - m2->x;
  double dy = m1->y - m2->y;
  double dz = m1->z - m2->z;
  return sqrt( dx*dx + dy*dy + dz*dz );
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
extern	void	step_masses(MASS_node *mh,
			double iax, double iay, double iaz,
			double time,
			double dx, double dy, double dz);

//-----------------------------------------------------------------------
// Helper functions to build structures.

extern	int	make_string(int num,
			MASS_node **mh, SPRING_node **sh, HINGE_node **hh,
			double mass, double damping,
			double rest, double springk,
			double hingek, double sbreak = 1e100, double radius = 0.0);

typedef	MASS_node *CUBE[8]; 
extern	int	make_cube(MASS_node **mh, SPRING_node **sh, HINGE_node **hh,
			double x, double y, double z,
			double mass, double damping,
			double rest, double sprink,
			double hingek, CUBE c);

extern	int	make_pyramid(SPRING_node **sh, MASS_node *apex,
			MASS_node *b1, MASS_node *b2,
			MASS_node *b3, MASS_node *b4,
			double rest, double springk);

typedef	MASS_node *CAPPED_CUBE[6];
extern	int	make_capped_cube(MASS_node **mh, SPRING_node **sh, 
			double x, double y, double z,
			double mass, double damping,
			double rest, double springk,
			MASS_node *minusx, MASS_node *minusy, MASS_node *minusz,
			CAPPED_CUBE cc, CUBE c);

// Read in a description of a single structure from a file; this includes
// reading in named masses, springs to connect masses, and hinges to connect
// triples of masses.  The format of the data in the file is:
//    structure {
//      mass_damping DAMPING  {default 5}
//      mass_radius RADIUS {default 0}
//      mass NAME MASS X Y Z                  (or)
//      mass NAME MASS X Y Z DAMPING          (or)
//      mass NAME MASS X Y Z DAMPING RADIUS
//      (however many masses, each with a unique name
//      spring_constant_over_length VALUE   {Defaults to 1.0}
//      rest_length_fraction VALUE   {Defaults to 1.0, 0.5 means extended to twice rest length}
//      spring NAME1 NAME2  {Computes based on spring_constant_over_length and rest_length_fraction and mass distance} (or)
//      spring NAME1 NAME2 REST K                 (or)
//      spring NAME1 NAME2 REST K BREAKING_FORCE
//      (however many springs)
//      general_strain_force_curve NAME {
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
extern  bool    parse_structure_from_file(FILE *f, MASS_node **mh, SPRING_node **sh, GENERAL_SPRING_node **gsh, HINGE_node **hh);

//-----------------------------------------------------------------------
// Functions to move structures.

extern	void  translate_masses(MASS_node *mlist, double x, double y, double z);

#endif
