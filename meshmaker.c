//-------------------------------------------------------------------------
// This program checks the rigidity of various structures by building them
// as mass-spring systems and seeing if they stay rigid under various
// forces and gravity.

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// Global variables describing default values.
double	g_mass_damping = 1.0;
double	g_mass_radius = 2.0;
double	g_mass = 1.0;
double	g_edge_mass_radius = 5.0;
double	g_edge_mass = 1.0e6;
double	g_spring_constant_over_length = 100.0;
double	g_rest_length_fraction = 0.9;

// Number of potential cross-link locations in X and Y
// Location to store an array of NX x NY crosslink values.
unsigned g_NX = 20;
unsigned g_NY = 16;
int *g_crosslinks = NULL;
double	g_fraction = 1.0;

// Functions to specify the X and Y coordinates given a mesh
// index.
double X(int i)
{
  return (i-10)*15;
}

double Y(int j)
{
  return (j-8)*15;
}

void Usage(const char *s)
{
  fprintf(stderr,"Usage: %s [-frac F] [outfile_name]\n", s);
  fprintf(stderr,"       outfile_name: (default stdout)\n");
  fprintf(stderr,"       -frac: Fraction of potential links made (default 1.0)\n");
  exit(-1);
}

int main(unsigned argc, const char *argv[])
{
  // File name to write, defaults to empty
  const char  *outfile_name = "";

  //-------------------------------------------------------------
  // Command-line parsing
  unsigned  i;
  unsigned  real_params = 0;  //< Number of non-flag parameters
  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i],"-help")) {
      Usage(argv[0]);
    } else if (!strcmp(argv[i],"-frac")) {
      if (++i >= argc) { Usage(argv[0]); }
      g_fraction = atof(argv[i]);
    } else if (argv[i][0] == '-') {
      Usage(argv[0]);
    } else {
      switch (real_params) {
      case 0:
        outfile_name = argv[i];
        break;
      default:
        Usage(argv[0]);
      }
    }
  }

  //-------------------------------------------------------------
  // Configure the output file.  Use stdout if the name is blank, or open
  // the file if it is not.

  FILE *f = stdout;
  if (strlen(outfile_name) > 0) {
	f = fopen(outfile_name, "w");
  }
  if (f == NULL) {
    char s[1024];
    sprintf(s, "Could not open output file '%s' for writing", outfile_name);
    perror(s);
    return(-3);
  }

  //-------------------------------------------------------------
  // Generate an NxM array to keep track of which entries have crosslinks.
  // Fill in its entries with 0 or 1 for whether it has a link or not.
  // The border entries always have crosslinks.
  g_crosslinks = malloc(g_NX * g_NY * sizeof(int));
  if (g_crosslinks == NULL) {
	fprintf(stderr,"Out of memory\n");
	return -4;
  }
  unsigned j;
  for (i = 0; i < g_NX; i++) {
    for (j = 0; j < g_NY; j++) {
      if ( (i == 0) || (j == 0) || (i == g_NX-1) || (j == g_NY-1) ) {
	g_crosslinks[i + j*g_NX] = 1;
      } else {
	if ( (random() % 1000) <= (g_fraction * 1000) ) {
		g_crosslinks[i + j*g_NX] = 1;
	} else {
		g_crosslinks[i + j*g_NX] = 0;
	}
      }
    }
  }

  //-------------------------------------------------------------
  // Generate the output describing the mesh structure.

  // Open the structure and set default parameters.
  fprintf(f, "structure {\n");
  fprintf(f, "  mass_damping %lf\n", g_mass_damping);
  fprintf(f, "  spring_constant_over_length %lf\n", g_spring_constant_over_length);
  fprintf(f, "  rest_length_fraction %lf\n", g_rest_length_fraction);
  fprintf(f, "\n");

  // Put in the masses.  All masses are named by their coordinates,
  // format XxY, where X and Y are the coordinates.
  // We only put in a crosslinked mass if it has a nonzero value in
  // its entry.  All of the edge masses are present and have large
  // masses.
  for (i = 0; i < g_NX; i++) {
    for (j = 0; j < g_NY; j++) {
      if ( (i == 0) || (j == 0) || (i == g_NX-1) || (j == g_NY-1) ) {
	fprintf(f, "  mass_radius %lf\n", g_edge_mass_radius);
	fprintf(f, "  mass\t%dx%d\t%lf\t%lf\t%lf\t%lf\n", i, j, g_edge_mass,
			X(i), Y(j), 0.0);
      } else if (g_crosslinks[i + j*g_NX] == 1) {
	fprintf(f, "  mass_radius %lf\n", g_mass_radius);
	fprintf(f, "  mass\t%dx%d\t%lf\t%lf\t%lf\t%lf\n", i, j, g_mass,
			X(i), Y(j), 0.0);
      }
    }
  }
  fprintf(f, "\n");

  // Put in the springs.  We add springs from the left to the right and from
  // top to bottom.  If we're at a node that has no entry, we don't add a spring.
  // If we're at a node that has an entry, we look to the right (or down) until
  // we find another node that does have an entry.  When we find it, we attach
  // a spring between the nodes.  We don't need to worry about whether we are
  // an edge node for this.
  unsigned k;
  for (i = 0; i < g_NX; i++) {
    for (j = 0; j < g_NY; j++) {
      // Only add links for nodes that exist.
      if (g_crosslinks[i + j*g_NX] == 1) {

	// Find and add the spring going in X
	// There is guaranteed to be one.
	for (k = i+1; k < g_NX; k++) {
		if (g_crosslinks[k + j*g_NX] == 1) {
			fprintf(f, "  spring %dx%d %dx%d\n", i,j, k,j);
			break;
		}
	}

	// Find and add the spring going in Y
	// There is guaranteed to be one.
	for (k = j+1; k < g_NY; k++) {
		if (g_crosslinks[i + k*g_NX] == 1) {
			fprintf(f, "  spring %dx%d %dx%d\n", i,j, i,k);
			break;
		}
	}

      }
    }
  }
  fprintf(f, "\n");

  // End of structure.
  fprintf(f, "\n");
  fprintf(f, "}\n");

  //-------------------------------------------------------------
  // Close the output file and return.
  free(g_crosslinks);
  fclose(f);
  return 0;
}

