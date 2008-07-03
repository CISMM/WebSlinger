/*	This program will draw a specified line on the screen.  The origin
  of the coordinate system is taken to be the center of the screen.  The line
  will not be clipped if it goes beyond the boundaries of the screen, but it
  will be scissored by the display routine. */

#include <stdio.h>
#include <graphics.h>
#include <math.h>

#define	PI	(3.14159)

main()
{
	void	*(drawfunc());	/* What drawing routine to use */
	int	x0,y0;		/* The first point on the line */
	int	x1,y1;		/* The last point on the line */
	float	rho,theta;
	static	int	graphdriver = DETECT;
	static	int	graphmode = DETECT;
	extern	void	Put_Point();
	short	loop;

/*	printf("Enter the coordinates of the first point (x y): ");
	scanf("%d %d", &x0,&y0);
	printf("Enter the coordinates of the second point (x y): ");
	scanf("%d %d", &x1, &y1);
*/
	initgraph(&graphdriver,&graphmode,"\\tc");
	graphdefaults();
	cleardevice();

	rho = 150.0;
	for (theta = 0.0; theta <= 2*PI; theta += PI/100) {
		x1 = (int)(rho*cos(theta));
		y1 = (int)(rho*sin(theta));
		Draw_Line(0,0, x1,y1, WHITE);
		Draw_Line(x1,y1, 0,0, RED);
	}

/*	for (loop = 0; loop < 200; loop++) {
		Draw_Line(0,loop, 150,loop,WHITE);
		Draw_Line(150,loop, 0,loop, RED);
	}
*/
	if (getch()) {}
	closegraph();

} /* end of main */


/*	This routine will draw an arbitrary line on a screen.  It is passed
  the endpoints of the line and its color.  It determines the correct
  parameters to call the midpoint() routine with and then calls it. */

Draw_Line(x0,y0, x1,y1, color)
	int	x0,y0;		/* The first end point */
	int	x1,y1;		/* The second end point */
	int	color;		/* The color for the line */
{
	int	cx, cy;		/* Change in x and y from p0 to p1 */
	int	a,b;		/* Used in finding changes in d */
	extern	void	Put_Point();
	extern	void	Put_Backwards();

	/* Determine which octant the line is in and find the parameters
	   to pass to the midpoint() routine.  First find the quadrant
	   and then find the octant within the quadrant. */

	cx = (x1 - x0);
	cy = (y1 - y0);
	a = cy;
	b = -cx;
	if ( (cx > 0) && (cy >= 0) ) {	/* first quadrant */
		if (cx > cy) { /* first octant */
			midpoint(x0,y0,x1, 1,0,1, 2*a+b,2*a,2*(a+b),
                                Put_Point,color);
		}
		else { /* second octant */
			midpoint(y0,x0,y1, 1,0,1, -(2*b+a),-(2*b),-(2*(b+a)),
				Put_Backwards,color);
		}
	}
	else if ( (cx <= 0) && (cy > 0) ) { /* second quadrant */
		if (-cx < cy) { /* third octant */
			midpoint(y0,x0,y1, 1,-1,0, -(2*b-a),-(2*(b-a)),-(2*b),
				Put_Backwards,color);
		}
		else { /* fourth octant */
			midpoint(x0,y0,x1, -1,1,0, -(2*a-b),-(2*(a-b)),-(2*a),
				Put_Point,color);
		}
	}
	else if ( (cx < 0) && (cy <= 0) ) { /* third quadrant */
		if (-cx > -cy) { /* fifth octant */
			midpoint(x0,y0,x1, -1,-1,0, (2*a+b),(2*(a+b)),(2*a),
				Put_Point,color);
		}
		else { /* sixth octant */
			midpoint(y0,x0,y1, -1,-1,0, -(2*b+a),-(2*(b+a)),-(2*b),
				Put_Backwards,color);
		}
	}
	else if ( (cx >= 0) && (cy < 0) ) { /* fourth quadrant */
		if (cx < -cy) { /* seventh quadrant */
			midpoint(y0,x0,y1, -1,0,1, -2*b+a,-2*b,2*(-b+a),
				Put_Backwards,color);
		}
		else { /* eight quadrant */
			midpoint(x0,y0,x1, 1,0,-1, -2*a+b,-2*a,2*(-a+b),
				Put_Point,color);
		}
	}
	else { /* Just a point */
                Put_Point(x0,y0, WHITE);
	}

} /* end of Draw_Line() */


/*	This routine will draw a line on the screen.  It uses the midline
  scan conversion algorithm to do this.  The values of X and Y are general
  in that the variables can refer to either the horizontal or vertical
  direction.  The ambiguity is resolved by which function is passed to this
  routine to draw the points. */

midpoint(x0, y0, xm, dx,dy1,dy2, d,d1,d2, plot, color)
	int	x0,y0;		/* The starting point */
	int	xm;		/* The maximum value for x */
	int	dx;		/* The change to be applied to x each step */
	int	dy1,dy2;	/* delta E and delta NE in the first octant */
	int	d;		/* Value of the decision variable at start */
	int	d1,d2;		/* Change in d in E and NE */
	void	*(plot());	/* The function to use to draw the point */
	int	color;		/* Color for the line */
{
	int	x,y;	/* Tracks along x and y */

	x = x0; y = y0;
	plot(x,y, color);

	while (x != xm) {
		if (d<=0) {	/* East in the 1st octant */
			x+= dx;
			y+= dy1;
			d+= d1;
		}
		else {
			x+= dx;
			y+= dy2;
			d+= d2;
		}
		plot(x,y, color);
	}

} /* end of midpoint() */


/*	This routine will plot a point on the screen in the color that is
  specified.  The origin is the center of the screen.  X goes to the right
  and Y goes up.
*/

void	Put_Point(x,y, c)
	int	x,y;
	int	c;
{
	putpixel( getmaxx()/2 + x, getmaxy()/2 - y, c);

} /* end of Put_Point() */


/*	This routine draw the point at (x,y), but the parameters are
  reversed so that y is passed first. */

void	Put_Backwards(y,x, c)
	int	x,y;
	int	c;
{
	putpixel( getmaxx()/2 + x, getmaxy()/2 - y, c);

} /* end of Put_Backwards() */
