#ifndef	GRAPHICS_H
#define GRAPHICS_H

#ifdef _WIN32
#include <windows.h>
#endif

#include "massmesh.h"

extern int getmaxx(void);
extern int getmaxy(void);

extern int draw_masses(const MASS_node *mn);
extern int draw_springs(const SPRING_node *sn);
extern int draw_general_springs(const GENERAL_SPRING_node *sn);

typedef void (*Display_Func)(void);

extern int init_graphics(const char *name, Display_Func df);
extern int setup_graphics_frame(void);

extern double mouse_x_in_world(int mouse_raw_x);
extern double mouse_y_in_world(int mouse_raw_y);

#endif
