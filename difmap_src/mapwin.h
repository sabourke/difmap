/*
 * Define a struct to record rectangular windows pertanining to maps
 * and also a container for a linked list of such windows.
 */

#ifndef mapwin_h
#define mapwin_h

#include "mapmem.h"

typedef struct subwin {
  float xmin, xmax;      /* Window X-axis limits (radians) */
  float ymin, ymax;      /* Window Y-axis limits (radians) */
  struct subwin *next;
} Subwin;

typedef struct {
  int nwin;         /* The number of windows in the linked list */
  Subwin *head;     /* Head of the linked list */
  Subwin *tail;     /* Tail of the linked list */
} Mapwin;

typedef struct {    /* For recording pixel limits of a window */
  int xa, xb;
  int ya, yb;
} Winran;

Mapwin *new_Mapwin(void);          /* Create a new Mapwin container */
Mapwin *del_Mapwin(Mapwin *mwin);  /* Delete container and constituents */

/*
 * Add a window to the linked list in a given Mapwin container.
 */
Subwin *add_win(Mapwin *mwin, float xmin, float xmax, float ymin, float ymax);

/*
 * See if a given coordinate is in one or more windows.
 */
int inmapwin(Mapwin *mwin, float xpos, float ypos);

Subwin *del_win(Subwin *win);     /* Delete a Subwin structure */

Subwin *rem_win(Mapwin *mwin, Subwin *prev, Subwin *win);

int win_pix(Subwin *win, MapBeam *mb, int ixmin, int ixmax, int iymin,
	    int iymax, Winran *wr);

void winstats(MapBeam *mb, int domap, Subwin *win, int doall);

void shiftwin(Mapwin *mwin, float east, float north);

int peakwin(MapBeam *mb, Mapwin *mw, float size, int doabs);

int rwins(Mapwin *mwin, char *filename, float east, float north);

int wwins(Mapwin *mwin, char *filename, float east, float north, int do_old);

#endif

