#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "logio.h"
#include "mapwin.h"
#include "units.h"
#include "vlbmath.h"
#include "vlbconst.h"
#include "ellips.h"

#define WIN_INC 20
typedef struct winm {
  Subwin wins[WIN_INC];
  struct winm *next;
} Winmem;

static Winmem *new_Winmem(void);
static Subwin *new_win(void);

static struct {
  Subwin *free;
  int nused;
  Winmem *memhead;
  Winmem *memtail;
} winmem = {0,0,0,0};

/*.......................................................................
 * Allocate and initialize a new rectangular window memory block and
 * append to the linked list of such in winmem.
 *
 * Inputs: NONE
 * Side efffects:
 *    modifies winmem in this file.
 * Output:
 *   return  Winmem *  The newly allocated block, or NULL on error.
 */
static Winmem *new_Winmem(void)
{
  Winmem *newmem;   /* A temporary pointer to the new memory block */
  int i;
/*
 * Attempt to allocate the new block.
 */
  newmem = (Winmem *) malloc(sizeof(Winmem));
  if(newmem == NULL) {
    lprintf(stderr, "Insufficient memory for new block of windows\n");
    return newmem;
  };
/*
 * This function should only be called when the free-list is exhausted, so
 * the first element of the new block will now become the head of the free
 * list.
 */
  winmem.free = &newmem->wins[0];
/*
 * Link this to all the other elements of the block to complete
 * the free-list.
 */
  for(i=0; i<WIN_INC-1; i++)
    newmem->wins[i].next = &newmem->wins[i+1];
/*
 * Terminate the free-list at the last element of newmem.
 */
  newmem->wins[WIN_INC-1].next = NULL;
/*
 * Install the new block as the tail of the list of such blocks
 * in 'winmem'. If this is the first entry in the list then
 * also initialize winmem.memhead.
 */
  if(winmem.memhead == NULL) {
    winmem.memhead = winmem.memtail = newmem;
    newmem->next = NULL;
  }
  else {
    winmem.memtail->next = newmem;
    winmem.memtail = newmem;
    newmem->next = NULL;
  };
/*
 * Return the new block.
 */
  return newmem;
}

/*.......................................................................
 * Get a new rectangular window structure from the free-list. Also
 * initialise it so that all members are zero, including the 'next' link.
 *
 * Output:
 *   return  Subwin *  The new window, or NULL on error.
 */
static Subwin *new_win(void)
{
  static Subwin *newwin;   /* The new sub-window */
/*
 * Get a new memory block of windows if the free-list is empty.
 */
  if(winmem.free == NULL && new_Winmem() == NULL) {
    lprintf(stderr, "Insufficient memory for a new window\n");
    return NULL;
  };
/*
 * Increment the count of subwindows used.
 */
  winmem.nused++;
/*
 * Extract the new subwindow from the free-list and fix up the
 * the free-list.
 */
  newwin = winmem.free;
  winmem.free = winmem.free->next;
/*
 * Initialise the sub-window's members before returning it.
 */
  newwin->next=NULL;
  newwin->xmin = 0.0f;
  newwin->xmax = 0.0f;
  newwin->ymin = 0.0f;
  newwin->ymax = 0.0f;
  return newwin;
}

/*.......................................................................
 * Place a sub-window back on the free-list. NOTE. The caller must
 * relink around this component (use rem_win()) before sending it to
 * this routine - otherwise the caller's window list will thereafter dive
 * into the free-list via the modified ->next member of this component.
 * eg. win = del_win(rem_win(mwin,prev,win));  Where prev->next=win .
 *
 * Input:
 *   win  Subwin *  The unlinked window to be returned to the free-list.
 * Output:
 *   Subwin * 0     Allways.
 */
Subwin *del_win(Subwin *win)
{
  win->next = winmem.free;
  winmem.free = win;
  winmem.nused--;
  return NULL;
}

/*.......................................................................
 * Dynamically allocate an empty sub-window container class
 * and return its pointer on success, or NULL on failure.
 *
 * Output:
 *   return   Mapwin *   The initialised window-list container or NULL
 *                       on error.
 */
Mapwin *new_Mapwin(void)
{
  Mapwin *mwin;
/*
 * Allocate the container class.
 */
  mwin = (Mapwin *) malloc(sizeof(Mapwin));
  if(mwin == NULL) {
    lprintf(stderr, "Insufficient memory to allocate new Mapwin\n");
    return mwin;
  };
/*
 * Initially the container is empty.
 */
  mwin->nwin = 0;
  mwin->head  = 0;
  mwin->tail = 0;
/*
 * return the initialized container.
 */
  return mwin;
}

/*.......................................................................
 * Delete a Mapwin container and its contents. This involves free()ing
 * all memory allocated and returning components to the free-list.
 *
 * Input:
 *   mwin  Mapwin *   The Mapwin to be deleted.
 * Output:
 *   return Mapwin * The deleted Mapwin (ie. a NULL Mapwin pointer!).
 */
Mapwin *del_Mapwin(Mapwin *mwin)
{
  Subwin *win;   /* The current window component being free'd. */
  Subwin *next;  /* The window following 'win' */
/*
 * Do nothing if the window-list container has already been deleted.
 */
  if(mwin == NULL)
    return mwin;
/*
 * Return all sub-windows to the free-list.
 */
  for(win=mwin->head; win != NULL; win=next) {
    next = win->next;
    del_win(win);
  };
/*
 * Finally, free the container itself.
 */
  free(mwin);
  return (Mapwin *) 0;
}

/*.......................................................................
 * Append a new subwindow to a given window container. Also, if
 * xmin is > xmax then the two values will be swapped. The same goes for
 * ymin and ymax. Coordinates are in radians.
 *
 * Input:
 *   mwin  Mapwin *   The Mapwin container in which the component is to
 *                    be added.
 *   float   xmin     X-coord of left edge of the window.
 *   float   xmax     X-coord of right edge of the window.
 *   float   ymin     Y-coord of bottom edge of the window.
 *   float   ymax     Y-coord of top edge of the window.
 * Output:
 *   return Subwin *  The subwindow added or NULL on error.
 */
Subwin *add_win(Mapwin *mwin, float xmin, float xmax, float ymin, float ymax)
{
  static Subwin *win;  /* The new component */
/*
 * Get a new component from the free list.
 */
  win = new_win();
  if(win == NULL)
    return win;
/*
 * Fill in the window bounds attributes.
 */
  if(xmin < xmax) {
    win->xmax = xmax;
    win->xmin = xmin;
  } else {
    win->xmax = xmin;
    win->xmin = xmax;
  };
  if(ymin < ymax) {
    win->ymax = ymax;
    win->ymin = ymin;
  } else {
    win->ymax = ymin;
    win->ymin = ymax;
  };
/*
 * Insert the new subwindow and return it.
 */
  if(mwin->head == NULL) {
    mwin->head=mwin->tail=win;
  }
  else {
    mwin->tail->next = win;
    mwin->tail = win;
  };
  mwin->nwin++;
  return win;
}

/*.......................................................................
 * Return 1 if the pixel number is within any of the sub-windows of the
 * Mapwin list. Return 0 otherwise.
 *
 * Input:
 *  mwin Mapwin *  The container of the windows to be searched.
 *  xpos float     The X-coordinate of the pixel.
 *  ypos float     The Y-coordinate of the pixel.
 * Output:
 *  return int     1 - if the pixel is in one or more windows.
 *                 0 - if the pixel is outside all windows.
 */
int inmapwin(Mapwin *mwin, float xpos, float ypos)
{
  static Subwin *win; /* The window being checked */
  for(win=mwin->head; win != NULL; win=win->next)
    if(xpos >= win->xmin && xpos <= win->xmax &&
       ypos >= win->ymin && ypos <= win->ymax)
      return 1;
  return 0;
}

/*.......................................................................
 * Shift all the windows in a window list by a given distance.
 *
 * Input/Output:
 *  mwin  Mapwin *  The list of windows to be shifted (can be NULL). 
 * Input:
 *  east  float     Window eastward shift (radians).
 *  north float     Window northward shift (radians).
 */
void shiftwin(Mapwin *mwin, float east, float north)
{
  Subwin *win; /* The window being shifted */
/*
 * No change if window list is NULL.
 */
  if(mwin==NULL)
    return;
  for(win=mwin->head; win != NULL; win=win->next) {
    win->xmin += east;
    win->xmax += east;
    win->ymin += north;
    win->ymax += north;
  };
  return;
}

/*.......................................................................
 * Remove a window from a window list and return it. In order to unlink
 * the window from the list, the function needs to know what the previous
 * window in the list is. If you don't know this set prev=NULL and a
 * linear search will be made through the list for the previous window.
 * Needless to say, if you are iterating over the window list it is
 * faster to keep track of the previous unremoved window than to get
 * this function to search for it. You must set prev=NULL for the first
 * window since there is no preceding window in the list.
 *
 * Input:
 *  mwin   Mapwin *  The list containing the window to be removed.
 *  prev   Subwin *  The window in the list preceding 'win', or
 *                   NULL either if unknown or win is the head of
 *                   the list.
 *  win    Subwin *  The window to be removed.
 * Output:
 *  return Subwin *  The extracted window. Apply del_win() to this
 *                   if it is no longer required.
 */
Subwin *rem_win(Mapwin *mwin, Subwin *prev, Subwin *win)
{
  static Subwin *tmpwin;
  if(win==NULL)
    return win;
/*
 * If prev is NULL search for 'win' and the window preceding it.
 */
  if(prev==NULL) {
    for(tmpwin=mwin->head; tmpwin!=win && tmpwin!=NULL; tmpwin=tmpwin->next)
      prev = tmpwin;
    if(tmpwin==NULL) {
      lprintf(stderr, "rem_win: programmer error: window not found\n");
      return win;
    };
  };
/*
 * If 'prev' is still NULL then 'win' is the first window of the list.
 */
  if(prev==NULL)
    mwin->head = win->next;
  else
    prev->next = win->next;
/*
 * Fix up the tail if win was the tail of the list.
 */
  if(win->next==NULL)
    mwin->tail = prev;
/*
 * Record the shortening of the list.
 */
  mwin->nwin--;
  return win;
}

/*.......................................................................
 * Take a window and a map-beam container and return the pixel limits
 * of the window within a given area of the map or beam grid.
 * The returned limits are adjusted to lie within the proscribed area
 * unless all window edges lie outside the area, in which case -1 is
 * returned.
 *
 * Input:
 *  win    Subwin  *  A single Subwin window element.
 *  mb     MapBeam *  The Map within which the window is to be defined.
 *  ixmin      int    Element number of left edge of grid area.
 *  ixmax      int    Element number of right edge of grid area.
 *  iymin      int    Element number of lower edge of grid area.
 *  iymax      int    Element number of upper edge of grid area.
 * Output:
 *  wr      Winran *  The Window range element to be filled.
 *  return int         0 - Window was OK.
 *                     1 - Window was wholly outside the CLEAN area.
 */
int win_pix(Subwin *win, MapBeam *mb, int ixmin, int ixmax, int iymin,
	    int iymax, Winran *wr)
{
  static int xa;    /* Working copy of min element number along X */
  static int xb;    /* Working copy of max element number along X */
  static int ya;    /* Working copy of min element number along Y */
  static int yb;    /* Working copy of max element number along Y */
  static int xcent; /* X-axis pixel of origin of map (nx/2) */
  static int ycent; /* Y-axis pixel of origin of map (ny/2) */
  static float wxa; /* Float pixel dist from map center to window X min */
  static float wxb; /* Float pixel dist from map center to window X max */
  static float wya; /* Float pixel dist from map center to window Y min */
  static float wyb; /* Float pixel dist from map center to window Y max */
/*
 * Determine the central pixel.
 */
  xcent = mb->nx/2;
  ycent = mb->ny/2;
/*
 * Convert the window coordinates to numbers of pixels wrt the centre
 * of the map.
 */
  wxa = win->xmin/mb->xinc;
  wxb = win->xmax/mb->xinc;
  wya = win->ymin/mb->yinc;
  wyb = win->ymax/mb->yinc;
/*
 * Convert the ranges from radians to map pixel element numbers in such
 * a way that all pixels whose centres are enclosed by the window
 * are included.
 */
  xa = xcent + (int) (wxa + ((wxa<0)?0.0:1.0));
  xb = xcent + (int) (wxb - ((wxb<0)?1.0:0.0));
  ya = ycent + (int) (wya + ((wya<0)?0.0:1.0));
  yb = ycent + (int) (wyb - ((wyb<0)?1.0:0.0));
/*
 * add_win() ensures that win->xmax >= win->xmin and win->ymax >= win->ymin.
 * Thus if xb < xa then it must be because of the overlap in the above
 * equations when win->xmin and win->xmax are less than half a pixel apart.
 * In this case set both to the nearest pixel centre to the mean position.
 */
  if(xa > xb)
    xa = xb = xcent + fnint((wxa+wxb)/2.0f);
  if(ya > yb)
    ya = yb = ycent + fnint((wya+wyb)/2.0f);
/*
 * Is the window completely outside the map in X or Y?
 */
  if( ( (xa < ixmin && xb < ixmin) || (xa > ixmax && xb > ixmax) ) ||
      ( (ya < iymin && yb < iymin) || (ya > iymax && yb > iymax) ) ) {
    return 1;
  };
/*
 * Enforce bounds.
 */
  if(xa < ixmin)
    xa = ixmin;
  if(xb > ixmax)
    xb = ixmax;
  if(ya < iymin)
    ya = iymin;
  if(yb > iymax)
    yb = iymax;
/*
 * Set up values for return.
 */
  wr->xa = xa;
  wr->xb = xb;
  wr->ya = ya;
  wr->yb = yb;
  return 0;
}

/*.......................................................................
 * Determine and report statistics of the map pixels within a given
 * window. The statistics found are, rms, mean, total, max and min flux
 * within the window.
 *
 * Input:
 *  mb   MapBeam *  The container of the plotted image.
 *  domap    int    If true get stats from the map, otherwise from the
 *                  beam.
 *  win   Subwin *  The delimiting window.
 *  doall    int    0 - Only generate stats for 'win'.
 *                  1 - Generate stats for all windows in the list
 *                      starting at 'win'.
 */
void winstats(MapBeam *mb, int domap, Subwin *win, int doall)
{
  Winran wr;          /* Contains pixel limits of window 'win' */
  int npts=0;         /* Number of points used */
  float *fptr;        /* Pointer into map array */
  float sum_sqr=0.0f; /* Sum of squared flux */
  float flux=0.0f;    /* Total flux */
  float fmin=0.0f;    /* Min pixel flux in window */
  float fmax=0.0f;    /* Max pixel flux in window */
  float beam_area;    /* Beam area in pixels */
  int xwid,ywid;      /* Number of pixels in window along each dimension */
  int ix,iy;          /* Pixel coords wrt bottom of window */
  int xskip;   /* The number of floats between a right and next left edge */
/*
 * Check map container.
 */
  if(mb==0) {
    lprintf(stderr, "winstats: syserror - NULL map intercepted\n");
    return;
  };
/*
 * Look at each window in turn.
 */
  while(win != NULL) {
    if(win_pix(win, mb, 0, mb->nx-1, 0, mb->ny-1, &wr)==0) {
      xwid = wr.xb - wr.xa + 1;
      ywid = wr.yb - wr.ya + 1;
/*
 * Determine the number of pixels between the right edge of the current
 * window and the left edge of the same window on the next row.
 */
      xskip = mb->nx - xwid;
/*
 * Pointer to corner of current window with lowest address.
 */
      fptr = (domap?mb->map:mb->beam) + wr.xa + wr.ya * mb->nx;
      if(npts==0)
	fmin = fmax = *fptr;
/*
 * Iterate over the window area.
 */
      for(iy=0; iy < ywid; iy++,fptr += xskip) {
	for(ix=0; ix < xwid; ix++,fptr++) {
	  npts++;
	  if(*fptr < fmin)
	    fmin = *fptr;
	  if(*fptr > fmax)
	    fmax = *fptr;
	  flux += *fptr;
	  sum_sqr += *fptr * *fptr;
	};
      };
/*
 * Next window.
 */
      if(doall)
	win=win->next;
      else
	win=NULL;
    };
  };
/*
 * Report the results.
 */
  if(npts!=0) {
/*
 * If the map has been restored then we have an estimate of the beam
 * size so the total flux can be given in sensible units.
 */
    if(mb->ncmp) {
      beam_area = pi/(4.0*log(2.0))*mb->bmaj*mb->bmin /	(mb->xinc*mb->yinc);
      lprintf(stdout, "Total flux=%g Jy\n", flux/beam_area);
    };
    lprintf(stdout, "Mean=%g  rms=%g  min=%g  max=%g Jy/beam\n",
	   flux/npts, sqrt(sum_sqr/npts), fmin, fmax);
  };
  return;
}

/*.......................................................................
 * Save windows to a file.
 *
 * Input:
 *  mwin    Mapwin *  The container of the list of windows to be saved.
 *  filename  char *  The name of the new file to hold the windows,
 *                    or NULL or "" to send output to stdout.
 *  east     float    Any eastward X-axis offset to be removed before writing.
 *  north    float    Any northward Y-axis offset to removed before writing.
 *  do_old     int    If true write windows in old format, ie.
 *                    LRTB = n,n,n,n,n...
 * Output:
 *  return     int    0 - OK.
 *                    1 - I/O error.
 */
int wwins(Mapwin *mwin, char *filename, float east, float north, int do_old)
{
  FILE *fp;    /* Pointer to new file */
  Subwin *win; /* The window being written */
  int ierr;    /* Stores error status */
/*
 * No windows to save?
 */
  if(mwin==NULL || mwin->nwin==0)
    return 0;
/*
 * Use stdout instead of a normal file?
 */
  if(filename==NULL || filename[0] == '\0') {
    filename = "(stdout)";
    fp = stdout;
  } else {
/*
 * Attempt to open the new file.
 */
    fp = fopen(filename, "w");
    if(fp==NULL) {
      lprintf(stderr, "wwins: Couldn\'t create window file: %s\n", filename);
      return 1;
    };
  };
/*
 * Write textual intro.
 */
  lprintf(fp, "! CLEAN windows written by wwins in difmap.\n");
  lprintf(fp, "! Windows are specified as xmin xmax ymin ymax (mas).\n");
/*
 * Old LRTB window format?
 */
  if(do_old)
    lprintf(fp, "LRTB = ");
/*
 * Write out each member of the window list.
 */
  for(win=mwin->head; win != NULL && !ferror(fp); win=win->next) {
/*
 * Get the un-shifted window coordinates - in milli-arc-sec.
 */
    float xmax = (win->xmax - east) * rtomas;
    float xmin = (win->xmin - east) * rtomas;
    float ymax = (win->ymax - north) * rtomas;
    float ymin = (win->ymin - north) * rtomas;
/*
 * Write the window in the required format.
 */
    if(do_old)
      lprintf(fp, "%#15.9g, %#15.9g, %#15.9g, %#15.9g%s", -xmax,-xmin,ymax,ymin,
	      win->next ? ",\n       " : "\n");
    else
      lprintf(fp, "%#15.9g %#15.9g %#15.9g %#15.9g\n", xmin, xmax, ymin, ymax);
  };
/*
 * File write error or close error?
 */
  ierr = ferror(fp);
  if((fp!=stdout && fclose(fp)==EOF) || ierr) {
    lprintf(stderr, "I/O Error writing windows to file: %s\n", filename);
    return 1;
  };
  lprintf(stdout, "wwins: Wrote %d windows to %s\n", mwin->nwin, filename);
  return 0;
}

/*.......................................................................
 * Read windows from a file (previously saved with wwins() or one
 * written in the same format).
 *
 * Input:
 *  mwin    Mapwin *  The container of the list of windows to append to.
 *  filename  char *  The name of the file containing the windows.
 *  east     float    Any eastward X-axis offset to be added after reading.
 *  north    float    Any northward Y-axis offset to be added after reading.
 * Output:
 *  return     int    0 - OK.
 *                    1 - I/O error.
 *                    2 - mwin is NULL.
 */
int rwins(Mapwin *mwin, char *filename, float east, float north)
{
  FILE *fp;   /* Pointer to window file */
#define MAXLEN 80
  char linbuf[MAXLEN]; /* Input line buffer */
  char *cptr; /* Pointer into 'linbuf' */
  int nread;  /* Number of numbers read from a line */
  int nline;  /* The sequential number of the latest line in the file */
  float xa,xb;/* RA limits read from file */
  float ya,yb;/* DEC limits read from file */
/*
 * Bad window container?
 */
  if(mwin==NULL) {
    lprintf(stderr, "rwins: Illegal NULL container intercepted\n");
    return 2;
  };
/*
 * Attempt to open the window file.
 */
  fp = fopen(filename, "r");
  if(fp==NULL) {
    lprintf(stderr, "rwins: Couldn\'t open window file: %s\n", filename);
    return 1;
  };
/*
 * Read the windows. Each window occupies a single line and consists of
 * 4 numbers separated by spaces. While reading, report incomplete window
 * lines and assume that lines where the first non-white-space is not
 * - + or a digit constitute comment lines.
 */
  for(nline=0; fgets(linbuf, MAXLEN, fp) != NULL; nline++) {
/*
 * Skip white space.
 */
    cptr = &linbuf[0];
    while(isspace( (int) *cptr))
      cptr++;
/*
 * Attempt to read 4 space separated numbers from the input line.
 */
    nread = sscanf(linbuf, "%f%f%f%f", &xa, &xb, &ya, &yb);
    if(nread > 0 && nread < 4)
      lprintf(stderr, "Ignoring incomplete window on line: %d\n", nline);
    else {
/*
 * Apply coordinate transformations.
 */
      xa = xa * mastor + east;
      xb = xb * mastor + east;
      ya = ya * mastor + north;
      yb = yb * mastor + north;
/*
 * Append the new window to the window container.
 */
      if(nread == 4 && add_win(mwin, xa, xb, ya, yb) == NULL) {
	lprintf(stderr, "Failed to include windows on line %d and beyond\n",
		nline+1);
	return 1;
      };
    };
  };
/*
 * File write error?
 */
  if(ferror(fp) || fclose(fp)==EOF) {
    lprintf(stderr, "I/O Error reading windows from file: %s\n", filename);
    return 1;
  };
  lprintf(stdout, "rwins: Read %d windows from %s\n", mwin->nwin, filename);
  return 0;
}

/*.......................................................................
 * If the position of the peak absolute flux recorded with a map, is
 * not already enclosed by any of the given CLEAN windows, append a
 * new CLEAN window centered on the peak position. The size of the
 * window used is given by the rectangular area of the CLEAN beam.
 *
 * Input:
 *  mb   MapBeam *  The map/beam descriptor of the map whose peak is
 *                  to be windowed. The map must be in a ready state.
 *  mw    Mapwin *  The container of the list of CLEAN windows
 *                  associated with the map.
 *  size   float    The relative size of the clean box wrt the
 *                  fwhm aspect of the beam. (>0).
 *  doabs    int    If true search out the peak absolute component.
 *                  If false, search out the peak positive component.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int peakwin(MapBeam *mb, Mapwin *mw, float size, int doabs)
{
  float xpos,ypos; /* X,Y coordinates of the max absolute flux in the map */
/*
 * Check arguments.
 */
  if(mb==NULL || mw==NULL) {
    lprintf(stderr, "peakwin: NULL %s descriptor intercepted.\n",
	    mb==NULL ? "MapBeam":"Mapwin");
    return 1;
  };
/*
 * Are the map and beam up to date? The map is required to be up to date
 * because the recorded details of the peak pixel will ptherwise be out of date.
 * The beam is required to be up to date because we use the beam size estimated
 * by 'invert' when it last inverted the beam.
 */
  if(mb->domap || mb->dobeam) {
    lprintf(stderr, "peakwin: The map and/or beam is out of date.\n");
    return 1;
  };
/*
 * Get the position of the max absolute flux in the map.
 */
  if(doabs && fabs(mb->minpix.value) > fabs(mb->maxpix.value)) {
    xpos = mb->minpix.xpos;
    ypos = mb->minpix.ypos;
  } else {
    xpos = mb->maxpix.xpos;
    ypos = mb->maxpix.ypos;
  };
/*
 * If the position is not within any window prepare to window it.
 */
  if(!inmapwin(mw, xpos, ypos)) {
    Ellipse el;
/*
 * Get a description of the eliptical beam aspect in relation to the
 * current map coordinate system.
 */
    el_define(&el, mb->e_bmin, mb->e_bmaj, mb->e_bpa, 0.0f, 0.0f);
/*
 * Require that size be +ve.
 */
    if(size < 0.0f)
      size = -size;
/*
 * Add a window equal to the size of the beam area.
 */
    if(add_win(mw, xpos - size*el.xwid/2.0f, xpos + size*el.xwid/2.0,
	           ypos - size*el.ywid/2.0f, ypos + size*el.ywid/2.0)==NULL)
      return 1;
/*
 * Report the successful addition of a new window.
 */
    lprintf(stdout, "Added new window around map position (%g, %g).\n",
	    radtoxy(xpos), radtoxy(ypos));
  };
  return 0;
}
