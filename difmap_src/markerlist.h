#ifndef markerlist_h
#define markerlist_h

#include "freelist.h"
#include "enumpar.h"

/*
 * Enumerate the plotting symbols.
 * Most of these are pgplot symbols.
 */
typedef enum {
  MK_UNKNOWN = -10,
  MK_ARROW = -9,
  MK_FILLED_DIAMOND = -4,
  MK_FILLED_TRIANGLE = -3,
  MK_OPEN_SQUARE = 0,
  MK_DOT = 1,    
  MK_PLUS = 2,
  MK_ASTERISK = 3,
  MK_CROSS = 5,
  MK_OPEN_TRIANGLE = 7, 
  MK_CIRCLE_PLUS = 8,
  MK_CIRCLE_DOT = 9,
  MK_KNOTTED_HANKY = 10,
  MK_OPEN_DIAMOND = 11,
  MK_OPEN_STAR = 12,
  MK_MALTESE_CROSS = 14,
  MK_STAR_OF_DAVID = 15,
  MK_FILLED_SQUARE = 16,
  MK_FILLED_CIRCLE = 17,
  MK_FILLED_STAR = 18,
  MK_CIRCLE1 = 20,
  MK_CIRCLE2 = 21,
  MK_CIRCLE3 = 22,
  MK_CIRCLE4 = 23
} MarkerSymbol;

/*
 * The definition of each marker is recorded in a list node of the
 * following form.
 */
typedef struct MarkerNode MarkerNode;
struct MarkerNode {
  double ra,dec;   /* The position of the marker (radians) */
  MarkerSymbol sym;/* The pgplot symbol to use to represent the object */
  int color;       /* The color to plot the symbol with */
  float size;      /* The PGPLOT character height to use */
  char *text;      /* The anotation string (or NULL if not wanted) */
  float just;      /* The justification of the text */
  float xpos;      /* The horizontal location of the justification point */
                   /*  of the text relative to the symbol (characters). */
  float ypos;      /* The vertical location of the justification point */
                   /*  of the text relative to the symbol (characters). */
  MarkerNode *next;/* The next node in the list */
};

/*
 * Encapsulate a list of markers.
 */
typedef struct {
  Enumtab *symtab;      /* The symbol table of marker symbols */
  FreeList *marker_mem; /* The freelist of marker list nodes */
  MarkerNode *head;     /* The head of the list of markers */
  MarkerNode *tail;     /* The tail of the list of markers */
} MarkerList;

MarkerList *new_MarkerList(void);
int clr_MarkerList(MarkerList *markers);
MarkerList *del_MarkerList(MarkerList *markers);

MarkerNode *add_MarkerNode(MarkerList *markers, double ra, double dec,
			   MarkerSymbol sym, int color, float size,
			   char *text, float just, float xpos, float ypos);
MarkerNode *closest_MarkerNode(MarkerList *markers, double ra, double dec);
MarkerNode *del_MarkerNode(MarkerList *markers, MarkerNode *marker);
MarkerSymbol lookup_marker_symbol(MarkerList *markers, char *name);
char *lookup_marker_name(MarkerList *markers, MarkerSymbol sym); 

#endif
