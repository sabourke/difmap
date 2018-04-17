#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "markerlist.h"
#include "vlbconst.h"
#include "logio.h"

/*
 * Set the MarkerNode allocation blocking factor.
 */
#define MARKER_BLK_FACT 50

/*
 * Enumerate the known symbols.
 */
static Enumpar marker_symbols[] = {
  {"arrow",           MK_ARROW},
  {"filled_diamond",  MK_FILLED_DIAMOND},
  {"filled_triangle", MK_FILLED_TRIANGLE},
  {"open_square",     MK_OPEN_SQUARE},
  {"dot",             MK_DOT},
  {"plus",            MK_PLUS},
  {"asterisk",        MK_ASTERISK},
  {"cross",           MK_CROSS},
  {"open_triangle",   MK_OPEN_TRIANGLE },
  {"circle_plus",     MK_CIRCLE_PLUS},
  {"circle_dot",      MK_CIRCLE_DOT},
  {"knotted_hanky",   MK_KNOTTED_HANKY},
  {"open_diamond",    MK_OPEN_DIAMOND},
  {"open_star",       MK_OPEN_STAR},
  {"maltese_cross",   MK_MALTESE_CROSS},
  {"star_of_david",   MK_STAR_OF_DAVID},
  {"filled_square",   MK_FILLED_SQUARE},
  {"filled_circle",   MK_FILLED_CIRCLE},
  {"filled_star",     MK_FILLED_STAR},
  {"circle1",         MK_CIRCLE1},
  {"circle2",         MK_CIRCLE2},
  {"circle3",         MK_CIRCLE3},
  {"circle4",         MK_CIRCLE4},
};
static int num_marker_sym = sizeof(marker_symbols)/sizeof(marker_symbols[0]);

/*.......................................................................
 * Create a new MarkerList object.
 *
 * Output:
 *  return  MarkerList *  The new object, or NULL on error.
 */
MarkerList *new_MarkerList(void)
{
  MarkerList *markers;  /* The object to be returned */
/*
 * Allocate the container.
 */
  markers = malloc(sizeof(MarkerList));
  if(!markers) {
    lprintf(stderr, "new_MarkerList: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_MarkerList().
 */
  markers->symtab = NULL;
  markers->marker_mem = NULL;
  markers->head = NULL;
  markers->tail = NULL;
/*
 * Create a symbol table of plot symbols.
 */
  markers->symtab = new_Enumtab(marker_symbols, num_marker_sym,"Marker symbol");
  if(!markers->symtab)
    return del_MarkerList(markers);
/*
 * Allocate a freelist from which to allocate marker list node.
 */
  markers->marker_mem = new_FreeList("new_MarkerList", sizeof(MarkerNode),
				 MARKER_BLK_FACT);
  if(!markers->marker_mem)
    return del_MarkerList(markers);
  return markers;
}

/*.......................................................................
 * Delete a MarkerList object.
 *
 * Input:
 *  markers  MarkerList *  The object to be deleted.
 * Output:
 *  return   MarkerList *  The deleted object (always NULL).
 */
MarkerList *del_MarkerList(MarkerList *markers)
{
  if(markers) {
/*
 * Delete the list of markers.
 */
    clr_MarkerList(markers);
/*
 * Discard the freelist.
 */
    markers->marker_mem = del_FreeList("del_MarkerList", markers->marker_mem,1);
/*
 * Delete the plot-symbol symbol table.
 */
    markers->symtab = del_Enumtab(markers->symtab);
/*
 * Delete the container.
 */
    free(markers);
  };
  return NULL;
}

/*.......................................................................
 * Prepend a marker to a given list of markers.
 *
 * Input:
 *  markers MarkerList *  The list to add to.
 *  ra,dec      double    The location of the marker (radians).
 *  sym   MarkerSymbol    The PGPLOT symbol to use when plotting the symbol.
 *  color          int    The PGPLOT color index to use when plotting the
 *                        marker.
 *  size         float    The character size to use relative to the
 *                        normal size of 1.0.
 *  text          char *  The anotation to place next to the marker (or NULL
 *                        if not needed).
 *  just         float    The justification of the text, with 0 meaning
 *                        left justify, 1 meaning right justify, and 0.5
 *                        meaning center the text.
 *  xpos         float    The x-offset of the justification point of the
 *                        text from the marker (characters).
 *  ypos         float    The y-offset of the vertical center of the text
 *                        from the the marker (characters).
 * Output:
 *  return  MarkerNode *  The list node of the newly added marker, or NULL on
 *                        error.
 */
MarkerNode *add_MarkerNode(MarkerList *markers, double ra, double dec,
			   MarkerSymbol sym, int color, float size,
			   char *text, float just, float xpos, float ypos)
{
  MarkerNode *marker;   /* The new marker list node */
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "add_MarkerNode: NULL arguments.\n");
    return NULL;
  };
  if(ra < 0.0 || ra > twopi) {
    lprintf(stderr, "Marker Right Ascension out of range.\n");
    return NULL;
  };
  if(dec < -halfpi || dec > halfpi) {
    lprintf(stderr, "Marker Declination out of range.\n");
    return NULL;
  };
  if(color < 0 || color > 15) {
    lprintf(stderr, "Marker PGPLOT color %d outside supported 0-15 range.\n",
	    color);
    return NULL;
  };
  if(size <= 0) {
    lprintf(stderr, "Marker character size '%f' out of range.\n", size);
    return NULL;
  };
/*
 * Get a new marker node.
 */
  marker = new_FreeListNode("add_MarkerNode", markers->marker_mem);
  if(!marker)
    return NULL;
/*
 * Record the marker attributes.
 */
  marker->ra = ra;
  marker->dec = dec;
  marker->sym = sym;
  marker->color = color;
  marker->text = NULL;
  marker->size = size;
  marker->just = just;
  marker->xpos = xpos;
  marker->ypos = ypos;
  marker->next = NULL;
/*
 * Attempt to allocate a copy of the anotation text.
 */
  if(text && *text) {
    marker->text = malloc(strlen(text) + 1);
    if(!marker->text) {
      lprintf(stderr,"Insufficient memory to record marker annotation text.\n");
      return del_FreeListNode("add_MarkerNode", markers->marker_mem, marker);
    };
    strcpy(marker->text, text);
  };
/*
 * Append the node to the list of markers.
 */
  if(markers->head)
    markers->tail->next = marker;
  else
    markers->head = marker;
  markers->tail = marker;
  return marker;
}

/*.......................................................................
 * Delete a given marker node from a list of markers.
 *
 * Input:
 *  markers  MarkerList *   The list of markers.
 *  marker   MarkerNode *   The marker node to be deleted.
 * Output:
 *  return   MarkerNode *   Always NULL.
 */
MarkerNode *del_MarkerNode(MarkerList *markers, MarkerNode *marker)
{
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "del_MarkerNode: NULL marker list.\n");
    return NULL;
  };
/*
 * If the marker hasn't already been deleted, locate it and delete it
 * from the list.
 */
  if(marker) {
    MarkerNode *prev; /* The node in the list that precedes the target marker */
    MarkerNode *node; /* The node of the list being compared to the target */
/*
 * Search for the marker in the list.
 */
    for(prev=NULL, node=markers->head; node && node!=marker;
	prev=node, node=node->next)
      ;
/*
 * If the marker wasn't found, complain.
 */
    if(!node) {
      lprintf(stderr, "del_MarkerNode: Unknown marker.\n");
      return NULL;
    };
/*
 * Relink around the redundant marker node.
 */
    if(prev)
      prev->next = marker->next;
    else
      markers->head = marker->next;
    if(!marker->next)
      markers->tail = prev;
    marker->next = NULL;
/*
 * Delete the text anotation string.
 */
    if(marker->text) {
      free(marker->text);
      marker->text = NULL;
    };
/*
 * Return the marker node to the freelist.
 */
    marker = del_FreeListNode("del_MarkerNode", markers->marker_mem, marker);
  };
  return NULL;
}

/*.......................................................................
 * Find the closest marker to a given RA, Dec position.
 *
 * Input:
 *  markers MarkerList *  The list of markers to be searched.
 *  ra,dec      double    The location to search for (radians).
 * Output:
 *  return  MarkerNode *  The nearest marker, or NULL if the list is empty.
 */
MarkerNode *closest_MarkerNode(MarkerList *markers, double ra, double dec)
{
  MarkerNode *marker;    /* The marker node being examined */
  MarkerNode *best=NULL; /* The closest marker so far */
  double minrsqr=0.0;    /* The minimum distance squared found so far */
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "closest_MarkerNode: NULL marker list.\n");
    return NULL;
  };
/*
 * Find the closest marker.
 */
  for(marker=markers->head; marker; marker=marker->next) {
    double dra = ra - marker->ra;
    double ddec = dec - marker->dec;
    double rsqr = dra * dra + ddec * ddec;
    if(rsqr < minrsqr || !best) {
      best = marker;
      minrsqr = rsqr;
    };
  };
/*
 * Return the best matching marker.
 */
  return best;
}

/*.......................................................................
 * Clear a list of markers.
 *
 * Input:
 *  markers  MarkerList *  The list of markers to be searched.
 * Output:
 *  return          int    0 - OK.
 *                         1 - Error.
 */
int clr_MarkerList(MarkerList *markers)
{
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "clr_MarkerNode: NULL marker list.\n");
    return 1;
  };
/*
 * Clear the list.
 */
  while(markers->head)
    del_MarkerNode(markers, markers->head);
  return 0;
}

/*.......................................................................
 * Lookup a marker plot symbol by name.
 *
 * Input:
 *  markers  MarkerList *  The resource object of the list of markers.
 *  name           char *  The name of the symbol.
 * Output:
 *  return MarkerSymbol    The named plot symbol, or MK_UNKNOWN if not
 *                         found.
 */
MarkerSymbol lookup_marker_symbol(MarkerList *markers, char *name)
{
  Enumpar *epar;   /* The enumeration entry of the symbol */
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "lookup_marker_symbol: NULL marker list.\n");
    return 1;
  };
/*
 * Lookup the symbol.
 */
  epar = find_enum(markers->symtab, name);
  return epar ? epar->id : MK_UNKNOWN;
}

/*.......................................................................
 * Return the name of a given marker symbol number.
 *
 * Input:
 *  markers   MarkerList *  The resource object of the list of markers.
 *  sym     MarkerSymbol    The symbol code to lookup.
 * Output:
 *  return          char *  The name of the marker, or "unknown".
 */
char *lookup_marker_name(MarkerList *markers, MarkerSymbol sym)
{
/*
 * Check the arguments.
 */
  if(!markers) {
    lprintf(stderr, "lookup_marker_symbol: NULL marker list.\n");
    return "unknown";
  };
  return name_enum(markers->symtab, sym, "unknown");
}

