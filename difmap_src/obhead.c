#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "logio.h"
#include "obs.h"

static char *headcop(const char *string);
static void headdel(char **name);

/*.......................................................................
 * Initialize the container of descriptive FITS keyword values in a
 * given observation with malloc'd copies of strings where provided, or
 * NULL otherwise. NB. This function over-writes whatever was in the
 * container, so if you want to change the contents of an initialized
 * container, call clr_Obhead first to free existing strings.
 *
 * Input:
 *  origin    char *  Organization that created the FITS file.
 *  date_obs  char *  Date of observation as DD/MM/YY.
 *  telescop  char *  Telescope used in taking the observation.
 *  instrume  char *  Instrument used in taking the observation.
 *  observer  char *  The name of the observer.
 *  bunit     char *  The units that the data are recorded with.
 *  equinox double    The equinox of the data - or 0.0 if not known.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int ini_Obhead(Observation *ob, char *origin, char *date_obs, char *telescop,
	       char *instrume, char *observer, char *bunit, double equinox)
{
/*
 * Check the Observation descriptor.
 */
  if(!ob_ready(ob, OB_ALLOC, "ini_Obhead"))
    return 1;
/*
 * Copy each string.
 */
  ob->misc.origin   = headcop(origin);
  ob->misc.date_obs = headcop(date_obs);
  ob->misc.telescop = headcop(telescop);
  ob->misc.instrume = headcop(instrume);
  ob->misc.observer = headcop(observer);
  ob->misc.bunit    = headcop(bunit);
  ob->misc.equinox  = equinox;
  return 0;
}

/*.......................................................................
 * Private function of ini_Obhead() used to return a dynamically allocated
 * copy of a string.
 *
 * Input:
 *  string    char *   The string to be copied, or NULL.
 * Output:
 *  return    char *   The copied string, or NULL if string was NULL or
 *                     on memory allocation errors.
 */
static char *headcop(const char *string)
{
  char *sptr;    /* Pointer to the copied string */
/*
 * Nothing to copy?
 */
  if(string==NULL)
    return NULL;
/*
 * Skip white-space.
 */
  while(isspace((int)*string))
    string++;
/*
 * Empty string?
 */
  if(*string == '\0')
    return NULL;
/*
 * Allocate sufficient space for a copy of string[].
 */
  sptr = malloc(strlen(string)+1);
  if(sptr)
/*
 * Copy the string.
 */
    strcpy(sptr, string);
  else
    lprintf(stderr,"headcop: Insufficient memory to copy string: %s\n", string);
  return sptr;
}

/*.......................................................................
 * Delete the contents of an Obhead container. ini_Obhead() must have
 * been called before this function on the given Observation descriptor.
 *
 * Input:
 *  ob  Observation *  The descriptor of the observation.
 */
void clr_Obhead(Observation *ob)
{
/*
 * Anything to delete?
 */
  if(ob) {
    Obhead *h = &ob->misc;
    headdel(&h->origin);
    headdel(&h->date_obs);
    headdel(&h->telescop);
    headdel(&h->instrume);
    headdel(&h->observer);
    headdel(&h->bunit);
  };
  return;
}

/*.......................................................................
 * Private function of clr_Obhead() used to free a given entry and assign
 * NULL to it.
 */
static void headdel(char **name)
{
  if(*name)
    free(*name);
  *name = NULL;
  return;
}
