#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "vlbconst.h"
#include "units.h"
#include "symtab.h"
#include "logio.h"

/* Create a units descriptor */

typedef struct {
  double conv;     /* Factor to multiply internal units to get user units */
  char *name;      /* The official name of the unit */
  char *tlabel;    /* A label to use in text */
  char *plabel;    /* A label to use in PGPLOT labels */
} Unittype;

/* The users chosen map units will be paired with appropriate UVW units
 * through a descriptor of the following type.
 */

typedef struct {
  Unittype map;    /* Map units descriptor */
  Unittype uvw;    /* UVW units descriptor */
} Skyunits;

/* List the attributes of each of the supported unit pairs.
 * NB. The first entry describes the default units.
 */

static Skyunits unit_table[] = {
  {
    {rtomas,  "mas",    "milli-arcsec",     "mas"},
    {1.0e-6,  "Mw",     "mega-wavelengths", "10\\u6 \\d\\gl"},
  },
  {
    {rtoas,  "arcsec", "arcsec",           "arcsec"},
    {1.0e-3, "kw",     "kilo-wavelengths", "10\\u3 \\d\\gl"},
  },
  {
    {rtoam,  "arcmin", "arcmin",           "arcmin"},
    {1.0e-3, "kw",     "kilo-wavelengths", "10\\u3 \\d\\gl"},
  },
};
static int nunit = sizeof(unit_table) / sizeof(unit_table[0]);

/*
 * The following file-scope variable points to the last descriptor
 * selected from unit_table[] via skyunits(). The initialization below
 * selects milli-arcseconds as the default unit.
 */

static Skyunits *sky_units = &unit_table[0];

/*.......................................................................
 * Stage a new map units descriptor for return by future calls to
 * get_Mapunits().
 *
 * Input:
 *  name    char * The name of the unit to stage a descriptor for.
 *                 The recognised unit names are:
 *                   "arcmin"  -  Arc minutes.
 *                   "arcsec"  -  Arc seconds.
 *                   "mas"     -  Milli-arcseconds.
 *                  The string comparison is caseless and ignores
 *                  leading and trailing white-space.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int skyunits(char *name)
{
  static Symtab *sym = NULL;   /* The units symbol table */
  Skyunits *sky;               /* The matching entry from the symbol table */
  int i;
/*
 * Bad name?
 */
  if(name==NULL) {
    lprintf(stderr, "skyunits: NULL name intercepted.\n");
    return 1;
  };
/*
 * Allocate and initialize the symbol table?
 */
  if(!sym) {
    if(!(sym = new_Symtab(nunit, "Map unit", amb_report, 0)))
      return 1;
    for(i=0; i<nunit; i++) {
      if(add_symbol(sym, unit_table[i].map.name, (void *) &unit_table[i], 0))
	return 1;
    };
  };
/*
 * Search for the given unit name in the symbol table.
 */
  sky = (Skyunits *) get_symbol(sym, name, 1);
  if(!sky)
    return 1;
  sky_units = sky;
  return 0;
}

/*.......................................................................
 * Convert the given argument from user map units to radians.
 *
 * Input:
 *  xy     double    The map coordinate in user units to convert to radians.
 * Output:
 *  return double    The converted coordinate, in radians.
 */
double xytorad(double xy)
{
  return xy / sky_units->map.conv;
}

/*.......................................................................
 * Convert the given argument from radians to user map units.
 *
 * Input:
 *  rad    double    The map coordinate in radians to convert to user units.
 * Output:
 *  return double    The converted coordinate, in user units.
 */
double radtoxy(double rad)
{
  return rad * sky_units->map.conv;
}

/*.......................................................................
 * Convert the given argument from user UVW units to wavelengths.
 *
 * Input:
 *  uv     double    The UVW coordinate in user units to convert to wavelengths.
 * Output:
 *  return double    The converted coordinate, in wavelengths.
 */
double uvtowav(double uv)
{
  return uv / sky_units->uvw.conv;
}

/*.......................................................................
 * Convert the given argument from wavelengths to user UVW units.
 *
 * Input:
 *  uv     double    The UVW coordinate in wavelengths to convert to user units.
 * Output:
 *  return double    The converted coordinate, in user units.
 */
double wavtouv(double wav)
{
  return wav * sky_units->uvw.conv;
}

/*.......................................................................
 * Return a string to be used to label the current user selected map units.
 *
 * Input:
 *  ltype   Ultype   The type of label to return.
 *                    U_NAME  -  The official name of the type.
 *                    U_TLAB  -  The label to use in text.
 *                    U_PLAB  -  The PGPLOT label to use in plots.
 * Output:
 *  return    char * The '\0' terminated label string.
 */
char *mapunits(Ultype ltype)
{
  switch(ltype) {
  case U_NAME:
    return sky_units->map.name;
  case U_TLAB:
    return sky_units->map.tlabel;
  case U_PLAB:
    return sky_units->map.plabel;
  default:
    return "Unknown";
  };
}

/*.......................................................................
 * Return a string to be used to label the current user selected UVW units.
 *
 * Input:
 *  ltype   Ultype   The type of label to return.
 *                    U_NAME  -  The official name of the type.
 *                    U_TLAB  -  The label to use in text.
 *                    U_PLAB  -  The PGPLOT label to use in plots.
 * Output:
 *  return    char * The '\0' terminated label string.
 */
char *uvwunits(Ultype ltype)
{
  switch(ltype) {
  case U_NAME:
    return sky_units->uvw.name;
  case U_TLAB:
    return sky_units->uvw.tlabel;
  case U_PLAB:
    return sky_units->uvw.plabel;
  default:
    return "Unknown";
  };
}

/*.......................................................................
 * Return the two character ordinal suffix of an integer (eg. "th" for 13)
 *
 * Input:
 *  n        int    The number for which the suffix applies.
 * Output:
 *  return  char *  The suffix.
 */
char *ordinal_suffix(int n)
{
/*
 * Handle special cases.
 */
  if(n>=11 && n<=13)
    return "th";
/*
 * General cases.
 */
  switch(n % 10) {
  case 1:
    return "st";
  case 2:
    return "nd";
  case 3:
    return "rd";
  default:
    return "th";
  };
}

