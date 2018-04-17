#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "logio.h"
#include "slalib.h"
#include "vlbconst.h"
#include "planet.h"


/*
 * Collect information about a given planet.
 */
typedef struct {
  const char *name;  /* The name of the planet */
  int slalib_id;     /* The slalib id of the planet */
  double radius;     /* The equatorial radius of the planet (m) */  
  double flattening; /* The flattening of the planet (a-b)/a */
} PlanetData;

/*
 * List the pertinent details of each planet.
 */
static const PlanetData planet_data[] = {
  {"Mercury", 1, 2439700.0,  0.0},
  {"Venus",   2, 6051900.0,  0.0},
  {"Mars",    4, 3397000.0,  0.0065},
  {"Jupiter", 5, 71492000.0, 0.06487},
  {"Saturn",  6, 60268000.0, 0.09796},
  {"Uranus",  7, 25559000.0, 0.02293},
  {"Neptune", 8, 24764000.0, 0.0171},
};

static const PlanetData *pln_lookup_planet(const char *name);

/*.......................................................................
 * Lookup a planet by name.
 *
 * Input:
 *  name   const char *  The name of the planet (case insensitive).
 * Output:
 *  return PlanetData *  The entry of the planet in the planet_data[]
 *                       array, or NULL if not found.
 */
static const PlanetData *pln_lookup_planet(const char *name)
{
  int i;
/*
 * Compare the specified name to each entry in the planet_data[] array.
 */
  for(i=0; i<sizeof(planet_data)/sizeof(planet_data[0]); i++) {
    const char *nptr = name;
    const char *pptr = planet_data[i].name;
/*
 * Perform a case-insensitive comparison of the specified planet name
 * with the name that appears in the current entry of the table.
 */
    while(*nptr && *pptr) {
/*
 * Convert both input characters to lower case.
 */
      char nc = islower((int)*nptr) ? *nptr : tolower((int)*nptr);
      char pc = islower((int)*pptr) ? *pptr : tolower((int)*pptr);
/*
 * If the characters differ, abort the current comparison.
 */
      if(nc != pc)
	break;
/*
 * Advance to the next character.
 */
      nptr++;
      pptr++;
    };
/*
 * If the names matched, return the corresponding table entry.
 */
    if(*nptr=='\0' && *pptr=='\0')
      return planet_data + i;
  };
/*
 * Planet not found.
 */
  lprintf(stderr, "Unable to find information on planet: %s.\n", name);
  return NULL;
}

/*.......................................................................
 * Return the approximate position and the angular equatorial diameter
 * of a given planet seen from the center of the Earth, at a given time.
 * This function is an adaptation of the slalib slaRdplan() function,
 * modified to return geocentric instead of topocentric results.
 *
 * Input:
 *  name     char *   The name of the planet to be looked up
 *                    (case insensitive). Note that the slalib slaPlanet()
 *                    function doesn't recognize Pluto, so this function
 *                    doesn't either.
 *  tt     double     The terrestrial time expressed as a Modified Julian
 *                    date (JD - 2400000.5). For most applications, UTC
 *                    can be used instead of TT.
 * Input/Output:
 *  ra     double *   The geocentric Right Ascension of the planet (radians),
 *                    will be assigned to *ra.
 *  dec    double *   The geocentric Declination of the planet (radians),
 *                    will be assigned to *dec.
 *  diam    float *   The angular equatorial diameter of the planet, as seen
 *                    from the center of the Earth (radians), will be
 *                    assigned to *diam.
 *  flat    float *   The flattening (a-b)/a will be assigned to *flat.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */ 
int planet_geometry(const char *name, double tt, double *ra, double *dec,
		    float *diam, float *flat)
{
  const PlanetData *data;   /* The characteristics of the planet */
  double vgm[6], v[6], rmat[3][3], vse[6], vsg[6], vsp[6], r, tl;
  int i, j;

/*
 * Lookup the characteristics of the planet.
 */
  data = pln_lookup_planet(name);
  if(!data)
    return 1;
/*
 * Get the geocentric distance and velocity of the Moon at the given date.
 */
  slaDmoon(tt, v);
/*
 * Nutation to true of date
 */
  slaNut(tt, rmat);
  slaDmxv(rmat, &v[0], &vgm[0]);
  slaDmxv(rmat, &v[3], &vgm[3]);
/*
 * Precession/nutation matrix, J2000 to date.
 */
  slaPrenut(2000.0, tt, rmat);
/*
 * Sun to Earth-Moon Barycentre (J2000)
 */
  slaPlanet(tt, 3, v, &j);
/*
 * Precession and nutation to date.
 */
  slaDmxv(rmat, &v[0], &vse[0]);
  slaDmxv(rmat, &v[3], &vse[3]);
/*
 * Sun to geocentre.
 */
  for(i = 0; i <= 5; i++) vsg[i] = vse[i] - 0.012150581 * vgm[i];
/*
 * Sun to Planet.
 */
  slaPlanet(tt, data->slalib_id, v, &j);
/*
 * Precession and nutation to date.
 */
  slaDmxv(rmat, &v[0], &vsp[0]);
  slaDmxv(rmat, &v[3], &vsp[3]);
/*
 * Geocentre to planet.
 */
  for(i = 0; i <= 5; i++) v[i] = vsp[i] - vsg[i];
/*
 * Geometric distance (m).
 */
  r = au_to_m * sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
/*
 * Light time (seconds).
 */
  tl = r / cvel;
/*
 * Correct position for planetary aberration.
 */
  for(i = 0; i <= 2; i++) v[i] -= tl * v[i+3];
/*
 * To RA,Dec.
 */
  slaDcc2s(v, ra, dec);
  *ra = slaDranrm(*ra);
/*
 * Angular diameter (radians).
 */
  *diam = 2.0 * asin(data->radius / r);
/*
 * Return the flattening of the planet.
 */
  *flat = data->flattening;
  return 0;
}
