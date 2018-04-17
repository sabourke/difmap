/* FITS Unknown extension-type HDU code */
#include <stdio.h>
#include <stdlib.h>

#include "sysfits.h"
#include "fits.h"
#include "utils.h"
#include "uhdu.h"
#include "fitkey.h"

static GETFN(get_uhdu);
static NEWFN(new_uhdu);
static DELFN(del_uhdu);
static SIZEFN(siz_uhdu);
static ADDFN(add_uhdu);
static COPYFN(cop_uhdu);
static ENDFN(end_uhdu);

Hdutab uhdufns={" ", " ", get_uhdu, new_uhdu, del_uhdu, siz_uhdu, add_uhdu,
		cop_uhdu, end_uhdu};

/*.......................................................................
 * Clear the derived part of an unknown-header descriptor to the
 * point at which it can be safely sent to del_Hdu().
 *
 * Input:
 *  hdu     Hdu *  The base-class pointer to a Uhdu descriptor.
 * Output:
 *  return  Hdu *  'hdu', or NULL if an error occurs and hdu is deleted.
 */
static NEWFN(new_uhdu)
{
  return hdu;
}

/*.......................................................................
 * Delete the derived parts of an Unknown-type HDU descriptor.
 *
 * Input:
 *  hdu      Hdu *  The base-class pointer to the Uhdu descriptor.
 */
static DELFN(del_uhdu)
{
  return;
}

/*.......................................................................
 * Read the Unkown-type HDU parts of a FITS header.
 *
 * Input:
 *  fits  Fits *  The descriptor of the FITS file.
 * Input/Output:
 *  hdu    Hdu *  The base-class descriptor to a Uhdu descriptor.
 * Output:
 *  return int    0 - OK.
 */
static GETFN(get_uhdu)
{
  return 0;
}

/*.......................................................................
 * Return the size of an Unkown-type HDU descriptor.
 *
 * Output:
 *  return  size_t  sizeof(Uhdu).
 */
static SIZEFN(siz_uhdu)
{
  return sizeof(Uhdu);
}

/*.......................................................................
 * Make it illegal to write an Unknown-type extension HDU into a FITS
 * file.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file to which the HDU is
 *                  being added.
 *  hdu      Hdu *  An initialized HDU descriptor for the new HDU.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ADDFN(add_uhdu)
{
  fprintf(stderr,"add_uhdu: Illegal attempt to add an unhandled extension type to a FITS file.\n");
  return 1;
}

/*.......................................................................
 * Copy the derived parts of a Uhdu descriptor into another Uhdu
 * descriptor. This function should be called only via copy_Hdu().
 *
 * Input:
 *  hdu    Hdu * The original hdu to be copied.
 * Output:
 *  return Hdu * The new copied version of 'hdu', or NULL on error.
 */
static COPYFN(cop_uhdu)
{
  return 0;
}

/*.......................................................................
 * Make it illegal to complete an Unknown-type extension HDU into a FITS
 * file.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file in which the HDU is
 *                  being completed.
 *  hdu      Hdu *  An initialized HDU descriptor for the new HDU.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static ADDFN(end_uhdu)
{
  fprintf(stderr,"end_uhdu: Illegal attempt to complete an unhandled extension type.\n");
  return 1;
}
