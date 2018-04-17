#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "sysfits.h"

/*
 * Get a mask with the most-significant-bit of a char set.
 */
#define CHR_SGN_MASK (1U<<CHAR_BIT-1U)

/*
 * The following are conversion/copying functions for use to convert
 * between native and FITS datatypes.
 */

#ifdef NEED_B2B4
/*.......................................................................
 * Orig: Two's complement big-endian 2-byte integers.
 * Dest: Two's complement big-endian 4-byte integers.
 * The 2 MSBs of the output are sign extended.
 */
void cp_b2b4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=2, dest+=4) {
    dest[0] = dest[1] = orig[0] & CHR_SGN_MASK ? ~0U : 0U;
    dest[2] = orig[0];
    dest[3] = orig[1];
  };
}
/*.......................................................................
 * Orig: Two's complement big-endian 4-byte integers.
 * Dest: Two's complement big-endian 2-byte integers.
 * The 2 MSBs of the input integers are dropped.
 */
void cp_b4b2(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=2) {
    dest[0] = orig[2];
    dest[1] = orig[3];
  };
}
#endif


#ifdef NEED_B8B4
/*.......................................................................
 * Orig: Two's complement big-endian 8-byte integers.
 * Dest: Two's complement big-endian 4-byte integers.
 * The 4 MSBs of the input integers are dropped.
 */
void cp_b8b4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=8, dest+=4) {
    dest[0] = orig[4];
    dest[1] = orig[5];
    dest[2] = orig[6];
    dest[3] = orig[7];
  };
}
/*.......................................................................
 * Orig: Two's complement big-endian 4-byte integers.
 * Dest: Two's complement big-endian 8-byte integers.
 * The 4 MSBs of the output are sign extended.
 */
void cp_b4b8(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=8) {
    dest[0] = dest[1] = dest[2] = dest[3] = orig[0] & CHR_SGN_MASK ? ~0U : 0U;
    dest[4] = orig[0];
    dest[5] = orig[1];
    dest[6] = orig[2];
    dest[7] = orig[3];
  };
}
#endif

#ifdef NEED_2R2
/*.......................................................................
 * Orig: 2-byte datatype.
 * Dest: 2-byte datatype with byte order reversed.
 */
void cp_2r2(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=2, dest+=2) {
    dest[0] = orig[1];
    dest[1] = orig[0];
  };
}
#endif

#ifdef NEED_4R4
/*.......................................................................
 * Orig: 4-byte datatype.
 * Dest: 4-byte datatype with byte order reversed.
 */
void cp_4r4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=4) {
    dest[0] = orig[3];
    dest[1] = orig[2];
    dest[2] = orig[1];
    dest[3] = orig[0];
  };
}
#endif

#ifdef NEED_8R8
/*.......................................................................
 * Orig: 8-byte datatype.
 * Dest: 8-byte datatype with byte order reversed.
 */
void cp_8r8(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=8, dest+=8) {
    dest[0] = orig[7];
    dest[1] = orig[6];
    dest[2] = orig[5];
    dest[3] = orig[4];
    dest[4] = orig[3];
    dest[5] = orig[2];
    dest[6] = orig[1];
    dest[7] = orig[0];
  };
}
#endif

#ifdef NEED_L2B4
/*.......................................................................
 * Orig: Two's complement little-endian 2-byte integers.
 * Dest: Two's complement big-endian 4-byte integers.
 * The 2 MSBs of the output integers are sign-extended.
 */
void cp_l2b4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=2, dest+=4) {
    dest[0] = dest[1] = orig[1] & CHR_SGN_MASK ? ~0U : 0U;
    dest[2] = orig[1];
    dest[3] = orig[0];
  };
}
/*.......................................................................
 * Orig: Two's complement big-endian 4-byte integers.
 * Dest: Two's complement little-endian 2-byte integers.
 * The 2 MSBs of the input integers are dropped.
 */
void cp_b4l2(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=2) {
    dest[0] = orig[3];
    dest[1] = orig[2];
  };
}
#endif

#ifdef NEED_L4B2
/*.......................................................................
 * Orig: Two's complement little-endian 4-byte integers.
 * Dest: Two's complement big-endian 2-byte integers.
 * The 2 MSBs of the input integers are dropped.
 */
void cp_l4b2(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=2) {
    dest[0] = orig[1];
    dest[1] = orig[0];
  };
}
/*.......................................................................
 * Orig: Two's complement big-endian 2-byte integers.
 * Dest: Two's complement little-endian 4-byte integers.
 * The 2 MSBs of the output integers are sign-extended.
 */
void cp_b2l4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=2, dest+=4) {
    dest[0] = orig[1];
    dest[1] = orig[0];
    dest[2] = dest[3] = orig[0] & CHR_SGN_MASK ? ~0U : 0U;
  };
}
#endif

#ifdef NEED_L8B4
/*.......................................................................
 * Orig: Two's complement little-endian 8-byte integers.
 * Dest: Two's complement big-endian 4-byte integers.
 * The 4 MSBs of the input integers are dropped.
 */
void cp_l8b4(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=8, dest+=4) {
    dest[0] = orig[3];
    dest[1] = orig[2];
    dest[2] = orig[1];
    dest[3] = orig[0];    
  };
}
/*.......................................................................
 * Orig: Two's complement big-endian 4-byte integers.
 * Dest: Two's complement little-endian 8-byte integers.
 * The 4 MSBs of the output integers are sign-extended.
 */
void cp_b4l8(unsigned char *dest, unsigned char *orig, size_t nitem)
{
  size_t i;
  for(i=0; i<nitem; i++, orig+=4, dest+=8) {
    dest[0] = orig[3];
    dest[1] = orig[2];
    dest[2] = orig[1];
    dest[3] = orig[0];
    dest[4] = dest[5] = dest[6] = dest[7] = orig[0] & CHR_SGN_MASK ? ~0U : 0U;
  };
}
#endif
