#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "cksum.h"

/*
 * Set the checksum divisor key.
 */
#define CKSUM_KEY 0x04c11db7UL  /* The key used by the ethernet protocol */

/*
 * Define a type of at least 32 bits.
 */
typedef unsigned long CrcType;

struct CheckSum {
  CrcType table[256];   /* A lookup table of byte-specific checksums */
};


/*.......................................................................
 * Create a new CheckSum object.
 *
 * Output:
 *  return  CheckSum *  The new object, or NULL on error.
 */
CheckSum *new_CheckSum(void)
{
  CheckSum *cs;  /* The object to be returned */
  CrcType mask;  /* The result of a checksum calculation on a byte */
  int i,j;
/*
 * Allocate the container.
 */
  cs = malloc(sizeof(CheckSum));
  if(!cs) {
    lprintf(stderr, "new_CheckSum: Insufficient memory.\n");
    return NULL;
  };
/*
 * Loop over all possible 8-bit bytes.
 */
  for(i=0; i<256; i++) {
/*
 * Place the byte in the most significant byte of least significant
 * 32 bits of the accumulation mask.
 */
    mask = i << 24;
/*
 * Loop over the 8 bits in a byte.
 */
    for(j=0; j<8; j++) {
/*
 * Get the current most significant bit in the 32 bit mask.
 */
      unsigned msb = mask & 0x80000000;
/*
 * Shift the mask left. Note that this shifts the bit that we just
 * stored above the last bit of the 32 bit number that we are working
 * on. Later, it will be masked off, so effectively we are disgarding
 * it here.
 */
      mask <<= 1;
/*
 * If the bit that we just dropped was on, exclusive OR the mask with
 * the division key (we are doing a weird form of long division, where
 * exclusive OR is used to do an addition, but without doing any carries).
 */
      if(msb)
	mask ^= CKSUM_KEY;
    };
/*
 * The mask now contains the remainder of the division, and thus the
 * checksum of the latest byte. Record this in the table, after masking
 * off any other bytes, just in case CrcType has more than 32 bits.
 */
    cs->table[i] = mask & 0xFFFFFFFFUL;
  };
  return cs;
}

/*.......................................................................
 * Delete a CheckSum object.
 *
 * Input:
 *  cs     CheckSum *  The object to be deleted.
 * Output:
 *  return CheckSum *  The deleted object (always NULL).
 */
CheckSum *del_CheckSum(CheckSum *cs)
{
  if(cs) {
    free(cs);
  };
  return NULL;
}

/*.......................................................................
 * Compute a 32-bit cyclic redundancy checksum of an array of bytes, using
 * the ethernet key.
 *
 * Input:
 *  cs          CheckSum *  A checksum resource object.
 *  obj             void *  A pointer to the first byte of the object to be
 *                          characterized.
 *  nbyte         size_t    The number of bytes in the object.
 * Output:
 *  return unsigned long    The checksum of the object, or 0 if cs or obj
 *                          is NULL, or nbyte < 1.
 */
unsigned long checksum_of_object(CheckSum *cs, void *obj, size_t nbyte)
{
  unsigned char *bptr;   /* A pointer to a byte in obj */
  CrcType sum = 0;       /* The checksum */
  size_t i;
/*
 * Check the arguments.
 */
  if(!cs || !obj) {
    lprintf(stderr, "checksum_of_object: NULL argument(s).\n");
    return 0;
  };
  if(nbyte < 1) {
    lprintf(stderr, "checksum_of_object: Zero sized object.\n");
    return 0;
  };
/*
 * Loop through the bytes of the object.
 */
  bptr = (unsigned char *) obj;
  for(i=0; i<nbyte; i++, bptr++) {
/*
 * Preserve a copy of the most significant byte of the 32 least significant
 * bits of the checksum accumulation mask.
 */
    unsigned msb = (sum >> 24) & 0xFFU;
/*
 * Shift the checksum of the byte into the least significant byte of
 * the checksum accumulation mask, dropping the most significant byte.
 */
    sum = (sum << 8) ^ cs->table[msb ^ (unsigned) *bptr];
  };
/*
 * Return the significant 32 bits of the checksum.
 */
  return sum & 0xFFFFFFFFUL;
}

