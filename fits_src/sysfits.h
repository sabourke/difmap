#ifndef sysfits_h
#define sysfits_h

/* System dependent fits definitions file */

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

/*
 * FITS datatypes are stored in an architecture independent binary format.
 * This means that on architectures where the native datatypes are not
 * represented in the same way as FITS datatypes, functions are needed
 * to convert between the two formats. The following macros are defined
 * with arguments:
 *
 *  unsigned char *to     A pointer to an array of 'nitem' elements of the
 *                        output datatype.
 *  unsigned char *from   A pointer to an array of 'nitem' elements of the
 *                        input datatype.
 *  size_t nitem          The number of items to copy from 'to' into 'from'.
 *
 * Each element in the 'from' array must be extracted, converted to
 * the equivalent output datatype and placed in the 'to' array.
 *
 * The following macros are required to convert from native to FITS datatypes:
 *  CHRTOFIT  -  Native char to 7-bit ASCII in an 8-bit two's compliment char.
 *  BYTTOFIT  -  Native char to 8-bit two's-complement integer.
 *  SHTTOFIT  -  Native short to 16-bit big-endian two's-complement integer.
 *  INTTOFIT  -  Native int to 32-bit big-endian two's-complement integer.
 *  LNGTOFIT  -  Native long to 32-bit big-endian two's-complement integer.
 *  FLTTOFIT  -  Native float to 32-bit big-endian IEEE 754 S format floats.
 *  DBLTOFIT  -  Native float to 64-bit big-endian IEEE 754 T format floats.
 * Macros to perform the opposite conversion are also required and are named
 * as above, but with the words that precede and follow the TO part of the
 * name reversed.
 */

/*
 * The following macros are taylored to machines with the following
 * native C datatypes:
 *
 *  datatype    size in bits    Byte order       Representation
 *  char             8                           Two's complement + ASCII.
 *  short           16          Big-endian       Two's complement.
 *  int             32          Big-endian       Two's complement.
 *  long            32          Big-endian       Two's complement.
 *  float           32          Big-endian       IEEE 754 S format.
 *  double          64          Big-endian       IEEE 754 T format.
 *
 * Note that the output FITS types are identical to the input, so no
 * conversion is needed, and the conversion operation is a direct copy
 * from the input array to the output array.
 */
#if defined(sparc) || defined(hpux) || defined(ibm_aix) || defined(apple_osx)
#define CHRTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOCHR(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define BYTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOBYT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define SHTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(short))
#define FITTOSHT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(short))
#define INTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(int))
#define FITTOINT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(int))
#define LNGTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(long))
#define FITTOLNG(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(long))
#define FLTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(float))
#define FITTOFLT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(float))
#define DBLTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(double))
#define FITTODBL(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(double))

/*
 * The following macros are taylored to machines with the following
 * native C datatypes:
 *
 *  datatype    size in bits    Byte order       Representation
 *  char             8                           Two's complement + ASCII.
 *  short           16          Big-endian       Two's complement.
 *  int             32          Big-endian       Two's complement.
 *  long            64          Big-endian       Two's complement.
 *  float           32          Big-endian       IEEE 754 S format.
 *  double          64          Big-endian       IEEE 754 T format.
 *
 * Note that the output FITS types are identical to the input, so no
 * conversion is needed, and the conversion operation is a direct copy
 * from the input array to the output array.
 */
#elif defined(ibm_aix_64)
#define CHRTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOCHR(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define BYTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOBYT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define SHTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(short))
#define FITTOSHT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(short))
#define INTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(int))
#define FITTOINT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(int))
#define NEED_B8B4 1
#define LNGTOFIT(to,from,nitem) cp_b8b4((to), (from), (nitem))
#define FITTOLNG(to,from,nitem) cp_b4b8((to), (from), (nitem))
#define FLTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(float))
#define FITTOFLT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(float))
#define DBLTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(double))
#define FITTODBL(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(double))

/*
 * The following macros are tailored to machines with the following
 * native C datatypes:
 *
 *  datatype    size in bits    Byte order       Representation
 *  char             8                           Two's complement + ASCII.
 *  short           16          Little-endian    Two's complement.
 *  int             32          Little-endian    Two's complement.
 *  long            64          Little-endian    Two's complement.
 *  float           32          Little-endian    IEEE 754 S format.
 *  double          64          Little-endian    IEEE 754 T format.
 *
 * Note that apart from (64-bit long <-> 32-bit FITS integer) the only
 * operation needed by the macros is byte-reversal while copying.
 */
#elif defined(__alpha__) || defined(ia64)
#define CHRTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOCHR(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define BYTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOBYT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define NEED_2R2 1
#define SHTTOFIT(to,from,nitem) cp_2r2((to), (from), (nitem))
#define FITTOSHT(to,from,nitem) cp_2r2((to), (from), (nitem))
#define NEED_4R4 1
#define INTTOFIT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FITTOINT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define NEED_L8B4 1
#define LNGTOFIT(to,from,nitem) cp_l8b4((to), (from), (nitem))
#define FITTOLNG(to,from,nitem) cp_b4l8((to), (from), (nitem))
#define FLTTOFIT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FITTOFLT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define NEED_8R8 1
#define DBLTOFIT(to,from,nitem) cp_8r8((to), (from), (nitem))
#define FITTODBL(to,from,nitem) cp_8r8((to), (from), (nitem))

/*
 * The following macros are tailored to machines with the following
 * native C datatypes:
 *
 *  datatype    size in bits    Byte order       Representation
 *  char             8                           Two's complement + ASCII.
 *  short           16          Little-endian    Two's complement.
 *  int             32          Little-endian    Two's complement.
 *  long            32          Little-endian    Two's complement.
 *  float           32          Little-endian    IEEE 754 S format.
 *  double          64          Little-endian    IEEE 754 T format.
 *
 * Note that apart from (64-bit long <-> 32-bit FITS integer) the only
 * operation needed by the macros is byte-reversal while copying.
 */
#elif defined(linux_i486_gcc) || defined(intel_osx)
#define CHRTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOCHR(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define BYTTOFIT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define FITTOBYT(to,from,nitem) memcpy((to), (from), (nitem)*sizeof(char))
#define NEED_2R2 1
#define SHTTOFIT(to,from,nitem) cp_2r2((to), (from), (nitem))
#define FITTOSHT(to,from,nitem) cp_2r2((to), (from), (nitem))
#define NEED_4R4 1
#define INTTOFIT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FITTOINT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define LNGTOFIT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FITTOLNG(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FLTTOFIT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define FITTOFLT(to,from,nitem) cp_4r4((to), (from), (nitem))
#define NEED_8R8 1
#define DBLTOFIT(to,from,nitem) cp_8r8((to), (from), (nitem))
#define FITTODBL(to,from,nitem) cp_8r8((to), (from), (nitem))

#else
#error "Architecture not represented in fits_src/sysfits.h"
#endif

/*
 * To check compilation of all conversion functions define
 * CHECK_ALL_CONV. This should be done using a compiler switch.
 * For example: cc -DCHECK_ALL_CONV
 */
#ifdef CHECK_ALL_CONV
#define NEED_B2B4
#define NEED_B8B4
#define NEED_2R2
#define NEED_4R4
#define NEED_8R8
#define NEED_L2B4
#define NEED_L4B2
#define NEED_L8B4
#endif

/* Big-endian two's complement integers: 2-byte to/from 4-byte */

#ifdef NEED_B2B4
void cp_b2b4(unsigned char *dest, unsigned char *orig, size_t nitem);
void cp_b4b2(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* Big-endian two's complement integers: Convert 8-byte to/from 4-byte */

#ifdef NEED_B8B4
void cp_b8b4(unsigned char *dest, unsigned char *orig, size_t nitem);
void cp_b4b8(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* 2-byte datatypes: Byte-reversal performed during copying */

#ifdef NEED_2R2
void cp_2r2(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* 4-byte datatypes: Byte-reversal performed during copying */

#ifdef NEED_4R4
void cp_4r4(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* 8-byte datatypes: Byte-reversal performed during copying */

#ifdef NEED_8R8
void cp_8r8(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* Two's complement integers: Little-endian 2-byte <-> big-endian 4-byte */

#ifdef NEED_L2B4
void cp_l2b4(unsigned char *dest, unsigned char *orig, size_t nitem);
void cp_b4l2(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* Two's complement integers: Little-endian 4-byte <-> big-endian 2-byte */

#ifdef NEED_L4B2
void cp_l4b2(unsigned char *dest, unsigned char *orig, size_t nitem);
void cp_b2l4(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

/* Two's complement integers: Little-endian 8-byte <-> big-endian 4-byte */

#ifdef NEED_L8B4
void cp_l8b4(unsigned char *dest, unsigned char *orig, size_t nitem);
void cp_b4l8(unsigned char *dest, unsigned char *orig, size_t nitem);
#endif

#endif
