#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>

#include "sysfits.h"
#include "recio.h"
#include "fits.h"
#include "utils.h"

static Fits *get_Fits(Fits *fits);
static void flagblank(Fittype type, int nobj, void *data, long blank,
		      Fitsflag *flags);
static void makeblank(Fittype type, int nobj, void *data, long blank,
		      Fitsflag *flags);
static void flagnan(Fittype type, int nobj, unsigned char *data,
		    Fitsflag *flags);
static void makenan(Fittype type, int nobj, unsigned char *data,
		    Fitsflag *flags);


/* Buffer used to send data to writedata and receive data from readdata */

static unsigned char fits_buff[FITSLEN];
static long readdata(Fits *fits, Hdu *hdu, long offset, size_t size,
		     long start, long nobj, int isdata);
static int writedata(Fits *fits, Hdu *hdu, long offset, size_t size,
		     long start, long nobj, int isdata);

/*
 * Temorary scratch buffer aligned for arrays of up to CNVBUF_LEN elements of
 * any supported type. This is separately used by get_data and put_data
 * as a type conversion buffer.
 */
enum {CNVBUF_LEN=200};
static union {
  char c; short s; int i; long l; float f; double d;
} cnvbuf[CNVBUF_LEN];

/*.......................................................................
 * Open either an existing or new FITS file and create a FITS descriptor
 * for the file.
 *
 * Input:
 *  name   const char *  The name of the FITS file.
 *  isold         int    If true then an existing file should be opened.
 *                       If false, a new FITS file will be created.
 *  readonly      int    When reading an exisiting file, this requests
 *                       that the file be opened without write access.
 *  pedantic      int    If true supply extra warnings about departures
 *                       from the FITS standard.
 *  aips          int    If true then substitute the pre-standard AIPS
 *                       EPOCH keyword for EQUINOX, and AIPS A3DTABLE
 *                       type name for the standardized BINTABLE name.
 * Output:
 *  return       Fits *  The descriptor of the FITS file, or NULL
 *                       on error.
 */
Fits *new_Fits(const char *name, int isold, int readonly, int pedantic,
	       int aips)
{
  Fits *fits;    /* Pointer to the new FITS descriptor */
/*
 * Check sanity of arguments.
 */
  if(readonly && !isold) {
    fprintf(stderr, "new_Fits: Can\'t create readonly FITS file.\n");
    return NULL;
  };
/*
 * Allocate memory for the new FITS descriptor.
 */
  fits = (Fits *) malloc(sizeof(Fits));
  if(fits==0) {
    fprintf(stderr, "new_Fits: Insufficient memory for new FITS file\n");
    return fits;
  };
/*
 * NULL all pointers in the FITS descriptor to ensure safe
 * deletion by del_Fits() if an error occurs.
 */
  fits->hdu = 0;
  fits->name = 0;
/*
 * Initialize other members.
 */
  fits->readonly = readonly;
  fits->pedantic = pedantic;
  fits->aips = aips;
  fits->modified = 0;
  fits->complete = 1;
  fits->pad = 0;
  fits->recnum = -1;
  fits->nullrec = 0;
/*
 * Record the file name.
 */
  fits->name = fitsstr(name);
  if(fits->name==NULL)
    return del_Fits(fits);
/*
 * Open the file.
 */
  fits->rec = new_Recio(name, (isold ? IS_OLD : IS_NEW), readonly, FITSLEN);
  if(fits->rec==NULL)
    return del_Fits(fits);
/*
 * If the FITS file exists then determine its structure.
 */
  return isold ? get_Fits(fits) : fits;
}

/*.......................................................................
 * Close a FITS file previously opened by new_Fits() and deallocate
 * memory in the Fits descriptor.
 *
 * Input:
 *  fits     Fits *  A FITS descriptor returned by new_Fits().
 * Output:
 *  return   Fits *  Allways NULL. Use like:  fits=del_Fits(fits);
 */
Fits *del_Fits(Fits *fits)
{
  Hdu *hdu;
  Hdu *next;
/*
 * Already deleted?
 */
  if(fits==NULL)
    return fits;
/*
 * Flush any pending data to the FITS file.
 */
  fits_flush(fits);
/*
 * Emit a warning if the last HDU written is incomplete?
 */
  if(fits->hdu) {
    for(hdu=fits->hdu; hdu->next != NULL; hdu = hdu->next);
    if(hdu->state != HDU_INFILE)
      fprintf(stderr, "Warning: Last HDU of FITS file is incomplete.\n");
  };
/*
 * Close the FITS file.
 */
  fits->rec = del_Recio(fits->rec);
/*
 * Free the memory used to store the file name.
 */
  if(fits->name)
    free(fits->name);
/*
 * Delete the linked list of Hdu descriptors.
 */
  for(hdu=fits->hdu; hdu != NULL; hdu=next) {
    next = hdu->next;
    del_Hdu(hdu);
  };
/*
 * Free the Fits descriptor.
 */
  free(fits);
  return NULL;
}

/*.......................................................................
 * Ascertain the structure of a FITS file by reading relevant parts of
 * each header, either up to the end of file or to the first non-standard
 * header. An error will be returned if the first header in the file is
 * not a PRIMARY header.
 *
 * Input:
 *  fits    Fits *  A Fits descriptor previously initialised by new_Fits().
 *                  Any previous FITS structure information will be deleted.
 * Output:
 *  return  Fits *  The modified Fits descriptor or NULL if an error
 *                  occurred and the corrupt descriptor was deleted.
 */
static Fits *get_Fits(Fits *fits)
{
  Hdu *hdu;     /* The latest Hdu read */
/*
 * Get the primary hdu.
 */
  fits->hdu = get_Hdu(fits, 0);
  if(fits->hdu==NULL)
    return del_Fits(fits);
/*
 * Read all the remaining HDUs in the FITS file up to the first non-standard
 * extension or the END of file.
 */
  for(hdu=fits->hdu; (hdu->next=get_Hdu(fits, hdu->headrec+hdu->nrec)); hdu=hdu->next);
/*
 * Check that the data-segment of the final HDU is fully there by
 * attempting to read its final record, unless the data segment is
 * supposed to be empty.
 */
  if(hdu->nrec - (hdu->datarec - hdu->headrec) > 0 && 
     fits_read(fits, hdu->headrec + hdu->nrec - 1, 0)) {
    fprintf(stderr, "FITS file shorter than expected.\n");
    return del_Fits(fits);
  };
  return fits;
}

/*.......................................................................
 * Read a given header line from a given HDU.
 *
 * Input:
 *  fits   Fits * The FITS file descriptor.
 *  hdu     Hdu * Base class descriptor of the HDU to be read from.
 *  lnum    int   0-relative line number within the HDU to be read.
 * Output:
 *  return char * Pointer to internal buffer containing the line read,
 *                or NULL on error.
 */
char *rheadline(Fits *fits, Hdu *hdu, int lnum)
{
  static char lbuff[81];  /* Input line buffer */
/*
 * Is the HDU in the fits file?
 */
  if(hdu->state==HDU_DESCR) {
    fprintf(stderr, "rheadline: The cited HDU is not in any FITS file\n");
    return NULL;
  };
/*
 * Reject the line if lnum lies outside the allocated header area.
 */
  if(lnum < 0 || (hdu->endline>=0 && lnum > hdu->endline)) {
    fprintf(stderr, "rheadline: Requested header line position is not in the header.\n");
    return NULL;
  };
/*
 * Read the new header line.
 */
  if(readdata(fits, hdu, lnum*80L, (size_t) 1, 0L, 80L, 0) < 80L) {
/*
 * If an error occurs on the first line then this simply indicates that
 * there is no further HDU in the file and is not an error.
 */
    if(lnum!=0) {
      if(rec_eof(fits->rec))
	fprintf(stderr, "rheadline: Premature end of file in header.\n");
      else
	fprintf(stderr, "rheadline: Unable to read line from FITS file.\n");
    };
    return NULL;
  };
/*
 * Copy and terminate the new line.
 */
  FITTOCHR((unsigned char *)lbuff, (unsigned char *)fits_buff, (size_t) 80);
  lbuff[80] = '\0';
/*
 * Update the record of the current line number to point at the next
 * line (ie where the file pointer now points).
 */
  if(hdu->endline<0 || lnum+1 <= hdu->endline)
    hdu->nextline = lnum+1;
  return &lbuff[0];
}

/*.......................................................................
 * Write a given header line to a given HDU.
 *
 * Input:
 *  fits   Fits * The FITS file descriptor.
 *  hdu     Hdu * Base class descriptor of the HDU to be written to.
 *  lnum    int   0-relative line number within the HDU to be written.
 *  line   char * String to be written. This must point to an array of
 *                at least 81 chars and the last element line[80] must
 *                be '\0'.
 * Output:
 *  return  int   0 - OK.
 *                1 - Error.
 */
int wheadline(Fits *fits, Hdu *hdu, int lnum, char *line)
{
/*
 * Was this file opened for writing?
 */
  if(fits->readonly) {
    fprintf(stderr, "wheadline: File not opened for writing.\n");
    return 1;
  };
/*
 * Is the HDU in the fits file?
 */
  if(hdu->state==HDU_DESCR) {
    fprintf(stderr, "wheadline: HDU is not in a FITS file - use add_Hdu().\n");
    return 1;
  };
/*
 * Is this the END line?
 */
  if(strncmp(line, "END     ", 8)==0) {
    if(hdu->endline!=lnum || lnum<0) {
      fprintf(stderr, "wheadline: Attempt to write misplaced END line.\n");
      return 1;
    };
/*
 * Reject the line if lnum lies outside the allocated header area.
 */
  } else if(lnum < 0 || lnum > hdu->endline) {
    fprintf(stderr, "wheadline: Out of bounds header line number - outside of the FITS header.\n");
    fprintf(stderr, "wheadline: Rejecting header line: \"%.30s...\"\n", line);
    return 1;
  } else if(lnum == hdu->endline) {
/*
 * Would this over-write the END line?
 */
    if(hdu->state==HDU_INFILE) {
      fprintf(stderr, "wheadline: No room for new header line before END line.\n");
      fprintf(stderr, "wheadline: Rejecting header line: \"%.30s...\"\n", line);
      return 1;
    };
/*
 * Move the end line to make way for the new line.
 */
    hdu->endline++;
  };
/*
 * Copy the header line into the FITS I/O buffer.
 */
  CHRTOFIT((unsigned char *)fits_buff, (unsigned char *)line, (size_t) 80);
/*
 * Write the header line.
 */
  if(writedata(fits, hdu, lnum*80L, (size_t) 1, 0L, 80L, 0)) {
    fprintf(stderr, "wheadline: Unable to write header line to FITS file.\n");
    return 1;
  };
/*
 * Update the record of the current line number to point at the next
 * line (ie where the file pointer now points).
 */
  hdu->wnxtline = lnum<hdu->endline ? lnum+1:hdu->endline;
  return 0;
}

/*.......................................................................
 * Make a dynamically allocated copy a (\0 terminated) string and return
 * either a pointer to the copy or NULL on memory allocation error.
 *
 * NB. Trailing spaces are not copied.
 *
 * Input:
 *  str    char *   Input string to be copied. Must be '\0' terminated.
 *                  If NULL, a NULL pointer will be returned without
 *                  soliciting an error message.
 * Output:
 *  return char *   Pointer to the dynamically allocated copy of 'str',
 *                  or NULL on error.
 */
char *fitsstr(const char *str)
{
  static size_t slen;       /* Length of string up to first trailing space */
  static const char *cptr;  /* Pointer into str */
  static const char *keep;  /* Place marker in str */
  static char *newstr;      /* Pointer to the copy of 'str' */
/*
 * No string to be allocated?
 */
  if(str==NULL)
    return NULL;
/*
 * Determine the length of the string up to the first trailing space, or
 * '\0' terminator.
 */
  for(keep=cptr=str; *cptr; cptr++) {
    if(*cptr!=' ')
      keep = cptr;
  };
/*
 * Record the required string length.
 */
  slen = (keep - str) + 1;
/*
 * Allocate the new string.
 */
  newstr = malloc(slen + 1);
  if(newstr==NULL) {
    fprintf(stderr, "fitsstr: Insufficient memory to copy string %s\n", str);
  } else {
    strncpy(newstr, str, slen);
    newstr[slen] = '\0';
  };
  return newstr;
}

/*.......................................................................
 * Return the size of a given FITS data-type in 8-bit bytes.
 *
 * Input:
 *  type   Fittype  The FITS data-type.
 * Output:
 *  return  size_t  The number of 8-bit bytes in the data-type, or 0
 *                  on error.
 */
size_t typesize(Fittype type)
{
  switch(type) {
  case DAT_NON: default:   /* void - unknown */
    fprintf(stderr, "typesize: Unrecognised data-type intercepted\n");
    return 0;
  case DAT_SHT:   /* (short) */
    return 2;
  case DAT_INT:   /* (int) */
  case DAT_LNG:   /* (long) */
    return 4;
  case DAT_FLT:   /* (float) */
    return 4;
  case DAT_DBL:   /* (double) */
    return 8;
  case DAT_CHR:   /* (char) */
    return 1;
  case DAT_BYT:   /* (unsigned char) representation of byte. */
    return 1;
  case DAT_BIT:   /* (unsigned char) bit array */
    return 1;
  case DAT_LOG:   /* (char) representation of FITS logical 'T' or 'F' */
    return 1;
  case DAT_SCMP:  /* (float)[2] representation of complex number */
    return 8;
  case DAT_DCMP:  /* (double)[2] representation of complex number */
    return 16;
  };
}

/*.......................................................................
 * Return the size of the equivalent machine dependant FITS data-type in
 * chars.
 *
 * Input:
 *  type   Fittype  The FITS data-type.
 * Output:
 *  return  size_t  The number of 8-bit bytes in the data-type, or 0
 *                  on error.
 */
size_t machsize(Fittype type)
{
  switch(type) {
  case DAT_NON:
    return 0;
  case DAT_SHT:   /* (short) */
    return sizeof(short);
  case DAT_INT:   /* (int) */
    return sizeof(int);
  case DAT_LNG:   /* (long) */
    return sizeof(long);
  case DAT_FLT:   /* (float) */
    return sizeof(float);
  case DAT_DBL:   /* (double) */
    return sizeof(double);
  case DAT_CHR:   /* (char) */
    return sizeof(char);
  case DAT_BYT:   /* (unsigned char) representation of byte. */
    return sizeof(unsigned char);
  case DAT_BIT:   /* (unsigned char) bit array */
    return sizeof(unsigned char);
  case DAT_LOG:   /* (char) representation of FITS logical 'T' or 'F' */
    return sizeof(char);
  case DAT_SCMP:  /* (float)[2] representation of complex number */
    return sizeof(float [2]);
  case DAT_DCMP:  /* (double)[2] representation of complex number */
    return sizeof(double [2]);
  default:   /* void - unknown */
    fprintf(stderr, "machsize: Unrecognised data-type intercepted\n");
    return 0;
  };
}

/*.......................................................................
 * Return the name of a given FITS data-type.
 *
 * Input:
 *  type   Fittype    The FITS data-type.
 * Output:
 *  return    char *  The name of the data-type or "(unknown)" on error.
 *                    NB. Pointers to internal static strings are returned.
 */
char *typename(Fittype type)
{
  switch(type) {
  case DAT_NON:   /* No value */
    return "no value";
  case DAT_SHT:   /* (short) */
    return "short";
  case DAT_INT:   /* (int) */
    return "int";
  case DAT_LNG:   /* (long) */
    return "long";
  case DAT_FLT:   /* (float) */
    return "float";
  case DAT_DBL:   /* (double) */
    return "double";
  case DAT_CHR:   /* (char) */
    return "char";
  case DAT_BYT:   /* (unsigned char) representation of byte. */
    return "byte";
  case DAT_BIT:   /* (unsigned char) bit array */
    return "bit";
  case DAT_LOG:   /* (char) representation of FITS logical 'T' or 'F' */
    return "logical";
  case DAT_SCMP:  /* (float)[2] representation of complex number */
    return "float-complex";
  case DAT_DCMP:  /* (double)[2] representation of complex number */
    return "double-complex";
  case DAT_COM:
    return "comment-string";
  case DAT_STR:
    return "string";
  default:   /* Unknown */
    fprintf(stderr, "typename: Unrecognised data-type intercepted\n");
    return "(unknown)";
  };
}

/*.......................................................................
 * Read an array of binary values from the data segment of an HDU.
 *
 * Input:
 *  fits      Fits *  The FITS file descriptor.
 *  hdu        Hdu *  The descriptor of the HDU to be read from.
 *  offset    long    The offset into the data-segment of the HDU to the
 *                    start of the array (measured in FITS 8-bit bytes).
 *  atype  Fittype    The data-type of the data in the FITS file.
 *  start     long    The start element (0-relative) in the array.
 *  nobj      long    The number of elements to be read.
 *  btype  Fittype    The data-type to convert to.
 *  zero    double    The zero offset to apply to arithmetic types.
 *  scale   double    The scale factor to apply to arithmetic types.
 *  os     Offscal *  If not NULL then os[] should be an array of
 *                    'nobj' different scale-factors and offsets to be
 *                    applied to the corresponding elements in the
 *                    input stream. If os==NULL then the above scalar
 *                    values of 'zero' and 'scale' will be used instead.
 *  flags Fitsflag *  If not NULL, then this should point to an array
 *                    of nobj elements to be used to record which elements
 *                    of data[] are blanked 1=blanked, 0=OK.
 *  blank     long    The value to assign to null integral fields.
 * Output:
 *  data      void *  An array of type 'btype' of at least 'nobj' elements.
 *                    The returned data in this array will have had any
 *                    system dependant conversions performed from the
 *                    FITS binary representation.
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int get_data(Fits *fits, Hdu *hdu, long offset, Fittype atype, long start,
	     long nobj, Fittype btype, double zero, double scale, Offscal *os,
	     Fitsflag *flags, long blank, void *data)
{
  size_t fsize;  /* Size of FITS data-type in 8-bit bytes */
  size_t msize;  /* Size of equivalent data-type on this machine (chars). */
  size_t bsize;  /* Size of an element of 'data[]' (8-bit bytes) */
  long nread;    /* Number of objects buffered from latest read */
  long ndone;    /* Number of objects so far read and processed */
  long nreq;     /* Number of elements requested in the next read */
  char *datptr;  /* Pointer into 'data[]' */
/*
 * Get the size of the data-type in FITS 8-bit bytes and size of the
 * equivalent machine dependant type, ans the size of an element of 'data[]'.
 */
  fsize = typesize(atype);
  msize = machsize(atype);
  bsize = machsize(btype);
  if(fsize==0 || msize==0 || bsize==0)
    return 1;
/*
 * Read data into the FITS buffer array.
 */
  nread = 0;
  datptr = (char *) data;
  for(ndone=0; ndone<nobj; ndone += nread, datptr += nread * bsize) {
/*
 * Read as many objects as will fit in the FITS buffer in fits_buff
 * and in the conversion buffer cnvbuf[].
 */
    nreq = nobj - ndone;   /* The number of elements remaining to be read */
    if(nreq > CNVBUF_LEN)
      nreq = CNVBUF_LEN;
    nread = readdata(fits, hdu, offset, fsize, start+ndone, nreq, 1);
    if(nread==0)
      return 1;
/*
 * Copy into the output array while performing any necessary conversions
 * between FITS types and the type in (void *data).
 */
    switch(atype) {
    case DAT_SHT:   /* (short) */
      FITTOSHT((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      if(flags)	flagblank(atype, nread, cnvbuf, blank, &flags[ndone]);
      break;
    case DAT_INT:   /* (int) */
      FITTOINT((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      if(flags)	flagblank(atype, nread, cnvbuf, blank, &flags[ndone]);
      break;
    case DAT_LNG:   /* (long) */
      FITTOLNG((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      if(flags)	flagblank(atype, nread, cnvbuf, blank, &flags[ndone]);
      break;
    case DAT_FLT:   /* (float) */
      if(flags) flagnan(atype, nread, fits_buff, &flags[ndone]);
      FITTOFLT((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      break;
    case DAT_DBL:   /* (double) */
      if(flags) flagnan(atype, nread, fits_buff, &flags[ndone]);
      FITTODBL((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      break;
    case DAT_CHR: case DAT_LOG:  /* (char) ASCII */
      FITTOCHR((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      break;
    case DAT_BYT: case DAT_BIT:  /* (unsigned char) representation of byte. */
      FITTOBYT((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      if(flags)	flagblank(atype, nread, cnvbuf, blank, &flags[ndone]);
      break;
    case DAT_SCMP:  /* (float)[2] representation of complex number */
      if(flags) flagnan(atype, nread, fits_buff, &flags[ndone]);
      FITTOFLT((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      break;
    case DAT_DCMP:  /* (double)[2] representation of complex number */
      if(flags) flagnan(atype, nread, fits_buff, &flags[ndone]);
      FITTODBL((unsigned char *)cnvbuf, (unsigned char *)fits_buff, nread);
      break;
    default: case DAT_NON:
      fprintf(stderr, "get_data: Don't know how to read type: %s\n",
	      typename(atype));
      return 1;
    };
/*
 * Copy the data into the output array. While doing this, apply
 * 'zero' and 'scale' or os[*].off and os[*].mul, and convert to the
 * requested output type. 
 */
    if(os ?
       arrconv(nread, atype, cnvbuf, &os[ndone], btype, datptr) :
       typeconv(nread, atype, cnvbuf, zero, scale, btype, datptr))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Write an array of binary values to the data segment of an HDU.
 *
 * Input:
 *  fits     Fits *  The FITS file descriptor.
 *  hdu       Hdu *  The descriptor of the HDU to be written to.
 *  offset   long    The offset into the data-segment of the HDU to the
 *                   start of the array (measured in FITS 8-bit bytes).
 *  atype Fittype    The data-type of the data in the FITS file.
 *  start    long    The start element (0-relative) in the array.
 *  nobj     long    The number of elements to be written.
 *  btype Fittype    The data-type to convert from.
 *  zero   double    The offset to apply to arithmetic types. Note that
 *                   this should have the opposite sign to the value
 *                   sent to a matching call to get_data().
 *  scale  double    The scale factor to apply to arithmetic types.
 *                   Note that this should be the reciprocal of the value
 *                   sent to a matching call to get_data().
 *  os    Offscal *  If not NULL then os[] should be an array of
 *                   'nobj' different scale-factors and offsets to be
 *                   applied to the corresponding elements in the
 *                   output stream. If os==NULL then the above scalar
 *                   values of 'zero' and 'scale' will be used instead.
 *  flags Fitsflag * If not NULL, then this should point to an array
 *                   of nobj elements specifying which elements
 *                   of data[] are to be blanked 1=blanked, 0=OK.
 *  blank    long    The value to assign to null integral fields.
 *  data     void *  An array of type 'btype' of at least 'nobj' elements
 *                   containing the data to be written.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
int put_data(Fits *fits, Hdu *hdu, long offset, Fittype atype, long start,
	     long nobj, Fittype btype, double zero, double scale, Offscal *os,
	     Fitsflag *flags, long blank, void *data)
{
  size_t fsize;  /* Size of FITS data-type in 8-bit bytes */
  size_t msize;  /* Size of equivalent data-type on this machine (chars). */
  size_t bsize;  /* Size of an element of 'data[]' (8-bit bytes) */
  long ndone;    /* Number of objects so far written and processed */
  long nnew;     /* Number of elements to be processed next */
  char *datptr;  /* Pointer into 'data[]' */
/*
 * Get the size of the data-type in FITS 8-bit bytes and size of the
 * equivalent machine dependant type, ans the size of an element of 'data[]'.
 */
  fsize = typesize(atype);
  msize = machsize(atype);
  bsize = machsize(btype);
  if(fsize==0 || msize==0 || bsize==0)
    return 1;
/*
 * Write data via the FITS buffer array.
 */
  datptr = (char *) data;
  for(ndone=0; ndone<nobj; ndone += nnew, datptr += nnew * bsize) {
/*
 * Determine the number of objects that will fit in the FITS buffer in
 * fits_buff and in the conversion buffer cnvbuf[].
 */
    nnew = nobj - ndone;   /* The number of elements remaining to be written */
    if(nnew > CNVBUF_LEN)
      nnew = CNVBUF_LEN;
/*
 * Copy the data into the output array. While doing this, apply
 * 'zero' and 'scale' or os[*].off and os[*].mul and convert to the
 * requested output type. 
 */
    if(os ?
       arrconv(nnew, btype, datptr, &os[ndone], atype, cnvbuf) :
       typeconv(nnew, btype, datptr, zero, scale, atype, cnvbuf))
      return 1;
/*
 * Copy into the fits buffer array while performing any necessary conversions
 * between FITS types and the machine type in (void *data).
 */
    switch(atype) {
    case DAT_SHT:   /* (short) */
      if(flags) makeblank(atype, nnew, cnvbuf, blank, &flags[ndone]);
      SHTTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      break;
    case DAT_INT:   /* (int) */
      if(flags)	makeblank(atype, nnew, cnvbuf, blank, &flags[ndone]);
      INTTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      break;
    case DAT_LNG:   /* (long) */
      if(flags)	makeblank(atype, nnew, cnvbuf, blank, &flags[ndone]);
      LNGTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      break;
    case DAT_FLT:   /* (float) */
      FLTTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      if(flags) makenan(atype, nnew, fits_buff, &flags[ndone]);
      break;
    case DAT_DBL:   /* (double) */
      DBLTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      if(flags) makenan(atype, nnew, fits_buff, &flags[ndone]);
      break;
    case DAT_CHR: case DAT_LOG:  /* (char) ASCII */
      CHRTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      break;
    case DAT_BYT: case DAT_BIT:  /* (unsigned char) representation of byte. */
      if(flags)	makeblank(atype, nnew, cnvbuf, blank, &flags[ndone]);
      BYTTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      break;
    case DAT_SCMP:  /* (float)[2] representation of complex number */
      FLTTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      if(flags) makenan(atype, nnew, fits_buff, &flags[ndone]);
      break;
    case DAT_DCMP:  /* (double)[2] representation of complex number */
      DBLTOFIT((unsigned char *)fits_buff, (unsigned char *)cnvbuf, nnew);
      if(flags) makenan(atype, nnew, fits_buff, &flags[ndone]);
      break;
    default: case DAT_NON:
      fprintf(stderr, "put_data: Don't know how to write type: %s\n",
	      typename(atype));
      return 1;
    };
/*
 * Write the used part of the fits buffer to the FITS file.
 */
    if(writedata(fits, hdu, offset, fsize, start+ndone, nnew, 1))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * Convert an array of one data-type to one of another.
 *
 * Input:
 *  ndata      long    The number of elements in 'adata[]' and 'bdata[]'.
 *  atype   Fittype    The data-type of the input array.
 *  adata      void *  The input array of dimension 'ndata' and type 'atype'.
 *  zero     double    Arithmetic data-types are offset by this value
 *                     during conversion.
 *  scal     double    Arithmetic data-types are scaled by this value
 *                     during conversion.
 *  btype   Fittype    The data-type of the output array (bdata).
 * Output:
 *  bdata      void *  The output array of dimension 'ndata' and type 'btype'.
 *  return      int    0 - OK.
 *                     1 - Error, eg. Unhandled conversion.
 */
int typeconv(long ndata, Fittype atype, void *adata, double zero, double scal,
	     Fittype btype, void *bdata)
{
  static int i;   /* Index into adata[] and bdata[] */
  int ierr = 0;   /* Error status */
  switch(atype) {
  case DAT_BYT:      /* Treat byte type as a small int */
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
       ((unsigned char *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((unsigned char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_SHT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((short *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_INT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((int *)adata)[i];
      break;
    case DAT_LOG:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((int *)adata)[i] ? 'T' : 'F';
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_LNG:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((long *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_FLT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((float *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_DBL:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = zero + scal * ((double *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_CHR:
    switch(btype) {
    case DAT_CHR:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_BIT:
    switch(btype) {
    case DAT_BIT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = ((unsigned char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_LOG:
    switch(btype) {
    case DAT_LOG:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((char *)adata)[i];
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = ((char *)adata)[i]=='T';
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_SCMP:
    switch(btype) {
    case DAT_SCMP:
      for(i=0;i<ndata;i++) {
	((float *)bdata)[2*i]   = zero + scal * ((float *)adata)[2*i];
	((float *)bdata)[2*i+1] = scal * ((float *)adata)[2*i+1];
      };
      break;
    case DAT_DCMP:
      for(i=0;i<ndata;i++) {
	((double *)bdata)[2*i]   = zero + scal * ((float *)adata)[2*i];
	((double *)bdata)[2*i+1] = scal * ((float *)adata)[2*i+1];
      };
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_DCMP:
    switch(btype) {
    case DAT_SCMP:
      for(i=0;i<ndata;i++) {
	((float *)bdata)[2*i]   = zero + scal * ((double *)adata)[2*i];
	((float *)bdata)[2*i+1] = scal * ((double *)adata)[2*i+1];
      };
      break;
    case DAT_DCMP:
      for(i=0;i<ndata;i++) {
	((double *)bdata)[2*i]   = zero + scal * ((double *)adata)[2*i];
	((double *)bdata)[2*i+1] = scal * ((double *)adata)[2*i+1];
      };
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_COM:
    switch(btype) {
    case DAT_STR: case DAT_COM:
      for(i=0;i<ndata;i++)
	((char **)bdata)[i] = ((char **)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_STR:
    switch(btype) {
    case DAT_COM: case DAT_STR:
      for(i=0;i<ndata;i++)
	((char **)bdata)[i] = ((char **)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  default:
    ierr = 1;
  };
/*
 * If ierr is set then the conversion was not one of the above legal
 * conversions.
 */
  if(ierr)
    fprintf(stderr, "typeconv: Unhandled conversion from (%s) to (%s)\n",
	    typename(atype), typename(btype));
  return ierr;
}

/*.......................................................................
 * Convert an array of one data-type to a scaled and offset array of
 * another a different data-type. This version of typeconv() allows for
 * a different scale-factor and offset for each element being converted.
 *
 * Input:
 *  ndata      long    The number of elements in 'adata[]' and 'bdata[]'.
 *  atype   Fittype    The data-type of the input array.
 *  adata      void *  The input array of dimension 'ndata' and type 'atype'.
 *  Offscal      os *  An array of ndata offsets and scale factors to apply
 *                     to corresponding elements of adata[] during conversion.
 *  btype   Fittype    The data-type of the output array (bdata).
 * Output:
 *  bdata      void *  The output array of dimension 'ndata' and type 'btype'.
 *  return      int    0 - OK.
 *                     1 - Error, eg. Unhandled conversion.
 */
int arrconv(long ndata, Fittype atype, void *adata, Offscal *os,
	    Fittype btype, void *bdata)
{
  static int i;   /* Index into adata[] and bdata[] */
  int ierr = 0;   /* Error status */
  switch(atype) {
  case DAT_BYT:      /* Treat byte type as a small int */
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
       ((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((unsigned char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_SHT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((short *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_INT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((int *)adata)[i];
      break;
    case DAT_LOG:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((int *)adata)[i] ? 'T' : 'F';
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_LNG:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((long *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_FLT:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((float *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_DBL:
    switch(btype) {
    case DAT_BYT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    case DAT_SHT:
      for(i=0;i<ndata;i++)
	((short *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    case DAT_LNG:
      for(i=0;i<ndata;i++)
	((long *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    case DAT_FLT:
      for(i=0;i<ndata;i++)
	((float *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    case DAT_DBL:
      for(i=0;i<ndata;i++)
	((double *)bdata)[i] = os[i].off + os[i].mul * ((double *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_CHR:
    switch(btype) {
    case DAT_CHR:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_BIT:
    switch(btype) {
    case DAT_BIT:
      for(i=0;i<ndata;i++)
	((unsigned char *)bdata)[i] = ((unsigned char *)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_LOG:
    switch(btype) {
    case DAT_LOG:
      for(i=0;i<ndata;i++)
	((char *)bdata)[i] = ((char *)adata)[i];
    case DAT_INT:
      for(i=0;i<ndata;i++)
	((int *)bdata)[i] = ((char *)adata)[i]=='T';
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_SCMP:
    switch(btype) {
    case DAT_SCMP:
      for(i=0;i<ndata;i++) {
	((float *)bdata)[2*i]   = os[i].off + os[i].mul * ((float *)adata)[2*i];
	((float *)bdata)[2*i+1] = os[i].mul * ((float *)adata)[2*i+1];
      };
      break;
    case DAT_DCMP:
      for(i=0;i<ndata;i++) {
	((double *)bdata)[2*i]   = os[i].off + os[i].mul * ((float *)adata)[2*i];
	((double *)bdata)[2*i+1] = os[i].mul * ((float *)adata)[2*i+1];
      };
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_DCMP:
    switch(btype) {
    case DAT_SCMP:
      for(i=0;i<ndata;i++) {
	((float *)bdata)[2*i]   = os[i].off + os[i].mul * ((double *)adata)[2*i];
	((float *)bdata)[2*i+1] = os[i].mul * ((double *)adata)[2*i+1];
      };
      break;
    case DAT_DCMP:
      for(i=0;i<ndata;i++) {
	((double *)bdata)[2*i]   = os[i].off + os[i].mul * ((double *)adata)[2*i];
	((double *)bdata)[2*i+1] = os[i].mul * ((double *)adata)[2*i+1];
      };
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_COM:
    switch(btype) {
    case DAT_STR: case DAT_COM:
      for(i=0;i<ndata;i++)
	((char **)bdata)[i] = ((char **)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  case DAT_STR:
    switch(btype) {
    case DAT_COM: case DAT_STR:
      for(i=0;i<ndata;i++)
	((char **)bdata)[i] = ((char **)adata)[i];
      break;
    default:
      ierr = 1;
    };
    break;
  default:
    ierr = 1;
  };
/*
 * If ierr is set then the conversion was not one of the above legal
 * conversions.
 */
  if(ierr)
    fprintf(stderr, "arrconv: Unhandled conversion from (%s) to (%s)\n",
	    typename(atype), typename(btype));
  return ierr;
}

/*.......................................................................
 * Given an array of flags and an equally sized array of values, set
 * the flag corresponding to each NaN value in the data array.
 *
 * Input:
 *  type       Fittype    The type of object in data[]. NB. Only the 4
 *                        floating-point and complex types are supported.
 *  nobj           int    The dimension of the 'flags' and the number of
 *                        objects in 'data'.
 *  data unsigned char *  A byte array containing 'nobj' objects of size
 *                        'size' bytes.
 * Output:
 *  flags     Fitsflag *  Send an array of 'nobj' flags. This will be
 *                        returned containing a 1 for each null value found.
 */
static void flagnan(Fittype type, int nobj, unsigned char *data,
		    Fitsflag *flags)
{
  int i;          /* Index into data[] */
  switch(type) {
  case DAT_FLT:
    for(i=0;i<nobj;i++)
      flags[i] = *data++ >= 0x7F && *data++ == 0xFF &&
	         *data++ == 0xFF && *data++ >= 0xFF;
    break;
  case DAT_DBL:
    for(i=0;i<nobj;i++)
      flags[i] = *data++ >= 0x7F && *data++ == 0xFF &&
	         *data++ == 0xFF && *data++ == 0xFF &&
	         *data++ == 0xFF && *data++ == 0xFF &&
	         *data++ == 0xFF && *data++ >= 0xFF;
    break;
  case DAT_SCMP:
    for(i=0;i<nobj;i++)
      flags[i] = (*data++ >= 0x7F && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ >= 0xFF)
	         ||
		 (*data++ >= 0x7F && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ >= 0xFF);
    break;
  case DAT_DCMP:
    for(i=0;i<nobj;i++)
      flags[i] = (*data++ >= 0x7F && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ >= 0xFF)
	         ||
		 (*data++ >= 0x7F && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ == 0xFF &&
	          *data++ == 0xFF && *data++ >= 0xFF);
    break;
  default:
    break;
  };
  return;
}

/*.......................................................................
 * Given an array of flags and an equally sized array of values, set
 * flagged values to IEEE NaNs.
 *
 * Input:
 *  type       Fittype    The type of object in data[]. NB. Only the 4
 *                        floating-point and complex types are supported.
 *  nobj           int    The dimension of the 'flags' and the number of
 *                        objects in 'data'.
 *  data unsigned char *  A byte array containing 'nobj' objects of size
 *                        'size' bytes.
 *  flags     Fitsflag *  Send an array of 'nobj' flags. This will be
 *                        returned containing a 1 for each null value found.
 */
static void makenan(Fittype type, int nobj, unsigned char *data,
		    Fitsflag *flags)
{
  int i;          /* Index into data[] */
  switch(type) {
  case DAT_FLT:
    for(i=0;i<nobj;i++) {
      if(flags[i]) {
	*data++ = 0x7F; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
      } else {
	data += 4;
      };
    };
    break;
  case DAT_DBL:
    for(i=0;i<nobj;i++) {
      if(flags[i]) {
	*data++ = 0x7F; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
	*data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
      } else {
	data += 8;
      };
    };
    break;
  case DAT_SCMP:
    for(i=0;i<nobj;i++) {
      if(flags[i]) {
	*data++ = 0x7F; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
	data += 4;
      } else {
	data += 8;
      };
    };
    break;
  case DAT_DCMP:   /* Only set the real part to NaN */
    for(i=0;i<nobj;i++) {
      if(flags[i]) {
	*data++ = 0x7F; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
	*data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF; *data++ = 0xFF;
	data += 8;
      } else {
	data += 16;
      };
    };
    break;
  default:
    break;
  };
  return;
}

/*.......................................................................
 * Given an array of flags and an equally sized array of values, set
 * the flag corresponding to each blank value in the data array.
 *
 * Input:
 *  type       Fittype    The type of object in data[]. NB. Only the 4
 *                        integral types are supported.
 *  nobj           int    The dimension of the 'flags' and the number of
 *                        objects in 'data'.
 *  data          void *  An array containing 'nobj' objects of size
 *                        'size' bytes.
 *  blank         long    The integer value of a binary-table TNULL keyword
 *                        or of an IMAGE or primary-header BLANK keyword
 *                        value.
 * Output:
 *  flags     Fitsflag *  Send an array of 'nobj' flags. This will be
 *                        returned containing a 1 for each null value found.
 */
static void flagblank(Fittype type, int nobj, void *data, long blank,
		      Fitsflag *flags)
{
  int i;          /* Index into data[] */
  switch(type) {
  case DAT_BYT:
    for(i=0;i<nobj;i++)
      flags[i] = ((unsigned char *)data)[i] == blank;
    break;
  case DAT_SHT:
    for(i=0;i<nobj;i++)
      flags[i] = ((short *)data)[i] == blank;
    break;
  case DAT_INT:
    for(i=0;i<nobj;i++)
      flags[i] = ((int *)data)[i] == blank;
    break;
  case DAT_LNG:
    for(i=0;i<nobj;i++)
      flags[i] = ((long *)data)[i] == blank;
    break;
  default:
    break;
  };
  return;
}

/*.......................................................................
 * Given an array of flags and an equally sized array of integral values,
 * set flagged values to a given integer used to denote blank.
 *
 * Input:
 *  type       Fittype    The type of object in data[]. NB. Only the 4
 *                        integral types are supported.
 *  nobj           int    The dimension of the 'flags' and the number of
 *                        objects in 'data'.
 *  data          void *  An array containing 'nobj' objects of size
 *                        'size' bytes.
 *  blank         long    The number to assign as a blank value.
 *  flags     Fitsflag *  Send an array of 'nobj' flags. This will be
 *                        returned containing a 1 for each null value found.
 */
static void makeblank(Fittype type, int nobj, void *data, long blank,
		      Fitsflag *flags)
{
  int i;          /* Index into data[] */
  switch(type) {
  case DAT_BYT:
    for(i=0; i<nobj; i++)
      if(flags[i]) ((unsigned char *)data)[i] = blank;
    break;
  case DAT_SHT:
    for(i=0; i<nobj; i++)
      if(flags[i]) ((short *)data)[i] = blank;
    break;
  case DAT_INT:
    for(i=0; i<nobj; i++)
      if(flags[i]) ((int *)data)[i] = blank;
    break;
  case DAT_LNG:
    for(i=0; i<nobj; i++)
      if(flags[i]) ((long *)data)[i] = blank;
    break;
  default:
    break;
  };
  return;
}

/*.......................................................................
 * Compare two strings.
 *
 * Input:
 *  sa     char *   The first of the two strings to be compared.
 *  sb     char *   The string to compare against.
 *  fixlen  int     If > 0, then this defines the max number of
 *                  characters to be compared. This enables searching
 *                  using prefixes.
 * Output:
 *  return  int     1  sa == sb.
 *                  0  sa != sb.
 */
int matchstr(const char *sa, const char *sb, int fixlen)
{
  return (fixlen>0 ? strncmp(sa, sb, fixlen) : strcmp(sa, sb))==0;
}

/*.......................................................................
 * Read a given number of FITS 8-bit bytes into fits_buff and return the
 * number of objects read.
 *
 * Input:
 *  fits    Fits *   The fits descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to read data from.
 *  offset  long     The offset into the respetive segment of the HDU to
 *                   the start of the array (measured in FITS 8-bit bytes).
 *  size  size_t     The size of a data-type element in FITS 8-bit bytes.
 *                   Use the return value from typesize().
 *  start   long     The start element (0-relative) in the array.
 *  nobj    long     The number of elements to be read.
 *                   If 0 then only the seek will be performed.
 *  isdata   int     If 1, then offset is wrt the start of the data
 *                   segment. If 0 then offset is wrt header segment.
 * Output:
 *  return  long     The number of objects succesfully read and placed
 *                   in fits_buff. NB. This may be less than nobj and
 *                   more than one call may be required to read all
 *                   elements of an array.
 */
static long readdata(Fits *fits, Hdu *hdu, long offset, size_t size,
		     long start, long nobj, int isdata)
{
  static size_t nread;    /* Number of chars in buffer */
  static size_t ndata;    /* Number of chars to read on a single pass */
  static size_t nbytes;   /* Total number of bytes to be read */
  static size_t recnum;   /* The required record number */
  static long recoff;     /* Record offset */
  static long bytoff;     /* Byte offset into record */
/*
 * Is the given HDU ready for use?
 */
  if(hdu->state==HDU_DESCR) {
    fprintf(stderr, "readdata: HDU is not in a FITS file - use add_Hdu().\n");
    return 0L;
  };
  if(isdata && hdu->state<HDU_DATA) {
    fprintf(stderr, "readdata: HDU not ready for data access - use end_header().\n");
    return 0L;
  };
/*
 * Determine how many FITS 8-bit bytes are to be read.
 */
  nbytes = (size * nobj * 8) / CHAR_BIT;
/*
 * Calculate the number of complete records + number of bytes corresponding
 * to the start of the array.
 */
  bytoff = offset+start*size;
  recoff = bytoff / FITSLEN;
  bytoff -= recoff * FITSLEN;
/*
 * Work out the start record number.
 */
  recnum = ((isdata) ? hdu->datarec : hdu->headrec) + recoff;
/*
 * Read data until all objects have been read or the output buffer is full.
 */
  for(nread=0; nread < nbytes && nread<FITSLEN; recnum++) {
/*
 * If this record isn't currently in the fits->buff buffer then flush
 * the previous record and read the new fits record into fits->buff.
 */
    if(fits->recnum != recnum && fits_read(fits, recnum, 1))
      break;
/*
 * How many bytes do we still require?
 */
    ndata = nbytes - nread;
/*
 * Apply buffer size limits on the number of points to be read on this
 * pass.
 */
    if(bytoff+ndata > FITSLEN)
      ndata = FITSLEN - bytoff;
/*
 * Copy the required data into the output buffer.
 */
    memcpy(&fits_buff[nread], &fits->buff[bytoff], ndata);
    nread += ndata;
    bytoff = 0L;
  };
/*
 * How many complete objects were read?
 */
  return (nread * CHAR_BIT)/(8*size);
}

/*.......................................................................
 * Write a given number of FITS 8-bit bytes from fits_buff into a FITS
 * file.
 *
 * Input:
 *  fits    Fits *   The fits descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to read data from.
 *  offset  long     The offset into the respetive segment of the HDU to
 *                   the start of the array (measured in FITS 8-bit bytes).
 *  size  size_t     The size of a data-type element in FITS 8-bit bytes.
 *                   Use the return value from typesize().
 *  start   long     The start element (0-relative) in the array.
 *  nobj    long     The number of elements to be read.
 *                   If 0 then only the seek will be performed.
 *  isdata   int     If 1, then offset is wrt the start of the data
 *                   segment. If 0 then offset is wrt header segment.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static int writedata(Fits *fits, Hdu *hdu, long offset, size_t size,
		     long start, long nobj, int isdata)
{
  static size_t nwrit;    /* Number of chars so far written */
  static size_t ndata;    /* Number of chars to read on a single pass */
  static size_t nbytes;   /* Total number of chars to be read */
  static long recnum;     /* The required record number */
  static long recoff;     /* Record offset */
  static size_t bytoff;   /* Byte offset into record */
/*
 * Was this file opened for writing?
 */
  if(fits->readonly) {
    fprintf(stderr, "writedata: File not opened for writing\n");
    return 1;
  };
/*
 * Is the given HDU ready for use?
 */
  if(hdu->state==HDU_DESCR) {
    fprintf(stderr, "writedata: HDU is not in a FITS file - use add_Hdu().\n");
    return 1;
  };
  if(isdata && hdu->state<HDU_DATA) {
    fprintf(stderr, "writedata: HDU not ready for data access - use end_header().\n");
    return 1;
  };
/*
 * Determine how many FITS 8-bit bytes are to be written.
 */
  nbytes = (size * nobj * 8) / CHAR_BIT;
/*
 * Calculate the number of complete records + number of bytes corresponding
 * to the start of the array.
 */
  bytoff = offset+start*size;
  recoff = bytoff / FITSLEN;
  bytoff -= recoff * FITSLEN;
/*
 * Work out the start record number.
 */
  recnum = ((isdata) ? hdu->datarec : hdu->headrec) + recoff;
/*
 * Write data until all objects have been written.
 */
  for(nwrit=0; nwrit < nbytes && nwrit<FITSLEN; recnum++) {
/*
 * If this record isn't currently in the fits->buff buffer then flush
 * the previous record and read the new fits record into fits->buff.
 * Note that an actual read is not performed when writing to un-written
 * records (in this case fits->complete is 0 and recnum >= fits->nullrec).
 */
    if(fits->recnum != recnum && fits_read(fits, recnum, 1))
      break;
/*
 * How many bytes are still to be written?
 */
    ndata = nbytes - nwrit;
/*
 * Apply buffer size limits on the number of points to be written on this
 * pass.
 */
    if(bytoff+ndata > FITSLEN)
      ndata = FITSLEN - bytoff;
/*
 * Copy the required data into the FITS I/O buffer.
 */
    memcpy(&fits->buff[bytoff], &fits_buff[nwrit], ndata);
    nwrit += ndata;
    bytoff = 0L;
    fits->modified = 1;
  };
/*
 * How many complete objects were written?
 */
  return (nwrit < nbytes);
}

/*.......................................................................
 * Flush a FITS I/O buffer to a FITS file.
 *
 * Input:
 *  fits    Fits *  The FITS file descriptor.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int fits_flush(Fits *fits)
{
/*
 * Don't write to the file if no modifications have been made to the
 * I/O buffer.
 */
  if(fits->modified) {
/*
 * Was this file opened for writing?
 */
    if(fits->readonly) {
      fprintf(stderr, "fits_flush: File not opened for writing\n");
      return 1;
    };
/*
 * If there is a gap between the end of the file and the record to be
 * written, fill the gap with padded records.
 */
    if(fits->recnum > fits->nullrec && fits_pad(fits, fits->recnum))
      return 1;
/*
 * Move the file pointer to the start of the fits record.
 */
    if(rec_seek(fits->rec, fits->recnum, 0L))
      return 1;
/*
 * Write the record to the FITS file.
 */
    if(rec_write(fits->rec, FITSLEN, sizeof(char), fits->buff)<FITSLEN)
      return 1;
/*
 * Record changed FITS status.
 */
    fits->modified = 0;
/*
 * If the newly written record is beyond the last recorded record in the file,
 * update fits->nullrec to record this.
 */
    if(fits->recnum >= fits->nullrec)
      fits->nullrec = fits->recnum+1;
  };
  return 0;
}

/*.......................................................................
 * Read a new record from a FITS file after flushing the previous one.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file.
 *  recnum  long    The record number to be read.
 *  doreport int    If true report read errors on stderr. Otherwise
 *                  keep silent about read errors.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int fits_read(Fits *fits, long recnum, int doreport)
{
/*
 * First flush the previous record if it has been modified.
 */
  if(fits->modified && fits_flush(fits))
    return 1;
/*
 * If the last HDU is marked as complete, or recnum precedes the
 * current estimate of the end of the file, attempt to read the requested
 * record. Otherwise, if recnum is at or beyond fits->nullrec, don't
 * do an actual read, since there is nothing there to be read. Simply
 * return the buffer padded with the current padding character fits->pad.
 * 
 */
  if(fits->complete || recnum < fits->nullrec) {
/*
 * Position the file at the start of the record.
 */
    if(rec_seek(fits->rec, recnum, 0L))
      return 1;
/*
 * Read the record.
 */
    if(rec_read(fits->rec, FITSLEN, sizeof(char), fits->buff)<FITSLEN) {
      fits->recnum = -1;
      if(doreport)
	fprintf(stderr, "fits_read: Error reading from file: %s\n", fits->name);
      return 1;
    };
/*
 * Update the current estimate of the end-of-file record.
 */
    if(recnum >= fits->nullrec)
      fits->nullrec = recnum+1;
  } else {
    memset(fits->buff, fits->pad, FITSLEN);
  };
/*
 * Record changed FITS status.
 */
  fits->recnum = recnum;
  fits->modified = 0;
  return 0;
}

/*.......................................................................
 * Fill the gap between fits->nullrec-1 and recnum with padded records.
 * Records fits->nullrec up to recnum-1 are padded with the current
 * padding character fits->pad.
 *
 * Input:
 *  fits    Fits *   The fits descriptor.
 *  recnum  long     The last record to be padded.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int fits_pad(Fits *fits, long recnum)
{
/*
 * Have the required records already been written?
 */
  if(recnum <= fits->nullrec)
    return 0;
/*
 * Move to the start of the first record to be padded.
 */
  if(rec_seek(fits->rec, fits->nullrec, 0L))
    return 1;
/*
 * Pad a temporary I/O buffer with the current padding character.
 */
  memset(fits_buff, fits->pad, FITSLEN);
/*
 * Write the contents of fits_buff[] to records in the range
 * fits->nullrec to recnum, inclusive.
 */
  for( ; fits->nullrec < recnum; fits->nullrec++) {
    if(rec_write(fits->rec, FITSLEN, sizeof(char), fits_buff) < FITSLEN)
      return 1;
  };
  return 0;
}
