#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <errno.h>

#include "recio.h"
#include "scrfil.h"

/* Cater for pre-ANSI-C C libraries such as SUNs bundled libraries */

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif

static int rec_open(Recio *rio, const char *name);
static int rec_bad(Recio *rio, const char *name);

/*.......................................................................
 * Open a binary file, and return a descriptor to be used in future record
 * oriented I/O on this file.
 *
 * Input:
 *  name      char *   The name of the file.
 *  status Fileuse     Enumeration used to specify the file disposition:
 *                       IS_OLD  -  Existing file.
 *                       IS_NEW  -  New file.
 *                       IS_SCR  -  Scratch file - the file will be
 *                                  deleted when del_Recio() is called.
 *  readonly   int     If true, the file will be opened only for reading.
 *                     The file must exist in this case.
 *  reclen    long     The size of one logical record, measured in
 *                     chars. This is used solely as an aid to
 *                     specifying file positions to rec_seek.
 * Output:
 *  return   Recio *   The record-I/O descriptor.
 */
Recio *new_Recio(const char *name, Fileuse status, int readonly, long reclen)
{
  Recio *rio;   /* The new descriptor */
/*
 * Sanity check the arguments.
 */
  if(name==NULL) {
    fprintf(stderr, "new_Recio: NULL file name intercepted.\n");
    return NULL;
  };
  if(reclen <= 0L) {
    fprintf(stderr, "new_Recio: Illegal negative record length specified.\n");
    return NULL;
  };
/*
 * Alocate a container for the descriptor.
 */
  rio = malloc(sizeof(*rio));
  if(rio==NULL) {
    fprintf(stderr, "new_Recio: Insufficient memory to open file.\n");
    return rio;
  };
/*
 * Initialize the members of the container, at least up to the point
 * at which it can safely be passed del_Recio().
 */
  rio->fp = NULL;
  rio->name = NULL;
  rio->readonly = readonly;
  rio->status = status;
  rio->lastio = REC_SK;
  rio->reclen = reclen;
/*
 * Record the max seekable record offset, measured in record lengths.
 */
  rio->reclim = LONG_MAX / rio->reclen;
/*
 * Open the file.
 */
  if(rec_open(rio, name))
    return del_Recio(rio);
/*
 * Set and record the current file position.
 */
  rewind(rio->fp);
  rio->recnum = rio->recoff = 0L;
/*
 * Return the successfully initialized descriptor.
 */
  return rio;
}

/*.......................................................................
 * Delete a Recio descriptor - this includes closing the file if open.
 *
 * Input:
 *  rio    Recio *  The descriptor to be deleted. This must have
 *                  originally been created by new_Recio().
 * Output:
 *  return Recio *  Always NULL. Use like:  rio = del_Recio(rio);
 */
Recio *del_Recio(Recio *rio)
{
  if(rio) {
/*
 * Close the file if it has been opened. Ignore close errors if the
 * file is a scratch file.
 */
    if(rio->fp && fclose(rio->fp)==EOF && rio->status!=IS_SCR)
      fprintf(stderr, "del_Recio: Error closing file: %s\n", rio->name);
/*
 * If the file was a scratch file, and deleting an open file only
 * removes its directory entry, remove it.
 */
#if defined(HIDE_SCRATCH_FILES) && HIDE_SCRATCH_FILES==1
    if(rio->status == IS_SCR && rio->name)
      remove(rio->name);
#endif
/*
 * Delete the copy of the file name.
 */
    if(rio->name)
      free(rio->name);
/*
 * Delete the container.
 */
    free(rio);
  };
  return NULL;
}

/*.......................................................................
 * Return a given portion of a Recio file.
 *
 * Input:
 *  rio    Recio *  The Recio descriptor of the file.
 *  nobj  size_t    The number of objects of size 'size' to be read.
 *  size  size_t    The size of the objects to be read, in chars.
 * Input/Output:
 *  buff    void *  A buffer of sufficient size to contain the
 *                  requested number of objects.
 * Output:
 *  return   int    The number of objects read, which may only differ
 *                  from nobj at the end of the file. On error -1 is
 *                  returned.
 */
int rec_read(Recio *rio, size_t nobj, size_t size, void *buff)
{
  long nnew;    /* The number of complete objects read in the latest read */
  long nreq;    /* The number of objects still to be read */
/*
 * Check args.
 */
  if(rec_bad(rio, "rec_read"))
    return -1;
  if(buff==NULL) {
    fprintf(stderr, "rec_read: NULL buffer received.\n");
    return -1;
  };
/*
 * Continue no further if an I/O error previously occured.
 */
  if(ferror(rio->fp))
    return -1;
/*
 * If the last I/O operation on the file was a write, then
 * use a file positioning command before continuing (The ANSI C standard
 * requires this when switching I/O direction).
 */
  if(rio->lastio==REC_WR && fseek(rio->fp, 0L, SEEK_CUR)) {
    fprintf(stderr, "rec_read: Error positioning file: %s\n", rio->name);
    return -1;
  };
/*
 * Ready to read.
 */
  rio->lastio = REC_RD;
/*
 * Read directly into the supplied buffer.
 */
  for(nreq=nobj; nreq>0; buff = (char *)buff + nnew) {
#ifdef EINTR
    errno = 0;
#endif
    nnew = fread(buff, size, nreq, rio->fp);
/*
 * Accumulate a count of the number of objects still to be read.
 */
    nreq -= nnew;
/*
 * I/O error?
 */
    if(ferror(rio->fp)) {
#ifdef EINTR
      if(errno==EINTR) {   /* Interrupted system call - prepare to retry */
	clearerr(rio->fp);
      } else
#endif
      {
	fprintf(stderr, "rec_read: Error reading from file: %s\n", rio->name);
	return -1;
      };
    }
/*
 * End of file?
 */
    else if(feof(rio->fp)) {
      break;
    };
  };
/*
 * Update the record of the file position.
 */
  {
    long recoff = (nobj-nreq) * size;
    long numoff = recoff / rio->reclen;
    rio->recnum += numoff;
    rio->recoff += recoff - numoff * rio->reclen;
  };
/*
 * Return a count of the number of complete objects read.
 */
  return nobj - nreq;
}

/*.......................................................................
 * Return a given portion of a Recio file.
 *
 * Input:
 *  rio    Recio *  The Recio descriptor of the file.
 *  nobj  size_t    The number of objects of size 'size' to be written.
 *  size  size_t    The size of the objects to be written, in chars.
 * Input/Output:
 *  buff    void *  The buffer containing the data to be written.
 * Output:
 *  return   int    The number of objects written, which is less than
 *                  the number specified only on error.
 */
int rec_write(Recio *rio, size_t nobj, size_t size, void *buff)
{
  long nnew;    /* The number of complete objects written in the latest write */
  long nreq;    /* The number of objects still to be written */
/*
 * Check args.
 */
  if(rec_bad(rio, "rec_write"))
    return 0;
  if(buff==NULL) {
    fprintf(stderr, "rec_write: NULL buffer received.\n");
    return 0;
  };
/*
 * If the last I/O operation on the file was a read, then
 * use a file positioning command before continuing (The ANSI C standard
 * requires this when switching I/O direction).
 */
  if(rio->lastio==REC_RD && fseek(rio->fp, 0L, SEEK_CUR)) {
    fprintf(stderr, "rec_write: Error positioning file: %s\n", rio->name);
    return 0;
  };
/*
 * Ready to write.
 */
  rio->lastio = REC_WR;
/*
 * Write directly from the supplied buffer.
 */
  for(nreq=nobj; nreq>0; buff = (char *)buff + nnew) {
#ifdef EINTR
    errno = 0;
#endif
    nnew = fwrite(buff, size, nreq, rio->fp);
/*
 * Update the count of the number of objects still to be written.
 */
    nreq -= nnew;
/*
 * I/O error?
 */
    if(ferror(rio->fp)) {
#ifdef EINTR
      if(errno==EINTR) {   /* Interrupted system call - prepare to retry */
	clearerr(rio->fp);
      } else
#endif
      {
	fprintf(stderr, "rec_write: Error writing to file: %s\n", rio->name);
	return -1;
      };
    };
  };
/*
 * Update the record of the file position.
 */
  {
    long recoff = (nobj - nreq) * size;
    long numoff = recoff / rio->reclen;
    rio->recnum += numoff;
    rio->recoff += recoff - numoff * rio->reclen;
  };
/*
 * Return a count of the number of complete objects written.
 */
  return nobj - nreq;
}

/*.......................................................................
 * Seek to a position specified by a logical record number and offset.
 *
 * Input:
 *  recnum   long   The 0-relative index of the record to position
 *                  within.
 *  recoff   long   The 0-relative offset (chars) within record recnum
 *                  to position at.
 * Output:
 *  return    int   0 - OK.
 *                  1 - Error (On error the file is rewound).
 */
int rec_seek(Recio *rio, long recnum, long recoff)
{
  long recdif;   /* The number of whole records that we are away from */
                 /*  the desired position. */
  long offdif;   /* The partial record not accounted for by recdif, */
                 /*  recorded in bytes. */
/*
 * Check args.
 */
  if(rec_bad(rio, "rec_seek"))
    return 1;
  if(recnum < 0 || recoff < 0) {
    fprintf(stderr, "rec_seek: Negative record %s.\n",
	    recnum<0 ? "number" : "offset");
    return 1;
  };
/*
 * How far off are we in whole records and residual bytes?
 */
  recdif = recnum - rio->recnum;
  offdif = recoff - rio->recoff;
  while(offdif >= rio->reclen) {
    recdif++;
    offdif -= rio->reclen;
  };
  while(offdif <= -rio->reclen) {
    recdif--;
    offdif += rio->reclen;
  };
/*
 * Are we already correctly positioned?
 */
  if(recdif==0 && offdif==0)
    return 0;
/*
 * First move by offdif bytes.
 */
  if(offdif != 0 && fseek(rio->fp, offdif, SEEK_CUR)) {
    fprintf(stderr, "rec_seek: Error positioning file: %s (%s)\n",
	    rio->name, strerror(errno));
    rec_rewind(rio);
    return 1;
  };
/*
 * Now move the required number of records, moving in steps of at
 * most rio->reclim records to avoid overflowing the (long int) offset.
 */
  while(recdif != 0) {
    long newoff = recdif > rio->reclim ? rio->reclim :
      (recdif < -rio->reclim ? -rio->reclim : recdif);
    if(fseek(rio->fp, newoff * rio->reclen, SEEK_CUR)) {
      fprintf(stderr, "rec_seek: Error positioning file: %s (%s)\n",
	      rio->name, strerror(errno));
      rec_rewind(rio);
      return 1;
    };
/*
 * Record the IO operation type.
 */
    rio->lastio = REC_SK;
    recdif -= newoff;
  };
/*
 * The file is now positioned.
 */
  rio->recnum = recnum;
  rio->recoff = recoff;
  return 0;
}

/*.......................................................................
 * Make a malloc'd copy of a given file name and record it in
 * rio->name, then open the file. If the file is to be a scratch file,
 * try not to overwrite any other files, by repeatedly postfixing a
 * number to the file name until failure to open the corresponding
 * file indicates that no other readable scratch file of that name exists.
 * Record the file pointer in rio->fp.
 *
 * Input:
 *  rio      Recio    The record-I/O descriptor.
 *  name      char *  The name of the file to be opened.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int rec_open(Recio *rio, const char *name)
{
  char *mode;        /* The mode string to be sent to fopen() */
/*
 * Valid arguments?
 */
  if(rio==NULL) {
    fprintf(stderr, "rec_open: Intercepted NULL Recio descriptor.\n");
    return 1;
  };
  if(name==NULL) {
    fprintf(stderr, "rec_open: Intercepted NULL file name argument.\n");
    return 1;
  };
/*
 * Sanity check the requested file disposition.
 */
  if(rio->readonly && rio->status != IS_OLD) {
    fprintf(stderr, "rec_open: Readonly requested for non-existent file.\n");
    return 1;
  };
/*
 * Select the appropriate fopen() mode string.
 */
  if(rio->readonly)
    mode = "rb";
  else if(rio->status==IS_OLD)
    mode = "r+b";
  else
    mode = "w+b";
/*
 * If this is to be a scratch file then get a unique name for the file,
 * formed by the given basis name in 'name' plus a
 * numeric postfix if other files of the same name exist.
 * Otherwise simply make a copy of 'name'.
 */
  if(rio->status==IS_SCR) {
    rio->name = scrname(name);
  } else {
    rio->name = (char *) malloc(strlen(name) + 1);
    if(rio->name)
      strcpy(rio->name, name);
  };
  if(rio->name==NULL) {
    fprintf(stderr, "rec_open: Insufficient memory to record file name.\n");
    return 1;
  };
/*
 * Attempt to open the file.
 */
  rio->fp = fopen(rio->name, mode);
  if(rio->fp==NULL) {
    fprintf(stderr, "rec_open: Unable to open file: %s\n", rio->name);
    return 1;
  };
/*
 * If this is a scratch file, remove its directory entry, so that it
 * gets deleted automatically when the program exits for whatever reason.
 */
  if(rio->status==IS_SCR)
    remove(rio->name);
/*
 * Turn off stdio buffering to avoid the overhead of the standard I/O
 * library reading ahead large amounts on each of the many fseek()
 * calls that this facility performs.
 */
  if(setvbuf(rio->fp, NULL, _IONBF, 0) != 0) {
    fprintf(stderr, "rec_open: setvbuf error (%s)\n", strerror(errno));
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Given a buffer of length nbuff, rewind the file, fill the buffer
 * with nbuff '\0's, write the buffer npad times to the file and
 * finally rewind the file.
 *
 * Input:
 *  rio    Recio *  The record I/O descriptor.
 *  buff    void *  A buffer of length 'nbuff' used to pad the file
 *                  npad times.
 *  nbuff   long    The length of buff[] in chars.
 *  npad    long    The number of times to write buff[] to the file.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int rec_pad(Recio *rio, void *buff, long nbuff, long npad)
{
  long ipad;   /* The number of buffers written */
/*
 * Check arguments.
 */
  if(rec_bad(rio, "rec_pad"))
    return 1;
  if(buff==NULL) {
    fprintf(stderr, "rec_pad: Intercepted NULL buffer argument.\n");
    return 1;
  };
/*
 * Rewind to the start of the file.
 */
  rec_rewind(rio);
/*
 * Record the IO operation type.
 */
  rio->lastio = REC_SK;
/*
 * Clear the buffer.
 */
  if(nbuff>0)
    memset(buff, 0, nbuff);
/*
 * Nothing to be written?
 */
  if(npad<=0 || nbuff<=0)
    return 0;
/*
 * File not opened for write?
 */
  if(rio->readonly) {
    fprintf(stderr, "rec_pad: File %s not opened for writing.\n", rio->name);
    return 1;
  };
/*
 * Record the pending IO operation type.
 */
  rio->lastio = REC_WR;
/*
 * Write it to the file npad times.
 */
  for(ipad=0; ipad<npad; ipad++) {
    if(fwrite(buff, 1, nbuff, rio->fp) < nbuff) {
      fprintf(stderr, "rec_pad: Error writing to file: %s\n", rio->name);
      return 1;
    };
  };
/*
 * Ready the file for subsequent I/O.
 */
  rec_rewind(rio);
/*
 * Record the IO operation type.
 */
  rio->lastio = REC_SK;
  return 0;
}

/*.......................................................................
 * Rewind the record I/O file. This also clears any file error status.
 *
 * Input:
 *  rio    Recio *  The descriptor of the record-I/O file.
 */
void rec_rewind(Recio *rio)
{
/*
 * Note that we mustn't call rec_bad() here because it would return 1
 * if the file error flag was set, and one of the things that rec_rewind
 * promises to do (through the advertised behavior of rewind()), is
 * reset any error conditions.
 */
  if(rio && rio->fp) {
    rewind(rio->fp);
    rio->recnum = rio->recoff = 0L;
  };
  return;
}

/*.......................................................................
 * Emit an error message and return 1 if the Recio descriptor given
 * is invalid. If evidence of a previous I/O error is detected, no
 * error message will be emitted.
 *
 * Input:
 *  rio   Recio *  The descriptor of a record I/O file.
 *  name   char *  The name of the calling function.
 * Output:
 *  return  int    0 - OK.
 *                 1 - rio is un-useable.
 */
static int rec_bad(Recio *rio, const char *name)
{
  if(rio==NULL) {
    fprintf(stderr, "%s: NULL Recio descriptor intercepted.\n", name);
    return 1;
  };
  return ferror(rio->fp);
}

/*.......................................................................
 * Ensure that the file is up to date by flushing pending I/O.
 *
 * Input:
 *  rio   Recio *  The descriptor of a record I/O file.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int rec_flush(Recio *rio)
{
/*
 * Check args.
 */
  if(rec_bad(rio, "rec_flush"))
    return 0;
/*
 * Have the file flushed.
 */
  return fflush(rio->fp)==EOF;
}

/*.......................................................................
 * Return true if EOF was reached on the last read.
 *
 * Input:
 *  rio   Recio *  The descriptor of a record I/O file.
 * Output:
 *  return  int    0 - OK.
 *                 1 - EOF flag is set.
 */
int rec_eof(Recio *rio)
{
/*
 * Check args.
 */
  if(rec_bad(rio, "rec_eof"))
    return 0;
/*
 * Have the file flushed.
 */
  return feof(rio->fp);
}

/*.......................................................................
 * Return true if the error flag is set on a recio file.
 *
 * Input:
 *  rio   Recio *  The descriptor of a record I/O file.
 * Output:
 *  return  int    0 - OK.
 *                 1 - EOF flag is set.
 */
int rec_error(Recio *rio)
{
  return rec_bad(rio, "rec_error");
}

/*.......................................................................
 * Return the current file position.
 *
 * Input:
 *  rio   Recio *  The descriptor of a record I/O file.
 * Input/Output:
 *  recnum long *  If not NULL, *recnum will be assigned with the 0-relative
 *                 index of the record in which the file pointer points.
 *  recoff long *  If not NULL, *recoff will be assigned with byte-offset
 *                 of the file pointer in record 'recnum'.
 * Output:
 *  return         0 - OK.
 *                 1 - Error.
 */
int rec_tell(Recio *rio, long *recnum, long *recoff)
{
  if(rec_bad(rio, "rec_tell"))
    return 1;
  if(recnum)
    *recnum = rio->recnum;
  if(recoff)
    *recoff = rio->recoff;
  return 0;
}
