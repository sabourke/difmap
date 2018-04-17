#ifndef recio_h
#define recio_h

/* Define an enumeration of possible file types */

typedef enum {IS_OLD, IS_NEW, IS_SCR} Fileuse;

/* Define an enumeration describing the nature of the last I/O operation */

typedef enum {REC_RD, REC_WR, REC_SK} Lastio;

typedef struct Recio {
  FILE *fp;            /* File descriptor of opened file */
  char *name;          /* Name of file or NULL if scratch file */
  int readonly;        /* True if file has been opened without write access */
  Fileuse status;      /* The file type */
  Lastio lastio;       /* The nature of the last I/O operation */
  long reclen;         /* Length of record (char's) */
  long reclim;         /* The max atomically seekable record offset */
  long recnum;         /* The record within which the file pointer lies */
  long recoff;         /* The offset of the file pointer into record recnum */
} Recio;

/* Open a new record addressable binary file */

Recio *new_Recio(const char *name, Fileuse status, int readonly, long reclen);

/* Close a file previously opened by new_Recio() */

Recio *del_Recio(Recio *rio);

/* Position a Recio file for read or write */

int rec_seek(Recio *rio, long recnum, long recoff);

/* Read from the current file position */

int rec_read(Recio *rio, size_t nobj, size_t size, void *buff);

/* Write at the current file position */

int rec_write(Recio *rio, size_t nobj, size_t size, void *buff);

/* Rewind the file */

void rec_rewind(Recio *rio);

/* Flush pending I/O */

int rec_flush(Recio *rio);

/* feof() equivalent to test for EOF */

int rec_eof(Recio *rio);

/* ferror() equivalent to test for I/O errors. */

int rec_error(Recio *rio);

/* Get the current file position */

int rec_tell(Recio *rio, long *recnum, long *recoff);

#endif
