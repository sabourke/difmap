#include <stdio.h>
#include <time.h>
#include <stdarg.h>
#include <string.h>
#include <stdlib.h>

#include "logio.h"
#include "scrfil.h"

/* Encapsulate local log file parameters to avoid consflicts with extern's */

static struct {
  FILE *stream;            /* The log file FILE pointer */
  char *name;              /* Malloc'd copy of file name */
  int waseol;              /* True if last write reached end of line */
  FILE *last_stream;       /* The previous stream written to. */
  int dolog;               /* True if last_stream was logged */
  int logset;              /* Set of streams to be logged */
} lg = {
  NULL, NULL, 1, NULL, 0, 0
};

static void logerror(void);
static LOGFN(logwrt);
static FILE *closelog(int waserr);

/*.......................................................................
 * Close any existing log file and if the file name given is non-NULL
 * open a new log file of that name.
 *
 * Input:
 *  name     char *   The name of the new log file. If NULL don't open
 *                    a new log file.
 * Output:
 *  return   FILE *   The file pointer of the log file, or NULL if the
 *                    open operation failed.
 */
FILE *logfile(const char *name)
{
  time_t t=time(NULL);   /* The current time */
/*
 * The next thing to be written to the log file will be at the start
 * of a new line.
 */
  lg.waseol = 1;
  lg.last_stream = NULL;
/*
 * Close any existing log file.
 */
  closelog(0);
/*
 * No new log file to be opened?
 */
  if(name == NULL)
    return closelog(0);
/*
 * Find an unambiguous file name for the log file.
 */
  lg.name = scrname(name);
  if(lg.name==NULL)
    return closelog(1);
/*
 * Open a new log file of the name given and specify line-buffering
 * so that the file will be kept up to date.
 */
  lg.stream = fopen(lg.name, "w");
  if(lg.stream==NULL || setvbuf(lg.stream, NULL, _IOLBF, BUFSIZ)) {
    fprintf(stderr, "logfile: Error opening file: %s\n", lg.name);
    return closelog(1);
  };
/*
 * Assign the default set of streams to be logged.
 */
  log_streams((LOG_IN|LOG_OUT|LOG_ERR), LOG_REP);
/*
 * Write logfile startup message.
 */
  lprintf(stdout, "Started logfile: %s on %s", lg.name,
	  (t == -1) ? "(date unavailable)\n" : ctime(&t));
/*
 * Return the pointer to the new log file stream.
 */
  return lg.stream;
}

/*.......................................................................
 * Close a log file and free the malloc'd copy of the file name.
 *
 * Input:
 *  waserr  int  If true, the stream is being closed due to an I/O error,
 *               so don't log anything to it before closing it.
 */
static FILE *closelog(int waserr)
{
/*
 * Write a closedown message?
 */
  if(lg.stream && lg.name) {
    time_t t=time(NULL);   /* The current time */
    if(waserr) {
      fprintf(stdout, "Log file %s closed on %s", lg.name,
	      (t == -1) ? "(date unavailable)\n" : ctime(&t));
    } else {
      lprintf(stdout, "Log file %s closed on %s", lg.name,
	      (t == -1) ? "(date unavailable)\n" : ctime(&t));
    };
  };
/*
 * Release the memory taken by the copy of the file name.
 */
  if(lg.name) {
    free(lg.name);
    lg.name = NULL;
  };
/*
 * Close the log file.
 */
  if(lg.stream && fclose(lg.stream) == EOF) {
    fprintf(stderr, "Warning: Error while closing log file\n");
  };
  return (lg.stream=NULL);
}

/*.......................................................................
 * Provide an equivalent function to fprintf() that writes both to the
 * designated stream and to the log file opened via logfile(). All
 * messages are written to the log prefixed with the sphere comment
 * character unless the cited stream is stdin. The log file is then
 * useable as a command file.
 *
 * Input:
 *  stream  FILE *   The normal stream for the message.
 *                    By default only stdin,stdout and stderr are logged.
 *                    stdin is logged as a command. All other logged streams
 *                    are logged as comment lines. See log_streams().
 *  format  char *   A standard printf() format string.
 *  ...              Variable argument list as in printf().
 * Output:
 *  return   int     The number of characters written on stream, or negative
 *                   on error.
 */
int lprintf(FILE *stream, const char *format, ...)
{
  static va_list ap;   /* Variable argument list */
  int nchar=0;         /* Number of characters written to output stream */
/*
 * Let vlprint() do the hard work.
 */
  va_start(ap, format);
  nchar = vlprintf(stream, format, ap);
  va_end(ap);
  return nchar;
}

/*.......................................................................
 * The equivalent of fputs() for log-file IO.
 *
 * Input:
 *  s   const char *  The string to be written.
 *  stream    FILE *  The stream to be written to.
 * Output:
 *  return     int    0   - OK.
 *                    EOF - Error.
 */
int lputs(const char *s, FILE *stream)
{
  return logtofile(stream, s, strlen(s)) != -1 ? 0 : EOF;
}

/*.......................................................................
 * The equivalent of fputc() for log-file IO.
 *
 * Input:
 *  c           int    The character to be written.
 *  stream     FILE *  The stream to be written to [or stdin - see lprintf()].
 * Output:
 *  return  int    The character written, or EOF on error on stream stream.
 */
int lputc(int c, FILE *stream)
{
  static char buff[2]={'\0', '\0'};
/*
 * Write the character to the output stream and to the log file.
 */
  buff[0] = c;
  return logtofile(stream, buff, 1)==1 ? c : EOF;
}

/*.......................................................................
 * Change the set of streams that get logged.
 *
 * Input:
 *  logset      int   A set of streams - This should be a union of one or
 *                    more of the following bit enumerators:
 *                     LOG_IN  - Log stdin.
 *                     LOG_OUT - Log stdout.
 *                     LOG_ERR - Log stderr.
 *                     LOG_ALL - Log all streams.
 *                    Thus (LOG_IN | LOG_ERR) means - stdin and stderr.
 *  oper        int   The operation to which the bit set refers.
 *                     LOG_SET   - Add the streams to the set.
 *                     LOG_CLR   - Remove the streams from the set.
 *                     LOG_REP   - Replace the set with that given.
 * Output:
 *  return      int   A copy of the modified bit set.
 */
int log_streams(int logset, Logoper oper)
{
  switch (oper) {
  case LOG_SET:          /* Add selected streams to the set */
    lg.logset |= logset;
    break;
  case LOG_CLR:          /* Remove selected streams from the set */
    lg.logset &= ~logset;
    break;
  case LOG_REP:          /* Replace the set with that given */
    lg.logset = logset;
    break;
  };
/*
 * Enforce a re-examination of whether the next stream to be logged
 * should be logged.
 */
  lg.last_stream = NULL;
  lg.dolog = 0;
  return lg.logset;
}

/*.......................................................................
 * After a log-file write error, write an error message and close the
 * log file.
 */
static void logerror(void)
{
  fprintf(stderr, "Error writing to log file - closing log file\n");
  clearerr(lg.stream);
  closelog(1);
  return;
}

/*.......................................................................
 * Provide the equivalent of vfprintf() for log I/O.
 *
 * Input:
 *  stream  FILE *   The normal stream for the message.
 *                    By default only stdin,stdout and stderr are logged.
 *                    stdin is logged as a command. All other logged streams
 *                    are logged as comment lines. See log_streams().
 *  format  char *   A standard printf() format string.
 *  arg  va_list     Variable argument list as in vfprintf().
 * Output:
 *  return   int     The number of characters written on stream, or negative
 *                   on error.
 */
int vlprintf(FILE *stream, const char *format, va_list ap)
{
  return lprint(logwrt, stream, format, ap);
}

/*.......................................................................
 * Output function for LOG I/O functions.
 *
 * Input:
 *  out     void *   The pointer to the stream to be written to.
 *                   [(FILE *) cast to (void *)].
 *  buff    char *   The buffer containing the text to be written.
 *  n        int     The number of characters to write from buff[].
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
static LOGFN(logwrt)
{
  return logtofile((FILE *) out, buff, n) != n;
}

/*.......................................................................
 * Write a NULL terminated line to the given stream and optionally
 * to a log file. All log file output must go through this function.
 *
 * Input:
 *  stream  FILE *   The normal stream for the message.
 *                    By default only stdin,stdout and stderr are logged.
 *                    stdin is logged as a command. All other logged streams
 *                    are logged as comment lines. See log_streams().
 *  message char *   A '\0' terminated message to be logged.
 *  nchar    int     The number of characters to write from message[].
 * Output:
 *  return   int     The number of characters written - or -1 on error.
 */
int logtofile(FILE *stream, const char *message, int nchar)
{
/*
 * Invalid stream?
 */
  if(stream==NULL) {
    fprintf(stderr, "logtofile: NULL stream intercepted\n");
    return -1;
  };
/*
 * Nothing to write?
 */
  if(nchar<1)
    return 0;
/*
 * Write the message to the given stream?
 */
  if(stream!=stdin && stream!=lg.stream) {
    if(fwrite(message, sizeof(char), (size_t) nchar, stream)!=nchar)
      return -1;
  };
/*
 * Is there a log file to be logged to?
 */
  if(lg.stream==NULL)
    return nchar;
/*
 * If the stream being logged has changed then determine whether
 * the new stream should be logged and ensure that the new line starts on
 * a new line.
 */
  if(stream!=lg.last_stream) {
    int logset = lg.logset;
    if((stream==stdout && logset & LOG_OUT) ||
       (stream==stdin  && logset & LOG_IN)  ||
       (stream==stderr && logset & LOG_ERR) ||
       stream==lg.stream || (logset & LOG_ALL)) {
/*
 * Record the fact that this stream was logged so that 
 * if the next call to this function uses the same stream, we don't
 * have to redo the above checks.
 */
      lg.dolog = 1;
/*
 * Ensure that the change of stream is accompanied by the start of a
 * new line in the log file.
 */
      if(!lg.waseol) {
	lg.waseol = 1;
	fputc('\n', lg.stream);
      };
    } else {
      lg.dolog = 0;   /* This stream is not to be logged */
    };
  };
/*
 * Keep a record of the latest stream.
 */
  lg.last_stream = stream;
/*
 * Logging required?
 */
  if(lg.dolog) {
    int waserr=0;      /* Error status while writing to the log file */
    const char *start; /* Pointer to start of line in message[] */
    const char *endp;  /* Pointer to end of line in message[] */
    int nc;            /* Number of characters in latest message line */
    int n=0;           /* Total number of characters used */
/*
 * Parse each line out of the buffer so that a comment character can be
 * prefixed to the start of each new comment line.
 */
    for(start=endp=message;  !waserr && n<nchar && *start!='\0';  start=endp) {
/*
 * If starting a new comment line, prefix it with the sphere comment
 * character and a space.
 */
      if(lg.waseol && stream!=stdin)
	waserr = fwrite("! ", sizeof(char), 2, lg.stream) != 2;
/*
 * Find the end of the next message line.
 */
      for(; n<nchar && *endp!='\0' && *endp!='\n' && *endp!='\r'; endp++, n++);
/*
 * End of line?
 */
      lg.waseol = *endp == '\n' || *endp == '\r';
      if(lg.waseol && n<nchar)
	endp++,n++;      /* Include the newline or carriage return */
/*
 * Write the new text to the log file.
 */
      nc = (endp-start);
      waserr = nc>0 && fwrite(start, sizeof(char), nc, lg.stream) != nc;
    };
/*
 * If an error occurred while writing the log file, report the error,
 * clear it, and close the file.
 */
    if(waserr)
      logerror();
  };
  return nchar;
}
