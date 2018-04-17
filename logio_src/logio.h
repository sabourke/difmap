#ifndef logio_h
#define logio_h

#include <stdarg.h>

/* Sphere message logging interface  */

FILE *logfile(const char *name); /* Start a new LOG file is name != NULL */
int logtofile(FILE *stream, const char *message, int nchar);

/* Functional equivalents to stdio output functions */

int lprintf(FILE *fp, const char *format, ...);             /* fprintf() */
int vlprintf(FILE *stream, const char *format, va_list ap); /* vfprintf() */
int lputs(const char *s, FILE *stream);                     /* fputs() */
int lputc(int c, FILE *fp);                                 /* fputc() */

/* Get and/or change the set of streams to be logged */

typedef enum {LOG_IN=1, LOG_OUT=2, LOG_ERR=4, LOG_ALL=8} Logset;
typedef enum {LOG_SET, LOG_CLR, LOG_REP} Logoper;
int log_streams(int logset, Logoper oper);

/*
 * Define the type of output function handed to lprint().
 * The 'out' argument is a pointer to a function specific repository
 * of extra info (such as a FILE * or buffer pointer) required by the
 * particular function sent. The function is required to output 'n'
 * characters from buff[] and return the number output, or -1 on error.
 */
#define LOGFN(fn) int (fn)(void *out, char *buff, int n)
typedef LOGFN(*Logfn);

/*
 * General printing function.
 */
int lprint(Logfn output, void *out, const char *format, va_list ap);

#endif
