#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>

#include "logio.h"

/*
 * As a format specifier is parsed, the attributes that are found
 * are recorded in an object of the following type.
 */
typedef struct {
  int left;        /* Left adjust field */
  int sign;        /* Prefix number with sign */
  int space;       /* Prefix number with space if no sign is to be printed */
  int zero;        /* Zero-pad to the field width with leading zeros */
  int alt;         /* Alternatate form of output */     
  int min;         /* Minimum field width */
  int prec;        /* Precision */
  int qual;        /* Qualifier character */
  int conv;        /* Conversion specification character */
} Options;

/*
 * Specify the size of the buffer used to contain the formatted output
 * associated with a single format directive. ANSI says that a conforming
 * program should not exceed 509 chars per format item. Checks will be made
 * to ensure that the buffer size is not exceeded.
 */

enum {LPBUFF_MAX=512};

static char *get_flags(Options *opts, char *form);
static int space_pad(Logfn output, void *out, int npad);
static int lperror(const char *format, ...);

/*.......................................................................
 * Function responsible for parsing the format string and sending the
 * lprintf output to a user provided output function.
 *
 * Input:
 *  output  Logfn    Pointer to function that takes the output buffer.
 *  out      void *  Pointer to output stream, output buffer or to other
 *                   information required by the particular output function
 *                   'output()'. This is the first argument to output().
 *  format   char *  The printf style format string.
 *  ap    va_list    The Argument list associated with 'format'.
 */
int lprint(Logfn output, void *out, const char *format, va_list ap)
{
  static Options def_opts = {0,0,0,0,0,0,-1,0,0};  /* Default options */
  static char buff[LPBUFF_MAX+10];                 /* Output buffer */
  static Options opts; /* Printing options */
  char *sptr;          /* Pointer to start of format item */
  char *eptr;          /* Pointer to one-past end of format item */
  int c;        /* The latest unprocessed character in the format string */
  int ntotal=0; /* The total number of characters so far written */
  int n;        /* The latest number of characters to write */
  int toobig;   /* True if the formatted output would overflow the buffer */
/*
 * Sanity check.
 */
  if(format==NULL)
    return lperror("lprint: NULL printf format string intercepted\n");
/*
 * Parse the format string.
 */
  for(sptr=eptr=(char *)format; *sptr; sptr=eptr) {
/*
 * Find the next % character.
 */
    do c = *eptr++; while(c && c!='%');
/*
 * Output everything up to the %.
 */
    n = (eptr-sptr) - 1;
    if(n>0 && output(out, sptr, n))
      return -1;
/*
 * Have we fallen off the end of the format string?
 */
    if(c=='\0')
      eptr--;
/*
 * Increment the count of the total number of characters so far written.
 */
    ntotal += n;
/*
 * Handle the new format directive.
 */
    if(c=='%') {
      char *bufptr = &buff[0]; /* The pointer to the converted output */
      int dopad = 0;           /* True when space padding is required */
      int prefix_len;          /* The length of the % and flags part of */
                               /*  the format directive. */
      char *suffix;            /* The part of the directive which follows */
                               /*  the precision */
/*
 * Keep a record of the start of the format directive (including the %).
 */
      sptr = eptr-1;
/*
 * Start with the default conversion flags.
 */
      opts = def_opts;
/*
 * Get the conversion flags.
 */
      eptr = get_flags(&opts, eptr);
      prefix_len = eptr - sptr;
/*
 * Get the minimum field width. If the next character is an asterix
 * then get the width from the next argument. Otherwise, if the next
 * character is a digit then read the int from the format string.
 */
      if(*eptr == '*') {
	eptr++;
	opts.min = va_arg(ap, int);
      } else if(isdigit((int) *eptr)) {
	opts.min = (int) strtol(eptr, &eptr, 10);
      } else {
	opts.min = 0;
      };
/*
 * If the next character is a decimal point, then get the precision that follows it. This can either
 * be specified as '*' to specify that the next argument contains the precision, or as a literal number
 * in the format string. If omitted, a value of -1 selects the default precision.
 */
      if(*eptr=='.') {
	eptr++;
	if(*eptr == '*') {
	  eptr++;
	  opts.prec = va_arg(ap, int);
	} else if(isdigit((int) *eptr)) {
	  opts.prec = (int) strtol(eptr, &eptr, 10);
	} else {
	  opts.prec = -1;
	};
      };
      suffix = eptr;
/*
 * Get any length modifier from the format string.
 */
      switch(*eptr) {
      case 'h': case 'l': case 'L':
	opts.qual = *eptr++;
      };
/*
 * Get the conversion specifier character.
 */
      opts.conv = c = *eptr++;
/*
 * Limit the size of the formatted output for all but string values.
 * String values can be passed directly without having to copy them
 * into the constraining print buffer.
 */
      if(opts.min > LPBUFF_MAX && opts.conv!='s') {
	toobig = 1;
      } else {
	enum {MAX_FORM=80};     /* Set the max width of any format directive */
	char subform[MAX_FORM+1];
/*
 * Estimate the lengths of the opts.min and opts.prec integers
 * rendered as strings.
 */
	int min_len = opts.min > 0 ? ceil(log10((double)opts.min)) : 0;
	int prec_len = opts.prec >= 0 ?
	  (opts.prec==0 ? 1 : ceil(log10((double)opts.prec))) : 0;
	toobig = 0;
/*
 * Write a version of the format directive with any *'s replaced by
 * the integers read from the argument list.
 */
	if(prefix_len + (eptr-suffix) + min_len + (prec_len ? prec_len+1:0)
	   >= MAX_FORM) {
	  return lperror("Format directive %s too large\n", sptr);
	};
/*
 * Start with the flags from the original format string.
 */
	strncpy(subform, sptr, prefix_len);
	subform[prefix_len] = '\0';
/*
 * If a minimum field width was specified, write it to the format.
 */
	if(opts.min > 0)
	  sprintf(subform + prefix_len, "%d", opts.min);
/*
 * If a precision was specified, write it to the format.
 */
	if(opts.prec >= 0)
	  sprintf(subform + strlen(subform), ".%d", opts.prec);
/*
 * Finally copy the original qualifiers and type specifier to the
 * format.
 */
	strncat(subform + strlen(subform), suffix, eptr - suffix);
/*
 * The precision implies a mimimum field width for some format types.
 * Check for potential buffer overflow and write the value into the buffer.
 */
	switch(c) {
	case 'd': case 'i':  /* Plain integer */
	case 'o':            /* Octal */
	case 'x': case 'X':  /* Hex */
	case 'u':            /* Unsigned */
	  toobig = opts.prec > LPBUFF_MAX-2;
	  if(!toobig) {
	    if(opts.qual=='l')
	      sprintf(buff, subform, va_arg(ap, long));
	    else
	      sprintf(buff, subform, va_arg(ap, int));
	  };
	  break;
	case 'c':            /* Output a single character: Padded below */
	  buff[0] = (unsigned char) va_arg(ap, int);
	  buff[1] = '\0';
	  dopad = 1;
	  break;
	case 's':           /* String: Handled below */
	  bufptr = va_arg(ap, char *);
	  dopad = 1;
	  break;
	case 'f':           /* Plain float */
	  toobig = opts.prec > LPBUFF_MAX-2;
	  if(!toobig) {
	    if(opts.qual=='L')
	      sprintf(buff, subform, va_arg(ap, long double));
	    else
	      sprintf(buff, subform, va_arg(ap, double));
	  };
	  break;
	case 'e': case 'E': /* Floating point with exponent */
	case 'g': case 'G': /* Floating point with or without exponent */
	  toobig = opts.prec > LPBUFF_MAX-7;
	  if(!toobig) {
	    if(opts.qual=='L')
	      sprintf(buff, subform, va_arg(ap, long double));
	    else
	      sprintf(buff, subform, va_arg(ap, double));
	  };
	  break;
	case 'p':           /* Pointer (machine specific representation) */
	  sprintf(buff, subform, va_arg(ap, void *));
	  break;
	case 'n':           /* Return count of chars output so far */
/*
 * Assign to the appropriately typed return argument.
 */
	  switch(opts.qual) {
	  case 'h':
	    *va_arg(ap, short *) = ntotal;
	    break;
	  case 'l':
	    *va_arg(ap, long *) = ntotal;
	    break;
	  default:
	    *va_arg(ap, int *) = ntotal;
	    break;
	  };
	  buff[0] = '\0';
	  break;
	case '%':           /* A plain % sign */
	  buff[0] = '%';
	  buff[1] = '\0';
	  break;
	default:
	  return lperror("lprint: Bad conversion character (%c) in format\n",
			 c);
	};
      };
/*
 * Was the formatted output too big?
 */
      if(toobig) {
	return lperror("lprint: Format \"%.10s...\" too wide for buffer\n",
		       sptr);
      };
/*
 * Determine the number of characters to be output from the buffer.
 */
      switch(opts.conv) {
      case 'n':
	n = 0;   /* Nothing to output */
	break;
      case 's':  /* Argument string */
/*
 * Work out the amount of the string to be written.
 */
	if(opts.prec < 0)
	  n = strlen(bufptr);
	else {
	  char *cptr=bufptr;
	  for(n=0; n<opts.prec && *cptr; n++,cptr++);
	};
	break;
      default:
	n = strlen(bufptr);
	break;
      };
/*
 * If a field width greater than the length of the string has been specified
 * and right-adjustment padding is required - pad to the required width
 * with spaces.
 */
      if(dopad && n<opts.min && !opts.left) {
	int npad=opts.min - n;
	if(space_pad(output, out, npad))
	  return -1;
	ntotal += npad;
      };
/*
 * Output the buffer.
 */
      if(output(out, bufptr, n))
	return -1;
      ntotal += n;
/*
 * If a field width greater than the length of the string has been specified
 * and left-adjustment padding is required - pad to the required width with
 * spaces.
 */
      if(dopad && n<opts.min && opts.left) {
	int npad=opts.min - n;
	if(space_pad(output, out, npad))
	  return -1;
	ntotal += npad;
      };
    };
  };
/*
 * Return a count of the number of characters output.
 */
  return ntotal;
}

/*.......................................................................
 * Register the format options that immediately follow a % format
 * directive.
 *
 * Input:
 *  opts     Options *   The descriptor to register the options in.
 *  form        char *   The pointer into the format string.
 * Output:
 *  return      char *   The pointer to the next un-processed character
 *                       in the format string.
 */
static char *get_flags(Options *opts, char *form)
{
/*
 * Identify all flags.
 */
  for(;;form++) {
    switch(*form) {
    case '-':
      opts->left = 1;
      break;
    case '+':
      opts->sign = 1;
      break;
    case ' ':
      opts->space = 1;
      break;
    case '0':
      opts->zero = 1;
      break;
    case '#':
      opts->alt = 1;
      break;
    default:
      return form;  /* Return the pointer to the unprocessed char */
    };
  };
}

/*.......................................................................
 * Output a given number of space padding characters.
 *
 * Input:
 *  output  Logfn    Pointer to function that takes the output buffer.
 *  out      void *  Pointer to output stream, output buffer or to other
 *                   information required by the particular output function
 *                   'output()'. This is the first argument to output().
 *  npad      int    The number of spaces to be written.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int space_pad(Logfn output, void *out, int npad)
{
  static char spaces[]="                                            ";
  static const int nmax = sizeof(spaces)/sizeof(char)-1;
  int nsent;  /* Total number of characters sent */
  int nnew;   /* Number of characters to send next */
/*
 * Loop until the required number of spaces have been output, or an error
 * occurs.
 */
  for(nsent=0; npad > nsent; nsent += nnew) { 
    nnew = npad - nsent;
    if(nnew > nmax) nnew = nmax;
    if(output(out, spaces, nnew))
      return 1;
  };
  return 0;
}

/*.......................................................................
 * lprint error function.
 *
 * Input:
 *  format   const char *  The printf format string to use.
 *  ...                    Variable argument list.
 * Output:
 *  return          int    -1.
 */
static int lperror(const char *format, ...)
{
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr, format, ap);
  va_end(ap);
  return -1;
}
