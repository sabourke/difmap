#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "vlbutil.h"
#include "logio.h"

/*.......................................................................
 * Given an input string and a length, place a '\0' terminator after the
 * last printable character (except space) to terminate the string there.
 * If the last character up to the length given is a non-space character
 * then that is overwritten with the '\0' terminator. Thus you should
 * send length of the string including space for the terminator.
 *
 * Input/Output:
 *   instr  char *   The input string of length >= slen.
 * Input:
 *   slen   int      The max length of the new string including the
 *                   character for the '\0' terminator.
 * Output:
 *   return int      The length of the reduced version of 'instr'
 *                   excluding the '\0' terminator.
 */
int termstr(char *instr, int slen)
{
  static int i;
  for(i=slen-1; i>=0 && !isgraph((int)instr[i]); i--);
/*
 * We have gone one character too far.
 */
  if(i<slen-1)
    i++;
  instr[i] = '\0';
  return i;
}

/*.......................................................................
 * Undo the effect of termstr(), by reaplacing the '\0' and following
 * characters in the string by blanks.
 *
 * Input/Output:
 *   instr  char *   The input string of length == slen.
 * Input:
 *   slen   int      The number of characters in the character array.
 */
void fillstr(char *instr, int slen)
{
  int i;
/*
 * Locate the '\0' terminator.
 */
  for(i=0; i<slen && instr[i]!='\0'; i++);
/*
 * If found before end of character array then fill the rest of the array
 * with blanks.
 */
  if(i < slen) {
    for(; i<slen; i++)
      instr[i]=' ';
  };
  return;
}

/*.......................................................................
 * Given a char array of ncmax elements, copy up to ncmax-1 chars of an
 * input string into it and terminate it with a '\0' after the last
 * printable (not including spaces) character in the resulting copy.
 *
 * Input:
 *  ostr   char *   The output string array of size ncmax.
 *  istr   char *   The string to be copied.
 *  ncmax   int     The number of char's available in ostr[], including
 *                  room for a '\0'. ie. If ostr is declared char ostr[n]
 *                  send ncmax=n.
 * Output:
 *  return char *   ostr.
 */
char *termcpy(char *ostr, const char *istr, int ncmax)
{
/*
 * Copy up to ncmax-1 chars into ostr[].
 */
  strncpy(ostr, istr, ncmax-1);
/*
 * Terminate the resulting string after the last non-space char.
 */
  termstr(ostr, ncmax);
  return ostr;
}

/*.......................................................................
 * Make a copy of a string with leading and trailing white-space removed.
 *
 * Input:
 *  ostr    char *  The output string, of declared dimension [nco].
 *  nco      int    The declared dimension of ostr[]. No more than nco-1
 *                  characters will be copied into ostr[], before
 *                  terminating the resulting copy with a '\0'.
 *  istr    char *  The string to be copied.
 *  nci      int    The max number of char's to use from istr[].
 * Output:
 *  return  char *  The copied string.
 */
char *stripcpy(char *ostr, int nco, const char *istr, int nci)
{
  int slen;   /* Length of tail of istr[] up to last white-space char */
  char *cptr; /* Pointer into istr[] */
  int ci;     /* Index of last used char in istr[] */
  int co;     /* Index of last used char in ostr[] */
/*
 * Bad output string.
 */
  if(ostr==NULL || nco<=0) {
    lprintf(stderr, "stripcpy: No output string provided.\n");
    return NULL;
  };
/*
 * Nothing to copy?
 */
  if(istr==NULL || nci<=0) {
    *ostr = '\0';
    return ostr;
  };
/*
 * Skip leading white-space in the input string.
 */
  for(ci=0; ci<nci && *istr && isspace((int)*istr); istr++,ci++);
/*
 * Locate the last printable character (except space, in the input
 * string.
 */
  for(slen=0,cptr = (char *)istr; ci<nci && *cptr; cptr++,ci++) {
    if(isgraph((int)*cptr))
      slen = (cptr - istr) + 1;
  };
/*
 * If copying slen characters would overflow ostr[] limit it.
 */
  if(slen > nco-1)
    slen = nco - 1;
/*
 * Copy slen characters from istr to ostr.
 */
  for(co=0, cptr=ostr; co<slen; co++)
    *cptr++ = *istr++;
/*
 * Terminate the result.
 */
  *cptr = '\0';
/*
 * Return the copied string.
 */
  return ostr;
}

/*.......................................................................
 * Remove both leading and trailing white-space from a string and
 * terminate the result.
 *
 * Input:
 *  istr    char *  The '\0' terminated string to be stripped.
 *  nci      int    The max number of char's to use from istr[].
 * Output:
 *  return  char *  The modified string.
 */
char *stripstr(char *istr, int nci)
{
  char *orig = istr;  /* Pointer to char to be copied in istr[] */
  char *dest = istr;  /* Pointer to destination of copy in istr[] */
  int ci;             /* The index of the last used char used in istr[] */
  int slen;           /* Length of copied string */
/*
 * NULL string?
 */
  if(istr==NULL || nci<=0) {
    lprintf(stderr, "stripstr: NULL input string.\n");
    return NULL;
  };
/*
 * Skip leading white-space.
 */
  for(ci=0; ci<nci && *orig && isspace((int)*orig); orig++,ci++);
/*
 * Copy the remaining characters, and maintain a record of the length
 * of the copied string including the last printable, no-space character.
 */
  for(slen=0; ci<nci && *orig; ci++) {
    if(isgraph((int)*orig))
      slen = (dest - istr) + 1;
    *dest++ = *orig++;
  };
/*
 * Terminate the string after the last printable, non-space character.
 */
  istr[slen] = '\0';
/*
 * Return the stripped string.
 */
  return istr;
}

/*.......................................................................
 * Write a string to a given stream using lprintf(), in a form acceptible
 * as a difmap string argument. This requires that special characters
 * be escaped.
 *
 * Input:
 *  fp     FILE *  The stream to write to.
 *  fname  char *  The name of the file to which fp writes, or NULL
 *                 if I/O errors shouldn't be reported.
 *  string char *  The string to write.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int write_string_arg(FILE *fp, char *fname, char *string)
{
  char *sptr;     /* A pointer into string[] */
  int waserr = 0; /* True after the first I/O error */
/*
 * Check the arguments.
 */
  if(!fp || !string) {
    lprintf(stderr, "write_string_arg: NULL argument(s).\n");
    return 1;
  };
/*
 * Initiate the string with speech marks.
 */
  waserr = lputc('\"', fp) == EOF;
/*
 * Examine each character of the string, escaping special characters and
 * writing the result to *fp.
 */
  for(sptr=string; *sptr && !waserr; sptr++) {
    char *cstr = NULL; /* An escape sequence to be printed instead of c */
    int c = *sptr;
    switch(c) {
    case '\a':
      cstr = "\\a";
      break;
    case '\b':
      cstr = "\\b";
      break;
    case '\\':
      cstr = "\\\\";
      break;
    case '\n':
      cstr = "\\n";
      break;
    case '\t':
      cstr = "\\t";
      break;
    case '\r':
      cstr = "\\r";
      break;
    case '\"':
      cstr = "\\\"";
      break;
    case '%':
      if(cstr[1]=='%')  /* Difmap interpretts %% as a preprocessor directive.*/
	cstr = "%\\";
      break;
    default:
      cstr = NULL;
      break;
    };
/*
 * If the character was replaced by an escape sequence, write the
 * escape sequence.
 */
    if(cstr) {
      waserr = lputs(cstr, fp) == EOF;
    } else {
/*
 * Only display printable characters.
 */
      waserr = isprint(c) && lputc(c, fp) == EOF;
    };
  };
/*
 * Terminate the string with speech marks.
 */
  waserr = waserr || lputc('\"', fp) == EOF;
/*
 * I/O error?
 */
  if(waserr && fname) {
    lprintf(stderr, "Error writing to file \"%s\".\n", fname);
    return 1;
  };
  return 0;
}
