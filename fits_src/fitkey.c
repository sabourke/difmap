#include <ctype.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "sysfits.h"
#include "fits.h"
#include "utils.h"
#include "fitkey.h"

static int rfitscom(char *hline, char *string, char **endp);
static int rfitsstr(char *hline, char *string, char **endp);
static int rfitsnum(char *hline, Keyval *value, Fitkey *key, char **endp);

static int chkkey(char *hline, Fitkey *key);
static int wdkey(double dval, char *string, int width);
static int wskey(char *sval, char *string, int minend, int maxend);

static void knamcpy(char *dest, const char *orig);
static int keymatch(const char *name, Fitkey *key);

/*.......................................................................
 * Get a named keyword-value pair from a FITS header.
 * A search for the keyword may be invoked. If so and the key is found,
 * a descriptor will be returned for it.
 * If an error occurs or no match with the keyword is found,
 * the header line pointer will be left as it was on entry to this
 * function. Otherwise, it will be advanced to the line following the
 * line on which the keyword was read from.
 *
 * Input:
 *  fits       Fits *  The FITS file descriptor.
 *  hdu         Hdu *  The base-class pointer to the required HDU.
 *  match      char *  The key name to be sought. Use NULL if any keyword
 *                     will do.
 *  type     Fitkey    The type of value expected. (Ignored if match==NULL).
 *  doseek Seektype    Specifies how the function should search for the
 *                     keyword.
 *                      NO_SEEK   - Only look at the next header line.
 *                      EOH_SEEK  - Search forwards for the keyword -
 *                                  stopping at the end of header (EOH).
 *                      LOOP_SEEK - As EOH_SEEK but if still not found
 *                                  by the EOH search from the start of
 *                                  the header, stopping before the
 *                                  initial line of the first pass.
 * Output:
 *  key     Fitkey *   Pointer to key descriptor to be filled.
 *  return  Keystat    Return status:
 *                      KEY_FOUND   -  Keyword read successfully.
 *                      KEY_EOH     -  Read past end of header.
 *                      KEY_UNKNOWN -  Search for keyword in list failed.
 *                      KEY_BAD     -  Error reading file.
 */
Keystat get_key(Fits *fits, Hdu *hdu, char *match, Fittype type,
		Seektype doseek, Fitkey *key)
{
  static Fitkey kval;  /* Used to compile temporary descriptor for key */
  Fitkey *kptr;        /* Pointer to key-word-value match key */
/*
 * Fill the fields of the input key descriptor.
 */
  if(match) {
    knamcpy(kval.name, match);
    kval.extn = 0;
    kval.keyid = 0;
    kval.type = type;
    kval.value = NULL;
    kval.comment =NULL;
    kptr = &kval;
  } else {
    kptr = NULL;
  };
  return next_key(fits, hdu, kptr, 1, doseek, key);
}

/*.......................................................................
 * Get the next keyword value pair from a FITS header.
 * A search for keywords from a given list may be invoked. If so and a
 * key in the list is found, a modified copy of its key descriptor will
 * be returned. If no list is provided a descriptor of the next keyword
 * in the header will be returned.
 * If an error occurs or no match with the list of keywords is found,
 * the header line pointer will be left as it was on entry to this
 * function. Otherwise, it will be advanced to the line following the
 * line on which the keyword was read from.
 *
 * Input:
 *  fits       Fits *  The FITS descriptor.
 *  hdu         Hdu *  The base-class pointer to the required HDU.
 *  keys     Fitkey *  Array of 'nkey' FITS key descriptors - send NULL
 *                     if any key will do.
 *  nkey        int    The number of entries in the 'keys' array.
 *  doseek Seektype    Specifies how the function should search for the
 *                     keyword.
 *                      NO_SEEK   - Only look at the next header line.
 *                      EOH_SEEK  - Search forwards for the keyword -
 *                                  stopping at the end of header (EOH).
 *                      LOOP_SEEK - As EOH_SEEK but if still not found
 *                                  by the EOH search from the start of
 *                                  the header, stopping before the
 *                                  initial line of the first pass.
 * Output:
 *  key      Fitkey *  Pointer to existing keyword-value pair descriptor.
 *                     A copy of the matched keys[] descriptor is normally
 *                     returned, with the value and comment fields pointing
 *                     to static internal buffers. If keys==NULL the key
 *                     attributes are determined from the header line.
 *  return  Keystat    Return status:
 *                      KEY_FOUND   -  Keyword read successfully.
 *                      KEY_EOH     -  Read past end of header.
 *                      KEY_UNKNOWN -  Search for keyword in list failed.
 *                      KEY_BAD     -  Error reading file.
 */
Keystat next_key(Fits *fits, Hdu *hdu, Fitkey *keys, int nkey, Seektype doseek,
	     Fitkey *key)
{
  Keystat kstat;       /* Return status */
  int initline;        /* Line from which search was initiated */
/*
 * Search forward for the keyword.
 */
  initline = hdu->nextline;
  while((kstat = read_key(fits, hdu, keys, nkey, key))==KEY_UNKNOWN && doseek);
/*
 * At end of file and keyword not found?
 * If requested, search again from the beginning of the header up to the
 * start line of the above search.
 */
  if(kstat==KEY_EOH && doseek==2) {
    hdu->nextline = 0;
    while((kstat=read_key(fits, hdu, keys, nkey, key))==KEY_UNKNOWN &&
	  hdu->nextline < initline);
  };
/*
 * If there was no match - restore the initial header line.
 */
  if(kstat!=KEY_FOUND)
    hdu->nextline = initline;
  return kstat;
}

/*.......................................................................
 * Read the next valid FITS keyword value pair from the FITS file.
 *
 * Input:
 *  fits      Fits *  The FITS file descriptor.
 *  hdu        Hdu *  The base-class pointer to the required HDU.
 *  keys    Fitkey *  An array of 'nkey' recognised keyword-value
 *                    descriptors to select keywords from. Send NULL
 *                    if any keyword will do.
 *  nkey       int    The number of keyword descriptors in 'keys[]'.
 * Output:
 *  key     Fitkey *  Pointer to existing keyword-value pair descriptor.
 *                    A copy of the matched keys[] descriptor is normally
 *                    returned, with the value and comment fields pointing
 *                    to internal values. If keys==NULL or nkey<1 the
 *                    keyword-value attributes are determined from the
 *                    header line. NOTE - the key->value and key->comment
 *                    fields point to static internal buffers.
 *  return Keystat    Error status - See fitkey.h for details.
 */
Keystat read_key(Fits *fits, Hdu *hdu, Fitkey *keys, int nkey, Fitkey *key)
{
  Fitkey *kmatch;/* Descriptor of matching keyword */
  char kword[9]; /* Keyword name read */
  char *hline;   /* Pointer to input line buffer */
  char *nextc;   /* Pointer to start of value */
  char *cptr;    /* General purpose pointer into string */
  int waseoh=0;  /* True if the new line is the end of header line */
  int ikey;      /* Index of matching key in keys[] */
  int domatch=0; /* True if a list of keys to match against was given */
  static Keyval value;          /* Keyword value union */
/*
 * Past end of header - NB this shouldn't be possible.
 */
  if(hdu->endline>=0 && hdu->nextline > hdu->endline)
    hdu->nextline = hdu->endline;
/*
 * Do we have a list of keys to match against?
 */
  domatch = keys!=NULL && nkey>0;
/*
 * Read the next line from the FITS header.
 */
  hline = rheadline(fits, hdu, hdu->nextline);
  if(hline==NULL)
    return KEY_BAD;
/*
 * Get a terminated copy of the keyword, including its extension
 * number (we don't yet know which of any digits that follow the
 * keyword name, are part of the name, and which are part of an extension).
 */
  knamcpy(kword, hline);
/*
 * Check that the keyword is a valid FITS keyword.
 * The FITS standard specifies that all characters in keywords must be
 * either upper-case alphabetical, underscore, hyphen or digits.
 */
  for(cptr = kword; *cptr; cptr++) {
    int c = *cptr;
    if(islower(c) && c != '_' && c != '-' && !isdigit(c)) {
      fprintf(stderr,"Illegal characters found in FITS keyword \'%.8s\'\n",
	      hline);
      return KEY_BAD;
    };
  };
/*
 * Trap the END card of the header.
 */
  if(strcmp(kword, "END") == 0) {
    waseoh = 1;
/*
 * If this is the first time that the END keyword has been seen record
 * its position.
 */
    if(hdu->endline<0)
      hdu->endline = --hdu->nextline;
  };
/*
 * Compare the keyword against the list of keywords supplied.
 */
  if(domatch) {
    for(ikey=0; ikey<nkey && !keymatch(kword, &keys[ikey]); ikey++);
    if(ikey>=nkey)
      return waseoh ? KEY_EOH : KEY_UNKNOWN;
    kmatch = &keys[ikey];
    *key = *kmatch;        /* Copy descriptor of located keyword */
  } else {
    kmatch = NULL;
    knamcpy(key->name, kword);
    key->keyid = 0;
    key->extn = 0;
  };
/*
 * No value or comment yet.
 */
  key->value = NULL;
  key->type = DAT_NON;
  key->comment = NULL;
/*
 * Determine the keyword type.
 *
 * A keyword with a value?
 */
  if(kword[0]!='\0' && hline[8]=='=' && hline[9]==' ') {
/*
 * Skip blanks up to the value.
 */
    nextc = &hline[10];
    while(*nextc == ' ')
      nextc++;
/*
 * String keyword?
 */
    if(*nextc == '\'') {
      key->type = DAT_STR;
      if(rfitsstr(hline, value.sval, &nextc))
	return KEY_BAD;
      key->value = (void *) value.sval;
    }
/*
 * Logical keyword?
 */
    else if(*nextc == 'T' || *nextc == 'F') {
      key->type = DAT_LOG;
      value.bval = *nextc++;
      key->value = (void *) &value.bval;
    }
/*
 * Numeric keyword?
 */
    else if(isdigit((int) *nextc) ||
	    *nextc=='.' || *nextc=='-' || *nextc=='+' || *nextc=='(') {
      if(rfitsnum(hline, &value, key, &nextc))
	return KEY_BAD;
    }
/*
 * Unknown value type?
 */
    else {
      fprintf(stderr,"Unable to determine type of value in header line:\n%s\n",
	      hline);
      return KEY_BAD;
    };
  }
/*
 * A history or comment line?
 */
  else {             /* HISTORY, COMMENT or keyword with no value */
    nextc = &hline[8];
/*
 * History or comment line?
 */
    if((kmatch && kmatch->type==DAT_COM) ||
       strcmp(kword, "HISTORY")==0 ||
       strcmp(kword, "COMMENT")==0 ||
       *kword=='\0') {
      key->type = DAT_COM;
      if(rfitscom(hline, value.sval, &nextc))
	return KEY_BAD;
      key->value = (void *) value.sval;
    } else {
      key->type = DAT_NON;
    };
  };
/*
 * Get the comment.
 */
  while(*nextc != '\0' && *nextc++!='/');
  key->comment = nextc;
/*
 * Does the value-type read agree with that expected for the given keyword?
 */
  if(kmatch && kmatch->type != key->type) {
/*
 * Try to convert the value to the expected type.
 */
    if(typeconv(1L, key->type, &value, 0.0f, 1.0f, kmatch->type, &value)) {
      fprintf(stderr,
	"read_key: Implicit type conversion failed for header line:\n%s\n",
	hline);
      return KEY_BAD;
    } else if(fits->pedantic) { 
      fprintf(stderr,
	"Warning: Implicit conversion -> (%s) applied on header line:\n%s\n",
	typename(kmatch->type), hline);
    };
    key->type = kmatch->type;
  };
  return KEY_FOUND;
}

/*.......................................................................
 * Read a value of a numeric FITS header parameter. (FORTRAN free format).
 *
 * Input:
 *  hline   char *  The input header line.
 * Output:
 *  value Keyval *  Send a pointer to an existing Keyval union of
 *                  possible fits keyword value types.
 *  key   Fitkey *  The type and value fields will be appropriately assigned.
 *  endp    char ** *endp will point to the next unprocessed character
 *                  in 'hline'.
 *  return   int    0 - OK.
 */
static int rfitsnum(char *hline, Keyval *value, Fitkey *key, char **endp)
{
  char *cptr;      /* Pointer into 'hline' */
  char *start;     /* Pointer to first char in number */
  char *keep;      /* Temporary storage of cptr value */
  int isflt=0;     /* True if number is a floating point number */
  double dval[2];  /* Decimal values of up to two numbers */
  int ndigit;      /* Number of digits read */
  int inum;        /* The index of the number being read */
  Fittype kt[2];   /* The type of numbers read */
/*
 * FITS header values start in or after column 11.
 */
  cptr = &hline[10];
/*
 * Read up to two numbers.
 */
  for(inum=0; inum<2; inum++) {
/*
 * The new number is not a floating point number until proved otherwise!
 */
/*
 * Skip leading spaces.
 */
    while(*cptr == ' ')
      cptr++;
/*
 * Record the start position. just in case the number turns out not to
 * be a number.
 */
    *endp = cptr;
/*
 * A '(' may preceed the first number, and a ',' preceed the second.
 */
    if( (inum==0 && *cptr == '(') || (inum==1 && *cptr==',') ) {
      cptr++;
      while(*cptr == ' ')    /* Skip spaces after '(' or ',' */
	cptr++;
    };
/*
 * Record the position of the first character in the number itself.
 */
    start = cptr;
/*
 * Skip the sign.
 */
    if(*cptr=='+' || *cptr=='-')
      cptr++;
/*
 * Skip the mantissa.
 */
    ndigit = 0;
    while(isdigit((int)*cptr)) /* Integer before decimal point */
      cptr++,ndigit++;
    if(*cptr == '.') {     /* Decimal point */
      isflt = 1;
      cptr++;
    };
    while(isdigit((int)*cptr)) /* Integer after decimal point */
      cptr++,ndigit++;
/*
 * Not a number?
 */
    if(ndigit==0) {
      cptr = *endp;
      break;
    };
/*
 * Is there an exponent?
 */
    switch(*cptr) {
    case 'D': case 'E': case 'd': case 'e':
      keep = cptr++;
      ndigit = 0;
      if(*cptr=='-' || *cptr=='+') /* Skip the sign */
	cptr++;
      while(isdigit((int)*cptr))  /* Skip the exponent integer */
	cptr++,ndigit++;
/*
 * Was there a valid exponent?
 */
      if(ndigit>0) {
	*keep = 'E';         /* Ensure that FORTRAN 'D' format is translated */
	isflt = 1;           /* The exponent implies a floating point type */
      } else {
	cptr = keep;
      };
      break;
    };
/*
 * Record the data type and read the number.
 */
    kt[inum] = (isflt) ? DAT_DBL:DAT_INT;
    dval[inum] = strtod(start, NULL);
  };
/*
 * How many numbers were read and what does this imply?
 *
 * Scalar number?
 */
  if(inum==1) {
    key->type = kt[0];
    switch(key->type) {
    case DAT_DBL:
      value->dval = dval[0];
      key->value = (void *) &value->dval;
      break;
    case DAT_INT:
      value->ival = dval[0];
      key->value = (void *) &value->ival;
      break;
    default:
      fprintf(stderr, "rfitsnum: Unhandled type (%s) in switch statement\n",
	      typename(key->type));
      return 1;
    };
/*
 * Complex number?
 */
  } else if(inum==2) {
    key->type = DAT_SCMP;
    value->cval[0] = dval[0];
    value->cval[1] = dval[1];
    key->value = (void *) &value->cval;
/*
 * No numbers read or inconsistent complex number.
 */
  } else {
    fprintf(stderr, "Garbled numeric value on header line:\n%s\n", hline);
    return 1;
  };
/*
 * Skip any trailing spaces and closing ')' character.
 */
  while(*cptr == ' ')
    cptr++;
  if(*cptr == ')') {
    while(*cptr == ' ')
      cptr++;
  };
  *endp = cptr;
  return 0;
}

/*.......................................................................
 * Read the value of a FITS string header parameter. (FORTRAN free format).
 *
 * Input:
 *  hline   char *  The input header line.
 * Output:
 *  string  char *  Send a pointer to an existing char array of at least
 *                  70 elements.
 *  endp    char ** *endp will point to the next unprocessed character
 *                  in 'hline'.
 *  return   int    0 - OK.
 */
static int rfitsstr(char *hline, char *string, char **endp)
{
  char *cptr;     /* Pointer into 'hline[]' */
  char *sptr;     /* Pointer into 'string[]' */
/*
 * FITS header values start in or after column 11.
 */
  cptr = &hline[10];
/*
 * Skip spaces.
 */
  while(*cptr == ' ')
    cptr++;
/*
 * Record the start position, just in case the string turns out to be
 * badly formed.
 */
  *endp = cptr;
/*
 * The first character must be an open quote.
 */
  if(*cptr != '\'') {
    fprintf(stderr, "Expected a \'string-value\' value on header line:\n%s\n",
	    hline);
    return 1;
  };
/*
 * Copy the string into 'string[]' stopping at the first single quote.
 * Note that two adjacent quote characters denote a single quote
 * in the string and not the terminating quote.
 */
  cptr++;
  sptr = &string[0];
  while(*cptr && (*cptr!= '\'' || *(++cptr)=='\''))
    *sptr++ = *cptr++;
/*
 * Terminate the string.
 */
  *sptr = '\0';
  *endp = cptr;  /* This points at the next character after the end quote */
  return 0;
}

/*.......................................................................
 * Read the value of a FITS HISTORY or COMMENT line.
 *
 * Input:
 *  hline   char *  The input header line.
 * Output:
 *  string  char *  Send a pointer to an existing char array of at least
 *                  73 elements.
 *  endp    char ** *endp will point to the next unprocessed character
 *                  in 'hline'.
 *  return   int    0 - OK.
 */
static int rfitscom(char *hline, char *string, char **endp)
{
  char *cptr;     /* Pointer into 'hline[]' */
  char *sptr;     /* Pointer into 'string[]' */
  char *eptr;     /* Pointer to last non-white-space character in output */
/*
 * FITS comment values start in column 9.
 */
  cptr = &hline[8];
/*
 * Copy the string (including leading spaces) into 'string[]'.
 */
  sptr = eptr = &string[0];
  while(*cptr) {
/*
 * Check for last non-white-space character.
 */
    if(!isspace((int)*cptr))
      eptr = sptr;
/*
 * Copy character.
 */
    *sptr++ = *cptr++;
  };
/*
 * Terminate the string.
 */
  *(eptr+1) = '\0';
  *endp = cptr;
  return 0;
}

/*.......................................................................
 * Write a FITS keyword into a given buffer. In the process check for illegal
 * characters and keywords that are too long, convert to upper case,
 * add any numeric postfix and add the '=' for value'd keywords.
 *
 * Input:
 *  hline   char *   Pointer to the 80 char header line buffer to be
 *                   written into.
 *  key   Fitkey *   The descriptor of the keyword-value pair.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Illegal characters detected.
 */
static int chkkey(char *hline, Fitkey *key)
{
  char *kptr = &key->name[0];   /* Pointer into 'key->name' */
  char *hptr = &hline[0];       /* Pointer into 'hline' */
  int slen;    /* Length of composed keyword. */
  int ipos;    /* Position in output string */
/*
 * Blank keyword - this is legal for comments.
 */
  if(key->name[0] == '\0')
    return 0;
/*
 * The first character must be an alphabetical character.
 */
  if(!isalpha((int) *kptr)) {
    fprintf(stderr,"Illegal character at start of FITS keyword \'%s\'\n",kptr);
    return 1;
  };
/*
 * All the other characters must be either alphabetical, underscore
 * or numbers.
 */
  for(ipos=0; ipos<8 && *kptr != '\0'; ipos++) {
    if(!isalnum((int) *kptr) && *kptr != '_' && *kptr != '-') {
      fprintf(stderr, "Illegal character %c in FITS keyword \'%s\'\n",
	      *kptr, key->name);
      return 1;
    };
    *(hptr++) = toupper((int) *(kptr++));
  };
/*
 * If a postfix integer was sent, write it.
 */
  if(key->extn > 0)
    sprintf(hptr,"%d", key->extn);
  else
    *hptr = '\0';
/*
 * Check that the keyword was written in <= 8 characters.
 */
  slen = strlen(hline);
  if(slen > 8) {
    fprintf(stderr, "Illegal FITS keyword > 8 characters: %s\n", hline);
    return 1;
  };
/*
 * Remove the '\0' terminator.
 */
  hline[slen] = ' ';
/*
 * If this is a keyword of a keyword/value pair, place an '=' in
 * column 9.
 */
  if(key->type != DAT_NON && key->type != DAT_COM)
    hline[8] = '=';
  return 0;
}

/*.......................................................................
 * Write a keyword-value pair + comment to a FITS header.
 *
 * Input:
 *  fits   Fits  *   The FITS file descriptor.
 *  hdu     Hdu  *   The descriptor of the HDU to be written to.
 *  key  Fitkey  *   The descriptor of the keyword-value pair to be
 *                   written.
 * Output:
 *  return  int      0 - OK.
 *                   1 - Error.
 */
int putkey(Fits *fits, Hdu *hdu, Fitkey *key)
{
  static char hline[81];      /* Header line to be written */
  static const int compos=40; /* Default column to start comment in */
  int toobig=0;               /* True if the value exceeds the field width */
  int endpos;                 /* Index of char following the value field */
  char *cptr;
  int slen;
/*
 * Sanity checks.
 */
  if(fits==NULL || hdu==NULL || key==NULL) {
    fprintf(stderr, "putkey: NULL argument intercepted\n");
    return 1;
  };
/*
 * Fill the header line output string with spaces and terminate.
 */
  memset(hline, ' ', (size_t) 80);
  hline[80] = '\0';
/*
 * Write the keyword (if there is one), its extension number and if
 * the keyword takes a value, write '=' in column 9 of hline[].
 */
  if(chkkey(hline, key))
    return 1;
/*
 * Blank keywords are only legal as comment lines.
 */
  if(key->name[0]=='\0' && key->type!=DAT_COM) {
    fprintf(stderr, "putkey: Blank keyword but not a comment?\n");
    return 1;
  };
/*
 * Write the keyword value.
 */
  switch(key->type) {
  case DAT_NON:
    endpos = compos-1;
    break;
  case DAT_INT: /* Integer - right-justified in columns 11->30 */
    sprintf(&hline[10], "%20d", KEYINT(*key));
    toobig = hline[30]!='\0';
    hline[30] = ' ';
    endpos = 30;
    break;
  case DAT_FLT: /* Real keyword right-justified in columns 11->30 */
    toobig = wdkey((double) KEYFLT(*key), &hline[10], 20);
    endpos = 30;
    break;
  case DAT_DBL: /* Real keyword right-justified in columns 11->30 */
    toobig = wdkey(KEYDBL(*key), &hline[10], 20);
    endpos = 30;
    break;
  case DAT_LOG: /* Logical T or F character in column 30 */
    hline[29] = KEYBOOL(*key)=='T' ? 'T' : 'F';
    endpos = 30;
    break;
  case DAT_SCMP: /* (Complex) Two reals in columns 11->30 and 31->50 */
    toobig = wdkey((double) KEYCMP(*key)[0], &hline[10], 20) ||
             wdkey((double) KEYCMP(*key)[1], &hline[30], 20);
    endpos = 50;
    break;
  case DAT_STR:
    endpos = 10 + wskey(KEYSTR(*key), &hline[10], 10, 70);
    toobig = endpos==10;
    break;
  case DAT_COM:  /* Comment within columns 9->80 */
    cptr = KEYSTR(*key);
    for(endpos=8; endpos<80 && *cptr; endpos++)
      hline[endpos] = *cptr++;
    break;
  default:
    fprintf(stderr, "putkey: Unhandled keyword-value type (%s)\n",
	    typename(key->type));
    return 1;
  };
/*
 * Did the value fit into its FITS proscribed field width.
 */
  if(toobig) {
    fprintf(stderr, "Value too big for field width of key %s\n", key->name);
    return 1;
  };
/*
 * Work out how long the header line has to be to hold the comment.
 */
  slen = endpos + 3 + (key->comment ? strlen(key->comment) : 0);
/*
 * Write the comment if there is room.
 */
  if(key->comment && key->comment[0]!='\0') {
    if(slen > 80) {
      fprintf(stderr,
	    "Warning: Insufficient room to write comment on header line:\n%s\n",
	    hline);
    } else {
      if(endpos < compos && slen+1+compos-endpos <= 80)
	endpos = compos;
      sprintf(&hline[endpos+1], "/%s", key->comment ? key->comment:"");
/*
 * Strip the extra '\0' placed in hline by sprintf().
 */
      slen = strlen(hline);
      if(slen < 80)
	hline[slen]=' ';
    };
  };
/*
 * Write the header line to the FITS file.
 */
  return wheadline(fits, hdu, hdu->wnxtline, hline);
}

/*.......................................................................
 * Write a double keyword value into a given number of characters, using
 * either normal or exponential notation as appropriate. NB. sprintf()
 * %G format can not be used here because it is allowed to suppress the
 * decimal point.
 *
 * Input:
 *  dval   double   The value to be written.
 *  string   char * The string to be written into. This should contain
 *                  sufficient characters to write any reasonably sized
 *                  float ie possibly > 'width'.
 *  width     int   The width of the field to write the number - right
 *                  justified - into.
 * Output:
 *  return    int   0 - OK.
 *                  1 - Number will not fit in the given field width.
 */
static int wdkey(double dval, char *string, int width)
{
  size_t slen;
/*
 * Write the number right-justified in the appropriate notation.
 */
  sprintf(string, "%#*.*G", width, width-6, dval);
/*
 * Find the length of the resulting string.
 */
  slen = strlen(string);
  if(slen > width)
    return 1;
/*
 * Remove the '\0' inserted by sprintf().
 */
  string[slen] = ' ';
  return 0;
}

/*.......................................................................
 * Write the string value of a string keyword into a given field width.
 * This includes placing quotes in the appropriate places and escaping
 * literal quotes in the string value.
 *
 * Input:
 *  sval   char *  The string value to be written.
 *  string char *  The string to be written into - this must contain
 *                 at least maxend+1 characters.
 *  minend  int    The mimimum field width (including quotes).
 *  maxend  int    The maximum field width (including quotes).
 * Output:
 *  return  int    The width of the field, or 0 if an error occured.
 */
static int wskey(char *sval, char *string, int minend, int maxend)
{
  int apos=0;  /* Index into input string 'sval' */
  int bpos=0;  /* Index into output string 'string' */
/*
 * Copy 'sval' to 'string', surrounded by single quotes.
 * Where "'" is found in the string, escape it with another quote
 * as required by the standard.
 */
  string[bpos++] = '\'';
  for( ; sval[apos]!='\0' && bpos<maxend; apos++, bpos++) {
    string[bpos] = sval[apos];
    if(sval[apos]=='\'' && ++bpos<maxend)    /* Escape a literal quote */
      string[bpos] = '\'';
  };
/*
 * Check for overflow of the field width.
 */
  if(bpos>=maxend)
    return 0;
/*
 * Place the closing quote in the appropriate position.
 */
  if(bpos < minend-1)
    bpos = minend-1;
  string[bpos++] = '\'';
  return bpos;
}

/*.......................................................................
 * Write an integer valued keyword-value header line to a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value    int     The value to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wintkey(Fits *fits, Hdu *hdu, char *name, int extn, int value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_INT;
  key.value = &value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write a floating point valued keyword-value header line to a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value double     The value to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wfltkey(Fits *fits, Hdu *hdu, char *name, int extn, double value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_DBL;
  key.value = &value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write a logical valued keyword-value header line to a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value   char     The value to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wlogkey(Fits *fits, Hdu *hdu, char *name, int extn, char value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_LOG;
  key.value = &value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write a single precision-complex valued keyword-value header line to
 * a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value  float *   Pointer to the array of two floats to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wcmpkey(Fits *fits, Hdu *hdu, char *name, int extn, float *value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_SCMP;
  key.value = value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write an string valued keyword-value header line to a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value   char *   Pointer to a terminated string array to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wstrkey(Fits *fits, Hdu *hdu, char *name, int extn, char *value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_STR;
  key.value = value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write a comment keyword-value header line to a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  value   char *   Pointer to a terminated string array to be written.
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wcomkey(Fits *fits, Hdu *hdu, char *name, int extn, char *value,
	    char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_COM;
  key.value = value;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Write the header line of a keyword that takes no value, into a FITS file.
 *
 * Input:
 *  fits    Fits *   The FITS file descriptor.
 *  hdu      Hdu *   The descriptor of the HDU to be written to.
 *  name    char *   The keyword name.
 *  extn     int     Extension number of keyword (send 0 if not required).
 *  comment char *   The comment to follow the value, or NULL if not required.
 * Output:
 *  return   int     0 - OK.
 *                   1 - Error.
 */
int wvoidkey(Fits *fits, Hdu *hdu, char *name, int extn, char *comment)
{
  Fitkey key;
  knamcpy(key.name, name);
  key.extn = extn;
  key.type = DAT_NON;
  key.value = NULL;
  key.comment = comment;
  return putkey(fits, hdu, &key);
}

/*.......................................................................
 * Public function used to specify the next line to be read or written
 * in the header of a given HDU.
 * If the requested line number is out-of-bounds the closest bound is
 * substituted.
 *
 * Input:
 *  hdu    Hdu *   The descriptor of the HDU.
 *  iline  int     The 0-relative line number required. 
 * Output:
 *  return int     The original line-number.
 */
int new_hline(Hdu *hdu, int iline)
{
  int saveline;
/*
 * Sanity checks.
 */
  if(hdu==NULL || hdu->state==HDU_DESCR) {
    fprintf(stderr, "set_line: Bad HDU descriptor recieved.\n");
    return 0;
  };
/*
 * Save the current line number.
 */
  saveline = hdu->nextline;
/*
 * Enforce bounds on the given line number.
 */
  if(iline > hdu->endline)
    iline = hdu->endline;
  else if(iline < 0)
    iline = 0;
/*
 * Install the new nextline value and return the original.
 */
  hdu->nextline = iline;
  hdu->wnxtline = iline;
  return saveline;
}

/*.......................................................................
 * Copy up to 8 keyword characters from a given string and terminate
 * the copy with a '\0'.
 *
 * Input:
 *  dest   char *   The destination array of at least 9 elements.
 *  orig   char *   The string to be copied, terminated with ' '
 *                  and/or '\0' if it contains fewer than 8 characters.
 */
static void knamcpy(char *dest, const char *orig)
{
  int c;  /* The source character to be checked before copying */
  int i;  /* Count of the number of characters written */
/*
 * Copy up to 8 characters from orig[], stopping if ' ' or '\0' is
 * encounterred.
 */
  for(i=0; i<8 && (c = *orig)!='\0' && c!=' '; i++)
    *dest++ = *orig++;
/*
 * Terminate the output string.
 */
  *dest = '\0';
  return ;
}

/*.......................................................................
 * Given a raw name from a FITS header line (including any extension
 * number), and a keyword descriptor check if the name matches the name
 * in the keyword descriptor without regard for trailing numeric digits
 * in the raw name.
 * If it does match, record the trailing digits as the keyword's
 * extension number in key->extn. If there are no trailing digits,
 * place 0 in key->extn.
 *
 * Input:
 *  name   char *  The '\0' terminated name to be matched.
 *  key  Fitkey *  Pointer to the keyword descriptor to be matched against.
 * Output:
 *  return  int    0 - No match.
 *                 1 - Keyword name matches.
 */
static int keymatch(const char *name, Fitkey *key)
{
  const char *nptr;  /* Pointer into name[] */
  const char *kptr;  /* Pointer into key->name[] */
  const char *eptr;  /* Pointer into the keyword extension digits */
/*
 * Find the position at which the given name and the keyword name
 * differ.
 */
  nptr = name;
  kptr = key->name;
  while(*nptr && *kptr && *kptr == *nptr)
    kptr++,nptr++;
/*
 * Reached end of name before reaching end of key->name[]?
 */
  if(*nptr != *kptr && (*nptr=='\0' || *kptr!='\0'))
    return 0;
/*
 * We must be at the end of key->name[], with extra non-matching characters
 * on the end of name[]. If the extra characters are ALL numeric
 * then the number will be treated as an extension number and the
 * keywords regarded as matching.
 */
  eptr = nptr;
  while(*eptr && isdigit((int)*eptr))
    eptr++;
/*
 * Not a valid extension?
 */
  if(*eptr)
    return 0;
/*
 * Record the extension number.
 */
  key->extn = *nptr=='\0' ? 0 : atoi(nptr);
/*
 * Match succeeded.
 */
  return 1;
}
