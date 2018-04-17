#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>

#include "sphere.h"
#include "utils.h"

#define MAXIO 20

typedef struct {
  FILE *fptr;
  char *filename;
  char is_text;
  char is_read;
} Filetable;

#ifndef SEEK_SET
#define SEEK_SET 0
#endif

#define MAXLINE 132
static char io_buff[MAXLINE];

Filetable fat[MAXIO];

FILE *fptr;
int file_lun;
int is_text;

/*.......................................................................
  Initialise the file allocation table such that all file pointers are
  NULL.
*/
void fat_init(void)
{
        int i;
        for(i=0; i<MAXIO; i++)
	  fat[i].fptr = NULL;
	fat[0].is_text = 1;
	return;
}

/*.......................................................................
  Given a user logical unit number, return the file pointer to the
  associated file in the file allocation table. If there is no open file
  associated with that pointer or it is read when write was required
  (or vice versa) then report an error and return. If neither read nor
  write are specifically required then want_read should be 2.
*/
FILE *check_lun(int lun, int want_read, int *is_text)
{
/*
  stdin and stdout are not actually in the file allocation table.
  Handle them separately.
*/
        if(lun == 0) {
	  *is_text = 1;
	  if(want_read)
	    return stdin;
	  else
	    return stdout;
	};
/*
  Check that the requested lun has an open file assigned to it.
*/
        if(lun >= MAXIO || lun < 1 || fat[lun].fptr == NULL) {
	  lprintf(stderr, "No file assigned with lun: %d\n",lun);
	  return NULL;
	};
/*
  Report its text/binary type.
*/
	*is_text = fat[lun].is_text;
/*
  Check that the file is open for data transfer in the required
  direction.
*/
	if(fat[lun].is_read != want_read && want_read != 2) {
	  if(want_read)
	    lprintf(stderr, "lun %d is write-only.\n",lun);
	  else
	    lprintf(stderr, "lun %d is read-only.\n",lun);
	  return NULL;
	};
/*
  No errors - return the required file pointer.
*/
	return fat[lun].fptr;

}
/*.......................................................................
  Utility routine to open a new user file and insert it in the file
  allocation table. It returns the FAT number if succesfull or -1
  otherwise.
*/
int file_open(char read_write_append, char is_text, char *filename)
{
        char mode[4];
	int i;
/*
  Find a free slot in the file allocation table.
*/
        for(i=1;i<MAXIO;i++) {
	  if(fat[i].fptr == NULL) {
/*
  Allocate memory for a copy of the filename, then copy it into the fat entry.
*/
	    if( (fat[i].filename=stralloc(strlen(filename))) == NULL)
	      return -1;
	    strcpy(fat[i].filename, filename);
/*
  Create the mode string.
*/
	    switch(read_write_append) {
	    case 0:
	      mode[0]='r';
	      fat[i].is_read=1;
	      break;
	    case 1:
	      mode[0]='w';
	      fat[i].is_read=0;
	      break;
	    case 2:
	      mode[0]='a';
	      fat[i].is_read=0;
	      break;
	    };
	    mode[1] = (is_text) ? '\0':'b';
	    mode[2] = '\0';
	    fat[i].is_text = is_text;
/*
  Atempt to open the file.
*/
	    if( (fat[i].fptr=fopen(filename, mode)) == NULL) {
	      lprintf(stderr, "Unable to open file: '%s'\n", filename);
	      free(fat[i].filename);
	      return -1;
	    };
	    return i;
	  };
	};
/*
  No free slots in the file allocation table - report error.
*/
	lprintf(stderr, "Sorry - no free slots are available in file allocation table.\n");
	return -1;
}

/*.......................................................................
  Utility routine to close a user file and release its FAT slot. When
  called with a FAT number that has not yet been allocated, nothing
  is done.
*/
int file_close(int lun)
{
/*
  Check that the passed lun has an open file assigned to it.
*/
        if(check_lun(lun, 2, &is_text) == NULL)
	  return -1;
/*
  Make it illegal to close stdin and stdout.
*/
	if(lun == 0) {
	  lprintf(stderr, "Illegal attempt to close stdin and stdout - lun: 0\n");
	  return -1;
	};
	fclose(fat[lun].fptr);
	fat[lun].fptr = NULL;
	free(fat[lun].filename);
	return 0;
}

/*.......................................................................
  Utility routine to rewind a user file, given its lun in the file
  allocation table.
*/
int file_rewind(int lun)
{
/*
  Check that the passed lun has an open (read-only) file assigned to it.
*/
        if(check_lun(lun, 1, &is_text) == NULL)
	  return -1;
	if(lun == 0) {
	  lprintf(stderr, "Illegal attempt to rewind stdin - lun: 0\n");
	  return -1;
	};
	rewind(fat[lun].fptr);
	return 0;
}

/*.......................................................................
  Utility routine to return the status of the EOF flag for a given user
  file. It returns 0 if the flag is not set or 1 if either the flag is
  set, or the lun referenced has no open file allocated to it.
*/
int file_check_eof(int lun)
{
/*
  Check that the passed lun has an open file assigned to it.
*/
        if((fptr=check_lun(lun, 1, &is_text)) == NULL)
	  return -1;
	return(feof(fptr) != 0);
}

/*.......................................................................
  Utility routine to check the error indicator for a given user file,
  report any system error message and then clear the error indicator.
  It returns -1 if the error indicator was set.
*/
int file_error(int lun)
{
/*
  Check that the passed lun has an open file assigned to it.
*/
        if((fptr=check_lun(lun, 2, &is_text)) == NULL)
	  return -1;
/*
  Return normally if no error indicator is set.
*/
	if(ferror(fptr) == 0)
	  return 0;
/*
  Otherwise report the error, clear it and return true.
*/
	perror("File i/o error");
	clearerr(fptr);
	return 1;
}

/*.......................................................................
  Utility routine to list open user files and there attributes.
*/
void file_cat(void)
{
        int i;
        lprintf(stdout, "Catalogue of user files:\n");
	lprintf(stdout, "LUN=0: standard input and output - normally the terminal.\n");
	for(i=1;i<MAXIO;i++) {
	  if(fat[i].fptr != NULL)
	    lprintf(stdout, "LUN=%d: NAME='%s': %s file, %s-only\n",i, fat[i].filename, (fat[i].is_text)?"Text":"Binary", (fat[i].is_read)?"Read":"Write");
	};
	return;
}

/*.......................................................................
  Utility routine to search a user file - specified via its FAT lun -
  for a passed string and position the file pointer immediately after it.
  If the string is found, 1 is returned. If it isn't then 0 is returned
  and if an error occurs -1 is returned.
  The caller may specify that the file position be returned to the start
  of the match by setting leave_start to true (1).
*/
int file_search(FILE *fptr, char string[], size_t slen, int leave_start)
{
	  int c,i;
	  long fpos;
	  if(fptr == stdin) {
	    lprintf(stderr, "Illegal search on stdin - lun: 0\n");
	    return -1;
	  };
/*
  Abort search if zero length string was sent.
*/
	  if(slen < 1) return 1;
/*
  Search until EOF is reached or the string is matched.
*/
	  do {
/*
  Check each character to see if it matches the first character of the
  search string.
*/
	    if( (c=fgetc(fptr)) == *string) {
/*
  Match succesful if the search string is only one character long.
*/
	      if(slen == 1) {
		if(leave_start) ungetc(c,fptr);
		return 1;
	      };
/*
  Keep a record of the current file position such that if the following characters
  don't match, the search can start again from the current position.
*/
	      fpos=ftell(fptr);
/*
  See if the following characters match the required string.
*/
	      for(i=1;;i++) {
		if( (c=fgetc(fptr)) != string[i]) {
/*
  Match failed for the current run, reset the file position to the character
  after the start character of this run so that the search can continue from
  there.
*/
		  fseek(fptr, fpos, SEEK_SET);
		  break;
		};
/*
  On a succesfull match return. First however check where to leave the file
  position. Return the file position to the start of the search match
  if this function was called with leave_start true.
*/
		if(i == slen-1) {
		  if(leave_start)
		    fseek(fptr, fpos-1, SEEK_SET);
		  return 1;
		};
	      };
	    };
	  } while (c != EOF);
/*
  No match before end of file.
*/
	  lprintf(stderr, "Reached end of file before finding '%.*s'\n", slen, string);
	  return 0;
}

/*.......................................................................
  Given a file pointer, a format string, a specification of the number
  of arguments and an array of those arguments (as pointers todescriptors),
  emulate fprintf.
*/
int user_printf(FILE *fptr, char fmt[], int nargs, Descriptor *args[])
{
        static int c;
	static char end_flags, was_point, end_all, was_error;
	static int arg;
	static char io_buff[MAXLINE], *cptr, *str_ptr, *fmt_ptr;
/*
  First argument - args[0].
*/
	arg=0;
/*
  Parse the format string.
*/
	for(str_ptr=cptr=fmt; (c = *cptr) != '\0'; cptr++) {
	  switch (c) {
/*
  Format item specifier?
*/
	  case '%':
	    cptr++;
/*
  %% escapes the first % so ignore it.
*/
	    if(*cptr == '%')
	      break;
/*
  Keep a record of the start of the specifier just in case.
*/
	    fmt_ptr=cptr;
/*
  Clear the specifier parse flags.
*/
	    end_flags = was_point = end_all = was_error = 0;
/*
  A valid specifier consists of '%<flags><number><.><number>type'.
  Where the <> denote optional arguments.
*/
	    for(; (c=*cptr) != '\0'; cptr++) {
	      switch (c) {
	      case '-': case '+': case ' ': case '0': case '#':
		was_error = (end_flags && c != '0');
		break;
	      case '1': case '2': case '3': case '4': case '5':
	      case '6': case '7': case '8': case '9':
		end_flags = 1;
		break;
	      case '.':
		was_error = was_point;
		was_point = end_flags = 1;
		break;
	      case 's': case 'e': case 'f': case 'g': case 'd': case 'i':
		end_all = 1;
		break;
	      default:
		was_error = 1;
	      };
/*
  Report syntax errors in format specifier.
*/
	      if(was_error) {
		lprintf(stderr, "fprintf(,'%s',...): Illegal format specifier: %s\n",fmt, fmt_ptr);
		return -1;
	      };
/*
  When the format specifier has been syntactically checked, check that
  the specifier type matches the user argument type. Then send the part
  of the format string starting at the end of the last format specifier
  and ending at the end of the current specifier, to fprintf together
  with the pertinent argument.
*/
	      if(end_all) {
/*
  First take a copy of the required section of the format string.
*/
		strncpy(io_buff,str_ptr, cptr-str_ptr+1);
		io_buff[cptr-str_ptr+1]='\0';
/*
  Check that there is an argument associated with the current specifier.
*/
		if(arg+1 > nargs) {
		  lprintf(stderr, "fprintf(,'%s',...): More specifiers than arguments?\n",fmt);
		  return -1;
		};
		switch (c) {
/*
  Now check coincidence of argument and specifier types before writing the
  required output.
*/
		case 's':
		  switch (args[arg]->atyp) {
		  case 'c':
		    lprintf(fptr, io_buff, *STRPTR(args[arg]));
		    break;
		  case 'l':
		    lprintf(fptr, io_buff, (*LOGPTR(args[arg])) ? "TRUE":"FALSE");
		    break;
		  default:
		    was_error = 1;
		    break;
		  };
		  break;
		case 'e': case 'f': case 'g':
		  was_error = (args[arg]->atyp != 'f');
		  if(!was_error)
		    lprintf(fptr, io_buff, *FLTPTR(args[arg]));
		  break;
		case 'd': case 'i':
		  was_error = (args[arg]->atyp != 'i');
		  if(!was_error)
		    lprintf(fptr, io_buff, *INTPTR(args[arg]));
		  break;
		};
/*
  Report specifier-argument mismatches.
*/
		if(was_error) {
		  lprintf(stderr, "fprintf(,'%s',...): argument specifier %d does not match its argument\n",fmt, arg+1);
		  return -1;
		};
/*
  Continue parsing the format string for the next argument.
*/
		arg++;
		str_ptr=cptr+1;
		break;
	      };
	    };
/*
  Check if the end of the format line was reached before the latest
  specifier was completed.
*/
	    if(c == '\0') {
	      lprintf(stderr, "printf(,'%s',...): Final format specifier incomplete.\n",fmt);
	      return -1;
	    };
	    break;
/*
  Check for Carriage return and line feed requests.
*/
	  case '\\':
/*
  If a valid escape sequence is encounterred - first copy any unwritten preceding
  string into io_buff[] - write the pre-escape string + the escape sequence, then
  advance cptr and str_ptr past the escape sequence.
*/
	    switch (*(cptr+1)) {
	    case 'n': case 'r': case 't':
	      strncpy(io_buff,str_ptr, cptr-str_ptr);
	      io_buff[cptr-str_ptr]='\0';
	      switch (*(cptr+1)) {
	      case 'n':
		lprintf(fptr, "%s\n",io_buff);
		break;
	      case 'r':
		lprintf(fptr, "%s\r",io_buff);
		break;
	      };
	      str_ptr = (++cptr) + 1;
	    };
	    break;
	  default:
	    break;
	  };
	};
/*
  The final character of the format string has been parsed, write
  any remaining un-written output.
*/
	if(cptr-str_ptr > 0) lprintf(fptr, str_ptr);
	return 0;
}

/*.......................................................................
  Query the user. First prompt the user using the string sent as the only
  argument with '? (y/n): ' appended. Return 1 if the user answers with
  'y' or presses carriage return without any string. Otherwise if the
  user enters 'n' return 0. If anything else is typed re-prompt
  with 'Please type y (or press return) for yes or n for no: '. On
  error -1 is returned.
*/
int ask_user(char prompt[])
{
        static char *ctst;
/*
  Prompt.
*/
        lprintf(stdout, "%s? (y/n): ",prompt);
	for (;;) {
	  ctst = fgets(io_buff, MAXLINE, stdin);
	  if(ctst == NULL) {
	    clearerr(stdin);
	    lprintf(stderr, "Aborted due to read error on stdin\n");
	    return -1;
	  };
/*
  Check what the user typed.
*/
	  if(io_buff[1]=='\n') {
	    switch (io_buff[0]) {
	    case 'y': case 'Y':
	      return 1;
	    case 'n': case 'N':
	      return 0;
	    };
	  }
/*
  If nothing was typed before new-line the default is true.
*/
	  else if(io_buff[0]=='\n')
	    return 1;
/*
  Re-prompt.
*/
	  lprintf(stdout, "Please answer y or n.  %s? (y/n): ",prompt); 
	};
}

/*.......................................................................
 * Query the user for a string. The user will be presented with the given
 * prompt, and the provided default value. If the user presses return
 * without entering anything then the default string will be substituted.
 *
 * Input:
 *  prompt     char *   The string the prompt the user with.
 *  defstr     char *   The default string to be returned if the user
 *                      doesn't enter anything. If this is NULL then
 *                      the user will be re-prompted for input.
 * Output:
 *  return     char *   The answer string or NULL on error. Note that
 *                      this is malloced memory which must be free'd
 *                      by the caller when no longer required.
 */
char *prompt_user(char *prompt, char *defstr)
{
  char *retval;          /* A malloc'd copy of the answer string */
  char *answer = NULL;   /* The user's answer */
/*
 * Ask the use for input and read the user's reply.
 * This may take more than one attempt.
 */
  while(!answer) {
    char *nl; /* A pointer to the newline in io_buff[] */
/*
 * Prompt for input.
 */
    if(!defstr || *defstr == '\0')
      lprintf(stdout, "%s: ", prompt);
    else
      lprintf(stdout, "%s (%s): ", prompt, defstr);
/*
 * Read the reply.
 */
    if(!fgets(io_buff, MAXLINE, stdin)) {
      clearerr(stdin);
      lprintf(stderr, "prompt_user: Aborted due to read error on stdin.\n");
      return NULL;
    };
/*
 * Find the newline.
 */
    nl = strchr(io_buff, '\n');
/*
 * String too long?
 */
    if(!nl) {
      lprintf(stderr, "prompt_user: String too long.\n");
      continue;
    };
/*
 * Remove the newline character.
 */
    *nl = '\0';
/*
 * If the string is empty and there is a default string, substitute
 * the default.
 */
    if(io_buff[0] != '\0')
      answer = io_buff;
    else if(defstr)
      answer = defstr;
    else
      answer = NULL;
  };
/*
 * Allocate storage for a copy of the string.
 */
  retval = stralloc(strlen(answer));
  if(retval)
    strcpy(retval, answer);
  return retval;
}

/*.......................................................................
  On input skip the next field in the file connected to FILE *fptr.
  A field is delimited either by white space or by matched quotes (or
  one quote and a carriage return character. Any preceding and following
  white space is also skipped. This function is only for use with text
  files. It is not recommended for use on stdin either. If the end of
  file mark is intercepted before the start of the required field, -1 is
  returned. If the end of file mark terminates the required field then
  0 is returned. Otherwise the terminating character of the field is
  returned.
*/
int skip_field(FILE *fptr)
{
        static int c;
/*
  Find the first non-white space character.
*/
	while((c=fgetc(fptr)) != EOF && (isspace(c) || c == ',') );
	if(c==EOF)
	  return -1;
/*
  If that character is a ' then skip everything up to the next ' or carriage return.
*/
	if(c=='\"')
	  while( (c=fgetc(fptr)) != EOF && c != '\n' && c != '\"');
/*
  Otherwise read up to the next white space.
*/
	else
	  while( !isspace(c=fgetc(fptr)) && c != EOF);
/*
  Check for end of file.
*/
	if(c == EOF)
	  return 0;
	return c;
}

/*.......................................................................
  Read a 1-D array of unknown size from a user file. The return array
  and its size are required. Reading of the array will be terminated
  when carriage return or end of file characters are intercepted, when
  non-numeric characters are intercepted (with the exception of tabs,
  commas ans spaces) or when the array has been filled up to its
  declared size. The return value of the function is either the size
  of the array read or, on error, -1. This function is only for use on
  text files.
*/
int input_array(FILE *fptr, float *array, size_t num_el)
{
        int c,i,num_read;
/*
  Skip tabs, commas and spaces.
*/
	for(i=0;i<num_el;) {
	  c=fgetc(fptr);
	  switch (c) {
/*
  Skip tabs, spaces and commas.
*/
	  case '\t': case ' ': case ',':
	    break;
/*
  Check for the start of a new number.
*/
	  case '-': case '+': case '.': case '0': case '1': case '2':
	  case '3': case '4': case '5': case '6': case '7': case '8':
	  case '9':
	    ungetc(c,fptr);
/*
  Read the number.
*/
	    if((num_read=fscanf(fptr, "%f", array+i)) == EOF || num_read == 0) {
	      lprintf(stderr, "Read error while reading a 1-D array.\n");
	      return -1;
	    };
/*
  Step the element pointer to the next empty slot of the return array.
*/
	    i++;
	    break;
/*
  A non-numeric, non-separator character has been intercepted - terminate
  the read operation and return the number of elements filled.
*/
	  default:
	    if(c != '\n') ungetc(c,fptr);
	    return i;
	    break;
	  };
	};
/*
  The maximum number have been read.
*/
	return i;
}

/*.......................................................................
  Given a file pointer, a format string, a specification of the number
  of arguments and an array of those arguments (as pointers to descriptors),
  handle formatted input.
*/
int fmt_read(FILE *fptr, char fmt[], int nargs, Descriptor *args[])
{
        static int c;
	static int arg, num, buf_pos, in_space, was_error;
	static char *cptr, *tmp_ptr;
/*
  First argument - args[0].
*/
	arg=0;
/*
  Parse the format string.
*/
	for(cptr=fmt; *cptr != '\0'; cptr++) {
/*
  Skip spaces tabs and commas in the format string.
*/
	  in_space = 1;
	  while(in_space) {
	    switch (*cptr) {
	    case ',': case ' ': case '\t':
	      cptr++;
	      break;
	    case '\0':
	      return 0;
	      break;
	    default:
	      in_space=0;
	    };
	  };
/*
  The next characters should consist of an optional integer followed by format
  type-specifier letter.
*/
	  if(isdigit((int)*cptr))
	    num = strtol(cptr, &cptr, 10);
	  else
	    num=0;
/*
  Obey the specifier type.
*/
	  switch(*cptr) {
/*
  Read a number from the file.
*/
	  case 'f': case 'i':
	    if(arg+1 > nargs) {
	      lprintf(stderr, "fprintf(,'%s',...): More specifiers than arguments?\n",fmt);
	      return -1;
	    };
	    switch (*cptr) {
	    case 'f':
	      was_error = args[arg]->atyp != 'f';
	      break;
	    case 'i':
	      was_error = args[arg]->atyp != 'i';
	      break;
	    };
	    if(was_error) {
	      lprintf(stderr, "fprintf(,'%s',...): argument specifier %d does not match its argument\n",fmt, arg+1);
	      return -1;
	    };
/*
  Skip white space preceding the number.
*/
	    while( (c=fgetc(fptr)) != EOF && (isspace(c) || c == ','));
	    if(c == EOF) {
	      lprintf(stderr, "Reached end of file while searching for start of argument %d of:\n\t%s\n", arg+1, fmt);
	      return -1;
	    };
	    ungetc(c,fptr);
/*
  Read the number.
*/
	    switch (*cptr) {
	    case 'f':
/*
  If the optional width specifier was provided, use it.
*/
	      if(num != 0) {
		sprintf(io_buff, "%%%df",num);
		c=fscanf(fptr, io_buff, FLTPTR(args[arg]));
	      }
	      else c=fscanf(fptr, "%f",FLTPTR(args[arg]));
	      break;
	    case 'i':
	      if(num != 0) {
		sprintf(io_buff, "%%%di",num);
		c=fscanf(fptr, io_buff, INTPTR(args[arg]));
	      }
	      else c=fscanf(fptr, "%i",INTPTR(args[arg]));
	      break;
	    };
/*
  Check that the read was successful.
*/
	    if(c == 0) {
	      lprintf(stderr, "Error in reading argument %d of '%s'\n", arg+1, fmt);
	      return -1;
	    };
/*
  The current argument has been assigned - prepare for the next one - if any.
*/
	    arg++;
/*
  Skip any following white-space up to the next record or up to any new-line
  character. This is here essentially to make sure that if the last record on
  a given line is simply followed by spaces, then the line will have been completely
  read - this is important for reads from stdin where the carriage return would
  be interpretted by the next read or the compiler.
*/
	    for(;;) {
	      c=fgetc(fptr);
	      if(c=='\n' || c==EOF)
		break;
	      if(!isspace(c) && c != ',') {
		ungetc(c,fptr);
		break;
	      };
	    };
	    break;
/*
  Read a string from the file.
*/
	  case 's':
	    if(arg+1 > nargs) {
	      lprintf(stderr, "fprintf(,'%s',...): More specifiers than arguments?\n",fmt);
	      return -1;
	    };
	    if(args[arg]->atyp != 'c') {
	      lprintf(stderr, "fprintf(,'%s',...): argument specifier %d does not match its argument\n",fmt, arg+1);
	      return -1;
	    };
/*
  Skip white space in the file.
*/
	    while( (c=fgetc(fptr)) != EOF && (isspace(c) || c == ','));
	    if(c == EOF) {
	      lprintf(stderr, "Reached end of file while searching for start of argument %d of:\n\t%s\n", arg+1, fmt);
	      return -1;
	    };
/*
  The string will be read into io_buff first and then copied to the return string.
*/
	    io_buff[0] = c;
	    buf_pos = 1;
/*
  The optional number in the format specifier is the minimum number of characters to be
  read. If this is specified simply read the requested number but stop at line breaks and
  EOF.
*/
	    if(num != 0) {
	      if(num > MAXLINE-1) num = MAXLINE-1;
	      while(buf_pos < num && (c=fgetc(fptr)) != EOF && c != '\n')
		io_buff[buf_pos++] = c;
	    }
/*
  If the initial character is a ' then read everything up to the next ' or carriage return.
*/
	    else if(c == '\"') {
	      buf_pos = 0;
	      while(buf_pos <MAXLINE-1 && (c=fgetc(fptr)) != EOF && c != '\n' && c != '\"')
		io_buff[buf_pos++] = c;
	      c=fgetc(fptr);
	    }
/*
  Otherwise read up to the next white space.
*/
	    else {
	      while(buf_pos < MAXLINE-1 && !isspace(c=fgetc(fptr)) && c != ',' && c != EOF)
		io_buff[buf_pos++] = c;
	    };
/*
  Terminate the string.
*/
	    io_buff[buf_pos] = '\0';
/*
  Free the string pointed to by the current argument.
*/
	    char_free(STRPTR(args[arg]));
/*
  Allocate memory for a copy of the new string.
*/
	    if( (tmp_ptr=stralloc(buf_pos)) == NULL)
	      return -1;
/*
  Now copy the string into it.
*/
	    strcpy(tmp_ptr, io_buff);
	    STRPTR(args[arg])[0] = tmp_ptr;
	    arg++;
/*
  Skip any following white-space up to the next record or up to any new-line
  character. This is here essentially to make sure that if the last record on
  a given line is simply followed by spaces, then the line will have been completely
  read - this is important for reads from stdin where the carriage return would
  be interpretted by the next read or by the compiler.
*/
	    for(;;) {
	      if(c=='\n' || c==EOF)
		break;
	      if(!isspace(c) && c != ',') {
		ungetc(c,fptr);
		break;
	      };
	      c=fgetc(fptr);
	    };
	    break;
/*
  Skip num fields in the file.
*/
	  case 'F':
	    if(num == 0) num=1;
	    while(skip_field(fptr) != -1 && --num > 0);
	    break;
/*
  Skip num lines in the file.
*/
	  case 'L':
	    for(;;)
	      if( (c=fgetc(fptr)) == EOF || (c=='\n' && --num <= 0) )
		break;
	    break;
/*
  Skip num characters in the file.
*/
	  case 'C':
	    do {
	      if(fgetc(fptr) == EOF)
		break;
	    } while(--num > 0);
	    break;
/*
  The following characters delimited by a matching } should be searched for
  num times. Copy the string into io_buff first.
*/
	  case '{':
	    buf_pos = 0;
	    cptr++;
	    while( (c = *cptr) != '}' && c != '\0') {
	      cptr++;
	      io_buff[buf_pos++] = c;
	    };
	    if(c == '\0') {
	      lprintf(stderr, "Unmatched '{' in read-format: %s\n", fmt);
	      return -1;
	    };
/*
  Perform the search num times.
*/
	    do {
	      if(buf_pos > 0 && file_search(fptr, io_buff, buf_pos, 0) <= 0)
		return -1;
	    } while(--num > 0);
	    break;
	  case '\0':
	    return 0;
	    break;
/*
  Unsupported specifier type - report error.
*/
	  default:
	    lprintf(stderr, "Un-recognised read-format specifier %c\n",*cptr);
	    return -1;
	    break;
	  };
	};
	return 0;
}
