#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <errno.h>

#if USE_TECLA == 1
#include "libtecla.h"
static GetLine *gl;
static CPL_MATCH_FN(tecla_match_fn);
#endif

#include "sphere.h"
#include "lex.h"
#include "table.h"
#include "ops.h"
#include "utils.h"

/* Global handle on command line details */

Comline comline = {'\0', NULL, NULL, 0};

enum {
  MAX_NAME=256,     /* The maximum size of a symbol name */
  MAX_LEV=6,        /* The maximum number of command levels */
  MAX_ARGS=40,      /* The maximum number of pre-processor arguments */
  MAX_MACRO=50,     /* Max number of macros */
  MAXFNAME = 132    /* Max char's per file name */
};

/*
 * This struct holds information about the current command level.
 */
static struct {
  FILE *unit;             /* Command file pointer or NULL if not a file level */
  char inbuff[MAX_LINE];  /* Line input buffer */
  char ppbuff[MAX_LINE];  /* Command sub-line after pre-processing */
  char script[MAXFNAME+1];/* The file-name of the script, or '\0' if stdin */
  int line_no;            /* Number of lines read at this command level */
  char *in_ptr;           /* Pointer into ppbuff[] */
  char was_eol;           /* True if inbuff[] has been consumed */
  Args args[MAX_ARGS];    /* Array of arguments to this command level */
  int nargs;              /* Number of arguments in args[] */
} com[MAX_LEV];

typedef struct {
  int c;                 /* The latest character read */
  int was_escape;        /* True if the character was an escape character */
} LexChar;

static void parse_compound_char(LexChar *lc);

static int comlev=0;     /* Current command level */

static Table *macro_table[MAX_MACRO];
static int num_macro=0;
static char namebuf[MAX_NAME];

static int get_literal(void);
static int com_close(void);
static size_t get_name(char **next);
static int pp_args(const char *argstr, Args *args, int *num_arg, int maxargs);
static char nxtchar(char **last, char **next);
static int add_macro(char inbuff[]);
static char *lex_pre_process(char *inbuff, char *outbuff, int nmax,
			     Args *args, int narg);

/*.......................................................................
 * Start a new command file.
 *
 * Input:
 *  filestr   char *   The string containing the file name and any
 *                     arguments for the file.
 * Output:
 *  return     int     0 - OK.
 *                    -1 - Error.
 */
int com_open(const char *filestr)
{
  char fname[MAXFNAME+1];   /* The name of the script file */
  FILE *fp;                 /* File descriptor of new command file */
  int i;
/*
 * Skip white-space preceding the file name.
 */
  while(isspace((int)*filestr))
    filestr++;
/*
 * Copy the file name into fname[].
 */
  for(i=0; i<MAXFNAME && *filestr && !isspace((int)*filestr); i++)
    fname[i] = *filestr++;
  fname[i] = '\0';
/*
 * File name too long?
 */
  if(i>=MAXFNAME) {
    lprintf(stderr, "com_open: Command file name too long.\n");
    return -1;
  };
/*
 * Attempt to open the command file.
 */
  if(i==0) {
    fp = stdin;
    lprintf(stdout,
	    "Starting new interactive shell. Use the EOF character to exit.\n");
  } else {
/*
 * Log the attempt to run the script.
 */
    lprintf(stdin, "![@%s%s]\n", fname, filestr);
/*
 * Attempt to open the script file.
 */
    if((fp=fopen(fname,"r")) == NULL) {
      lprintf(stderr,"Unable to open command file: %s\n",fname);
      return -1;
    };
  };
/*
 * Create a new command level for the new file.
 */
  return push_command(fp, NULL, fname, filestr);
}

/*.......................................................................
 * Initialize the file buffers.
 *
 * Input:
 *  bootenv    char *   If not NULL then this specifies the name of the
 *                      environment variable that will be queried for the
 *                      name of a startup script.
 */
int com_init(const char *bootenv)
{
  char *boot_ptr;  /* Pointer to value of the bootenv environment variable */
/*
 * If needed, create the terminal input resource object.
 */
#if USE_TECLA == 1
  gl = new_GetLine(MAX_LINE, MAX_LINE * 10);
  if(!gl)
    return 1;
  if(gl_customize_completion(gl, NULL, tecla_match_fn))
    return 1;
#endif
/*
 * Initialize the command line file pointer to stdin and signal
 * that the stdin buffer is currently empty by signalling
 * end of line.
 */
  com[0].unit = stdin;
  com[0].was_eol = 1;
  com[0].script[0] = '\0';
/*
 * If the enviroment variable BOOT is set, try to open the file
 * named within it as an initialisation file to be obeyed.
 */
  if(bootenv && (boot_ptr=getenv(bootenv)) != NULL) {
    if(com_open(boot_ptr)==0)
      lprintf(stdout, "Obeying initialization file: %s\n",boot_ptr);
    else
      lprintf(stderr,"The file name was taken from environment variable %s\n", bootenv);
  };
  return 0;
}

/*.......................................................................
 * This routine is used for the lexical analysis of expressions. On each
 * succesive call it attempts to identify the object at the next position
 * in the input string as a numeric or string constant, an
 * operator symbol, a function name, or a user-variable name. In addition,
 * if instructed to via the passed value of optyp being 'l', then it
 * will return a literal string without cross checking
 * against the symbol tables. Similarly, if optyp is 'n' then it will
 * return the next identifier without cross checking against the symbol
 * table.
 *  If a succesfull parse of the new word or operator sequence is made
 * then the return value will either be a pointer to the matched function,
 * variable or operator table entry, or a pointer to a new table entry,
 * with a NULL name entry but an active Descriptor type entry describing
 * the constant quantity parsed. The latter case handles string, numeric
 * and literal constants.
 *  In order that different strings can be parsed alternately, the
 * string to be parsed is sent to the routine packaged within a
 * Comline structure, which in addition to a pointer to the
 * beginning of the string, also contains a pointer to the an arbitrary
 * location in the string. On entry to this routine, this pointer should
 * point to the next character to be parsed. On exit it will be left at
 * the character following the last character identified.
 *  If the analyser fails to identify the next character sequence, or
 * the first character is a null terminator string '\0', a NULL pointer
 * will be returned. NB. The string must be terminated by a '\0'.
 */
Table *lex_expr(char optyp)
{
  static int is_start;
  static size_t slen;
  static char c;
  static short int shint;
  static float fnum;
  static double num;
  static Table ret_table, *ttst;
/*
 * Set a flag to indicate if the next token comes from the start of
 * a new line.
 */
  is_start = (comline.next == com[comlev].ppbuff);
/*
 * Find the next non-white-space character.
 */
  while( isspace( (int) (c = *(comline.next++)) ) != 0);
/*
 * Make the string pointer point at the new character.
 */
  comline.next--;
/*
 * If optyp = 'l' then a literal string is required.
 */
  switch(optyp) {
  case 'l':
    if( get_literal() == -1)
      return NULL;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
    return store_const('c', (void *) namebuf);
/*
 * If optyp has been sent with the value 'n', then the caller explicitly
 * requires an identifier which probably does not exist in any of the symbol
 * tables - so don't check the symbol against the symbol tables.
 */
  case 'n':
/*
 * Identifiers must start with a letter - check this.
 */
    if(isalpha((int) c) == 0) {
      lex_err(++comline.next);
      lprintf(stderr,"Illegal non-letter first character of an identifier.\n");
      return NULL;
    };
/*
 * Get the identifier name and convert it to lower case.
 */
    if( (slen = get_name(&comline.next)) == -1)
      return NULL;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
/*
 * Return the identifier name in the name field of the returned Table
 * structure.
 */
    if( (ret_table.name=stralloc(slen)) == NULL)
      return NULL;
    strcpy(ret_table.name, namebuf);
    return &ret_table;
    break;
  };
/*
 * If the new character is alphabetic then copy both it and subsequent
 * alphanumeric characters into namebuf.
 */
  if(isalpha((int) c) != 0) {
    if( (slen = get_name(&comline.next)) == -1)
      return NULL;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
/*
 * Search for an unambiguous or exact match to one of the names in the
 * function and operator symbol tables.
 */
    return match_name(namebuf);
  }
/*
 * See if the next few characters constitute a numeric constant.
 */
  else if(isdigit((int) c) != 0) {
    num = strtod(comline.next,&comline.next);
    if(fabs(num) >= FLT_MAX) {
      lex_err(comline.last);
      lprintf(stderr,"Number too big to be read\n");
      return NULL;
    };
    fnum = (float) num;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
    return store_const('f', (void *) &fnum);
  }
/*
 * See if the following characters represent a string constant.
 */
  else if(c == '\"') {
    if( get_literal() == -1)
      return NULL;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
    return store_const('c', (void *) namebuf);
  }
/*
 * Check for #num, the pseudo variable that, during array expression
 * evaluations holds the the current index along the axis num.
 */
  else if(c == '#') {
    if(!isdigit((int) (c = *(++comline.next))) || (shint=c-'0') >= 3) {
      lex_err(comline.next);
      lprintf(stderr, "Format: #N where N is an integer between 0 and the dimension required.\n");
      return NULL;
    };
/*
 * Allocate the necessary table entry for the # instruction. The
 * type field will hold a short int denoting the dimension required.
 */
    if((ttst=table_alloc(HASH, NULL)) == NULL)
      return NULL;
    TABICODE(ttst) = shint;
    comline.next++;
    comline.nxtc = nxtchar(&comline.last,&comline.next);
    return ttst;
  }
/*
 * Finally try and match the next characters or couple of characters
 * with an operator symbol.
 */
  else {
    ttst = find_ops(&comline.next, namebuf);
    if(ttst==NULL) {
/*
 * Report error in operator name.
 */
      lex_err(comline.last);
      lprintf(stderr,"Unable to identify \"%s\" with an operator.\n",namebuf);
      return NULL;
    } else {
      comline.nxtc = nxtchar(&comline.last,&comline.next);
      return ttst;
    };
  };
  return NULL;
}


/*.......................................................................
 * Split a string into arguments for the pre-processor. An argument is
 * anything delimited by commas or end of string. Commas inside string
 * constants and between matched parenthesis that don't lie in string
 * constants, are ignored. On error returns -1. Otherwise 0.
 */
static int pp_args(const char *argstr, Args *args, int *num_arg, int maxargs)
{
  static const char *start_ptr;
  static char num_open, in_string;
/*
 * *num_arg holds a record of the current argument being scanned. in_string
 * is a flag which is 1 when the current character is within a string constant.
 * num_open is the current number of un-closed parenthesis - commas are ignored
 * unless num_open is zero.
 */
  *num_arg = in_string = num_open = 0;
/*
 * Skip white space to the start of the first argument.
 */
  while(isspace((int) *argstr) || *argstr==',')
    argstr++;
/*
 * Record the address of the start of the first argument.
 */
  start_ptr=argstr;
/*
 * Step through the string, searching for arguments, until the string
 * null terminator is encounterred.
 */
  for(; *argstr != '\0'; argstr++) {
/*
 * Only perform checks for brackets outside of string constants.
 */
    switch (in_string) {
    case 0:
      switch (*argstr) {
/*
 * Increment or decrement the recorded number of open brackets
 * when an open or close bracket is encounterred.
 */
      case '(':
	num_open++;
	break;
      case ')':
	if(--num_open < 0) {
	  lex_err(argstr);
	  lprintf(stderr, "Unmatched bracket.\n");
	  return -1;
	};
	break;
/*
 * When a comma is encounterred outside of any parentheses, then
 * this denotes the end of the current argument - record its
 * start pointer and length and increment arg for the next
 * argument.
 */
      case ',':
	if(num_open == 0) {
	  args[*num_arg].arg_ptr = (char *) start_ptr;
	  args[*num_arg].arg_len = argstr - start_ptr;
	  (*num_arg)++;
	  start_ptr = ++argstr;
/*
 * Trap too many arguments.
 */
	  if(*num_arg >= maxargs) {
	    lex_err(start_ptr);
	    lprintf(stderr, "Too many arguments\n");
	    return -1;
	  };
/*
 * Skip white space preceding the start of the next argument.
 */
	  while(isspace((int) *argstr))
	    argstr++;
	  argstr--;
	};
	break;
/*
 * Flag the start of a string constant.
 */
      case '\"':
	in_string = 1;
	break;
      };
/*
 * While in a string constant only check for the terminating
 * quote.
 */
    default:
      if(*argstr == '\"' && *(argstr-1) != '\\')
	in_string = 0;
    };
  };
/*
 * Record the final argument.
 */
  args[*num_arg].arg_ptr = (char *) start_ptr;
  args[*num_arg].arg_len = argstr - start_ptr;
/*
 * Check that all brackets have been closed.
 */
  if(num_open != 0) {
    lex_err(argstr);
    lprintf(stderr, "Un-matched parentheses\n");
    return -1;
  };
  if(args[*num_arg].arg_len != 0) (*num_arg)++;
  return no_error;
}

/*.......................................................................
 * Copy the next word (delimited by non-alphanumeric characters) into
 * namebuf, from the command line (Comline) line. All the characters
 * in the copy are converted to lower case.  If the word takes
 * more than MAX_NAME characters then the routine returns -1. Otherwise
 * it returns the length of the word (not including the '\0'.
 */
static size_t get_name(char **next)
{
  static size_t i;
  static int c;
/*
 * Copy at most MAX_NAME alphanumeric (including '_') characters.
 */
  for(i=0; i < MAX_NAME; i++) {
    c = (int) **next;
    if(isalnum(c) != 0 || c == '_'){
      namebuf[i] = (char) tolower(c);
      (*next)++;
    }
    else {
      break;
    };
  };
/*
 * Terminate the string-copy.
 */
  namebuf[i] = '\0';
/*
 * Check to ensure that the name buffer hasn't overflowed.
 */
  if(i == MAX_NAME) {
    lex_err(*next);
    lprintf(stdout, "Name too long: \"%s...\"\n", namebuf);
    return -1;
  };
  return i;
}
  
/*.......................................................................
 * Copy the next literal string into namebuf, from the command line
 * (Comline) line. The string is either delimited by ',', '\0', or if the
 * first character is a "'" character by the next occurrence of a
 * non-escaped (\') "'" or a null string terminator. 
 * If the string takes more than MAX_NAME characters then the routine
 * returns -1. Otherwise it returns the length of the string (not including
 * the '\0'). Trailing spaces are not included in the returned string.
 */
static int get_literal(void)
{
        LexChar lc;
        size_t i;
	int nb=0;
/*
 * Copy at most MAX_NAME characters.
 * If the first character is a "'" copy all characters up to the end quote or
 * end of line.
 */
        if(*comline.next == '\"') {
	  comline.next++;
          for(i=0;i < MAX_NAME-1; i++) {
	    parse_compound_char(&lc);
/*
 * End of string?
 */
	    if(lc.c=='\0' || (lc.c=='\"' && !lc.was_escape)) {
	      namebuf[i]='\0';
	      return i;
	    };
/*
 * Record the latest character.
 */
	    namebuf[i] = lc.c;
	  };
        }
/*
 * For unquoted strings copy all characters up to the next
 * comma outside of matching parentheses, unmatched close
 * paren, or end-of-line.
 */
        else {
	  int end_of_string = 0;
          for(i=0; i < MAX_NAME-1; i++) {
	    parse_compound_char(&lc);
	    if(lc.c == '\0') {
	      end_of_string = 1;
	      break;
	    };
	    if(!lc.was_escape) {
	      if((lc.c == ',' && nb<=0) ||
		 ((lc.c==')' || lc.c==']' || lc.c=='}') && --nb < 0)) {
		comline.next--;
		end_of_string = 1;
		break;
	      };
	      if((lc.c=='(' || lc.c=='[' || lc.c=='{'))
		nb++;
	    };
	    namebuf[i] = lc.c;
          };
/*
 * Did we reach the end of the string before overflowing the output
 * buffer?
 */
	  if(end_of_string) {
/*
 * Terminate the string copy before the first trailing white-space.
 */
	    while(i>0 && isspace((int)namebuf[i-1]))
	      i--;
	    namebuf[i] = '\0';
	    return i;
	  };
        };
/*
 * Check to ensure that the name buffer hasn't overflowed.
 */
        lex_err(comline.next);
        lprintf(stdout, "String too long: \"%s...\"\n",namebuf);
        return -1;
}

/*.......................................................................
 * This is a private function of get_literal(). It returns the next
 * character from the input stream, with due regard for escapes.
 *
 * Input/Output:
 *  lc        LexChar *  The container that the results are returned in.
 */
static void parse_compound_char(LexChar *lc)
{
/*
 * An escape character?
 */
  switch(*comline.next) {
  case '\\':
    switch(*(++comline.next)) {
    case 'a':
      lc->c = '\a';
      break;
    case 'b':
      lc->c = '\b';
      break;
    case '\\':
      lc->c = '\\';
      break;
    case 'n':
      lc->c = '\n';
      break;
    case 't':
      lc->c = '\t';
      break;
    case 'r':
      lc->c = '\r';
      break;
    case '\"':
      lc->c = '\"';
      break;
    default:  /* Un-recognized escape - simply discard the '\' character */
      lc->c = *comline.next;
      break;
    };
    lc->was_escape = 1;
    comline.next++;
    break;
  case '\0':       /* Be careful not to advance comline.next past '\0' */
    lc->c = '\0';
    lc->was_escape = 0;
    break;
  default:
    lc->c = *(comline.next++);
    lc->was_escape = 0;
    break;
  };
}

/*.......................................................................
 * This routine should be called whenever a syntax error is
 * encounterred in the command during lexical analysis. It echoes the
 * line to stderr and then places a pointer directly underneath the
 * end of the word that caused the error condition. (Comline line) is
 * the command line in which the error was detected. (char *err_ptr)
 * is the pointer to one character past the offending character. If
 * the char *fmt argument is non NULL it will be used as a vfprintf
 * format string for the variable lenght argument list that follows
 * it.
 */
void lex_err(const char *err_ptr)
{
  char *cptr;    /* A pointer into com[comlev].ppbuf */
/*
 * Echo the command line that holds the erroneous name.
 */
  lprintf(stderr, "Syntax error in line:\n%s\n", com[comlev].ppbuff);
/*
 * Print a carat under the problematic character, noting that err_ptr
 * points one character beyond the character that upset the parser.
 */
  for(cptr=com[comlev].ppbuff; cptr < err_ptr-1; cptr++)
    lputc(' ',stderr);
  lputc('^',stderr);
  lputc('\n',stderr);
  return;
}


/*.......................................................................
 * Abort user input. Close all command files, returning control to stdin,
 * and find the next newline, ignoring anything left on the current line.
 */
void flush_input(void)
{
/*
 * Close all command files until control has been returned to the lowest
 * shell level.
 */
  while(comlev > 0) com_close();
/*
 * Mark the stdin buffer as empty.
 */
  com[comlev].was_eol = 1;
  com[comlev].script[0] = '\0';
  comline.nest_block=0;
  lprintf(stdout, "\n");
  return;
}

/*.......................................................................
 * This function is used by lex_expr to find and return the next
 * non-white-space character in the current input line. It leaves the
 * line pointer pointing at the new character. If called when the line
 * pointer points at a '\0' terminator, it returns '\0'.
 */
static char nxtchar(char **last, char **next)
{
  *last = *next;
  while( isspace( (int) **next ) != 0)
    (*next)++;
  return **next;
}
/*.......................................................................
 * Close a command file. If an attempt is made to close stdin (comlev = 0)
 * then signal an error by returning -1.
 */
static int com_close(void)
{
  if(comlev == 0) {
    lprintf(stderr,"<Aborted by %s>\n", errno==EINTR ? "Interrupt" : "EOF");
/*
 * Mark the stdin buffer as empty and clear the stdin error flag.
 */
    com[0].was_eol = 1;
    clearerr(stdin);
    return -1;
  };
/*
 * Close the current command file or macro expansion buffer.
 */
  if(com[comlev].unit == stdin) {
    lprintf(stdout, "\nTerminated interactive shell.\n");
    clearerr(stdin);
  } else if(com[comlev].unit != NULL) {
    lprintf(stdin, "![Exited script file: %s]\n", com[comlev].script);
    fclose(com[comlev].unit);
  };
  comlev--;
  return no_error;
}

/*.......................................................................
 * This routine handles command input. Each call extracts and
 * pre-processes the next command line. Since it is legal to have
 * multiple commands per line, such a line may actually be a portion of
 * the line that the user enterred. Such sub-lines are delimited by the
 * ';' command separator.  When a full input line has been finished
 * with, a new line is read via fgets(). In order that command files may
 * be nested, separate buffers and the information to maintain them are
 * kept. The current command file nest level is held in comlev, its
 * buffer is com[comlev].inbuff, the current position in that buffer is
 * maintained in char com[comlev].in_ptr, and if the end of the buffer
 * was reached on the preceding call then com[comlev].was_eol is 1. When
 * the end of the current file is reached, the file is closed and comlev
 * is decremented to the preceding command level.
 */
int newline(void)
{
  int bot,top;             /* Indexes of closest matches to a macro name */
  char *inptr;             /* Pointer into the input line */
  int num_eof=0;           /* Count of consecutive end-of-file's on stdin */
  enum {MAX_PROMPT=10};    /* The maximum length of a prompt, including '\0' */
  char prompt[MAX_PROMPT]; /* The prompt string to display before each line */
/*
 * If an end of line was detected on the last call then get a new line.
 */
  if(com[comlev].was_eol) {
/*
 * Read new lines until a non-empty line is located.
 */
    do {
      com[comlev].was_eol=0;
/*
 * Increment the line count for the current input stream.
 */
      com[comlev].line_no++;
/*
 * Reset the character pointer to the start of the buffer for the
 * current input stream.
 */
      inptr = com[comlev].in_ptr = com[comlev].inbuff;
/*
 * If the current input stream is stdin, prompt for the new line.
 * The prompt starts with an indication of the current command
 * block nest level (0 outside of any block).
 */
      if(com[comlev].unit == stdin)
	sprintf(prompt, "%d>", comline.nest_block);
      else
	prompt[0] = '\0';
/*
 * Get a new line.
 */
      if(lexgets(com[comlev].inbuff, MAX_LINE, com[comlev].unit, prompt)!=0) {
/*
 * End of file - close command file. If EOF detected on stdin increment
 * a count of how many times this happens. If it happens more than 10
 * times, exit the program - this probably indicates that stdin is not open.
 */
	if(com_close() == -1 && errno!=EINTR && ++num_eof>3) {
	  lprintf(stderr, "%d consecutive EOF's on stdin - 10 will exit program\n", num_eof);
	  if(num_eof > 9)
	    closedown(1, DO_QUIT);
	};
	inptr = com[comlev].in_ptr;
      }
      else {
/*
 * Reset the count of consecutive end-of-files received.
 */
	num_eof=0;
/*
 * Skip white space at the start of the string.
 */
	while(isspace((int) *inptr))
	  inptr++;
/*
 * Check the first character for line based operators and for the null
 * terminator. Line based operators are those operators that require
 * a line of text that has not been pre-processed or split into
 * sub-lines.
 */
	switch (*inptr) {
/*
 * When either a blank line or a comment line (signalled by an '!') is
 * detected, signal the requirement for a new line by setting the
 * end of line flag.
 */
	case '\0': case '!':
	  com[comlev].was_eol=1;
	  break;
/*
 * Check for macro declarations - these are signalled by a # symbol
 * at the start of the line, followed by the macro name. The rest of
 * the input line is the actual macro.
 */
	case '#':
	  com[comlev].was_eol=1;
	  if(add_macro(inptr) == -1)
	    return -1;
	  break;
/*
 * If the first character on the line is a '$' send the rest of the line
 * to the operating system to be executed.
 */
	case '$':
	  com[comlev].was_eol = 1;
	  system(inptr+1);
	  break;
/*
 * A normal line.
 */
	default:
	  com[comlev].in_ptr = inptr;
	  break;
	};
      };
    } while(com[comlev].was_eol);
  }
/*
 * If no new line had to be read, make inptr point at the next character
 * to be proccessed from the current input buffer.
 */
  else {
    inptr = com[comlev].in_ptr;
  };
/*
 * Pre-process the latest sub-line (delimited by ; or '\0') into
 * com[comlev].ppbuff.
 */
  inptr = lex_pre_process(inptr, com[comlev].ppbuff, MAX_LINE, com[comlev].args,
			  com[comlev].nargs);
  if(inptr==NULL)
    return -1;
/*
 * Find the start of the next sub-line - if there is one.
 */
  while(isspace((int) *inptr))
    inptr++;
  com[comlev].in_ptr = inptr;
/*
 * If there isn't one then set the end of line flag.
 */
  com[comlev].was_eol = (*inptr == '\0');
/*
 * Initialize all pointers into the pre-processed-line buffer.
 */
  comline.last = comline.next = com[comlev].ppbuff;
  comline.nxtc = nxtchar(&comline.last, &comline.next);
/*
 * If the first word of the pre-processed line is the name of a macro,
 * have it translated.
 */
  inptr = comline.next;
  if(isalpha((int)comline.nxtc) && get_name(&inptr) > 0 &&
     find_symbol(namebuf, macro_table, num_macro, &bot, &top)=='e') {
    if(push_command(NULL, TABSTR(macro_table[bot]), NULL, inptr)!=0)
      return -1;
    return newline();
  };
/*
 * Command file indirection?
 */
  if(comline.nxtc == '@') {
    if(com_open(++comline.next)!=0)
      return -1;
    return newline();
  };
/*
 * Log the new input sub-line - prefixed by 2 spaces per nest level.
 */
  lprintf(stdin, "%*s%s\n", 2*comline.nest_block, "", comline.next);
  return no_error;
}

/*.......................................................................
 * Pre-process a given sub-line into a given output buffer, using given
 * arguments. Return a pointer to the next un-processed character in the
 * command line, or NULL on error.
 *
 * Input:
 *  inbuff  char *  The pointer to the start of the sub-line.
 *                  Sub-lines are terminated either by ';' or '\0'.
 *  outbuff char *  Pointer to start of output buffer.
 *  nmax     int    Dimension of outbuff[] array.
 *  args    Args *  Array of 'narg' arguments.
 *  narg     int    The number of arguments in args[].
 * Output:
 *  return  char *  Pointer to sub-line terminator in inbuff[], or NULL
 *                  on error.
 */
static char *lex_pre_process(char *inbuff, char *outbuff, int nmax,
			     Args *args, int narg)
{
  int arg;          /* The index of the argument being inserted */
  int first_arg;    /* First argument to be substituted */
  int last_arg;     /* Last argument to be substituted */
  int do_comma;     /* True if a , is to be prepended to each argument */
  int do_fwd;       /* Increment the argument index unless 0 */
  char *outptr;     /* Pointer into outbuff[] */
  char *inptr;      /* Pointer into inbuff[] */
  int in_quotes=0;  /* True between speech marks */
  int was_eol=0;    /* True when the sub-line has been fully processed */
/*
 * Pre-process the string - stopping when a sub-line terminator is seen.
 */
  for(outptr=outbuff, inptr=inbuff; !was_eol; ) {
/*
 * Check that there is room in the output string for the next character.
 */
    if(outptr-outbuff >= nmax-2) {
      lprintf(stderr,
	"Buffer overflow while pre-proccessing:\n'%.76s...'\n", outbuff);
      return NULL;
    };
/*
 * Check for special characters.
 */
    switch(*inptr) {
/*
 * Paired speech marks delimit string constants. In order that the
 * ';' command separator can appear in a string it is neccesary to
 * detect when characters are in such a string. The in_quotes
 * flag is used to signal this and is toggled on when a start
 * " is found and off when the matching " is encounterred. In
 * order that " can appear as part of a string let it be escaped
 * as in C by \".
 */
    case '\"':
      in_quotes = !in_quotes || *(inptr-1)=='\\';
      *outptr++ = *inptr++;
      break;
/*
 * End of input line?
 */
    case '\0':
      was_eol = 1;
      break;
/*
 * End of command sub-line? Ignore the semi-colon if it
 * appears within a string constant.
 */
    case ';':
      was_eol = !in_quotes;
      inptr++;
      break;
/*
 * Pre-processor insert-arguments directive?
 */
    case '%':
/*
 * In order to separate pre-processor directives from user required
 * percent signs, a double % is required for directives within
 * string constants.
 */
      if(in_quotes) {
	if(*(inptr+1) != '%') {
	  *outptr++ = *inptr++;
	  break;
	} else {
	  inptr += 2;
	};
      } else {
	inptr++;
      };
/*
 * The insert-argument directive has syntax %<,><first_arg><.last_arg>,
 * where <> denotes optional. eg. %3 specifies argument 3; %1.4 specifies
 * arguments 1,2,3,4 and % specifies argument 1. The optional comma
 * specifies that each argument should be prefixed with a comma.
 */
      if(*inptr == ',') {
	do_comma = 1;
	inptr++;
      } else {
	do_comma=0;
      };
/*
 * first_arg specifier?
 */
      if(isdigit((int) *inptr)) {
	first_arg = (int) strtol(inptr, &inptr, 10) - 1;
/*
 * last_arg specifier?
 */
	if(*inptr == '.') {
	  if(isdigit( (int) *(inptr+1))) {
	    last_arg = (int) strtol(inptr+1, &inptr,10) - 1;
	  } else if(inptr[1] == '*') {
	    last_arg = narg - 1;
	    inptr += 2;
	  } else {
	    last_arg = first_arg;
	  };
	} else {
	  last_arg = first_arg;
	};
      }
/*
 * % followed by 'n' means substitute the number of arguments.
 */
      else if(*inptr == 'n') {
	first_arg = last_arg = -1;
	inptr++;
/*
 * % or %% with no following number. Default to argument 1.
 */
      } else {
	last_arg=first_arg=1;
      };
/*
 * Substitute arguments?
 */
      if(first_arg >= 0) {
/*
 * Determine whether to substitute arguments in forward or reverse
 * order ie. %1.3 vs %3.1
 */
	do_fwd = last_arg >= first_arg;
/*
 * Locate the index of the first argument to be used.
 */
	arg = (!do_fwd && first_arg>narg-1) ? narg-1 : first_arg;
/*
 * Copy the requested argument(s) into the output string.
 */
	for( ;; arg += (do_fwd ? 1 : -1)) {
/*
 * The check for completion of the argument substitutions depends
 * upon whether the arguments are required in forward or reverse order.
 */
	  if(do_fwd) {
	    if(arg>last_arg || arg>narg-1) break;
	  }
	  else {
	    if(arg<last_arg || arg<0) break;
	  };
/*
 * Check if the addition of the new argument would overflow the output string.
 */
	  if((outptr-outbuff) + args[arg].arg_len + (do_comma?1:0) >= nmax-1) {
	    lprintf(stderr,
		    "Buffer overflow while pre-proccessing:\n'%.76s...'\n",
		    outbuff);
	    return NULL;
	  };
/*
 * Precede the argument with a comma if so requested.
 */
	  if(do_comma) {
	    *outptr = ',';
	    outptr++;
	  };
/*
 * Copy the argument into the output string.
 */
	  strncpy(outptr, args[arg].arg_ptr, args[arg].arg_len);
	  outptr += args[arg].arg_len;
	};
/*
 * Substitute the count of the number of arguments.
 */
      } else {
	char buf[10];  /* A work buffer */
	int buf_len;   /* The length of the contents of buf[] */
/*
 * Write the number of arguments into a work buffer, and get its length.
 */
	sprintf(buf, "%d", narg);
	buf_len = strlen(buf);
/*
 * Check that there is room in the output buffer to append the contents
 * of the work buffer.
 */
	if((outptr-outbuff) + buf_len >= nmax-1) {
	  lprintf(stderr,
		  "Buffer overflow while pre-proccessing:\n'%.76s...'\n",
		  outbuff);
	  return NULL;
	};
	strncpy(outptr, buf, buf_len);
	outptr += buf_len;
      };
      break;
    default:
      *outptr++ = *inptr++;   /* Copy normal characters */
      break;
    };
  };
/*
 * Terminate the sub-line.
 */
  *outptr = '\0';
  return inptr;
}

/*.......................................................................
 * This function is called by newline() whenever the first character in
 * a physical line is the # character. This signals the start of a
 * pre-processor directive to deal with new macro definitions, deletions,
 * and listings. To add a new macro name or delete an existing one, the
 * user is required to type '#+new_name definition' or '#-old_name'
 * respectively. To obtain a listing of macro names with their definitions
 * the user is required to type '#?abbreviated_name'. Anything else is an
 * error.
 */
static int add_macro(char inbuff[])
{
  size_t slen;
  char retv, *cptr, *inptr;
  int top,bot,i,tab_pos;
  Table *ttst;
/*
 * Check that a recognised directive follows the #.
 */
  switch (inbuff[1]) {
  case '?': case '-': case '+':
    break;
  default:
    lprintf(stderr,"Unrecognised pre-processor directive '%.2s...'\n", inbuff);
    lprintf(stderr, " '#+name definition'\tto add a macro,\n '#-existing_name'\tto delete one,\n '#?abbreviated_name'\tto list macros matching the abbreviation.\n");
    return -1;
  };
/*
 * Get the name that follows the directive character.
 */
  inptr = inbuff+2;
  if( (slen=get_name(&inptr)) == -1)
    return -1;
/*
 * Search for the name that follows the directive in the macro symbol
 * table. (It doesn't matter if the symbol is already aliased to a
 * command name, function name, or variable name. It is intended that
 * one should be able to re-define command names with macro aliases.
 */
  if(slen != 0)
    retv = find_symbol(namebuf, macro_table, num_macro, &bot, &top);
  else {
    retv='n';
    bot=0;
    top=num_macro-1;
  };
/*
 * Obey the directive.
 */
  switch (inbuff[1]) {
/*
 * List a subset of macro's?
 */
  case '?':
    if(retv == 'n' && slen != 0)
      lprintf(stderr, "No macro begins with '%s'\n",namebuf);
    else {
      for(i=bot;i<=top;i++)
	lprintf(stdout, "#%s = %s\n", macro_table[i]->name,
		TABSTR(macro_table[i]));
    };
    break;
/*
 * Delete an existing macro?
 */
  case '-':
    switch (retv) {
    case 'n':
      lprintf(stderr, "#-%s matches no existing macro name.\n",namebuf);
      break;
    case 'a':
      lprintf(stderr, "#-%s is ambiguous - it could match:\n",namebuf);
      for(i=bot;i<=top;i++)
	lprintf(stderr, "#%s = %s\n", macro_table[i]->name,
		TABSTR(macro_table[i]));
      return -1;
      break;
    default:
/*
 * To delete the macro, first zap the name and type.str strings, then
 * zap the table entry and finally pull the macro entries above it,
 * down over the vacant slot.
 */
      free(macro_table[bot]->name);
      free(TABSTR(macro_table[bot]));
      free(macro_table[bot]);
      for(i=bot+1;i<num_macro;i++)
	macro_table[i-1]=macro_table[i];
      num_macro--;
      break;
    };
    break;
/*
 * Add a new macro definition?
 */
  case '+':
/*
 * If no name was sent, signal an error.
 */
    if(slen==0) {
      lprintf(stderr, "No macro symbol name given: '%.6s...'\n",inbuff);
      return -1;
    };
/*
 * Does the macro symbol already exist? If so the new definition will supersede
 * the old one. Keep the table entry and its name field intact, but free its
 * definition entry from the type.str field.
 */
    if(retv == 'e') {
      tab_pos=bot;
      free(TABSTR(macro_table[tab_pos]));
    }
/*
 * If a new table entry is required, shift the entries above the required
 * position out of the way to make way for it. Then allocate a new
 * table structure and its name field.
 */
    else {
/*
 * Attempt to allocate memory for the new table entry that
 * will hold the new macro definition. For the moment make it
 * point to a null string.
 */
      if( (ttst=table_alloc(0,namebuf)) == NULL) {
	lprintf(stderr, "Unable to define macro name: %s\n",namebuf);
	return -1;
      };
      TABITEM(ttst) = null_string;
/*
 * Make room for the new entry in the macro symbol table.
 */
      tab_pos = (retv == 'n') ? top : bot;
      if(up_shift( macro_table, &num_macro, MAX_MACRO, tab_pos) == -1) {
	lprintf(stderr, "Unable to add macro name: %s\n",namebuf);
	free(ttst->name);
	free(ttst);
	return -1;
      };
/*
 * Install the new entry.
 */
      macro_table[tab_pos]=ttst;
    };
/*
 * The definition to be aliased with the macro name should reside
 * in the remainder of the command line, following the macro name.
 * Skip leading white-space.
 */
    while(*inptr && isspace((int)*inptr))
      inptr++;
/*
 * Allocate memory for a copy of the macro definition.
 */
    if( (cptr=stralloc(strlen(inptr))) == NULL) {
      lprintf(stderr, "Insufficient memory for macro alias to: %s\n",namebuf);
      return -1;
    };
/*
 * Make the copy of the macro definition.
 */
    strcpy(cptr, inptr);
/*
 * Install the definition in the macro table.
 */
    TABITEM(macro_table[tab_pos]) = cptr;
    break;
  };
  return no_error;
}

/*.......................................................................
 * Read a new line from the current input stream. This includes
 * concatenating continuation lines. No \n characters will remain in the
 * buffer and the buffer will be '\0' terminated.
 *
 * Input:
 *  buff    char *   The buffer to read into.
 *  nmax     int     The dimension of 'buff'.
 *  stream  FILE *   The text stream to read from. If the stream is
 *                   NULL (used to indicate an internal macro file)
 *                   end of file is returned.
 *  prompt  char *   The prompt.
 * Output:
 *  return   int       0 - OK.
 *                     1 - End of file.
 *                    -1 - Error.
 */
int lexgets(char *buff, int nmax, FILE *stream, char *prompt)
{
  char *sptr;    /* Pointer to start of latest string from fgets() */
  char *eptr;    /* Pointer to end of latest string from fgets() */
  int eol=0;     /* True when end of line has been detected */
  int c;         /* A single character in buff or from the input stream */
/*
 * Clear the error indicator.
 */
  errno = 0;
/*
 * Read lines from 'stream' until all continuation lines have been
 * concatentated, the buffer is full, or end of file is reached.
 */
  for(sptr=buff; !eol; sptr=eptr) {
/*
 * Compute the remaining space left in the input buffer.
 */
    int nleft = nmax - (sptr-buff);
    char *tmpbuf;
/*
 * Read the input line.
 */
    if(stream==stdin) {
#if USE_TECLA == 1
      tmpbuf = gl_get_line(gl, prompt, NULL, -1);
#else
      fputs(prompt, stdout);
      tmpbuf = fgets(sptr, nleft, stream);
#endif
    } else if(stream) {
      tmpbuf = fgets(sptr, nleft, stream);
    } else {
      tmpbuf = NULL;
    };
/*
 * Continuation lines are introduced with a '?' prompt.
 */
    prompt = "?";
/*
 * Did we reach the end of the file?
 */
    if(!tmpbuf) {
      buff[0] = '\0';
      if(stream!=NULL)
	clearerr(stream);
      return 1;
    };
#if USE_TECLA == 1
/*
 * When using gl_get_line() to read an input line from stdin, a private
 * string is returned. Copy this into the stream buffer.
 */
    if(stream == stdin)
      strncpy(sptr, tmpbuf, nleft);
#endif
/*
 * Locate the end of the string.
 */
    for(eptr=sptr; *eptr; eptr++)
      ;
/*
 * If the buffer was not terminated with a newline character then
 * either the buffer was too short to receive the whole line, or
 * end of file was reached. If only white space or no further
 * characters separate the end of the read input from the end of 
 * file or newline, then the read command line is complete, otherwise
 * signal buffer overflow and abort.
 */
    if(*sptr=='\0' || (stream!=stdin && eptr[-1]!='\n')) {
      do {
	c=fgetc(stream);
      } while(c!='\n' && c!=EOF && isspace(c));
/*
 * More text before end of line?
 */
      if(c!='\n' && c!=EOF) {
	lprintf(stderr, "lexgets: Input line too long for input buffer\n");
/*
 * Skip to the real end of line.
 */
	do {
	  c=fgetc(stream);
	} while(c!='\n' && c!=EOF);
	return -1;
      };
    };
/*
 * Skip trailing spaces - keep eptr pointing at the position at which
 * the '\0' terminator is to be placed.
 */
    while(eptr>sptr && isspace((int)eptr[-1]))
      eptr--;
/*
 * If the last non-white-space character is the escape character '\\'
 * then a continuation line follows.
 */
    if(eptr>sptr && eptr[-1]=='\\') {
      eptr--;
    } else {
      eol = 1;  /* No continuation line expected */
    };
/*
 * Terminate the string.
 */
    *eptr = '\0';
  };
  return 0;
}

/*.......................................................................
 * Push a line onto the command-line stack. This involves creating a
 * new command input level. All existing lines are preserved for when
 * the pushed line has been executed.
 *
 * Input:
 *  fp        FILE *  If the subsequent commands are to be read from
 *                    a file, provide its file pointer here. Otherwise
 *                    send NULL. On error fp will be closed.
 *  comstr    char *  The command line to install - or NULL if the first
 *                    command is to come from a file.
 *  filename  char *  The name of the script file from which commands
 *                    will be read, or NULL if from stdin.
 *  argstr    char *  A string containing comma-delimited arguments
 *                    to this command level - or NULL if no arguments
 *                    are to be supplied. This string MUST NOT
 *                    be changed until this command level has been
 *                    exited.
 * Output:
 *  return     int    0 - OK.
 *                   -1 - Error.
 */
int push_command(FILE *fp, const char *comstr, const char *filename,
		 const char *argstr)
{
  int slen;    /* Length of comstr */
/*
 * First check that there is an input buffer available.
 */
  if(comlev >= MAX_LEV-1) {
    lprintf(stderr, "push_command: No more command buffers available.\n");
    if(fp && fp!=stdin) fclose(fp);
    return -1;
  };
/*
 * Check that there is enough room to record the script filename, if any.
 */
  if(filename && strlen(filename) > MAXFNAME) {
    lprintf(stderr, "push_command: Filename too long.\n");
    return -1;
  };
/*
 * Determine the length of the new command string.
 */
  slen = comstr ? strlen(comstr) : 0;
/*
 * Initialize the new command level.
 */
  comlev++;
  com[comlev].unit = fp;
  com[comlev].in_ptr=com[comlev].inbuff;
  com[comlev].was_eol = (slen==0);
  com[comlev].line_no = 0;
  if(filename)
    strcpy(com[comlev].script, filename);
  else
    com[comlev].script[0] = '\0';
/*
 * Check that the string length does not exceed the size of the command
 * buffers.
 */
  if(slen > MAX_LINE) {
    lprintf(stderr, "push_command: Command string too long:\n%s\n", comstr);
    com_close();
    return -1;
  };
/*
 * If an argument string was provided, split it into the
 * arguments for this command level.
 */
  if(argstr==NULL) {
    com[comlev].nargs = 0;
  } else {
    if(pp_args(argstr, com[comlev].args, &com[comlev].nargs, MAX_ARGS)
       == -1) {
      com_close();
      return -1;
    };
  };
/*
 * Copy the command line into the new buffer.
 */
  strcpy(com[comlev].inbuff, slen>0 ? comstr:"");
  return no_error;
}

/*.......................................................................
 * This function is intended for use with new_Pager(), but it can also be
 * used by other applications to pause output. It asks the
 * user whether they want to continue. The user can press 'q' or 'Q'
 * to quit, CR to continue, or a command line that will be staged for
 * compilation, before quitting.
 *
 * Output:
 *  return     int    0 - Continue paging.
 *                    1 - Stop paging.
 *                    2 - Run external pager if available.
 */
int pause_output(void)
{
  static char endline[MAX_LINE]; /* Command line typed by user */
/*
 * Prompt the user and await their request.
 */
  fprintf(stdout,"Press return to continue, Q [or command] to quit, or P to page.\n");
/*
 * Read reply.
 */
  if(lexgets(endline, MAX_LINE, stdin, "#")==0) {
    if(endline[0] == '\0') {
      return 0;  /* Continue listing */
    } else if(tolower(endline[0])=='q' && endline[1]=='\0') {   /* Quit? */
      return 1;
    } else if(tolower(endline[0])=='p' && endline[1]=='\0') {
      return 2;  /* Run external pager ? */
    } else {     /* Command enterred - stage it for compilation and quit */
      push_command(NULL, endline, NULL, NULL);
      return 1;
    };
  };
  return 1;  /* Quit */
}

#if USE_TECLA == 1
/*.......................................................................
 * When using the tecla library to read terminal input, the following
 * callback is used to implement tab completion.
 */
static CPL_MATCH_FN(tecla_match_fn)
{
  char *stab;        /* The symbol table to look up the word in */
  CplFileArgs cfa;   /* The control parameters of cpl_file_completions() */
  int word_start;    /* The start of the word to be completed */
  int word_len;      /* The length of the word to be completed */
  int c;             /* A character being tested */
  int i,j;
/*
 * Get the symbol table to look up command names in.
 */
  stab = data;
/*
 * Search backwards for the start of a symbol to be completed.
 */
  for(word_start = word_end - 1;
      word_start>=0 && (isalnum((c=line[word_start])) || c == '_');
      word_start--)
    ;
/*
 * We will have gone one character too far.
 */
  word_start++;
/*
 * How long is the word that is to be completed?
 */
  word_len = word_end - word_start;
/*
 * Copy the word into the name buffer if it will fit.
 */
  if(word_len > MAX_NAME)
    return 0;
  strncpy(namebuf, line + word_start, word_len);
  namebuf[word_len] = '\0';
/*
 * Search for the first non-space character that precedes the start of
 * the word.
 */
  for(i=word_start-1; i>=0 && isspace((int) line[i]); i--)
    ;
/*
 * If the user is entering the first word of a new logical line,
 * attempt to complete the word as a command name.
 */
  if(i<0 || line[i]==';') {
    int lolim, uplim;    /* The first and last ambiguous matching entries */
/*
 * Lookup the word.
 */
    switch(find_symbol(namebuf, main_table, num_main, &lolim, &uplim)) {
    case 'a': case 'f': case 'e':
/*
 * Report each of the ambiguous matches.
 */
      for(j=lolim; j<=uplim; j++) {
	Table *sym = main_table[j];
/*
 * Only report command and variable names.
 */
	switch(sym->class) {
	case FUNC:
	  if(TABFUNC(sym)->type[0]==' ' || TABFUNC(sym)->access[0]=='?') {
	    if(cpl_add_completion(cpl, line, word_start, word_end,
				  sym->name + word_len, "", " "))
	      return 1;
	  };
	  break;
	case VAR:
	  if(cpl_add_completion(cpl, line, word_start, word_end,
				sym->name + word_len, "=", " = "))
	    return 1;
	  break;
	};
      };
      break;
    default:
      return 0;
    };
/*
 * Complete words that aren't at the start of a logical line, as filenames.
 */
  } else {
/*
 * Search backwards for the start of a filename to be completed,
 * stopping when a space or comma is encountered.
 */
    for(word_start = word_end - 1;
	word_start>=0 && (c=line[word_start]) != ',' &&
	!isspace((int) line[word_start]);
	word_start--)
      ;
/*
 * We will have gone one character too far.
 */
    word_start++;
/*
 * If the first character appears to be an '@', but this the first
 * non-space character at the start of the line, then it is actually the
 * 'obey-file' operator and should not be included in the filename.
 */
    if(line[word_start] == '@') {
      for(i=word_start-1; i>=0 && isspace((int) line[i]); i--)
	;
      if(i<0 || line[i]==';')
	word_start++;
    };
/*
 * Configure the filename completion callback.
 */
    cpl_init_FileArgs(&cfa);
    cfa.file_start = word_start;
    if(cpl_file_completions(cpl, &cfa, line, word_end))
      return 1;
  };
  return 0;
}
#endif
