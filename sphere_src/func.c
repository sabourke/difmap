/*
  This file contains an assortment of user functions and user
  accessible variables concerned with these functions.
*/

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "sphere.h"
#include "helpdir.h"
#include "utils.h"
#include "help.h"

/*
  Declare variables that are to be aliased as user variables below. Only
  float, integer, char and logical variables are supported.
  NB. character strings must NOT be initialised here unless they are marked
  as R_ONLY parameters. This is to allow variable length strings where
  the previous string is often free'd first on the assumption that the
  memory for the string was allocated using malloc(), not by the
  compiler).
*/

extern float pi;
static char true = 1;
static char false = 0;
char debug = 0;
static int ii = 0;
static int jj = 0;
static float xx=0.0, yy=0.0;
static char wrap_print_output = 1;

static Descriptor genv_type[] = {
   {'l' , '0' ,R_ONLY ,1, {1,1,1}, &true},
   {'l' , '0' ,R_ONLY ,1, {1,1,1}, &false},
   {'l' , '0' ,NO_DEL ,1, {1,1,1}, &debug},
   {'i' , '0' ,NO_DEL ,1, {1,1,1}, &ii},
   {'i' , '0' ,NO_DEL ,1, {1,1,1}, &jj},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &xx},
   {'f' , '0' ,NO_DEL ,1, {1,1,1}, &yy},
   {'l' , '0' ,NO_DEL ,1, {1,1,1}, &wrap_print_output},
};

/*
  In the same order as the above array of types, define the array
  of user names for the arrays.
*/

static char *genv_name[] = {
   "true",
   "false",
   "debug",
   "i",
   "j",
   "x",
   "y",
   "wrap_print_output",
};

/*
  Declare the user functions here.
*/

static Template(exit_fn);
static Template(quit_fn);
static Template(type_fn);
static Template(print_fn);
static Template(error_fn);
static Template(date_fn);
static Template(getenv_fn);
static Template(query_fn);
static Template(system_fn);
static Template(index_fn);
static Template(len_fn);
static Template(flagdel_fn);
static Template(dim_fn);
static Template(strnum_fn);
static Template(near_fn);
static Template(newlog_fn);
static Template(makeidx_fn);
static Template(aprop_fn);
static Template(prompt_fn);

/*
  Declare the function types below.
*/

static Functype genf_type[] = {
   {NULL, DECLARE     , 0,3,    "ci",   "00",     "vv",   1 },
   {NULL, DECLARE     , 0,3,    "fi",   "00",     "vv",   1 },
   {NULL, DECLARE     , 0,3,    "li",   "00",     "vv",   1 },
   {NULL, DECLARE     , 0,3,    "ii",   "00",     "vv",   1 },
   {NULL, START_BLOCK , 1,1,    " l",   " 0",     " v",   1 },
   {NULL, START_BLOCK , 0,0,    " l",   " 0",     " v",   1 },
   {NULL, START_BLOCK , 0,0,    " l",   " 0",     " v",   1 },
   {NULL, START_BLOCK , 1,1,    " l",   " 0",     " v",   1 },
   {NULL, END_BLOCK   , 1,1,    " l",   " 0",     " v",   1 },
   {NULL, END_BLOCK   , 0,1,    " l",   " 0",     " v",   1 },
   {NULL, END_BLOCK   , 1,1,    " l" ,  " 0",     " v",   1 },
   {NULL, END_BLOCK   , 0,0,    " " ,   " ",      " ",    1 },
   {NULL, CONT_BLOCK  , 0,0,    " " ,   " ",      " ",    1 },
   {NULL, BRK_BLOCK   , 0,0,    " " ,   " ",      " ",    1 },
   {NULL, STOP_EXE    , 0,0,    " " ,   " ",      " ",    1 },
   {NULL, WHATVAR     , 0,0,    " " ,   " ",      " ",    1 },
   {NULL, HELP        , 0,0,    " " ,   " ",      " ",    1 },
   {aprop_fn, NORM    , 1,1,    " C" ,  " 0",     " v",   1 },
   {exit_fn,  NORM    , 0,0,    "  ",   "  ",     "  ",   1 },
   {quit_fn,  NORM    , 0,0,    "  ",   "  ",     "  ",   1 },
   {type_fn,   NORM   , 1,MAXARG, " *", " *",     " v",   1 },
   {print_fn,  NORM   , 1,MAXARG, " *", " *",     " v",   1 },
   {error_fn,  NORM   , 1,MAXARG, " *", " *",     " v",   1 },
   {date_fn,   NORM   , 0,0,    "c",    "0",      "v",    1 },
   {getenv_fn, NORM   , 1,1,    "cC",   "00",     "vv",   1 },
   {query_fn,  NORM   , 1,1,    "lc",   "00",     "vv",   0 },
   {system_fn, NORM   , 1,1,    " c",   " 0",     " v",   1 },
   {index_fn,  NORM   , 2,2,     "icc", "000",    "vvv",  0 },
   {len_fn,    NORM   , 1,1,     "ic",  "00",     "vv",   0 },
   {flagdel_fn,NORM   , 3,11," lififififif"," 10*0*0*0*0*"," vvNvNvNvNvN",1},
   {dim_fn,    NORM   , 2,2,     "ii*",   "00*",    "vvv", 1 },
   {strnum_fn, NORM   , 1,2,     "cni",   "000",    "vvv", 1 },
   {near_fn,   NORM   , 2,2,     "iff",   "010",    "vvv", 1 },
   {newlog_fn, NORM   , 0,1,     " C",    " 0",     " v",  1 },
   {makeidx_fn,NORM   , 0,MAXARG," C",    " 0",     " v",  1 },
   {prompt_fn, NORM   , 1,2,     "ccc",   "000",    "vvv", 1 },
};

/*
  In the same order as the above array of types, define the array
  of user names for the functions.
*/

static char *genf_name[] = {
   "string",
   "float",
   "logical",
   "integer",
   "while",
   "repeat",
   "do",
   "if",
   "elseif",
   "else",
   "until",
   "end",
   "continue",
   "break",
   "stop",
   "varlist",
   "help",
   "apropos",
   "exit",
   "quit",
   "type",
   "print",
   "error",
   "date",
   "getenv",
   "query",
   "system",
   "index",
   "len",
   "flagdel",
   "dim",
   "strnum",
   "nearest",
   "logfile",
   "makeindex",
   "prompt_user",
};

/*
  Record the above declarations etc for this module in a global
  structure for use when building the main symbol table.
*/

Module m_general = {
  "general",
  HELP_DIR,
  NULL, 0,
  genv_type, genv_name, COUNT(genv_name),
  genf_type, genf_name, COUNT(genf_name)
};

/*
  Local variables and functions.
*/

#define MAXSTR 132
static char work_string[MAXSTR];

/*.......................................................................
  A simple means by which to dump descriptors to the users terminal.
*/
static Template(type_fn)
{
        int i,j,k,adim[3];
        int np;
        float *fval;
	int *ival;
        char **cval,atyp,*lval;
/*
  Print out each descriptor one at a time.
*/
        for(np=0; np<npar; np++) {
/*
  Take local copies of relevant parts of the descriptor.
*/
          atyp = invals[np]->atyp;
          adim[0] = invals[np]->adim[0];
          adim[1] = invals[np]->adim[1];
          adim[2] = invals[np]->adim[2];
	  ival = INTPTR(invals[np]);
          fval = FLTPTR(invals[np]);
          lval = LOGPTR(invals[np]);
          cval = STRPTR(invals[np]);
/*
  With due allowance for the type of variable show the value of the variable
  associated with the current descriptor. First get a pointer to the
  start of the value.
*/
          switch (atyp){
          case 'f': case 'l': case 'c': case 'i':
/*
  Write each row of the first dimension of the variable as a single line,
  and leave a space between lines when the third dimension is incremented.
*/
            for(i = 0;i<adim[2];i++) {
              for(j = 0;j<adim[1];j++) {
                lprintf(stdout, "\t\t");
                for(k = 0;k<adim[0];k++) {
/*
  Write the latest value to the terminal.
*/
                  switch (atyp) {
                  case 'f':
                    lprintf(stdout, "%3.2g ",*(fval++));
                    break;
                  case 'i':
                    lprintf(stdout, "%d ",*(ival++));
                    break;
                  case 'l':
                    lprintf(stdout, "%s ", *(lval++) ? "TRUE" : "FALSE" );
                    break;
                  case 'c':
                    lprintf(stdout, "\"%s\" ", *(cval++) );
                    break;
                  };
/*
  Stop on user interrupt (ctrl-c etc..).
*/
		  if(no_error) return no_error;
                };
                lprintf(stdout, "\n");
              };
              lprintf(stdout, "\n");
            };
            break;
          default:
            lprintf(stderr, "Unknown variable type: %c",atyp);
            return -1;
          };
        };
        return no_error;
}

/*.......................................................................
  A test command to print the values of as many user expressions as
  given as arguments to the command.
*/
static Template(print_fn)
{
        int i;
        int np;
        int nvals;
        char atyp,nch=0;
/*
  Print out each descriptor one at a time.
*/
        for(np=0; np<npar; np++) {
/*
  Find out how many values there are in the current descriptor.
*/
          nvals = 1;
          for(i=0;i<3;i++)
            nvals *= invals[np]->adim[i];
/*
  Find out the variable type.
*/
          atyp = invals[np]->atyp;
/*
  With due allowance for the type of variable show the value of the variable
  associated with the current descriptor. First get a pointer to the
  start of the value.
*/
          for(i=0;i<nvals;i++) {
            switch (atyp){
            case 'f':
              nch += lprintf(stdout, "%g ", *(FLTPTR(invals[np])+i) );
              break;
            case 'i':
              nch += lprintf(stdout, "%d ", *(INTPTR(invals[np])+i) );
              break;
            case 'l':
              nch += lprintf(stdout, "%s ", *(LOGPTR(invals[np])+i) ? "TRUE":"FALSE" );
              break;
            case 'c':
              nch += lprintf(stdout, "%s ", *(STRPTR(invals[np])+i) );
              break;
            };
/*
 * If enabled, force a line break when the line gets too long.
 */
            if(nch > 60 && wrap_print_output) {
              nch=0;
              lprintf(stdout, "\n");
            };
/*
  Stop on user interrupts (ctrl-c etc..).
*/
	    if(no_error) return no_error;
          };
        };
        lprintf(stdout, "\n");
        return no_error;
}

/*.......................................................................
 * Print an error message by calling the user "print" command, then raise
 * an error condition by returning -1.
 */
static Template(error_fn)
{
  (void) print_fn(invals, npar, outvals);
  return -1;
}

/*.......................................................................
  Provide for a user request to exit the program.
*/
static Template(exit_fn)
{
  lprintf(stderr, "Exiting program\n");
  closedown(0, DO_EXIT);
  return no_error;
}

/*.......................................................................
  Provide for a user request to quit the program.
*/
static Template(quit_fn)
{
  lprintf(stderr, "Quitting program\n");
  closedown(0, DO_QUIT);
  return no_error;
}

/*.......................................................................
  Return the current date and time.
*/
static Template(date_fn)
{
        static time_t tp;
/*
  Get the current time.
*/
        if(time(&tp) == -1) {
          lprintf(stderr, "Sorry the date is not available on your machine.\n");
          return -1;
        };
/*
  Convert it into an ascii string representation of the current date.
*/
        strcpy(work_string,ctime(&tp));
/*
  ctime() places a \n at the end of the date string - zap it.
*/
        *strchr(work_string,'\n')='\0';
/*
  Allocate a one element character string array for the return value.
*/
        if( (VOIDPTR(outvals)=valof_alloc(1,'c')) == NULL )
          return -1;
/*
  We need to allocate a return string long enough to hold the date.
*/
        if( (*STRPTR(outvals)=stralloc(strlen(work_string))) == NULL)
          return -1;
/*
  Copy the date string.
*/
        strcpy(*STRPTR(outvals), work_string);
        return no_error;
}

/*.......................................................................
  Return the value of the environment variable named by the user.
*/
static Template(getenv_fn)
{
        char *cptr;
/*
  Ask for the environment variable.
*/
        if((cptr=getenv(*STRPTR(invals[0]))) == NULL) {
          lprintf(stderr, "getenv: Unable to get equivalence of: '%s'.\n", *STRPTR(invals[0]));
          return -1;
        };
/*
  Allocate a one element character string array for the return value.
*/
        if( (VOIDPTR(outvals)=valof_alloc(1,'c')) == NULL )
          return -1;
/*
  We need to allocate a return string long enough to hold the equivalence
  string.
*/
        if( (*STRPTR(outvals)=stralloc(strlen(cptr))) == NULL)
          return -1;
/*
  Copy the equivalence string to the return descriptor.
*/
        strcpy(*STRPTR(outvals), cptr);
        return no_error;
}

/*.......................................................................
  Query the user, given the prompt sent as the only argument and return
  logical true if the user answers with y or presses carriage return
  without any string. Otherwise if the user enters n return false. If
  Anything else is typed re-prompt. Append the string '(y/n)?' to the
  user's prompt.
*/
static Template(query_fn)
{
        int was_yes;
/*
  Read the response.
*/
	was_yes = ask_user(*STRPTR(invals[0]));
	if(was_yes == -1)
	  return -1;
	*LOGPTR(outvals) = was_yes;
	return no_error;
}

/*.......................................................................
 * Prompt the user for a string reply.
 *
 * Input:
 *  prompt     char *   The string the prompt the user with.
 *  defstr     char *   The optional default string to be returned if the
 *                      user doesn't enter anything. If this is omitted then
 *                      the user will be re-prompted for input.
 */
static Template(prompt_fn)
{
  char *prompt;       /* The requested prompt string */
  char *defstr=NULL;  /* The requested default reply string */
  char *answer;       /* The user's answer string */
/*
 * Get the arguments.
 */
  switch(npar) {
  default:
  case 2:
    defstr = *STRPTR(invals[1]);
  case 1:
    prompt = *STRPTR(invals[0]);
  };
/*
 * Acquire the user's response into a dynamically allocated string.
 */
  answer = prompt_user(prompt, defstr);
  if(!answer)
    return -1;
/*
 * Allocate a one element character string array for the return value.
 */
  if( (VOIDPTR(outvals)=valof_alloc(1,'c')) == NULL ) {
    free(answer);
    return -1;
  };
/*
 * Install the string into the return descriptor.
 */
  *STRPTR(outvals) = answer;
  return no_error;
}

/*.......................................................................
  Send the user's string argument to the operating system to be executed.
  Returns -1 if there is no comand processor to handle the command.
*/
static Template(system_fn)
{
	system(*STRPTR(invals[0]));
	return no_error;
}

/*.......................................................................
  Given two string arguments, return the index position of the second
  within the first or zero if not found.
*/
static Template(index_fn)
{
        char *cptr;
	cptr=strstr(*STRPTR(invals[0]), *STRPTR(invals[1]));
	if(cptr == NULL)
	  *INTPTR(outvals) = 0;
	else
	  *INTPTR(outvals) = (cptr - *STRPTR(invals[0])) + 1;
	return no_error;
}

/*.......................................................................
  Return the length of the string sent.
*/
static Template(len_fn)
{
        *INTPTR(outvals) = strlen(*STRPTR(invals[0]));
	return no_error;
}

/*.......................................................................
  Given an axis designation, an array of flags and a data
  array of the same size, delete the elements specified, by shufling the
  elements above them down one place over the deleted element.
*/
static Template(flagdel_fn)
{
        int npts, i, p[3], dim[3], axis, nflag, arg;
	int *ip, *jp, *kp, *idim, *jdim, *kdim, *flag_el;
	float *outptr, *inptr;
	char *flagptr;
/*
  Determine the number of elements in the flag array and get a pointer
  to the start of the flag array.
*/
	npts = invals[0]->adim[0];
	flagptr = LOGPTR(invals[0]);
/*
  Count the number of points to be flagged.
*/
	for(nflag=0,i=0; i<npts; i++) {
	  if(flagptr[i]) nflag++;
	};
/*
  If the deletions would leave the data array with no elements
  signal an error.
*/
	if(nflag == npts ) {
	  lprintf(stderr, "No un-flagged elements?\n");
	  return -1;
	};
/*
  No deletions requested?
*/
	if(nflag == 0) return no_error;
/*
  Check that each axis specification is legal and that there are
  the same number of elements along the axis as in the flag array.
*/
	for(arg=1; arg<npar; arg += 2) {
/*
  Get the latest axis designation.
*/
	  axis = *INTPTR(invals[arg]);
/*
  Check its legality.
*/
	  if(axis < 0 || axis > 2) {
	    lprintf(stderr, "remove(): Axis specification (%d) out of bounds.\n", axis);
	    return -1;
	  };
/*
  Check that the data array has the same number of points as the
  flag array along axis 'axis';
*/
	  if(npts != invals[arg+1]->adim[axis]) {
	    lprintf(stderr, "The flag and data arrays differ in size.\n");
	    return -1;
	  };
	};
/*
  To avoid inefficient array indexing, get pointers to each
  element of axis increment arrays.
*/
	ip = &p[0];
	jp = &p[1];
	kp = &p[2];
	idim = &dim[0];
	jdim = &dim[1];
	kdim = &dim[2];
/*
  Apply the deletions to each array in turn.
*/
	for(arg=1; arg<npar; arg += 2) {
/*
  Get the axis designation.
*/
	  axis = *INTPTR(invals[arg]);
/*
  Get two pointers to the start of the data array.
*/
	  outptr = inptr = FLTPTR(invals[arg+1]);
/*
  Determine the number of elements on each axis of the data array.
*/
	  for(i=0; i<3; i++)
	    dim[i]=invals[arg+1]->adim[i];
/*
  Get a pointer to the index element for the flag array.
*/
	  flag_el = &p[axis];
/*
  Go through each element of the data array. At each one check that
  if it is flagged. If not flagged copy it to *ouptr and increment
  outptr, otherwise don't increment outptr.
*/
	  for(*kp=0; *kp < *kdim; (*kp)++) {
	    for(*jp=0; *jp < *jdim; (*jp)++) {
	      for(*ip=0; *ip < *idim; (*ip)++) {
		if(!flagptr[*flag_el]) {
		  *outptr = *inptr;
		  outptr++;
		};
		inptr++;
	      };
	    };
	  };
/*
  Change the declared dimensions of the data array to reflect
  the deletions.
*/
	  invals[arg+1]->adim[axis] = npts-nflag;
	};
	return no_error;
}

/*.......................................................................
  A utility function to return the current number of elements along a
  given axis of a variable.
*/
static Template(dim_fn)
{
        int axis;
/*
  Find out the axis for which the dimension is required.
*/
	axis = *INTPTR(invals[0]);
	if(axis < 0 || axis > 2) {
	  lprintf(stderr, "No such axis: %d\n",axis);
	  return -1;
	};
/*
  Return the dimension.
*/
	*INTPTR(outvals) = invals[1]->adim[axis];
	return no_error;
}

/*.......................................................................
  Return the string representation of a numeric argument. The final
  optional argument specifies the maximum precision for floating point
  numbers.
*/
static Template(strnum_fn)
{
        int prec;
/*
  Get sprintf to write a float or int into the global work string
  declared at the top of this file. Don't test the return value
  of sprintf() - this is because some machines still have this
  function as returning a char * pointer instead the int specified
  in ANSI-C.
*/
	switch (invals[0]->atyp) {
	case 'f':
/*
  Get the required maximum-precision - default to 2 if not specified.
*/
	  prec = (npar>1) ? *INTPTR(invals[1]) : 4;
/*
  Enforce min/max precisions of 0 and 10.
*/
	  if(prec < 0)
	    prec=0;
	  else if(prec > 10)
	    prec=10;
/*
  Write a temporary version of the string, using the precision given.
*/
	  sprintf(work_string,"%.*g", prec, *FLTPTR(invals[0]));
	  break;
/*
  An int - ignore any precision argument.
*/
	case 'i':
	  sprintf(work_string,"%d", *INTPTR(invals[0]));
	  break;
	};
/*
  Allocate a one element character array of strings for the return value.
*/
        if( (VOIDPTR(outvals)=valof_alloc(1,'c')) == NULL )
          return -1;
/*
  Allocate a return string for a copy of the string representation
  of the number, currently in work_string and perform the copy.
*/
	if( (*STRPTR(outvals)=stralloc(strlen(work_string))) == NULL)
	  return -1;
	strcpy(*STRPTR(outvals), work_string);
        return no_error;
}

/*.......................................................................
  This user function takes a 1D array and a scalar number and returns the
  element which holds a value closest to that number.
*/
static Template(near_fn)
{
	float *fptr;	/* Pointer to data array */
	float fnum;	/* Number to be located */
	int npts;	/* Size of data array */
	float dmin;	/* Smallest difference between fnum and fptr[i] */
	float imin;	/* Element of smallest difference */
	float diff;	/* Current difference */
	int i;
/*
  The array comes first - get a copy of the pointer to its first element.
*/
	fptr = FLTPTR(invals[0]);
/*
  How big is the data array?
*/
	npts = invals[0]->adim[0];
/*
  Get the number to be located.
*/
	fnum = *FLTPTR(invals[1]);
/*
  Loop for the closest difference.
*/
	imin = 0;
	dmin = fabs(fnum - *fptr);
	for(i=1; i < npts; i++) {
	  fptr++;
	  diff = fabs(fnum - *fptr);
/*
  Check the latest difference against the current minimum difference.
*/
	  if(diff < dmin) {
	    dmin = diff;
	    imin = i;
	  };
	};
/*
  Return the element number with the closest match.
*/
	*INTPTR(outvals) = imin+1;
	return no_error;
}

/*.......................................................................
 * Close the previous log file and if a file name is provided, open a new
 * log file of that name.
 *
 * Input:
 *  name   char *  The name of the file. If omitted no log file is opened.
 */
static Template(newlog_fn)
{
  char *name;
/*
 * Get the name of the file, or use NULL if no name given.
 */
  name = npar>0 ? *STRPTR(invals[0]) : NULL;
/*
 * Open the log file.
 */
  return (logfile(name)==NULL && name!=NULL) ? -1 : no_error;
}

/*.......................................................................
 * Create a module index file.
 *
 * Input:
 *  name   char *  The module to be indexed.
 *  ...
 */
static Template(makeidx_fn)
{
  int waserr = 0;  /* True after an error */
  int arg;
/*
 * Index each module named.
 */
  if(npar < 1) {
    waserr = index_module(NULL);
  } else {
    for(arg=0; arg<npar; arg++) {
      if(index_module(*STRPTR(invals[arg])))
	waserr = 1;
    };
  };
  return waserr ? -1 : no_error;
}

/*.......................................................................
 * Search for the occurrence of a given keyword in help topic names
 * and one-line intro's.
 *
 * Input:
 *  key    char *    The key to look up.
 */
static Template(aprop_fn)
{
  apropos(*STRPTR(invals[0]));
  return no_error;
}

