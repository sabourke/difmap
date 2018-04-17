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
#include "utils.h"
#include "helpdir.h"
#include "matrix_blocks.h"

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
static int user_lun=0;

static Descriptor iov_type[] = {
   {'i' , '0' ,NO_DEL ,1, {1,1,1}, &user_lun}
};

/*
  In the same order as the above array of types, define the array
  of user names for the arrays.
*/

static char *iov_name[] = {
  "lun"
};

/*
  Declare the user functions here.
*/

static Template(infile_fn);
static Template(outfile_fn);
static Template(close_fn);
static Template(cat_fn);
static Template(eof_fn);
static Template(rewind_fn);
static Template(read_int);
static Template(read_float);
static Template(read_string);
static Template(read_array);
static Template(fprintf_fn);
static Template(fread_fn);
static Template(read_fn);
static Template(write_fn);
static Template(search_fn);
static Template(read_table_fn);

/*
  Declare the function types below.
*/

static Functype iof_type[] = {
   {infile_fn, NORM   , 1,2,    "iCl",  "000",    "vvv",  1 },
   {outfile_fn,NORM   , 1,3,    "iCll", "0000",   "vvvv", 1 },
   {close_fn,  NORM   , 1,1,    " i",   " 0",     " v",   1 },
   {cat_fn,    NORM   , 0,0,    " ",    " ",      " ",    1 },
   {eof_fn,    NORM   , 1,1,    "li",   "00",     "vv",   1 },
   {rewind_fn, NORM   , 1,1,    " i",   " 0",     " v",   1 },
   {read_int,  NORM   , 0,2,    "iic",  "000",    "vvv",  0 },
   {read_float,NORM   , 0,2,    "fic",  "000",    "vvv",  0 },
   {read_string,NORM  , 0,2,    "cic",  "000",    "vvv",  0 },
   {read_array,NORM   , 2,2,    "fif",  "100",    "vvv",  1 },
   {fprintf_fn,NORM   , 2,MAXARG," ic*"," 000",   " vvv", 1 },
   {fread_fn,  NORM   , 2,MAXARG," ic*"," 000",   " vvv", 1 },
   {read_fn,   NORM   , 2,MAXARG," if", " 03",    " vr",  1 },
   {write_fn,  NORM   , 2,MAXARG," if", " 03",    " vv",  1 },
   {search_fn, NORM   , 2,2,     "lic", "000",    "vvv",  1 },
   {read_table_fn, NORM, 2,3,   " Cfi", " 020",   " vNv", 1 }
};

/*
  In the same order as the above array of types, define the array
  of user names for the functions.
*/

static char *iof_name[] = {
   "infile",
   "outfile",
   "close",
   "catalogue",
   "eof",
   "rewind",
   "read_int",
   "read_float",
   "read_string",
   "read_array",
   "fprintf",
   "fread",
   "read",
   "write",
   "search",
   "read_table"
};

/*
  Record the above declarations etc for this module in a global
  structure for use when building the main symbol table.
*/

Module m_iolib = {
  "file_io",
  HELP_DIR,
  NULL, 0,
  iov_type, iov_name, COUNT(iov_name),
  iof_type, iof_name, COUNT(iof_name)
};

/*
  Local variables and functions.
*/

static FILE *file_ptr;
static int file_lun,is_text;

/*.......................................................................
  Open a user file for input. The first argument is the filename, the
  second and third - optional flags to determine whether to open the
  file for write/append and whether to open it asa binary file.
*/
static Template(infile_fn)
{
        char rwa,is_text;
/*
  Enforce default modes if non-given. The default is to open the file
  as a text file.
*/
	if(npar > 1)
	  is_text = !(*LOGPTR(invals[1]));
	else
	  is_text = 1;
	rwa = 0;
/*
  Attempt to open the file.
*/
	if( (*INTPTR(outvals)=file_open(rwa, is_text, *STRPTR(invals[0]))) == -1)
	  return -1;
	return no_error;
}

/*.......................................................................
  Open a user file for output. The first argument is the filename, the
  second and third - optional flags to determine whether to open the
  file for write/append and whether to open it asa binary file.
*/
static Template(outfile_fn)
{
  char *name = NULL; /* The name of the file */
  char rwa = 1;      /* The file access mode 0=read,1=write,2=append */
  char is_text = 1;  /* Whether to open the file as a text or binary file */
/*
 * Intepret command-line arguments.
 */
  switch(npar) {
  case 3:
    is_text = !(*LOGPTR(invals[2]));
  case 2:
    rwa = (*LOGPTR(invals[1])) ? 2:1;
  case 1:
    name = *STRPTR(invals[0]);
    break;
  default:
    lprintf(stderr, "outfile(): Wrong number of arguments.\n");
    return -1;
  };
/*
 * Attempt to open the file.
 */
  if((*INTPTR(outvals)=file_open(rwa, is_text, name)) == -1)
    return -1;
  return no_error;
}

/*.......................................................................
  Close a user file given its unit number.
*/
static Template(close_fn)
{
        return file_close(*INTPTR(invals[0]));
}

/*.......................................................................
  Write out a catalogue of open user-files to stdout.
*/
static Template(cat_fn)
{
        file_cat();
	return no_error;
}

/*.......................................................................
  If the end of file indicator is set for the lun sent then return true
  otherwise false.
*/
static Template(eof_fn)
{
        if( (*LOGPTR(outvals)=file_check_eof(*INTPTR(invals[0]))) == -1)
	  return -1;;
	return no_error;
}

/*.......................................................................
  Rewind the file referenced by the lun sent by the user.
*/
static Template(rewind_fn)
{
        return file_rewind(*INTPTR(invals[0]));
}

/*.......................................................................
  Read a single float from the file associated with the file lun sent
  by the user. If the user does not specify a lun then default to stdin
  ie. lun=0.  A second optional argument is a format-string.
*/
static Template(read_float)
{
	file_lun = (npar>0) ? *INTPTR(invals[0]) : 0;
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "read_float(): File %d is binary\n",file_lun);
	  return -1;
	};
/*
  Specify the return type for type checking by the reading function.
*/
	outvals->atyp = 'f';
/*
  If the user specified a format string then have it obeyed.
*/
	if(npar>1) {
	  if(fmt_read(file_ptr, *STRPTR(invals[1]), 1, &outvals) == -1)
	    return -1;
	}
/*
  No format string specified, so make contrive one.
*/
	else {
/*
  Read the number.
*/
	  if(fmt_read(file_ptr, "f", 1, &outvals) == -1)
	    return -1;
	};
	return no_error;
}


/*.......................................................................
  Read a single integer from the file associated with the file lun sent
  by the user. If the user does not specify a lun then default to stdin
  ie. lun=0.  A second optional argument is a format-string.
*/
static Template(read_int)
{
	file_lun = (npar>0) ? *INTPTR(invals[0]) : 0;
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "read_int(): File %d is binary\n",file_lun);
	  return -1;
	};
/*
  Specify the return type for type checking by the reading function.
*/
	outvals->atyp = 'i';
/*
  If the user specified a format string then have it obeyed.
*/
	if(npar>1) {
	  if(fmt_read(file_ptr, *STRPTR(invals[1]), 1, &outvals) == -1)
	    return -1;
/*
  No format string specified, so make contrive one.
*/
	}
	else {
/*
  Read the number.
*/
	  if(fmt_read(file_ptr, "i", 1, &outvals) == -1)
	    return -1;
	};
	return no_error;
}

/*.......................................................................
  Read a single string from the file associated with the file lun sent
  by the user. If the user does not specify a lun then default to stdin
  ie. lun=0. A second optional argument is a format-string.
*/
static Template(read_string)
{
	file_lun = (npar>0) ? *INTPTR(invals[0]) : 0;
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "read_string(): File %d is binary\n",file_lun);
	  return -1;
	};
/*
  The reading function will attempt to delete the old string of the return value
  before assigning a new one. Since the return value is to be placed in outvals
  where no string currently resides, assign it with a null string. It also performs
  checking on the return type so set that up as well.
*/
	*STRPTR(outvals) = NULL;
	outvals->atyp = 'c';
/*
  If the user specified a format string then have it obeyed.
*/
	if(npar>1) {
	  if(fmt_read(file_ptr, *STRPTR(invals[1]), 1, &outvals) == -1)
	    return -1;
	}
/*
  No format string specified, so contrive one.
*/
	else {
/*
  Read the string.
*/
	  if(fmt_read(file_ptr, "s", 1, &outvals) == -1)
	    return -1;
	};
	return no_error;
}

/*.......................................................................
  Read an array of indeterminate size from a user file. The user is
  required to supply the file unit number the maximum size of the
  resulting array.
*/
static Template(read_array)
{
        size_t num_el, num_read;
/*
  Check the file unit number and get the file pointer.
*/
        file_lun = *INTPTR(invals[0]);
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "read_array(): File %d is binary\n",file_lun);
	  return -1;
	};
/*
  Get the max number of elements required by the user and check it.
*/
	num_el = *FLTPTR(invals[1]);
	if(num_el < 1) {
	  lprintf(stderr, "Illegal max size (%d) for return array in read_array()\n", num_el);
	  return -1;
	};
/*
  Allocate sufficient memory for the return array.
*/
	if( (VOIDPTR(outvals) = valof_alloc(num_el, 'f')) == NULL) {
	  lprintf(stderr, "Memory allocation failure in read_array()\n");
	  return -1;
	};
/*
  Read the array.
*/
	num_read=input_array(file_ptr, FLTPTR(outvals), num_el);
/*
  On error zap the temporary array.
*/
	if(num_read == -1) {
	  valof_free(outvals);
	  return -1;
	};
/*
  Fill in the descriptor attributes.
*/
	outvals->num_el = num_el;
        outvals->adim[0] = num_read;
	return no_error;
}

/*.......................................................................
  A filter for fprintf() to provide the user with formatted io.
*/
static Template(fprintf_fn)
{
	file_lun = *INTPTR(invals[0]);
	if( (file_ptr=check_lun(file_lun, 0, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "fprint(): File %d is binary\n",file_lun);
	  return -1;
	};
	return user_printf(file_ptr, *STRPTR(invals[1]), npar-2, &invals[2]);
}


/*.......................................................................
  A filter for fprintf() to provide the user with formatted io.
*/
static Template(fread_fn)
{
	file_lun = *INTPTR(invals[0]);
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Make it illegal to read from a binary file.
*/
	if(!is_text) {
	  lprintf(stderr, "fread(): File %d is binary\n",file_lun);
	  return -1;
	};
	return fmt_read(file_ptr, *STRPTR(invals[1]), npar-2, &invals[2]);
}

/*.......................................................................
  Write an array of numbers to a given text or binary file.
*/
static Template(write_fn)
{
        int dim[3],i,j,k,ntot,nn;
	float *inptr;
/*
  Get the file lun from the user's first argument.
*/
	file_lun = *INTPTR(invals[0]);
/*
  Check that it belongs to an open file with the correct attributes.
*/
	if( (file_ptr=check_lun(file_lun, 0, &is_text)) == NULL)
	  return -1;
/*
  Write each array individually.
*/
	for(nn=1; nn<npar; nn++) {
/*
  Determine the dimensions of the input array and its total
  number of elements.
*/
	  for(ntot=1,i=0; i<3; i++)
	    ntot *= (dim[i] = invals[nn]->adim[i]);
/*
  If the file is a binary file the write process is trivial.
*/
	  if(!is_text) {
	    if(fwrite(FLTPTR(invals[nn]), sizeof(float), ntot, file_ptr) < ntot) {
/*
  If the error indicator for the file is set then use the system
  error message.
*/
	      if(ferror(file_ptr) == 0)
		lprintf(stderr, "write: File write error.\n");
	      else
		perror("write: write error");
	      return -1;
	    };
	  }
/*
  Text file.
*/
	  else {
	    inptr = FLTPTR(invals[nn]);
/*
  Write each row of the first dimension of the variable as a single line,
  and leave a space between lines when the third dimension is incremented.
*/
            for(i = 0; i<dim[2]; i++) {
              for(j = 0; j<dim[1]; j++) {
                for(k = 0; k<dim[0]; k++) {
		  lprintf(file_ptr,"%f",*(inptr++));
/*
  Stop on user interrupt (ctrl-c etc..).
*/
		  if(no_error) return no_error;
/*
  If the error indicator for the file is set then use the system
  error message and return abnormally.
*/
		  if(ferror(file_ptr) != 0) {
		    perror("write: write error");
		    return -1;
		  };
                };
                lprintf(file_ptr,"\n");
              };
              lprintf(file_ptr,"\n");
            };
	  };
	};
	return no_error;
}

/*.......................................................................
  Read an array of numbers from a given text or binary file.
*/
static Template(read_fn)
{
        int dim[3],i,ntot,nn;
	float *inptr;
/*
  Get the file lun from the user's first argument.
*/
	file_lun = *INTPTR(invals[0]);
/*
  Check that it belongs to an open file with the correct attributes.
*/
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Read each array individually.
*/
	for(nn=1; nn<npar; nn++) {
/*
  Determine the dimensions of the input array and its total
  number of elements.
*/
	  for(ntot=1,i=0; i<3; i++)
	    ntot *= (dim[i] = invals[nn]->adim[i]);
/*
  If the file is a binary file the write process is trivial.
*/
	  if(!is_text) {
	    if(fread(FLTPTR(invals[nn]), sizeof(float), ntot, file_ptr) < ntot) {
/*
  If the error indicator for the file is set then use the system
  error message.
*/
	      if(ferror(file_ptr) == 0)
		lprintf(stderr, "read: File read error.\n");
	      else
		perror("read: read error");
	      return -1;
	    };
	  }
/*
  Text file.
*/
	  else {
	    inptr = FLTPTR(invals[nn]);
/*
  Read the numbers one at a time.
*/
	    for(i=0; i<ntot; i++, inptr++) {
	      if(fscanf(file_ptr, "%f",inptr) < 1) {
/*
  If the error indicator for the file is set then use the system
  error message and return abnormally.
*/
		if(ferror(file_ptr) == 0)
		  lprintf(stderr, "read: File read error.\n");
		else
		  perror("read: read error");
		return -1;
	      };
	      if(no_error) return no_error;
            };
	  };
	};
	return no_error;
}

/*.......................................................................
  Search a file, referenced by its lun, for the string 2nd argument and
  return false or true depending on whether the string is found.
  Leave the file position at the start of the match.
*/
static Template(search_fn)
{
        int match;
/*
  Get the file lun from the user's first argument.
*/
	file_lun = *INTPTR(invals[0]);
/*
  Check that it belongs to an open file with the correct attributes.
*/
	if( (file_ptr=check_lun(file_lun, 1, &is_text)) == NULL)
	  return -1;
/*
  Search the file.
*/
	if((match=file_search(file_ptr, *STRPTR(invals[1]), strlen(*STRPTR(invals[1])), 1)) == -1)
	  return -1;
/*
  Return true if sound, false otherwise.
*/
	*LOGPTR(outvals) = match;
	return no_error;
}

/*.......................................................................
 * Read a text file that contains one or more columns of numbers, into a
 * given 2-D array, resizing the array as needed.
 *
 * Input:
 *  fname     char *  The name of the file.
 *  matrx    float *  The matrix into which to read the numbers.
 */
static Template(read_table_fn)
{
  char *fname;         /* The file name */
  Descriptor *matrx;   /* The return matrix argument */
  MatrixBlocks *mb;    /* A container of a list of blocks of numbers, that */
                       /*  represent the contents of a matrix read from a */
                       /*  file. */
  int nrow,ncol;       /* The numbers of rows and columns in the matrix */
  long dims[3];        /* The dimensions of the return array */
  int nskip = 0;       /* The number of header lines to skip */
/*
 * Get the arguments.
 */
  switch(npar) {
  case 3:
    nskip = *INTPTR(invals[2]);
  case 2:
    fname = *STRPTR(invals[0]);
    matrx = invals[1];
    break;
  default:
    lprintf(stderr, "Unexpected number of arguments.\n");
    return -1;
  };
/*
 * Read the file into a list of dynamically allocated blocks.
 */
  mb = new_MatrixBlocks(fname, nskip);
  if(!mb)
    return -1;
/*
 * Query the size of the matrix.
 */
  mb_matrix_size(mb, &nrow, &ncol);
/*
 * Resize the return matrix to the size of the matrix that was read.
 */
  dims[0] = ncol;
  dims[1] = nrow;
  dims[2] = 1;
  if(re_declare(matrx, dims)) {
    mb = del_MatrixBlocks(mb);
    return -1;
  };
/*
 * Copy the matrix values into the return matrix.
 */
  if(mb_to_float_array(mb, FLTPTR(matrx), nrow*ncol)) {
    mb = del_MatrixBlocks(mb);
    return -1;
  };
/*
 * Clean up.
 */
  mb = del_MatrixBlocks(mb);
/*
 * Keep the user informed.
 */
  lprintf(stderr, "Read a table of %d rows and %d columns from file: %s.\n",
	  nrow, ncol, fname);
  return no_error;
}
