#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "scrfil.h"

/*.......................................................................
 * Given the name for a file, either return a copy of the name if
 * no file currently exists with that name, or return a copy of the file
 * name postfixed with the lowest version number for which no file exists.
 * The returned name is malloc()'d and must subsequently be free()d by the
 * caller.
 *
 * Input:
 *  name    char *   The name of the file.
 * Output:
 *  return  char *   A malloc()d copy of 'name', postfixed with a version
 *                   number if their is a file of the original name.
 *                   Note that the caller is responsible for free()'ing
 *                   this after use. On error, NULL is returned.
 */
char *scrname(const char *name)
{
  const int MAX_VER=999; /* Max postfixed version number */
  size_t slen;       /* The length of the original file name */
  char *file;        /* The full name of the file */
  int ver;           /* File version number */
/*
 * Valid arguments?
 */
  if(name==NULL) {
    fprintf(stderr, "scrname: Intercepted NULL name argument\n");
    return NULL;
  };
/*
 * Allocate sufficient memory to compile a name with up to 4 postfix
 * characters added to the original file name.
 */
  slen = strlen(name);
  file = (char *) malloc(slen+5);
  if(file==NULL) {
    fprintf(stderr,
	    "scrname: Insufficient memory to compile unambiguous file name\n");
    return NULL;
  };
/*
 * Append incrementally higher version numbers until a unique file name
 * is found.
 */
  strcpy(file, name);
  for(ver=0; ver<MAX_VER; ver++) {
    if(ver>0) sprintf(&file[slen], "_%d", ver);
    if(!file_exists(file))
      return file;
  };
/*
 * Max version number exceeded.
 */
  fprintf(stderr, "scrname: Max version number exceeded\n");
  free(file);
  return NULL;
}
