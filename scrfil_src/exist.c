#include <stdio.h>

#include "scrfil.h"

/*.......................................................................
 * Return 1 if the given named file exists and is readable, 0 otherwise.
 *
 * Input:
 *  name   char *  The name of the file to check for.
 * Output:
 *  return  int    0 - File doesn't exist.
 *                 1 - File does exist.
 */
int file_exists(const char *name)
{
  FILE *fp = fopen(name, "r");
  if(fp)
    fclose(fp);
  return fp!=NULL;
}
