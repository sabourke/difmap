#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "scrfil.h"

/*
 * If no default editor was specified at compile time, either set the
 * default to the NULL string (no external editor available) or
 * substitute the editor appropriate to the operating system being used.
 */
#ifndef EDITOR
#ifdef unix
#define EDITOR "vi"
#else
#define EDITOR NULL
#endif
#endif

/*.......................................................................
 * Edit a named file.
 *
 * Input:
 *  name   char *  The name of the file to be edited. If this is NULL,
 *                 "" will be substituted.
 * Output:
 *  return  int    0 - OK.
 *                 1 - Error.
 */
int ed_file(const char *name)
{
  char *command;           /* The command name of the editor. */
/*
 * If no file name argument was provided, substitute an empty string.
 */
  if(name==NULL)
    name = "";
/*
 * See if the user has set the EDITOR environment variable.
 */
  command = getenv("EDITOR");
/*
 * Skip white-space and check for the empty string?
 */
  if(command) {
    while(isspace((int)*command))
      command++;
    if(*command == '\0')
      command = NULL;
  };
/*
 * No valid EDITOR variable? Then substitute the default editor if there is one.
 */
  if(command==NULL)
    command = EDITOR;
/*
 * Skip white-space and check for the empty string?
 */
  if(command) {
    while(isspace((int)*command))
      command++;
    if(*command == '\0')
      command = NULL;
  };
/*
 * If an external editor has been selected, compose a shell command to have
 * the scratch file edited, then have it executed.
 */
  if(command) {
/*
 * How long must the edit command line be?
 */
    size_t slen = strlen(command) + 1 + strlen(name) + 1;
/*
 * Allocate a string long enough to contain the editor invokation line.
 */
    char *edstr = (char *) malloc(sizeof(char) * slen);
    if(edstr==NULL) {
      fprintf(stderr,"ed_file: Insufficient memory to compose command line.\n");
      return 1;
    };
/*
 * Compose the command line.
 */
    sprintf(edstr, "%s %s", command, name);
/*
 * Have it executed.
 */
    system(edstr);
/*
 * Release the memory used for the command line.
 */
    free(edstr);
/*
 * No editor command?
 */
  } else {
    fprintf(stderr,
"ed_file: Default editor unknown. Before starting this program, create an\n");
    fprintf(stderr,
"ed_file: EDITOR environment variable containing your editor's command name.\n");
    return 1;
  };
/*
 * File edit succesfully completed.
 */
  return 0;
}
