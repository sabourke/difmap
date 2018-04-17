#include <stdio.h>
#include <stdlib.h>

#include "logio.h"
#include "sphere.h"

/* Define the node of a list of functions */

typedef struct Fnnode {
  EXITFN(*fn);            /* The function to be called */
  struct Fnnode *next;    /* Pointer to next function node in the list */
} Fnnode;

/* The head of the list of functions to call at exit from the program */

static Fnnode *head=NULL;

/*.......................................................................
 * Add a node to the list of functions to be called at exit from the
 * program.
 *
 * Input:
 *  fn     EXITFN(*fn)   Pointer to the function to be installed.
 * Output:
 *  return         int    0 - OK.
 *                       -1 - Error.
 */
int add_exit_fn(EXITFN(*fn))
{
  Fnnode *node=NULL;  /* The new function node */
/*
 * Sanity check.
 */
  if(fn==NULL) {
    lprintf(stderr, "add_exit_fn: NULL function pointer intercepted\n");
    return -1;
  };
/*
 * Allocate memory for the new function node.
 */
  node = (Fnnode *) malloc(sizeof(Fnnode));
  if(node==NULL) {
    lprintf(stderr,
	    "add_exit_fn: Insufficient memory to register new function\n");
    return -1;
  };
/*
 * Register the function.
 */
  node->fn = fn;
/*
 * Install the new node at the head of the list. This orderring ensures
 * that functions are called in the reverse order to their registration.
 */
  node->next = head==NULL ? NULL : head;
  head = node;
  return 0;
}

/*.......................................................................
 * This function should be called to exit the program. Before calling
 * exit, it calls all the functions in the exit-function list in the
 * reverse order to their registration.
 *
 * Input:
 *  status      int   This value will be passed on to exit().
 *  code   Exitcode   This will be passed on to all cleanup functions:
 *                     DO_EXIT  -  Do cleanup and any optional closedown
 *                                 tasks.
 *                     DO_QUIT  -  Do minimal cleanup.
 */
int closedown(int status, Exitcode code)
{
  Fnnode *node;   /* A node in the list of functions to be called */
  Fnnode *next;   /* The node following 'node' */
/*
 * Call cleanup functions.
 */
  for(node=head; node!=NULL; node=next) {
    next = node->next;  /* Record the next node so that node can be free'd */
    node->fn(code);     /* Call the cleanup function */
    free(node);         /* Delete the list node */
  };
  exit(status);
}
