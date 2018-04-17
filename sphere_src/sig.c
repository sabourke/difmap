#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>

#include "sphere.h"
#include "utils.h"
#include "sig.h"

#define MAXLINE 132
/*
  The following flag is set whenever an event is trapped. It tells
  the compile and run time systems to abort and cleanup. It is 0
  normally, but when an error occurs, it is set to -1. This enables
  it to be used instead of (return 0) from all functions. Thus
  Whenever an error occurs the next routine to return will return -1
  and cause a compile or run time abort.
*/
int no_error=0;

/*.......................................................................
  Initialize the siganl handlers for user interrupts and for arithmetic
  exceptions.
*/
void sig_init(void)
{
        signal(SIGINT, interrupt_handler);
        signal(SIGFPE, float_exception);
        return;
}

/*.......................................................................
  This is the handler that gets called whenever there is an interrupt
  request from the user. It asks the user for confirmation, before
  aborting the current compile or execution. If confirmation is not
  recived then this function will return normally and execution will
  continue from where the interrupt was set. The interrupt signal
  will be completely ignored if intercepted while the user is currently
  enterring new commands via stdin.
*/
void interrupt_handler(int sig)
{
  char *reply_str;
/*
  Query the user.
*/
  if(in_run_mode) {
    reply_str = query_user("Abort command (y/n) or quit program [without saving data] (q)? ");
    if(reply_str != NULL) {
      switch(*reply_str) {
      case 'y':
	no_error = -1;          /* Set abort flag */
	break;
      case 'q':
	lprintf(stdout, "Quitting program\n");
        closedown(1, DO_QUIT);  /* Exit program quietly */
	break;
      };
    };
  };
/*
  Make sure that the signal handler is set up correctly for next time.
*/
    sig_init();
    return;
}

/*.......................................................................
  Trap floating point exceptions to prevent them from crashing the prog.
  Report the error and set the interrupt flag.
*/
void float_exception(int sig)
{
/*
  See if there is a system error message for the exception.
*/
        switch(errno) {
        case EDOM:
          lprintf(stderr,"Domain error\n");
          break;
        case ERANGE:
          lprintf(stderr,"Floating point overflow.\n");
          break;
        default:
          lprintf(stderr,"Floating point exception.\n");
          break;
        };
        no_error=-1;
/*
  Make sure that the signal handler is set up correctly for next time.
*/
        sig_init();
        return;
}

/*......................................................................
  Query the user via stdin and return a pointer to the lower case string
  reply. The prompt for the request is passed as the sole argument to
  this function. If an error or end of file is received then NULL will
  be returned.
*/
char *query_user(const char request_str[])
{
        static char reply_str[MAXLINE];
        static char *ctst;
/*
  Repeat the query if the user simply presses return
  without enterring any text.
*/
        lprintf(stderr,"\n");
        for(;;) {
/*
  Output the query to the terminal.
*/
          lprintf(stderr, "%s", request_str);
/*
  Read an input line from stdin for the reply.
*/
          ctst = fgets(reply_str, MAXLINE, stdin);
          if(ctst == NULL) {
            clearerr(stdin);
            break;
          }
          else if(lowstr(reply_str) > 1)
            break;
        };
        return ctst;
}
