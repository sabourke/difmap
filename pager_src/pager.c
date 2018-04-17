#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "pager.h"

/*
 * The following macro places double quotes around the expansion of
 * a macro.
 */
#define QUOTE_VALUE(s) #s
#define QUOTE_MACRO(s) QUOTE_VALUE(s)

/*
 * If no default pager was specified on the compile either set the
 * default to the NULL string (no external pager available) or
 * substitute the pager appropriate to the operating system being used.
 */
#ifdef PAGER
#define PAGER_CMD QUOTE_MACRO(PAGER)
#else
#define PAGER_CMD "more"
#endif

typedef struct {         /* Send pager characteristics between functions */
  int (*queryfn)(void);  /* Internal pager query function */
  Pagetype type;         /* Required pager type */
} Ptype;

static void do_page(Pager *page, Ptype *ptype);
static int p_getline(char *reply, int nmax, FILE *fd);
static int p_query(void);
static int int_pager(Pager *page, Ptype *ptype);
static int ext_pager(Pager *page);

/*.......................................................................
 * Initialize a new text pager and return a dynamically allocated
 * descriptor to be used with calls to the pager method functions.
 * Text lines supplied via the page_file() and pprintf() method
 * functions will be written to a temporary file until end_Pager() is
 * called. At that point the temporary file will be presented to the
 * pager and then deleted. In principle any number of pagers can be
 * in use at once (until memory or file space runs out).
 *
 * Input:
 * Output:
 *  return Pager *        The new pager descriptor, or NULL on error.
 */
Pager *new_Pager(void)
{
  Pager *page=NULL; /* The pointer to the pager descriptor to be returned */
/*
 * Allocate a new pager descriptor.
 */
  page = (Pager *) malloc(sizeof(Pager));
  if(page==NULL) {
    fprintf(stderr, "new_Pager: Insufficient memory to initialize pager\n");
    return page;
  };
/*
 * Initialize pointers to NULL for the sake of del_Pager().
 */
  page->fd = 0;
  page->line_no = 0;
  page->headlen = -1;
/*
 * Open a temporary file to compile the text in.
 */
  page->fd = fopen(tmpnam(page->fname), "w");
  if(page->fd==NULL) {
    fprintf(stderr, "new_Pager: Unable to open scratch file\n");
    return del_Pager(page);
  };
/*
 * Return the initialized container.
 */
  return page;
}

/*.......................................................................
 * If the given Pager descriptor contains a valid pointer to a scratch
 * file, close the file, have its contents paged and then delete the file.
 * Finally release all memory allocated in the descriptor.
 *
 * Input:
 *  page  Pager *         A pager descriptor returned by new_Pager().
 *  dopage  int           If false, don't page the scratch file.
 *  queryfn int (*)(void) Optional pointer to a user query function.
 *                        This will be used if type==PAGE_INT or no
 *                        external pager is available to prompt the
 *                        user at the end of each page. Send NULL to
 *                        select the default function. This function
 *                        should return one of the following:
 *                         0 - Display next page.
 *                         1 - Quit pager.
 *                         2 - Switch to external pager if available.
 *  type     Pagetype     The type of paging wanted.
 *                         PAGE_INT  -  Use internal pager.
 *                         PAGE_EXT  -  Use external pager if defined.
 *                                      Otherwise use internal pager.
 *                         PAGE_OFF  -  No pager - simply echo the text
 *                                      to stdout.
 * Output:
 *  return  Pager *       Always NULL. Use like:
 *                        page = end_Pager(page,1,NULL,PAGE_EXT,-1);
 */
Pager *end_Pager(Pager *page, int dopage, int (*queryfn)(void), Pagetype type)
{
  Ptype ptype;   /* Used to send pager characteristics between functions */
/*
 * No container?
 */
  if(page==NULL)
    return page;
/*
 * Collect pager characteristics into container.
 */
  ptype.type = type;
  ptype.queryfn = queryfn ? queryfn : p_query;
/*
 * Is there a scratch file to be paged?
 */
  if(page->fd==NULL) {
    fprintf(stderr, "end_Pager: No scratch file to page.\n");
  } else {
/*
 * Close the file, have it paged if it contains anything and delete it.
 */
    if(fclose(page->fd)==EOF) {
      if(dopage)
	fprintf(stderr, "end_Pager: Error closing scratch file\n");
    } else if(dopage && page->line_no>0) {
      do_page(page, &ptype);
    };
/*
 * del_Pager() will try to close the file and delete it again if fd!=NULL.
 */
    page->fd = NULL;
/*
 * Delete the scratch file.
 */
    if(remove(page->fname)) {
      fprintf(stderr, "end_Pager: Error deleting scratch file: %s\n",
	      page->fname);
    };
  };
/*
 * Delete the container.
 */
  return del_Pager(page);
}

/*.......................................................................
 * Delete the pager scratch file and pager descriptor.
 *
 * Input:
 *  page   Pager *   The descriptor of the pager.
 * Output:
 *  return Pager *   Allways NULL. Use like:  page=del_Pager(page);
 */
Pager *del_Pager(Pager *page)
{
/*
 * Container already deleted?
 */
  if(page==NULL)
    return page;
/*
 * Is there a scratch file to be deleted?
 */
  if(page->fd) {
/*
 * Close the file - ignore errors since the state of the contents no
 * longer concern us.
 */
    fclose(page->fd);
/*
 * Delete it.
 */
    if(remove(page->fname)) {
      fprintf(stderr, "del_Pager: Error deleting scratch file: %s\n",
	      page->fname);
    };
  };
/*
 * Delete the container.
 */
  free(page);
  return NULL;
}

/*.......................................................................
 * Private function for the use of end_Pager(). This function takes
 * the already closed scratch file named in page->fname[] and has it
 * paged.
 *
 * Input:
 *  page  Pager *  The descriptor of the pager.
 *  ptype Ptype *  Pager characteristics.
 */
static void do_page(Pager *page, Ptype *ptype)
{
/*
 * Run the internal pager if reuested, or if the external pager returns
 * with an error.
 */
  if(ptype->type!=PAGE_EXT || (ptype->type==PAGE_EXT && ext_pager(page)))
    int_pager(page, ptype);
  return;
}

/*.......................................................................
 * Run the internal pager on the closed scratch file named in page->fname.
 *
 * Input:
 *  page    Pager *  The pager descriptor.
 *  ptype   Ptype *  Pager characteristics.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int int_pager(Pager *page, Ptype *ptype)
{
  FILE *scrfd;       /* Pointer to open scratch file */
  int reply=0;       /* Reply from user query function */
  char *cptr;        /* Pointer to reply from getenv() */
  int nlines=24;     /* The number of lines per page */
  int ncolumns=80;   /* The width of the page */
  int first=1;       /* True only for first page */
/*
 * Open the scratch file for read.
 */
  scrfd = fopen(page->fname, "r");
  if(scrfd==NULL) {
    fprintf(stderr, "int_pager: Error opening scratch file\n");
    return 1;
  };
/*
 * See if the environment has recorded the number of lines available in
 * the LINES environment variable. (This is normally the case on
 * SYSV machines).
 */
  cptr = getenv("LINES");
  if(cptr)
    nlines = atoi(cptr);
/*
 * We require at least 3 lines - two for the query and prompt and one
 * for a single line of text.
 */
  if(nlines < 3)
    nlines = 3;
/*
 * See if the environment has recorded the number of columns available in
 * the COLUMNS environment variable. (This is normally the case on
 * SYSV machines).
 */
  cptr = getenv("COLUMNS");
  if(cptr)
    ncolumns = atoi(cptr);
/*
 * Non-sensical value?
 */
  if(ncolumns < 2)
    ncolumns = 80;
/*
 * Write pages of 'nlines-2' lines, and prompt at the end of each.
 */
  do {
    int lnum=0;   /* Line number within current page */
/*
 * On the first page, advance the line count to within headlen of the
 * end of the page so that the pager prompt will be presented after
 * just headlen lines.
 */
    if(first && page->headlen>=0 && page->headlen <= nlines-2) {
      first = 0;
      lnum = nlines - 2 - page->headlen;
    };
/*
 * Print up to nlines-2 lines.
 */
    for(; lnum<nlines-2; lnum++) {
/*
 * Read the next line from the scratch file and echo it to stdout.
 */
      if(p_getline(page->buffer, PAGE_WIDTH, scrfd))
	break;
/*
 * Assume that the screen is 'ncolumns' characters across and truncate the
 * output line if it will not fit in this space.
 */
      fprintf(stdout, "%-.*s\n", ncolumns-1, page->buffer);
    };
  } while(!feof(scrfd) && !ferror(scrfd) &&
	  (ptype->type==PAGE_OFF || (reply=ptype->queryfn())==0));
/*
 * Close the scratch file.
 */
  fclose(scrfd);
/*
 * Run external pager?
 */
  return reply==2 ? ext_pager(page) : 0;
}

/*.......................................................................
 * Run the external pager on the closed scratch file named in page->fname.
 *
 * Input:
 *  page    Pager *  The pager descriptor.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error.
 */
static int ext_pager(Pager *page)
{
  char *command;        /* The pager command. */
/*
 * See if the user has set the PAGER environment variable.
 */
  command = getenv("PAGER");
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
 * No valid PAGER variable? Then substitute the default pager if there is one.
 */
  if(command==NULL)
    command = PAGER_CMD;
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
 * If an external pager has been selected, compose a shell command to have
 * the scratch file paged, then have it executed.
 */
  if(command) {
/*
 * Is there insufficient room in page->buffer to compose the full
 * command in?
 */
    if(strlen(command) + 1 + strlen(page->fname) + 1 > PAGE_WIDTH) {
      fprintf(stderr, "Pager command too long for buffer.\n");
      return 1;
    } else {
/*
 * Compose the full pager command.
 */
      sprintf(page->buffer, "%s %s", command, page->fname);
/*
 * Have it executed.
 */
      system(page->buffer);
    };
/*
 * External pager requested, but non available?
 */
  } else {
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Take either the pointer to an open file, or the name of a file to
 * read from. Read and discard the first nskip lines, then prepend the
 * supplied prefix to each line of the remaining lines and echo them
 * to into the pager scratch file.
 *
 * Input:
 *  page    Pager *  A pager descriptor from new_Pager().
 *  name     char *  The name of the file to page, or NULL if the fd
 *                   argument is to be used istead.
 *  fd       FILE *  If name==NULL, then this must be a non-NULL
 *                   pointer to a file opened for reading. Otherwise
 *                   send NULL.
 *  nskip     int    The number of lines to skip in the input file -
 *                   send 0 if non are to be skipped.
 *  prefix   char *  The string to prepend to each line in the scratch
 *                   file. Send "" or NULL if not required.
 * Output:
 *  return    int    0 - OK.
 *                   1 - Error (don't forget to call del_Pager() after
 *                       this happens).
 */
int page_file(Pager *page, char *name, FILE *fd, int nskip, char *prefix)
{
  int waserr=0;   /* Error status */
/*
 * Check arguments.
 */
  if(page==NULL) {
    fprintf(stderr, "page_file: NULL Pager descriptor received.\n");
    return 1;
  };
/*
 * At least one file specification must have been provided.
 */
  if(name==NULL && fd==NULL) {
    fprintf(stderr, "page_file: No file specification provided.\n");
    return 1;
  };
/*
 * Previous error writing the scratch file?
 */
  if(ferror(page->fd))
    return 1;
/*
 * Open the named file for read if given.
 */
  if(name) {
    fd = fopen(name, "r");
    if(fd==NULL) {
      fprintf(stderr, "page_file: Error opening file: %s\n", name);
      return 1;
    };
  };
/*
 * Skip the first nskip lines of the input file.
 */
  if(nskip > 0) {
    int nlines;
    for(nlines=0; nlines<nskip; nlines++) {
      if(p_getline(page->buffer, PAGE_WIDTH, fd))
	break;
    };
  };
/*
 * If there are any remaining lines in the input file, echo them to
 * the scratch file.
 */
  while(!waserr && !feof(fd) && !ferror(fd) && p_getline(page->buffer, PAGE_WIDTH, fd)==0) {
    waserr = pprintf(page, "%s%s\n", prefix ? prefix : "", page->buffer) < 0;
  };
/*
 * Error reading file?
 */
  if(ferror(fd) && !feof(fd)) {
    fprintf(stderr, "page_file: Error reading file.\n");
    waserr = 1;
  };
/*
 * Close the file if it was opened in this function.
 */
  if(name)
    fclose(fd);
  return waserr;
}

/*.......................................................................
 * Write text to the pager scratch file.
 *
 * Input:
 *  page    Pager *  A pager descriptor from new_Pager().
 *  format   char *  A printf() style format string. Newline
 *                   characters will not be inserted for you. For
 *                   the purposes of line counting, it will be assumed
 *                   that all newline ('\n') characters will be in the
 *                   format string.
 *  ...              Variable argument list matching format items in
 *                   'format'.
 * Output:
 *  return    int    The number of characters printed, or -ve on error.
 */
int pprintf(Pager *page, const char *format, ...)
{
  va_list args;   /* Variable argument list */
  char *cptr;     /* Pointer into format[] */
  int iret=0;     /* Return value */
/*
 * Check arguments.
 */
  if(page==NULL) {
    fprintf(stderr, "pprintf: NULL Pager descriptor received\n");
    return -1;
  };
/*
 * Nothing to write to the file?
 */
  if(format==NULL)
    return iret;
/*
 * File no longer valid?
 */
  if(ferror(page->fd))
    return -1;
/*
 * Count newline characters in the format string.
 */
  for(cptr=(char *)format; *cptr!='\0'; cptr++)
    if(*cptr=='\n') page->line_no++;
/*
 * Write the text line to the scratch file.
 */
  va_start(args, format);
  iret = vfprintf(page->fd, format, args);
  va_end(args);
/*
 * Error writing file?
 */
  if(iret < 0)
    fprintf(stderr, "Error writing to pager scratch file: %s\n", page->fname);
  return iret;
}

/*.......................................................................
 * Read a line of text from a text file pointed to by fd. If any '\n'
 * is found in the read string, it is replaced with '\0'. If no '\n'
 * is found then the line has clearly been truncated, and all
 * characters up to and including the next '\n' or EOF are consumed
 * such that the next read will start at the start of a new line. The
 * text read is left in page->buffer[].
 *
 * Input/Output:
 *  reply   char *  A buffer for at least nmax characters. This will
 *                  contain the line that was read except after an error.
 * Input:
 *  nmax     int    The number of elements in reply[].
 *  fd      FILE *  The pointer to a file opened for read.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int p_getline(char *reply, int nmax, FILE *fd)
{
  char *cptr;   /* Pointer into reply[] */
/*
 * Read as much of the next line as will fit in the buffer.
 */
  if(fgets(reply, nmax, fd) == NULL)
    return 1;
/*
 * Check for '\n'.
 */
  cptr = strchr(reply, '\n');
/*
 * If no '\n' was found, then the line has been truncated. Read to the
 * start of the next line so that the next fgets() starts at the right
 * place.
 */
  if(cptr==NULL) {
    int c;
    while((c=getc(fd)) != '\n' && c != EOF);
  }
/*
 * Remove trailing '\n' if present.
 */
  else {
    *cptr = '\0';
  };
  return 0;
}

/*.......................................................................
 * Default user query function to be used if the user doesn't supply one.
 *
 * Output:
 *  return     int   0 - Continue paging.
 *                   1 - Stop paging.
 */
static int p_query(void)
{
  enum {MAXREPLY=2};
  char reply[MAXREPLY+1];  /* Buffer to receive reply string */
/*
 * Prompt the user.
 */
  fprintf(stdout, "Press return for the next page, Q to quit, or P for external pager.\n#");
/*
 * Read user's selection.
 */
  if(p_getline(reply, MAXREPLY, stdin)==0) {
    char *cptr=reply;
/*
 * Skip white space.
 */
    while(isspace((int)*cptr))
      cptr++;
/*
 * Nothing typed?
 */
    if(*cptr=='\0')
      return 0;     /* Keep paging */
/*
 * External pager requested?
 */
    if(tolower(*cptr)=='p')
      return 2;
  };
  return 1;         /* Stop paging */
}

/*.......................................................................
 * With the internal pager this sets a premature end for the first
 * displayed page, thus giving the opportunity via the prompt to
 * switch to an external pager after displaying only a minimal header.
 *
 * Input:
 *  page    Pager *   The pager descriptor.
 */
void page_mark(Pager *page)
{
  page->headlen = page->line_no;
}
