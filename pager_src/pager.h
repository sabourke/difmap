#ifndef pager_h
#define pager_h

#include <stdarg.h>

/* Type of paging required */

typedef enum {
  PAGE_INT,     /* Use internal pager */
  PAGE_EXT,     /* Use external pager */
  PAGE_OFF      /* No paging required - echo direct to stdout */
} Pagetype;

/* Define the pager descriptor - the members are private and should only */
/* be accessed by the following pager method functions. */

enum{PAGE_WIDTH=132};

typedef struct {
  char buffer[PAGE_WIDTH+1]; /* I/O buffer */
  char fname[L_tmpnam];      /* Name of temporary file */
  FILE *fd;                  /* File pointer of output file */
  int line_no;               /* The number of lines in the scratch file */
  int headlen;               /* Number of lines in header page. */
} Pager;

/* Start a new pager */

Pager *new_Pager(void);

/* Add the contents of a file to the pager scratch file */

int page_file(Pager *page, char *name, FILE *fd, int nskip, char *prefix);

/* Write text to the pager scratch file */

int pprintf(Pager *page, const char *format, ...);

/* Mark the end of the optional header page */

void page_mark(Pager *page);

/* Have the pager scratch file paged and delete it and the pager descriptor */

Pager *end_Pager(Pager *page, int dopage, int (*queryfn)(void), Pagetype type);

/* Delete the pager scratch file and the pager descriptor without paging */

Pager *del_Pager(Pager *page);

#endif
