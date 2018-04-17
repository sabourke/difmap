#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "sphere.h"
#include "table.h"
#include "help.h"
#include "lex.h"
#include "logio.h"
#include "utils.h"
#include "pager.h"


/* Declare a container type for help file info */

typedef struct {
  enum {HLP_LEN=160} len;
  char args[HLP_LEN+1];    /* Arguments string from first line of help file */
  char intro[HLP_LEN+1];   /* One line intro from second line of help file */
  char *help_name;         /* Full name of help file */
  FILE *fp;                /* Pointer to help file */
} Helpfile;

static char *fgetl(char *reply, int nmax, FILE *fd);

static Helpfile *open_help(const char *help_dir, const char *topic);
static Helpfile *close_help(Helpfile *hfile);
static int whatisfunc(char *name, Functype *fdsc, Helpfile *hfile, Pager *page);

/* Declare identifiers for each type of index topic entry */

typedef enum {IDX_NONE, IDX_GENERAL, IDX_COMMAND} Topictype; 

/* Declare an index file descriptor */

typedef struct {
  enum {IDX_LEN=80} len;
  char *general_title;     /* Title used to introduce general help topics */
  char *command_title;     /* Title used to introduce command help topics */
  char *nohelp;            /* Intro string used when there is no help file */
  char topic[IDX_LEN+1];   /* The latest topic name read from the index */
  char intro[IDX_LEN+1];   /* The 1-line intro associated with 'topic' */
  Topictype type;          /* The type of topic */
  char *index_name;        /* Full name of index file */
  FILE *fp;                /* File pointer of opened index file */
} Indexfile;

static Indexfile *open_index(const char *help_dir, const char *module,
			     const char *mode);
static int read_index(Indexfile *ifile);
static Indexfile *close_index(Indexfile *ifile);
static void low_strcpy(char *orig, char *dest);

/*.......................................................................
 * Translate the enumerated storage type into a human readable type name.
 * A pointer to a static string is returned.
 *
 * Input:
 *  type   char   Variable type enumerator.
 * Output:
 *  return char * Pointer to static string containing the associated
 *                name for the type, or "(unknown)" if not recognised.
 */
char *type_string(char type)
{
  switch (type) {
  case 'c':
    return "string";
  case 'f':
    return "float";
  case 'i':
    return "integer";
  case 'l':
    return "logical";
  case 'n':
    return "number";
  case '*':
    return "(any type)";
  case 'C':
    return "literal";
  case ' ':
    return "void";
  };
  return "(unknown)";
}

/*.......................................................................
 * Translates an array dimension, into a descriptive string.
 *
 * Input:
 *  dim     char   The integral dimension to be named (scalars are 0,
 *                 1D array is 1 etc..).
 * Output:
 *  return  char * Pointer to a static string containing a description
 *                 of the dimension, or "(unknown)" if not recognised.
 */
char *dims_string(char dim)
{
  switch (dim) {
  case '0': case ' ':
    return "SCALAR";
    break;
  case '1':
    return "1D ARRAY";
    break;
  case '2':
    return "2D ARRAY";
    break;
  case '3':
    return "3D ARRAY";
    break;
  };
  return "(unknown)";
}

/*.......................................................................
 * Return a text description of an enumerated expression access type.
 *
 * Input:
 *  access   char    The enumerated access permission.
 * Output:
 *  return   char *  Pointer to a static string describing the access type,
 *                   or "(unknown)" if not recognised.
 */
char *access_string(char access)
{
  switch (access) {
  case 'N':
    return "ARRAY NAME";
  case 'r':
    return "POINTER";
  case 'v': case 'V': case '?': case ' ':
    return "VALUE";
  };
  return "(unknown)";
}

/*.......................................................................
 * Translates the declaration of the function with descriptor, fdsc,
 * and the two header lines of the corresponding help file, into a
 * human readable form. This is then echoed to the user's terminal.
 *
 * Input:
 *  name      char *  The name of the function.
 *  fdsc  Functype *  The function descriptor.
 *  hfile Helpfile *  The help file descriptor.
 *  page     Pager *  The pager descriptor.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
static int whatisfunc(char *name, Functype *fdsc, Helpfile *hfile, Pager *page)
{
  char work_line[HLP_LEN+1]; /* Work array */
  char *type; /* Pointer to type string. */
  int wlen;   /* Number of characters in work_line[] */
  int ierr=0; /* True after an error */
/*
 * Echo the full name of the function.
 */
  sprintf(work_line, "%.*s", HLP_LEN, name);
  wlen = strlen(work_line);
/*
 * If it is a function then follow it with an open bracket otherwise
 * follow it with a space as it is a command.
 */
  if(wlen < HLP_LEN)
    work_line[wlen++] = (*fdsc->type == ' ') ? ' ' : '(';
/*
 * The first line in the help file holds textual names for the arguments.
 */
  strncpy(&work_line[wlen],
	  hfile ? hfile->args : "(Help file not available)", HLP_LEN-wlen);
  work_line[HLP_LEN] = '\0';
  wlen = strlen(work_line);
/*
 * Close function argument delimiter brackets.
 */
  if(wlen < HLP_LEN && *fdsc->type!=' ')
    work_line[wlen++] = ')';
/*
 * Make sure that the string is terminated.
 */
  work_line[wlen] = '\0';
/*
 * First line complete - echo it.
 */
  ierr = ierr || pprintf(page, "%s\n", work_line) < 0;
/*
 * Next line - echo the 1 line description of the function given as the
 * second line of the help file.
 */
  strncpy(work_line, hfile ? hfile->intro : "(Help file not available)",
	  HLP_LEN);
  work_line[HLP_LEN] = '\0';
/*
 * Second line complete - echo it.
 */
  ierr = ierr || pprintf(page, "  %s\n", work_line) < 0;
/*
 * Reserved word command?
 */
  if(fdsc->sub_class != NORM) {
    ierr = ierr || pprintf(page, "  A special command.\n") < 0;
  } else {
/*
  Command or function?
*/
    if(*fdsc->type == ' ') {
      strcpy(work_line, "A command.");
    } else {
/*
 * Report the function type.
 */
      if(fdsc->once == 1)
	strncpy(work_line, "A non-elemental", HLP_LEN);
      else
	strncpy(work_line, "An elemental", HLP_LEN);
      wlen = strlen(work_line);
      if(fdsc->access[0] == '?') {
	strncpy(&work_line[wlen], " command or", HLP_LEN-wlen);
	wlen = strlen(work_line);
      };
      strncpy(&work_line[wlen], " function: returns ", HLP_LEN-wlen);
      wlen = strlen(work_line);
/*
  Determine the return type.
  And the return dimensional type.
*/
      type = type_string(*fdsc->type);
      if(*fdsc->dim != '0') {
	if(HLP_LEN-wlen > 10 + strlen(type))
	  sprintf(&work_line[wlen], "%cD %s array.", *fdsc->dim, type);
      } else {
	if(HLP_LEN-wlen > 7 + strlen(type))
	  sprintf(&work_line[wlen], "%s scalar.", type);
      };
    };
/*
 * Third line complete.
 */
    ierr = ierr || pprintf(page, "  %s\n", work_line) < 0;
/*
 * No arguments?
 */
    if(fdsc->nmax == 0) {
      ierr = ierr || pprintf(page, "     No arguments.\n") < 0;
    } else {
/*
 * The rest of the lines display the declaration of each argument, one per line.
 */
      char *args = fdsc->type;       /* Pointer to argument declarations */
      int narg = strlen(args)-1;     /* Number of declarations in args[] */
      int arg;                       /* Argument number */
/*
 * Loop for each declared argument.
 */
      for(arg=1; arg<=narg; arg++) {
/*
 * Introduce the the start of the argument list and the start of the
 * optional argument list.
 */
	if(arg-1==fdsc->nmin)
	  ierr = ierr || pprintf(page, "Optional args:\n") < 0;
	else if(arg==1)
	  ierr = ierr || pprintf(page, "Args:\n") < 0;
/*
 * Report the dimensional type unless scalar.
 */
	wlen = 0;
	if(fdsc->dim[arg] != '0') {
	  sprintf(work_line, "%cD ", fdsc->dim[arg]);
	  wlen = strlen(work_line);
	};
/*
 * Get the argument storage type.
 */
	type = type_string(fdsc->type[arg]);
/*
 * Report the storage type and its access class.
 */
	switch (fdsc->access[arg]) {
	case 'v':
	  sprintf(&work_line[wlen], "%s value",type);
	  break;
	case 'r':
	  sprintf(&work_line[wlen], "%s reference",type);
	  break;
	case 'N':
	  sprintf(&work_line[wlen], "%s variable_name",type);
	  break;
	};
/*
 * Line complete.
 */
	ierr = ierr || pprintf(page, "  %s\n", work_line) < 0;
      };
/*
 * If the declaration string is exhausted before the max number of arguments
 * has been reached, write ellipses.
 */
      if(arg < fdsc->nmax)
	ierr = ierr || pprintf(page, "  ...\n") < 0;
    };
  };
/*
 * Finally add a blank line.
 */
  ierr = ierr || pprintf(page, "\n") < 0;
/*
 * Mark this point as the end of the pager header.
 */
  page_mark(page);
  return ierr ? 1 : 0;
}

/*.......................................................................
 * List the declarations of variables that ambiguously, or exactly match
 * a give identifier name, or list all declarations if no name is supplied.
 *
 * Input:
 *  name    char *  The name to be matched with variable names. Send
 *                  NULL to have all variables listed.
 */
void whatisvar(char *name)
{
  char match_typ;   /* The type of match as returned by find_symbol() */
  int bot, top;     /* Indexes of first and last matches in symbol table */
  int count=0;      /* The number of matching variable names */
  int i,j;
/*
 * Search the main symbol table for the symbol name sent.
 */
  if(name != NULL)
    match_typ = find_symbol(name, main_table, num_main, &bot, &top);
  else {
    match_typ = 'a';   /* Pretend that everything matched */
    bot = 0;
    top = num_main-1;
  };
/*
 * No match?
 */
  if(match_typ == 'n') {
    lprintf(stderr, "No function or variable symbol matches %s\n",name);
    return;
  };
/*
 * Loop for each match.
 */
  for(i=bot; i<=top; i++) {
    Table *symbol = main_table[i];
/*
 * Single out variable names.
 */
    if(symbol->class == VAR) {
      Descriptor *vdsc = TABDESC(symbol);
/*
 * Record the number of matching variables.
 */
      count++;
/*
 * Report the access class.
 */
      switch (vdsc->access) {
      case R_ONLY:
	fprintf(stdout,"r--");
	break;
      case NO_DEL:
	fprintf(stdout,"rw-");
	break;
      case RWD:
	fprintf(stdout,"rwd");
	break;
      };
/*
 * Report the variable type.
 */
      fprintf(stdout, "  %s %s", type_string(vdsc->atyp), symbol->name);
/*
 * Now handle the dimensional type.
 */
      if(vdsc->dim != '0') {
	int dim = vdsc->dim - '0';
	fputc('(',stdout);
	for(j=0; j<dim; j++)
	  fprintf(stdout, "%ld%c", vdsc->adim[j], (j==dim-1) ? ')':',');
      };
      fputc('\n',stdout);
    };
  };
/*
 * If no variables matched the name then report this.
 */
  if(count == 0)
    lprintf(stderr, "No variable names match %s\n", name);
  return;
}

/*.......................................................................
 * Top-level help function. This function takes the name that it is sent
 * and matches it with symbols in the help and main symbol tables.
 * If an exact match is found, then hep information is displayed for that
 * symbol.
 *
 * Input:
 *  name    char *  The name of the symbol to look up.
 */
void help(char *name)
{
  int bot, top;        /* Indexes of first and last matching symbols in table */
  char match_typ;      /* Type of match in symbol table */
  Table *symbol;       /* Pointer to table entry matching 'name' */
/*
 * If no name was provided then list all the function module names.
 */
  if(name == NULL) {
    list_modules();
    return;
  };
/*
 * Search for the name in the symbol table.
 */
  match_typ = find_symbol(name, main_table, num_main, &bot, &top);
/*
 * Given no match, report this fact then list the special help symbols.
 */
  switch (match_typ) {
  case 'n':
    lprintf(stderr, "help: no such symbol '%s'\n",name);
    list_modules();
    return;
    break;
  case 'a':                   /* Ambiguous match */
    lex_err(comline.last);
    list_matches(bot,top,name);
    return;
    break;
  };
/*
 * Act upon the symbol type found.
 */
  symbol = main_table[bot];
  switch (symbol->class) {
  case VAR:
    whatisvar(name);
    return;
    break;
  case FUNC:
    help_function(symbol);
    break;
  case HELP_SYM:
    help_topic(symbol);
    break;
/*
 * A module name has been found - list all the functions in that module.
 */
  case MODULE_SYM:
    help_module(symbol);
    break;
  };
  return;
}

/*.......................................................................
 * List help information on a given function symbol.
 *
 * Input:
 *  symbol   Table *  The symbol table entry of the function.
 */
void help_function(Table *symbol)
{
  Helpfile *hfile;   /* Help file descriptor */
  Pager *page;       /* File pager */
  int waserr=0;      /* True after errors */
/*
 * Not a function?
 */
  if(symbol->class != FUNC) {
    lprintf(stderr, "help_function: Received non-function type\n");
    return;
  };
/*
 * Attempt to open the help file for the matched symbol.
 */
  hfile = open_help(TABSTR(TABFUNC(symbol)->help), symbol->name);
/*
 * Start a pager.
 */
  page = new_Pager();
  if(page) {
/*
 * Translate the function declaration and help-file header lines
 * into human readable form and echo these to the user's terminal.
 */
    whatisfunc(symbol->name, TABFUNC(symbol), hfile, page);
/*
 * Add the rest of the help file to the pager.
 */
    waserr = hfile && page_file(page, NULL, hfile->fp, 0, " ")!=0;
/*
 * End the pager.
 */
    page = end_Pager(page, !waserr, pause_output, PAGE_INT);
    fprintf(stdout, "Help listing for %s completed\n", symbol->name);
  };
/*
 * Close the help file if open.
 */
  close_help(hfile);
  return;
}

/*.......................................................................
 * Show the help file associated with a given help topic.
 *
 * Input:
 *  symbol   Table *  The symbol table entry of the help topic.
 */
void help_topic(Table *symbol)
{
  Helpfile *hfile;   /* Help file descriptor */
/*
 * Not a help topic?
 */
  if(symbol->class != HELP_SYM) {
    lprintf(stderr, "help_topic: Symbol \'%s\' is not a help topic.\n",
	    symbol->name);
    return;
  };
/*
 * Attempt to open the help file for the matched symbol.
 * The value associated with a help topic symbol is a pointer to the
 * symbol table entry of the module with which it is associated.
 * This in turn holds the help directory as its value.
 */
  hfile = open_help(TABSTR(TABTAB(symbol)), symbol->name);
/*
 * See if the user wants any further help information from the help file.
 */
  if(hfile) {
/*
 * Run a pager on the help file.
 */
    Pager *page = new_Pager();
    if(page) {
      int dopage = page_file(page, NULL, hfile->fp, 0, " ")==0;
      page = end_Pager(page, dopage, pause_output, PAGE_INT);
    };
    fprintf(stdout, "Help listing for %s completed\n", symbol->name);
/*
 * Close the help file.
 */
    close_help(hfile);
  };
  return;
}

/*.......................................................................
 * List functions from a given module, along with their one-line
 * descriptions.
 *
 * Input:
 *  symbol   Table *  The symbol table entry of the module.
 */
void help_module(Table *symbol)
{
  Pager *page;     /* The descriptor of the file pager */
  Helpfile *hfile; /* Help file descriptor */
/*
 * Check that the symbol sent is a module symbol.
 */
  if(symbol->class != MODULE_SYM) {
    lprintf(stderr, "help_module: Symbol received is not a module symbol.\n");
    return;
  };
/*
 * Start a file pager.
 */
  page = new_Pager();
  if(page) {
    int endit=0;           /* Used to signal when to stop paging */
/*
 * See if there is a help file for this module.
 */
    hfile = open_help(TABSTR(symbol), symbol->name);
    if(hfile) {
      endit = endit || page_file(page, NULL, hfile->fp, 0, "");
      close_help(hfile);
    } else {
/*
 * Supply a default header if no module help file was found.
 */
      endit = endit || pprintf(page, "\nHelp available for module: %s\n",
			       symbol->name) < 0;
    };
/*
 * Also list the index file for this module.
 */
    if(!endit) {
      Indexfile *ifile = open_index(TABSTR(symbol), symbol->name, "r");
      if(ifile) {
	endit = endit || page_file(page, NULL, ifile->fp, 0, "");
	close_index(ifile);
      };
    };
/*
 * Have the pager display the listing.
 */
    page = end_Pager(page, !endit, pause_output, PAGE_INT);
  };
  fprintf(stdout, "Listing completed.\n");
  return;
}

/*.......................................................................
 * List the names of all modules cited in the main symbol table.
 */
void list_modules(void)
{
  int i;
  fprintf(stdout, "List of function modules.\n");
/*
 * Search through the main symbol table for help symbols and list them.
 */
  for(i=0; i<num_main; i++) {
    if(main_table[i]->class == MODULE_SYM)
      fprintf(stdout, "\t%s\n", main_table[i]->name);
  };
  fprintf(stdout, "For more help, type HELP module_name or HELP function_name.\n");
  return;
}

/*.......................................................................
 * Look for the help file corresponding to a given function table entry,
 * open it and return its file pointer.
 *
 * Input:
 *  help_dir   char *  The directory path to look for the help file under.
 *                     This must be a full pathname (in UNIX this must
 *                     include the trailing /).
 *  topic      char *  The topic to file help on.
 * Output:
 *  return  Helpfile * Pointer to a dynamically allocated help file
 *                     container, or NULL on error. This contains the
 *                     initial two line header and the pointer to the
 *                     help file, opened for read at the first line
 *                     that follows the header.
 */
static Helpfile *open_help(const char *help_dir, const char *topic)
{
  static char ext[] = ".hlp";   /* Help file extension */
  Helpfile *hfile=NULL;         /* Help file descriptor */
/*
 * Allocate the help file container.
 */
  hfile = (Helpfile *) malloc(sizeof(Helpfile));
  if(hfile==NULL) {
    lprintf(stderr, "open_help: Insuffient memory for descriptor.\n");
    return hfile;
  };
/*
 * Intialize the contents of the container.
 */
  hfile->len = HLP_LEN;
  hfile->intro[0] = hfile->args[0] = '\0';
  hfile->help_name = NULL;
  hfile->fp = NULL;
/*
 * Compile the help file name.
 */
  hfile->help_name = (char *) malloc(strlen(help_dir) + strlen(topic) +
				     strlen(ext) + 1);
  if(hfile->help_name==NULL) {
    lprintf(stderr, "open_help: Insufficient memory to compile file name.\n");
    return close_help(hfile);
  };
  sprintf(hfile->help_name, "%s%s%s", help_dir, topic, ext);
/*
 * Attempt to open the file.
 */
  hfile->fp = fopen(hfile->help_name, "r");
  if(hfile->fp) {
    if(fgetl(hfile->args, hfile->len, hfile->fp) == NULL)
      strcpy(hfile->args, "(Can't read help file)");
    if(fgetl(hfile->intro, hfile->len, hfile->fp) == NULL)
      strcpy(hfile->intro, "(Can't read help file)");
  } else {
    return close_help(hfile);
  };
/*
 * Return the new container.
 */
  return hfile;
}

/*.......................................................................
 * Close a help file if open and delete its help file container.
 *
 * Input:
 *  hfile  Helpfile *  The static Helpfile descriptor returned by
 *                     open_help().
 * Output:
 *  return Helpfile *  Allways NULL. Use like  hfile=close_help(hfile);
 */
static Helpfile *close_help(Helpfile *hfile)
{
  if(hfile) {
    if(hfile->help_name)
      free(hfile->help_name);
    if(hfile->fp)
      fclose(hfile->fp);
    free(hfile);
  };
  return NULL;
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
 *  nmax     int    The number dimension of reply[].
 *  fd      FILE *  The pointer to a file opened for read.
 * Output:
 *  return  char *  Pointer to reply, or to NULL on error.
 */
static char *fgetl(char *reply, int nmax, FILE *fd)
{
  char *cptr;   /* Pointer into reply[] */
/*
 * Read as much of the next line as will fit in the buffer.
 */
  if(fgets(reply, nmax, fd) == NULL)
    return NULL;
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
  return &reply[0];
}

/*.......................................................................
 * Implement a facility like the UNIX apropos for finding functions by
 * keyword. Search for occurences of the word given by the user in function
 * names and the one line descriptions of help files. When found echo the
 * function name and one-line description to the terminal.
 *
 * Input:
 *  name   char *  The name to search for.
 */
void apropos(char *name)
{
  Pager *page=NULL; /* File pager descriptor */
  int tabpos;	    /* Current entry number in main symbol table */
  int endit=0;      /* True when user signals to finish listing */
  typedef struct {
    char *key;            /* Lower-case copy of the apropos key */
    char topic[IDX_LEN];  /* Lowe-case copy of help topic name */
    char intro[IDX_LEN];  /* Lower-case copy of a help intro entry */
  } Apropos;
  Apropos *ap;      /* Apropos state container */
/*
 * If the user didn't supply a word to be searched for, signal an
 * error.
 */
  if(name == NULL) {
    lprintf(stderr, "apropos: No search string given\n");
    return;
  };
/*
 * Allocate a container for the apropos search.
 */
  ap = (Apropos *) malloc(sizeof(Apropos));
  if(ap==NULL || (ap->key = (char *) malloc(strlen(name)+1)) == NULL) {
    lprintf(stderr, "apropos: Insufficient memory.\n");
    return;
  };
/*
 * Make a lower case copy of the apropos string to enable case-insensitive
 * searches.
 */
  low_strcpy(ap->key, name);
  ap->intro[0] = ap->topic[0] = '\0';
/*
 * Loop to find help modules.
 */
  for(tabpos=0; tabpos < num_main && !endit; tabpos++) {
    Table *symbol = main_table[tabpos];
/*
 * Ignore all but help modules.
 */
    if(symbol->class == MODULE_SYM) {
/*
 * Open the index file of the module for reading.
 */
      Indexfile *ifile = open_index(TABSTR(symbol), symbol->name, "r");
      if(ifile) {
	Topictype type=IDX_NONE;  /* Current topic type */
/*
 * Get each entry of the index file.
 */
	while(!endit && read_index(ifile) == 0) {
/*
 * Convert the topic and intro strings to lower case, to enable
 * case-insensitive comparisons.
 */
	  low_strcpy(ap->topic, ifile->topic);
	  low_strcpy(ap->intro, ifile->intro);
/*
 * See if the name sent exists in the topic name or 1-line intro.
 */
	  if(strstr(ap->topic, ap->key) || strstr(ap->intro, ap->key)) {
/*
 * Start pager only after first match.
 */
	    if(page==NULL) {
	      page = new_Pager();
	      endit = page==NULL;
	    };
/*
 * If the topic type changed, between the last entry and this one,
 * write a new topic-type title.
 */
	    if(type != ifile->type) {
	      char *title=NULL;
	      type = ifile->type;
	      switch(type) {
	      case IDX_GENERAL:
		title = "general help topics";
		break;
	      case IDX_COMMAND:
		title = "functions and commands";
		break;
	      default:
		title = "topics of unknown type";
		break;
	      };
	      endit = endit || pprintf(page,
			       "\nMatching %s in module: %s\n", title, symbol->name) < 0;
	    };
/*
 * Display the function or command name, and its one line description.
 */
	    endit = endit || pprintf(page, "%s\n", ifile->topic) < 0;
	    endit = endit || pprintf(page, "%s\n", ifile->intro) < 0;
	  };
	};
/*
 * Close the module index file.
 */
	ifile = close_index(ifile);
      };
    };
  };
/*
 * End the pager if used.
 */
  if(page)
    page = end_Pager(page, !endit, pause_output, PAGE_INT);
  else
    lprintf(stdout, "No commands match: %s\n", name);
/*
 * Delete the apropos container.
 */
  if(ap) {
    if(ap->key)
      free(ap->key);
    free(ap);
  };
  return;
}

/*.......................................................................
 * Open a module index file.
 *
 * Input:
 *  help_dir    char *  The directory path to look for the index file under.
 *                      This must be a full pathname (in UNIX this must
 *                      include the trailing /).
 *  module      char *  The module whose index is to be opened.
 *  mode        char *  The fopen() mode required.
 * Output:
 *  return Indexfile *  Pointer to a dynamically allocated index file
 *                      object, or NULL on error.
 */
static Indexfile *open_index(const char *help_dir, const char *module,
			     const char *mode)
{
  static char ext[] = ".idx";      /* Index file extension */
  Indexfile *ifile;                /* Return descriptor */
/*
 * Allocate the help file container.
 */
  ifile = (Indexfile *) malloc(sizeof(Indexfile));
  if(ifile==NULL) {
    lprintf(stderr, "open_index: Insuffient memory for descriptor.\n");
    return ifile;
  };
/*
 * Intialize ifile.
 */
  ifile->len = IDX_LEN;
  ifile->general_title = "General help topics:";
  ifile->command_title = "Functions and commands:";
  ifile->nohelp = "(No help file)";
  ifile->topic[0] = ifile->intro[0] = '\0';
  ifile->type = IDX_NONE;
  ifile->index_name = NULL;
  ifile->fp = NULL;
/*
 * Compile the index file name.
 */
  ifile->index_name = (char *) malloc(strlen(help_dir) + strlen(module) +
				      strlen(ext) + 1);
  if(ifile->index_name==NULL) {
    lprintf(stderr, "open_index: Insufficient memory to compile file name.\n");
    return close_index(ifile);
  };
  sprintf(ifile->index_name, "%s%s%s", help_dir, module, ext);
/*
 * Attempt to open the index file.
 */
  ifile->fp = fopen(ifile->index_name, mode);
  if(ifile->fp == NULL) {
    lprintf(stderr, "open_index: Unable to open %s.\n", ifile->index_name);
    lprintf(stderr, "open_index: Use the 'makeindex' command.\n");
    return close_index(ifile);
  };
/*
 * The first line of the index file is an informational message
 * specifying that the user should not edit it.
 *
 * When creating a new index file, write the informational message.
 */
  if(*mode=='w') {
    fputs("Update this file with the 'makeindex' command. Do not edit directly.\n", ifile->fp);
/*
 * When opening an index file for reading, skip the informational line.
 */
  } else if(*mode=='r') {
    int c;
    do c = fgetc(ifile->fp); while(c!='\n' && c!=EOF);
  };
  return ifile;
}

/*.......................................................................
 * Close an index file if open, and delete its container.
 *
 * Input:
 *  ifile  Indexfile *  The static Indexfile descriptor returned by
 *                      open_index().
 * Output:
 *  return Indexfile *  Allways NULL. Use like  ifile=close_index(ifile);
 */
static Indexfile *close_index(Indexfile *ifile)
{
  if(ifile) {
    if(ifile->index_name)
      free(ifile->index_name);
    if(ifile->fp)
      fclose(ifile->fp);
    free(ifile);
  };
  return NULL;
}

/*.......................................................................
 * Create an index file for a named module.
 *
 * Input:
 *  name      char *  The name of the module to index, or NULL to index
 *                    all modules.
 * Output:
 *  return     int    0 - OK.
 *                    1 - Error.
 */
int index_module(char *name)
{
  Table *symbol;    /* A symbol table entry */
  Helpfile *hfile;  /* Help file descriptor */
  Indexfile *ifile; /* Module index file */
  int bot, top;     /* Indexes of first and last matching symbols in table */
  int waserr = 0;   /* Becomes true after an error */
  int pos;          /* The index of the symbol table entry being checked */
  int i;
/*
 * A specific module name is given, look it up in the main symbol table.
 */
  if(name) {
/*
 * Search for the name in the symbol table.
 */
    char match_typ = find_symbol(name, main_table, num_main, &bot, &top);
    switch (match_typ) {
    case 'n':
      lprintf(stderr, "index_module: No module matches '%s'\n",name);
      return 1;
      break;
    case 'a':                   /* Ambiguous match */
      lex_err(comline.last);
      list_matches(bot,top,name);
      return 1;
      break;
    default:
      symbol = main_table[bot];
      if(symbol->class != MODULE_SYM) {
	lprintf(stderr, "index_module: Symbol %s does not name a module.\n",
		symbol->name);
	return 1;
      };
      break;
    };
/*
 * No module name given, so specify that all modules in the table be
 * indexed.
 */
  } else {
    bot = 0;
    top = num_main - 1;
  };
/*
 * Locate modules within the given table index range, and index them.
 */
  for(pos=bot; pos<=top && !waserr; pos++) {
    symbol = main_table[pos];
    if(symbol->class == MODULE_SYM) {
/*
 * Report progress.
 */
      lprintf(stdout, "Indexing module: %s\n", symbol->name);
/*
 * Open the output index file.
 */
      ifile = open_index(TABSTR(symbol), symbol->name, "w");
      if(ifile == NULL)
	return 1;
/*
 * List help topics.
 */
      waserr = waserr || fprintf(ifile->fp, "\n%s\n", ifile->general_title)<0 ||
                         fprintf(ifile->fp,   "-------------------\n") < 0;
/*
 * Loop through each symbol in the main table and single out help topics.
 */
      for(i=0; i<num_main && !waserr; i++) {
	Table *sym = main_table[i];
	switch (sym->class) {
	case HELP_SYM:
/*
 * If the latest function is in the named module, echo its name
 * and the one line description of its help file to the terminal.
 */
	  if(TABTAB(sym) == symbol) {
/*
 *  First echo the topic name name.
 */
	    waserr = waserr || fprintf(ifile->fp, " %s\n", sym->name) < 0;
/*
 * Open the help file.
 */
	    hfile = open_help(TABSTR(TABTAB(sym)), sym->name);
/*
 * Echo the function description to the pager.
 */
	    waserr = waserr || fprintf(ifile->fp, "   %s\n",
		     hfile==NULL ? ifile->nohelp : hfile->intro) < 0;
	    close_help(hfile);
	  };
	  break;
	default:   /* Not a help topic symbol */
	  break;
	};
      };
/*
 * Now list functions and commands.
 */
      waserr = waserr || fprintf(ifile->fp, "\n%s\n", ifile->command_title)<0 ||
                         fprintf(ifile->fp,   "----------------------\n") < 0;
/*
 * Loop through each symbol in the main table and single out functions.
 */
      for(i=0; i<num_main && !waserr; i++) {
	Table *fnsym = main_table[i];
	switch (fnsym->class) {
	case FUNC:
/*
 * If the latest function is in the named module echo its name
 * and the one line description of its help file to the terminal.
 */
	  if(TABFUNC(fnsym)->help == symbol) {
/*
 *  First echo the function name.
 */
	    waserr = waserr || fprintf(ifile->fp, " %s\n", fnsym->name) < 0;
/*
 * Open the help file.
 */
	    hfile = open_help(TABSTR(TABFUNC(fnsym)->help), fnsym->name);
/*
 * Echo the function description to the pager.
 */
	    waserr = waserr || fprintf(ifile->fp, "   %s\n",
                     hfile==NULL ? ifile->nohelp : hfile->intro) < 0;
	    close_help(hfile);
	  };
	  break;
	default:   /* Not a function symbol */
	  break;
	};
      };
/*
 * Close the completed index file.
 */
      ifile = close_index(ifile);
    };
  };
  return waserr;
}

/*.......................................................................
 * Read the next entry from an open index file.
 *
 * Input/Output:
 *  ifile   Indexfile *  The descriptor of the index file.
 *                       If an entry is found before EOF, then the topic
 *                       name and 1-line intro will be left in
 *                       ifile->topic[] and ifile->intro[].
 * Output:
 *  return        int    0 - OK.
 *                       1 - End of file reached.
 */
static int read_index(Indexfile *ifile)
{
/*
 * Initialize the return values.
 */
  ifile->intro[0] = ifile->topic[0] = '\0';
/*
 * No more entries?
 */
  if(feof(ifile->fp))
    return 1;
/*
 * Help topic names are preceded by exactly one space.
 * Search for the next line that matches this criterion.
 */
  do {
    if(fgetl(ifile->topic, ifile->len, ifile->fp)==NULL)
      return 1;
/*
 * Check for the title lines that indicate what part of the list we
 * are in.
 */
    if(strcmp(ifile->topic, ifile->general_title)==0)
      ifile->type = IDX_GENERAL;
    if(strcmp(ifile->topic, ifile->command_title)==0)
      ifile->type = IDX_COMMAND;
  } while( !(ifile->topic[0]==' ' && isalpha((int)ifile->topic[1])) );
/*
 * The next line contains the 1-line intro for the topic.
 */
  if(fgetl(ifile->intro, ifile->len, ifile->fp)==NULL)
    return 1;
/*
 * Check for the (No help file) string.
 */
  if(strstr(ifile->intro, ifile->nohelp))
    ifile->intro[0] = '\0';
  return 0;
}

/*.......................................................................
 * Make a lower-case copy of a string.
 *
 * Input:
 *  dest   char *  The string to copy to.
 *  orig   char *  The string to copy from.
 */
static void low_strcpy(char *dest, char *orig)
{
  int c;
  while( (c = *orig++) != '\0' )
    *dest++ = islower(c) ? c : tolower(c);
  *dest = '\0';
  return;
}
