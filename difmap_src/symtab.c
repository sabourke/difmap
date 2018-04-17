#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "symtab.h"
#include "logio.h"

/* Declare a type to record the results of a binary search */

typedef struct {
  int first;      /* Index of first ambiguous match */
  int last;       /* Index of last ambiguous match */
  int slot;       /* Where the symbol belongs */
  enum {
    NO_MATCH,     /* No match - 'slot' specifies where to add the symbol */
    BAD_MATCH,    /* No match - the symbol contains non-alphanumeric chars */
    AMBIG_MATCH,  /* Symbols 'first' to 'last' ambiguously match. 'slot=first'*/
    EXACT_MATCH   /* Symbol matches entry 'slot=first=last'. */
  } result;       /* The type of match */
} Match;

static Match *search_table(Symtab *tab, const char *name);
static int compare_prefix(const char *prefix, const char *name);
static int char_coll(int c);

/*.......................................................................
 * Create a new empty symbol table.
 *
 * Input:
 *  size      int   The number of entries to reserve. This need only be
 *                  a guess. The actual table will be resized as necessary.
 *  name     char * The generic name to refer to the symbols in the table
 *                  by. If NULL, "Symbol" will be substituted. The name
 *                  is used at the start of error messages (hence the
 *                  initial capital of the default name) and followed
 *                  the word, name.
 *  SYM_AMB(*amb)   Function to resolve or report ambiguous matches.
 *                  This function is called if no exact match is detected.
 *                  Its arguments are the symbol table, the symbol name
 *                  and the indexes of the first and last ambiguous symbols
 *                  in the table. It should return the table index of a
 *                  resolved match, or -1 if there is still no match. 
 *  SYM_DEL(*del)   Function to be used to delete the values of symbols.
 *                  If del!=0, this function is called by del_Symtab() to
 *                  the values of each symbol contained in the table.
 *                  It is also called by add_Symbol() when one symbol
 *                  supersedes an old one of the same name, to delete the
 *                  old symbol value. Note that the value sent should be
 *                  cast back from (void *) to the actual type of the value
 *                  before using it in the deletion function.
 * Output:
 *  return Symtab * The newly allocated symbol table, or NULL on error.
 */
Symtab *new_Symtab(int size, char *name, SYM_AMB(*sym_amb), SYM_DEL(*sym_del))
{
  Symtab *tab;   /* The new table to be returned */
/*
 * Set limits on size.
 */
  if(size <= 0)
    size = SYM_INC;
/*
 * Attempt to allocate the symbol table container.
 */
  tab = (Symtab *) malloc(sizeof(Symtab));
  if(tab == NULL) {
    lprintf(stderr, "new_Symtab: Insufficient memory to create table.\n");
    return tab;
  };
/*
 * Initialize the table at least to the point at which it can safely be sent
 * to del_Symtab().
 */
  tab->nsym = 0;
  tab->nmax = 0;
  tab->type = NULL;
  tab->symbols = NULL;
  tab->amb = sym_amb;
  tab->del = sym_del;
/*
 * Allocate a copy of the table name.
 */
  if(!name)
    name = "Symbol";
  tab->type = (char *) malloc(strlen(name) + 1);
  if(tab->type == NULL) {
    lprintf(stderr, "new_Symtab: Insufficient memory.\n");
    return del_Symtab(tab);
  };
  strcpy(tab->type, name);
/*
 * Allocate the initial array of symbol entries.
 */
  tab->symbols = (Symbol *) malloc(sizeof(Symbol) * size);
  if(tab->symbols == NULL) {
    lprintf(stderr, "new_Symtab: Insufficient memory to create %d entries.\n",
	    size);
    return del_Symtab(tab);
  };
/*
 * Record the number of symbol entries allocated.
 */
  tab->nmax = size;
/*
 * Return the table for use.
 */
  return tab;
}

/*.......................................................................
 * Delete a symbol table. Symbol values in the table are not deleted.
 *
 * Input:
 *  tab    Symtab *  The symbol table to be deleted.
 * Output:
 *  return Symtab *  Allways NULL.
 */
Symtab *del_Symtab(Symtab *tab)
{
  int i;
  if(tab) {
/*
 * Delete the table name.
 */
    if(tab->type)
      free(tab->type);
/*
 * Free each symbol name and value.
 */
    if(tab->symbols) {
      for(i=0; i<tab->nsym; i++) {
	Symbol *sym = &tab->symbols[i];
	free(sym->name);
	if(tab->del)
	  tab->del(sym->value);
      };
      free(tab->symbols);
    };
  };
  return NULL;
}

/*.......................................................................
 * Add a new symbol to a symbol table.
 *
 * Input:
 *  tab    Symtab *   The symbol table created by new_Symtab() to add the
 *                    symbol to.
 *  name     char *   The name of the new symbol. This is copied.
 *  value    void *   The symbol value to be stored with 'name' cast to
 *                    (void *). Note that the symbol is not copied, so you
 *                    should send a pointer to data that will remain in
 *                    scope until after del_Symtab() is called.
 *  replace   int     If the new symbol matches an existing one, replace
 *                    it with the new one unless replace==0.
 * Output:
 *  return    int     0 - OK.
 *                    1 - Error.
 */
int add_symbol(Symtab *tab, const char *name, void *value, int replace)
{
  Match *match; /* Descriptor of symbol table match */
  int i;
/*
 * See where the new symbol should be placed.
 */
  match = search_table(tab, name);
  if(match==NULL)
    return 1;
/*
 * Replace existing entry?
 */
  switch(match->result) {
  case EXACT_MATCH:
    {
      Symbol *sym = &tab->symbols[match->slot];
      if(replace) {
	if(tab->del)
	  tab->del(sym->value);  /* Delete the existing value if necessary */
	sym->value = value;
      } else {
	lprintf(stderr, "add_symbol: %s name \"%s\" already exists.\n",
		tab->type, name);
	return 1;
      };
    };
    break;
  case NO_MATCH:
  case AMBIG_MATCH:
    {
      char *name_copy;   /* Lower case copy of symbol name */
      int name_len;      /* Length of symbol name */
/*
 * If necessary, extend the symbol table to accomodate the new entries.
 */
      if(tab->nsym + 1 > tab->nmax) {
	size_t ntotal = tab->nsym + SYM_INC;
	Symbol *symbols = (Symbol *) realloc(tab->symbols,
					     ntotal * sizeof(Symbol));
	if(symbols == NULL) {
	  lprintf(stderr, "add_symbol: Insufficient memory to add \"%s\"\n",
		  name);
	  return 1;
	};
	tab->symbols = symbols;
	tab->nmax = ntotal;
      };
/*
 * Make a lower case copy of the symbol name.
 */
      name_len = strlen(name);
      name_copy = (char *) malloc(name_len + 1);
      if(name_copy == NULL) {
	lprintf(stderr,"add_symbol: Insufficient memory to add \"%s\"\n", name);
	return 1;
      };
      for(i=0; i<name_len; i++) {
	int c = name[i];
	name_copy[i] = isupper(c) ? tolower(c) : c;
      };
      name_copy[i] = '\0';
/*
 * Make way for the new entry by sifting up entries.
 */
      for(i=tab->nsym; i>match->slot; i--)
	tab->symbols[i] = tab->symbols[i-1];
/*
 * Store the new symbol.
 */
      tab->symbols[match->slot].name = name_copy;
      tab->symbols[match->slot].value = value;
      tab->nsym++;
    };
    break;
  case BAD_MATCH:
    lprintf(stderr,
	    "add_symbol: Illegal character in symbol name: \"%s\".\n", name);
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Remove a symbol from a symbol table.
 *
 * Input:
 *  tab    Symtab *  The symbol table to look up the symbol in.
 *  name     char *  The name of the symbol to be removed.
 * Output:
 *  return   void *  The value that corresponded to the symbol,
 *                   or NULL if not found.
 */
void *rem_symbol(Symtab *tab, const char *name)
{
  Match *match;    /* Descriptor of symbol table match */
  void *value=NULL;/* The value of the removed symbol */
  int slot = -1;   /* The index of the symbol table entry to be removed. */
/*
 * Locate the symbol to be removed.
 */
  match = search_table(tab, name);
  if(match==NULL)
    return NULL;
/*
 * Act on the result.
 */
  switch(match->result) {
  case NO_MATCH:
  case BAD_MATCH:
    slot = -1;
    break;
  case EXACT_MATCH:
    slot = match->slot;
    break;
  case AMBIG_MATCH:
    slot = tab->amb(tab, name, match->first, match->last, 0);
    break;
  };
/*
 * Remove the symbol if found.
 */
  if(slot >= 0) {
    Symbol *sym = &tab->symbols[slot];
    int i;
/*
 * Delete the symbol name.
 */
    free(sym->name);
/*
 * Preserve the symbol value for return.
 */
    value = sym->value;
/*
 * Reuse the vacated symbol table slot.
 */
    for(i=slot; i<tab->nsym-1; i++)
      tab->symbols[i] = tab->symbols[i+1];
    tab->nsym--;
  };
  return value;
}

/*.......................................................................
 * Perform a min-match lookup of a symbol in a given symbol table.
 *
 * Input:
 *  tab    Symtab *    A symbol table created by new_Symtab().
 *  name     char *    The name of the symbol to be located.
 *  report    int      If true, report failed/ambgiuous matches.
 * Output:
 *  return   void *    The value of the symbol, or NULL if not found.
 */
void *get_symbol(Symtab *tab, const char *name, int report)
{
  Match *match;    /* Descriptor of symbol table match */
/*
 * Attempt to match the given name with a symbol in the table.
 */
  match = search_table(tab, name);
  if(match==NULL)
    return NULL;
/*
 * Handle the type of match obtained.
 */
  switch(match->result) {
  case EXACT_MATCH:
    break;
  case NO_MATCH:
  case BAD_MATCH:
    if(report)
      lprintf(stderr, "%s name \"%s\" not recognised.\n", tab->type, name);
    match->slot = -1;
    return NULL;
  case AMBIG_MATCH:
    {
      SYM_AMB(*ambfn) = tab->amb ? tab->amb : amb_report;
      match->slot = ambfn(tab, name, match->first, match->last, report);
    };
    break;
  };
  return match->slot < 0 ? NULL : tab->symbols[match->slot].value;
}

/*.......................................................................
 * Perform a caseless min-match search of a symbol table for a given
 * symbol name.
 *
 * Input:
 *  tab    Symtab *  The symbol table to look up the name in.
 *  name     char *  The name of the symbol to be returned.
 * Output:
 *  return  Match *  A pointer to a static internal descriptor describing
 *                   the match, or NULL on error.
 */
static Match *search_table(Symtab *tab, const char *name)
{
  static Match match;           /* The return descriptor. */
  char name_copy[MAX_SYM_LEN+1];/* Buffer for lower case copy of 'name'. */
  int name_len;                 /* Length of 'name' including '\0' */
  int cmp;                      /* Comparison result of two strings */
  int i;
/*
 * Determine the length of 'name'.
 */
  name_len = strlen(name);
  if(name_len > MAX_SYM_LEN) {
    lprintf(stderr, "%s name \"%s\" too long.\n", tab->type, name);
    return NULL;
  };
/*
 * Make a lower-case copy of 'name'.
 */
  for(i=0; i<name_len; i++) {
    int c = name[i];
/*
 * Check for illegal characters.
 */
    if(!(i>0 ? isalnum(c) : isalpha(c)) && c!='_') {
      match.result = BAD_MATCH;
      return &match;
    };
    name_copy[i] = isupper(c) ? tolower(c) : c;
  };
  name_copy[name_len] = '\0';
/*
 * Perform a binary search of the table.
 */
  match.first = 0;
  match.last = tab->nsym-1;
  match.slot = 0;
  cmp = -1;
  while(cmp && match.first <= match.last) {
    match.slot = (match.first + match.last) / 2;
    cmp = compare_prefix(name_copy, tab->symbols[match.slot].name);
    if(cmp < 0)
      match.last  = match.slot - 1;
    else if(cmp > 0)
      match.first = match.slot + 1;
    else
      match.first = match.last = match.slot;
  };
/*
 * Not even an ambiguous match?
 */
  if(cmp != 0) {
    match.result = NO_MATCH;
    if(cmp > 0)
      match.slot++;
    match.first = match.last = match.slot;
/*
 * Exact match?
 */
  } else if(strcmp(name_copy, tab->symbols[match.slot].name) == 0) {
    match.result = EXACT_MATCH;
/*
 * The prefix matches.
 * Determine whether we have a full match or an ambiguous match.
 */
  } else {
    int slot;
/*
 * Search backwards through the table for the first ambiguous match.
 */
    for(slot=match.slot-1; slot>=0; slot--) {
      if(strncmp(tab->symbols[slot].name, name_copy, name_len) == 0)
	match.first = slot;
      else
	break;
    };
/*
 * Search forwards through the table for the last ambiguous match.
 */
    for(slot=match.slot+1; slot<tab->nsym; slot++) {
      if(strncmp(tab->symbols[slot].name, name_copy, name_len) == 0)
	match.last = slot;
      else
	break;
    };
/*
 * Is the symbol unambiguous?
 */
    if(match.first == match.last) {
      match.result = EXACT_MATCH;
/*
 * Ambiguous match.
 */
    } else {
      match.result = AMBIG_MATCH;
/*
 * Record where the symbol would go if it were to be inserted in the table.
 */
      match.slot = match.first;
    };
  };
/*
 * Return the results.
 */
  return &match;
}

/*.......................................................................
 * Compare a target symbol name prefix to a given table symbol name,
 * using the normal alphanumeric collating sequence to determine which
 * name appears earlier or later in the alphabet.
 *
 * Input:
 *  prefix   char *   The name to be compared.
 *  name     char *   Compare up to strlen(prefix) characters of 'prefix'
 *                    to 'name'.
 * Output:
 *  return    int     -1 - prefix < name up to the length of prefix..
 *                     0 - prefix = name up to the length of prefix.
 *                    +1 - prefix > name up to the length of prefix..
 */
static int compare_prefix(const char *prefix, const char *name)
{
  const char *pptr;   /* Pointer into 'prefix' */
  const char *nptr;   /* Pointer into 'name' */
/*
 * Locate the first character at which the two names differ.
 */
  pptr = prefix;
  nptr = name;
  for(pptr=prefix,nptr=name; *pptr && *nptr && *pptr == *nptr; nptr++,pptr++);
/*
 * Exact match to the given prefix.
 */
  if(*pptr=='\0')
    return 0;
/*
 * Name finished before prefix?
 */
  if(*nptr=='\0')
    return 1;
/*
 * Compare the first differing characters using an alphabetic collating
 * sequence to determine which string appears first in the alphabet.
 */
  return (char_coll(*pptr) - char_coll(*nptr)) > 0 ? 1 : -1;
}

/*.......................................................................
 * Return the position of a character in internally defined alphabet.
 * The alphabet is defined as [0..9] [_] [a..z]. Other characters are
 * not defined.
 *
 * Input:
 *  c        int    The character to be located.
 * Output:
 *  return   int    The location [>0], or 0 if not located.
 */
static int char_coll(int c)
{
  static char *alphabet = "0123456789_abcdefghijklmnopqrstuvwxyz";
  char *cptr = strchr(alphabet, c);
  return cptr ? (cptr-alphabet)+1 : 0;
}

/*.......................................................................
 * Utility ambiguity report function.
 *
 * Input:
 *  tab     Symtab *   The symbol table to report on.
 *  name      char *   The name of the symbol that is ambiguous.
 *  first      int     The index of the first ambiguous match in tab->symbols[].
 *  last       int     The index of the last ambiguous match in tab->symbols[].
 *  report     int     If true, report failed/ambiguous matches.
 * Output:
 *  return     int     The index of the resolved match in tab->symbols[],
 *                     or -1 if there is still no match.
 */
int amb_report(Symtab *tab, const char *name, int first, int last, int report)
{
  int slot;  /* The index of the symbol being reported */
  int llen=0;/* The length of the current line */
/*
 * Report all ambiguous matches.
 */
  lprintf(stderr, "%s name \"%s\" is ambiguous with:\n", tab->type, name);
  for(slot=first; slot<=last; slot++) {
    char *sym_name = tab->symbols[slot].name;
    int slen = strlen(sym_name);
/*
 * Terminate the current line if there is insufficient room on it for
 * the next symbol plus ", " separator.
 */
    if(llen>0 && llen+slen+2 > 80) {
      lprintf(stderr, "\n");
      llen = 0;
    };
/*
 * Append the latest symbol to the current line.
 */
    llen += lprintf(stderr, "%s%s", llen==0 ? "  ":", ", sym_name);
  };
/*
 * Terminate the last line.
 */
  lprintf(stderr, "\n");
/*
 * Ambiguity not resolved.
 */
  return -1;
}
