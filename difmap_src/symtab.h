#ifndef symtab_h
#define symtab_h

/*.......................................................................
 * Include file for construction and maintenance of min-match symbol tables.
 */

/* Function to report or resolve ambiguous matches */

#define SYM_AMB(fn) int (fn)(struct Symtab *tab, const char *name, \
			     int first, int last, int report)

/* Function to delete a given symbol value. */

#define SYM_DEL(fn) void *(fn)(void *value)

/* Set the table increment size for use with realloc() */

#define SYM_INC 10

/* Set the maximum allowable symbol name length */

#define MAX_SYM_LEN 127

/* Symbol table entry */

typedef struct {
  char *name;     /* Lower case copy of symbol name */
  void *value;    /* Symbol value */
} Symbol;

typedef struct Symtab {
  int nsym;       /* The current number of symbols in the table */
  int nmax;       /* The allocated size of the table */
  char *type;     /* A generic name for the symbols in the table */
  Symbol *symbols;/* Array of 'nsym' sorted symbols */
  SYM_AMB(*amb);  /* Function to resolve or report ambiguous/failed matches. */
  SYM_DEL(*del);  /* Function called to delete a symbol value */
} Symtab;

Symtab *new_Symtab(int size, char *type, SYM_AMB(*sym_amb), SYM_DEL(*sym_del));
Symtab *del_Symtab(Symtab *tab);
int add_symbol(Symtab *tab, const char *name, void *value, int replace);
void *rem_symbol(Symtab *tab, const char *name);
void *get_symbol(Symtab *tab, const char *name, int report);

/* Utility functions */

SYM_AMB(amb_report);  /* Report ambiguous matches and return -1 */

#endif
