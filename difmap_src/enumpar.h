#ifndef enumpar_h
#define enumpar_h

#include "symtab.h"

/* Declare a struct used to associate names with enumerators */

typedef struct {
  char *name;     /* Enumeration name */
  int id;         /* Enumeration ID */
} Enumpar;

/* Declare a type to contain a min-match enumeration symbol table */

typedef Symtab Enumtab;

/* Enumeration symbol table constructor. */

Enumtab *new_Enumtab(Enumpar *epar, int npar, char *type);

/* Enumeration symbol table destructor. */

Enumtab *del_Enumtab(Enumtab *etab);

/* Match a given name to a given Enumpar table entry */

Enumpar *find_enum(Enumtab *etab, const char *name);

/* Return the name of a given enumeration identifier */

char *name_enum(Enumtab *etab, int id, char *def);

#endif
