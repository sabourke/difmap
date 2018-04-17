#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "logio.h"
#include "enumpar.h"

/*.......................................................................
 * Create a min-match symbol table for a given enumeration list.
 *
 * Input:
 *  epar   Enumpar *   The array of 'npar' enumerators to be placed in
 *                     the table.
 *  npar      char *   The number of elements in 'epar'.
 *  type      char *   The type name of the enumerator.
 * Output:
 *  return Enumtab *   The min-match symbol table, or NULL on error.
 */
Enumtab *new_Enumtab(Enumpar *epar, int npar, char *type)
{
  Enumtab *etab;   /* The return table */
  int i;
/*
 * Allocate the symbol table.
 */
  etab = new_Symtab(npar, type, amb_report, (SYM_DEL(*)) 0);
  if(etab == NULL)
    return NULL;
/*
 * Install each of the enumerators in the symbol table.
 */
  for(i=0; i<npar; i++) {
    if(add_symbol(etab, epar[i].name, (void *) &epar[i], 0))
      return del_Enumtab(etab);
  };
  return etab;
}

/*.......................................................................
 * Delete an enumeration symbol table.
 *
 * Input:
 *  etab    Enumtab *   The table to be deleted.
 * Output:
 *  return  Enumtab *   Allways NULL.
 */
Enumtab *del_Enumtab(Enumtab *etab)
{
  if(etab) {
    del_Symtab(etab);
  };
  return NULL;
}

/*.......................................................................
 * Perform a min-match case-less search for a name in a table of
 * enumerations, and return the matching enumerator descriptor,
 * or NULL if not found.
 *
 * Input:
 *  etab   Enumtab *  The enumeration symbol table.
 *  name      char *  The enumeration name to be matched.
 * Output:
 *  return Enumpar *  The pointer to the matching enumeration descriptor,
 *                    or NULL if not found.
 */
Enumpar *find_enum(Enumtab *etab, const char *name)
{
/*
 * Check args.
 */
  if(etab==NULL || name==NULL) {
    lprintf(stderr, "find_enum: NULL %s intercepted.\n",
	    etab==NULL ? "Symbol table" : "name");
    return NULL;
  };
/*
 * Search for the given name.
 */
  return (Enumpar *) get_symbol(etab, name, 1);
}

/*.......................................................................
 * Return the symbol table name associated with a given enumeration
 * indentifier.
 *
 * Input:
 *  etab   Enumtab *  The enumeration symbol table.
 *  id         int    The enumeration identifier to be located.
 *  def       char *  The default string to return if the enumerator.
 *                    This can be anything, including NULL.
 * Output:
 *  return    char *  The name associated with the identifier, or NULL
 *                    if not found.
 */
char *name_enum(Enumtab *etab, int id, char *def)
{
  int i;
  if(!etab) {
    lprintf(stderr, "name_enum: NULL symbol table received.\n");
    return def;
  };
  for(i=0; i<etab->nsym; i++) {
    Symbol *sym = etab->symbols + i;
    Enumpar *epar = (Enumpar *) sym->value;
    if(epar->id == id)
      return epar->name;
  };
  return def;
}
