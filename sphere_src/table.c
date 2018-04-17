#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "sphere.h"
#include "table.h"
#include "lex.h"
#include "logio.h"

/*
  Declare a pointer to the main symbol table and an integer that will
  hold a record of its size.
*/

Table **main_table;
int main_max=0;
int num_main=0;

static int add_module(Module *module);
static int check_func(Functype *func, char *name);
static int check_var(Descriptor *dsc, char *name);

/*.......................................................................
 * Build symbol tables from an array of modules.
 *
 * Input:
 *  modules   Module ** The array of modules to be added (see func.h
 *                      for the contents of a Module.
 *  nmodule      int    The number of elements in the 'modules[]' array.
 * Output:
 *  return       int     0 - OK.
 *                      -1 - Abort.
 */
int module_init(Module **modules, int nmodule)
{
  int i;
/*
  Calculate the number of entries required in the main symbol table
  when all of the module symbol tables have been added. Leave extra
  room for user run-time variable declarations.
*/
  main_max = MAXVAR;
  for(i=0; i<nmodule; i++)
    main_max += modules[i]->v_num + modules[i]->f_num + modules[i]->h_num;
/*
  Allocate memory for the main symbol table.
*/
  if( (main_table = (Table **) calloc(main_max, sizeof(Table *))) == NULL) {
    lprintf(stderr, "Insufficient memory to allocate main symbol table\n");
    return -1;
  };
/*
  Concantenate all the independant module symbol tables into the main symbol
  table.
*/
  for(i=0; i<nmodule; i++) {
    if(add_module(modules[i]) == -1)
      return -1;
  };
  return 0;
}

/*.......................................................................
  Perform a binary search through the name fields of a table (stab[]),
  with (ntab) entries, for the character string sent in (*name).
  Of each tabled string, only the number of characters in (*name) will
  be compared with (*name).
    The return state is 'e' for an exact match, 'f' for
  a non-exact match (ie where the matched string has further characters
  after those specified in (*name).), 'a' for an ambiguous match (ie.
  (*name) matches the beginnings of a number of table entries) and
  'n' for no match with any table entry. For each of these cases
  (*lolim) and (*uplim) behave as:

   return state.
        'a'     *lolim = the lowest entry number of an istance where a
                         string matches (*name).
                *uplim = the highest entry number of an istance where a
                         string matches (*name).

        'f'     *lolim = *hilim = the entry number of the only matching string

        'e'     (as 'a' but note the lowest match is also the entry
                 of the exact match.) ie:
                *lolim = the entry of the exact match.
                *uplim = the highest entry number of an istance where a
                         string matches (*name).

        'n'     *lolim = The table entry number below the hypothetical
                         position where (*name) would reside if it existed.
                *uplim = *lolim +1
*/
char find_symbol(char *name, Table *stab[], int ntab, int *lolim, int *uplim)
{
        static int mid,high,low;
        static int test;
        static size_t nlen, slen;
/*
  Case convert name and find its length up to the, '\0', null terminator.
*/
        nlen = strlen(name);
/*
  Start with low and high equal to the bounds of the symbol
  table array.
*/
        low=0;
        high=ntab-1;
/*
  Binary search for nlen-1 characters of name in the symbol table.
  ie. all characters of name except the null terminator '\0'.
*/
        while (low <= high) {
          mid = (low+high) / 2;
          test = strncmp( stab[mid]->name, name, nlen);
          if (test > 0)
            high = mid - 1;
          else if (test < 0)
            low = mid + 1;
          else {
/*
  A match has been found. Check the neighbouring symbol table entries
  to see if any of them also match. In this case ambiguity will
  have to be signalled and lolim and uplim should point to the
  first and last instances of matches in the symbol table.
*/
            *lolim = *uplim = mid;
            for(*lolim = mid; (*lolim > 0) && (strncmp(stab[*lolim-1]->name,name,nlen) == 0) ; (*lolim)-- );
            for(*uplim = mid; (*uplim < ntab-1) && (strncmp(stab[*uplim+1]->name,name,nlen) == 0); (*uplim)++ );
/*
  The first match should be the shortest string - see if the first match
  is an exact match by comparing its length to (name).
*/
            slen = strlen(stab[*lolim]->name);
/*
  Exact match found.
*/
            if (slen == nlen)
              return 'e';
/*
  Ambiguous or not ambiguous?
*/
            else if (*lolim == *uplim)
              return 'f';
            else 
              return 'a';
          };
        };
/*
  No match found. lolim and uplim will be returned with the neighbouring
  entry numbers straddling the position at which name would reside if
  it existed.
*/
        *lolim = high;
        *uplim = low;
        return 'n';
}

/*.......................................................................
  This routine is used to shift all the members of a table
  (stab[tab_size]) lying at or above a given table position, tab_pos, up
  by one position to make way for the insertion of a new pointer entry
  at position tab_pos.
     It first checks that the current number of entries in the table
  (*ntab) is less than the total size of the table to make sure that
  nothing will fall off its end, and returns -1 if this is not so. If
  ok then *ntab is incremented by 1 to reflect the increase in the table
  size.
     Since this task requires only the swapping of one pointer per
  table entry shifted, it is a quick operation.
*/
int up_shift(Table *stab[], int *ntab, int tab_size, int tab_pos)
{
        int i;
/*
  Check that there is room to shift the table into.
*/
        if (*ntab >= tab_size ) {
          lprintf(stderr, "Symbol table full\n");
          return -1;
        };
/*
  Shift the table entries, by shifting their pointers.
*/
        for (i = (*ntab); i > tab_pos; i--)
          stab[i] = stab[i-1];
/*
  Increment the recorded number of variables in the table.
*/
        (*ntab)++;
        return no_error;
}


/*.......................................................................
  Take the symbol name in char *name and attempt to match it with an
  entry in the main symbol table. If the symbol name is
  ambiguous or not found, then errors will be reported and NULL returned.
  Otherwise the matched Table entry will be returned.
*/
Table *match_name(char *name)
{
        static int bot,top;
        static char retv;
/*
  Search for the symbol name in the main table.
*/
	retv = find_symbol(name, main_table, num_main, &bot, &top);
	switch (retv) {
	case 'e': case 'f':
	  return main_table[bot];
	  break;
/*
  If ambiguous list each name that could have matched.
*/
        case 'a':
          lex_err(comline.last);
	  list_matches(bot, top, name);
          break;
/*
  No match found at all.
*/
        case 'n':
          lex_err(comline.last);
          lprintf(stderr,"Unable to identify keyword \"%s\"\n",name);
          break;
        };
        return NULL;
}

/*.......................................................................
  This routine should be called when the user has erroneously typed an
  ambiguous symbol name. It reports the ambiguity as an error and lists
  all the possible matches. The extent of the ambiguity is passed to this
  routine via arguments bot, and top, which hold the number of the first
  last entries in the main symbol table that match. char *name is the
  name that the user typed.
*/
int list_matches(int bot, int top, char *name)
{
        int i;
	lprintf(stderr,"\"%s\" is ambiguous and could match any of:\n",name);
	for(i=bot; i <= top; i++) {
	  switch(main_table[i]->class) {
	  case FUNC:
	    lprintf(stderr,"Function: %s()\n",main_table[i]->name);
	    break;
	  case VAR:
	    lprintf(stderr,"Variable: %s\n",main_table[i]->name);
	    break;
	  case MODULE_SYM:
	    lprintf(stderr,"Module help topic: %s\n",main_table[i]->name);
	    break;
	  case HELP_SYM:
	    lprintf(stderr,"Help topic: %s\n",main_table[i]->name);
	    break;
	  };
	};
	return 0;
}

/*.......................................................................
  Given a symbol table, stab[] and its maximum and current sizes max_sym
  and num_sym, a new symbol name, sname, the object that it represents,
  (object) and its object class, class, insert a new symbol table entry
  at the appropriate position in stab[]. If there is a clash with an
  existing symbol name or if there is in-sufficient room in the symbol
  table for the new entry, NULL is returned. Otherwise a pointer to the
  new entry is returned.
  *num_sym will be incremented by one if the new symbol is successfully
  installed.
*/
Table *install_new_symbol(Table *stab[], int max_sym, int *num_sym, char *sname, void *object, int class)
{
        int lolim,uplim,tab_pos;
        char retv;
	Table *ttst;
/*
  First see if there is any room for the new symbol in stab[].
*/
	if(*num_sym >= max_sym) {
	  lprintf(stderr, "No room in symbol table for symbol '%s'\n", sname);
	  return NULL;
	};
/*
  Find the correct position for the new symbol in the symbol table.
*/
	retv = find_symbol(sname, stab, *num_sym, &lolim, &uplim);
/*
  Work out where the new symbol should be inserted between the existing
  members of the table.
*/
	switch (retv) {
/*
  A non-exact match was found with one or more entries of the table.
  In this case lolim is the table location of the first mimimal match
  and this must be superseded by the new symbol - after shifting
  those at and above lolim out of the way.
*/
	case 'f': case 'a':
	  tab_pos = lolim;
	  break;
/*
  No match was found - uplim and lolim point to the table entries between which
  the new entry should reside. We will want to put the new symbol at uplim,
  after shifting the entries at and above uplim out of the way.
*/
	case 'n':
	  tab_pos = uplim;
	  break;
	case 'e':
/*
  There is already a symbol with the same name.
*/
	  lprintf(stderr, "System: Multiple declaration of %s.\n",sname);
	  return NULL;
	  break;
	default:
	  lprintf(stderr, "Unknown return value \'%c\' from find_symbol()\n", retv);
	  return NULL;
	  break;
	};
/*
  Allocate memory for the new table entry.
*/
	if( (ttst=table_alloc(class, sname)) == NULL) {
	  lprintf(stderr, "Insufficient memory to install symbol '%s'\n", sname);
	  return NULL;
	};
/*
  Shift each variable in the table, from element (tab_pos) upwards,
  by one element to allow for the new entry to be inserted at
  (tab_pos). There shouldn't be any error returned by up_shift()
  since we already checked that there was room in the table. To be
  on the safe side leave the check here.
*/
	if(up_shift( stab, num_sym, max_sym, tab_pos) == -1) {
	  free(ttst);
	  return NULL;
	};
/*
  Hang the object from the type field of the new table entry.
*/
	TABITEM(ttst) = object;
/*
  Install the new entry in element (tab_pos) of stab[].
*/
	stab[tab_pos] = ttst;
        return ttst;
}

/*.......................................................................
  Install the functions and variables of a given module, plus the symbols
  that represent them, in the main symbol table. Then run the module
  initialization function if provided.
*/
static int add_module(Module *module)
{
        int i;
	Table *help_entry;
	Descriptor *dsc;
/*
  Install the module name in the symbol table under the MODULE_SYM
  table class, recording the name of the associated help directory
  in the type field. Keep a record of the table entry pointer that
  it gets allocated.
*/
	help_entry = install_new_symbol(main_table, main_max, &num_main,
		     module->name, module->help_dir, MODULE_SYM);
	if(help_entry == NULL)
	  return -1;
/*
 * Also install help topics for the current module. The type field
 * will contain the above module table entry so that we know which
 * module it came from and which help directory to use.
 */
	for(i=0; i<module->h_num; i++) {
	  if(install_new_symbol(main_table, main_max, &num_main,
		 module->h_name[i], help_entry, HELP_SYM) == NULL)
	    return -1;
	};
/*
  Check and install each function separately. Also record
  the pointer to the module help entry in each Functype descriptor.
*/
	for(i=0; i<module->f_num; i++) {
	  module->f_type[i].help = help_entry;
	  if(check_func(&module->f_type[i], module->f_name[i]) == -1)
	    return -1;
	  if(install_new_symbol(main_table, main_max, &num_main,
	     module->f_name[i], &module->f_type[i], FUNC) == NULL)
	    return -1;
	};
/*
  Do the same for the variables.
*/
	for(i=0; i<module->v_num; i++) {
	  dsc = &module->v_type[i];
/*
 * If the new variable type is marked as a descriptor substitute the
 * descriptor given in the valof field.
 */
	  if(dsc->atyp=='D')
	    dsc = (Descriptor *) VOIDPTR(dsc);
	  if(check_var(dsc, module->v_name[i]) == -1)
	    return -1;
	  if(install_new_symbol(main_table, main_max, &num_main,
             module->v_name[i], dsc, VAR) == NULL)
	    return -1;
	};
/*
 * Run initialization code for the module.
 */
	if(module->begin && module->begin() != 0)
	  return -1;
/*
 * Register the exit function for this module.
 */
	if(module->end && add_exit_fn(module->end) != 0)
	  return -1;
	return 0;
}

/*.......................................................................
  Check the validity of a function declaration. Return 0 if OK or -1
  otherwise.
*/
static int check_func(Functype *func, char *name)
{
        size_t slen,i;
/*
  Check the number of declared max and min number of arguments.
*/
        if(func->nmin > func->nmax || func->nmin < 0) {
	  lprintf(stderr,"Syserror: Function declaration of '%s' is invalid:\n\t nmin=%d, nmax=%d\?\n",name,func->nmin, func->nmax);
	  return -1;
	};
/*
  Check that there are equal numbers of type, dim and access, argument-declarations.
*/
	if((slen=strlen(func->type)) != strlen(func->dim) || slen != strlen(func->access)) {
	  lprintf(stderr,"Syserror: Function declaration of '%s' is invalid:\n\t differing numbers of argument declarators\n",name);
	  return -1;
	};
/*
  There must at least be a return declarator. If the function takes any arguments then
  there must also be at least one argument declarator.
*/
	if(slen == 0 || (func->nmin != 0 && slen == 1)) {
	  lprintf(stderr, "Too few argument declarators in function %s()\n",name);
	  return -1;
	};
/*
 * Also check that all function declarators are consistent on whether the
 * fuction has a return type or is a command with no return type.
 */
	if( (*func->type==' ' || *func->dim==' ' || *func->access==' ') &&
	    (*func->type!=' ' || *func->dim!=' ' || *func->access!=' ') ) {
	  lprintf(stderr, "Inconsistent return type declaration of function: %s()\n",name);
	  return -1;
	};
/*
  It is only necessary to check up to nmax declarators.
*/
	if(slen > func->nmax+1) slen = func->nmax+1;
/*
  Now check the type declarations.
*/
	for(i=0; i<slen; i++) {
	  switch (func->type[i]) {
	  case 'f': case 'l': case 'c': case 'i':
	    break;
/*
  Types that can't be used as function return types.
*/
	  case 'n': case 'C':
	    if(i==0) {
	      lprintf(stderr, "Illegal return type in function: %s()\n",name);
	      return -1;
	    };
	    break;
/*
  No declarator - only legal as a return type when the function is actually
  a command.
*/
	  case ' ':
	    if(i != 0) {
	      lprintf(stderr, "Null argument type declaration in function: %s()\n", name);
	      return -1;
	    };
	    break;
/*
  Wild-card legal for the return type if their is an argument to
  inherit the actual return type from.
*/
	  case '*':
	    if(i==0 && slen == 1) {
	      lprintf(stderr, "Wild-card return type without argument to copy from in function: %s()\n",name);
	      return -1;
	    };
	    break;
	  default:
	    lprintf(stderr,"Unrecognised type declarator in function: %s()\n",name);
	    return -1;
	  };
	};
/*
  Now check the dimension declarations.
*/
	for(i=0; i<slen; i++) {
	  switch (func->dim[i]) {
	  case '0':
	    break;
/*
  Non-scalar return is illegal for elemental functions.
*/
	  case '1': case '2': case '3':
	    if(i==0 && !func->once) {
	      lprintf(stderr, "Non-scalar declarator for elemental function: %s()\n",name);
	      return -1;
	    };
	    break;
/*
  No declarator - only legal as a return type - the function is actually a command.
*/
	  case ' ':
	    if(i != 0) {
	      lprintf(stderr, "Null argument dim declaration in function: %s()\n", name);
	      return -1;
	    };
	    break;
/*
  Wild-card legal in arguments. Illegal for return type except when the function
  is of the return once type with at least one argument whose dimenional type
  can be copied.
*/
	  case '*':
	    if(i==0) {
	      if(!func->once) {
		lprintf(stderr, "Non-scalar declarator for function: %s()\n",name);
		return -1;
	      }
	      else if(slen == 1) {
		lprintf(stderr, "Wild-card dim return type without argument to copy from in function: %s()\n",name);
		return -1;
	      };
	    };
	    break;
	  default:
	    lprintf(stderr,"Unrecognised dim declarator in function: %s()\n",name);
	    return -1;
	  };
	};
/*
  Now check the access declarations.
*/
	for(i=0; i<slen; i++) {
	  switch (func->access[i]) {
	  case 'v': case '?':
	    break;
	  case 'r': case 'N':
	    if(i==0) {
	      lprintf(stderr, "Illegal use of non-value return access declaration in function: %s()\n",name);
	      return -1;
	    };
	    break;
/*
  No declarator - only legal as a return type - the function is actually a
  command.
*/
	  case ' ':
	    if(i != 0) {
	      lprintf(stderr, "Null argument access declaration in function: %s()\n", name);
	      return -1;
	    };
	    break;
/*
  Wild-card return access declarators are illegal.
*/
	  case '*':
	    if(i==0) {
	      lprintf(stderr, "Illegal wild-card return access declarator in function: %s()\n",name);
	      return -1;
	    };
	    break;
	  default:
	    lprintf(stderr,"Unrecognised access declarator in function: %s()\n",name);
	    return -1;
	  };
	};
/*
  Declaration passed.
*/
	return 0;
}

/*.......................................................................
  Check the validity of a variable declaration. Return 0 if OK or -1
  otherwise.
*/
static int check_var(Descriptor *dsc, char *name)
{
        size_t i;
/*
 * Ensure that a valid descriptor was sent.
 */
	if(dsc==NULL) {
	  lprintf(stderr, "Variable %s has a NULL descriptor\n", name);
	  return -1;
	};
/*
  Now the type declaration.
*/
	switch (dsc->atyp) {
	case 'f': case 'l': case 'c': case 'i':
	  break;
	default:
	  lprintf(stderr,"Unrecognised type declarator in variable: %s\n",name);
	  return -1;
	};
/*
  Now check the dimension declaration.
*/
	switch (dsc->dim) {
	case '0': case '1': case '2': case '3':
	  break;
	default:
	  lprintf(stderr,"Unrecognised dim declarator in variable: %s()\n",name);
	  return -1;
	};
/*
  Check the number of elements.
*/
	i = dsc->adim[0] * dsc->adim[1] *dsc->adim[2];
	if(dsc->num_el <= 0 || i > dsc->num_el || i == 0) {
	  lprintf(stderr, "Invalid element number declarations of variable: %s\n",name);
	  return -1;
	};
/*
 * If the descriptor doesn't point to a variable allocate memory for
 * its value field.
 */
	if(VOIDPTR(dsc)==NULL &&
	  (VOIDPTR(dsc)=valof_alloc(dsc->num_el, dsc->atyp)) == NULL) {
	  lprintf(stderr, "Unable to allocate memory for variable: %s\n",name);
	  return -1;
	};
/*
  If the variable is a string variable and is not marked as a constant,
  make each element point to the null_string, recognised by this prog.
*/
	if(dsc->atyp == 'c' && dsc->access!=R_ONLY) {
	  for(i=0; i<dsc->num_el; i++)
	    STRPTR(dsc)[i] = null_string;
	};
/*
  Declaration passed.
*/
	return 0;
}
