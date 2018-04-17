#include <stdio.h>
#include <stdlib.h>

#include "obs.h"
#include "enumpar.h"
#include "logio.h"

static Enumtab *Stokes_table(void);

/*.......................................................................
 * Return a Stokes enumeration symbol table.
 *
 * Output:
 *  return   Enumtab *  The stokes name symbol table, or NULL on error.
 */
static Enumtab *Stokes_table(void)
{
  static Enumtab *etab = NULL;  /* The stokes name symbol table */
/*
 * List the assignment of names to Stokes enumerators.
 */
  static Enumpar spar[] = {
    {"I", SI}, {"Q", SQ}, {"U", SU},  {"V", SV}, 
    {"RR",RR}, {"LL",LL}, {"RL",RL},  {"LR",LR},
    {"XX",XX}, {"YY",YY}, {"XY",XY},  {"YX",YX},
    {"PI",PI_POL}
  };
/*
 * Allocate and initialize the symbol table?
 */
  if(!etab)
    etab = new_Enumtab(spar, sizeof(spar) / sizeof(spar[0]), "Polarization");
  return etab;
}

/*.......................................................................
 * Return the name of the polarization associated with a given Stokes
 * enumerator.
 *
 * Input:
 *  pol   Stokes   The enumerator to name.
 * Output:
 *  return  char * The name of the enumerator, or "" if not found.
 */
char *Stokes_name(Stokes pol)
{
  Enumtab *stab = Stokes_table();
  return stab ? name_enum(stab, (int)pol, "") : "";
}

/*.......................................................................
 * Lookup a Stokes enumerator by name.
 *
 * Input:
 *  name   char *   The name of the polarization.
 * Output:
 *  return Stokes   The stokes enumeration corresponding to name[], or
 *                  NO_POL if not found.
 */
Stokes Stokes_id(char *name)
{
  Enumtab *stab = Stokes_table();
  Enumpar *epar = stab ? find_enum(stab, name) : NULL;
  return epar ? (Stokes) epar->id : NO_POL;
}

