#ifndef table_h
#define table_h

/*
  All symbol tables have a name field and
  a field for a pointer to some kind of descriptor. The class field
  specifies which of the possible descriptor types is actually
  being pointed to in the type field. The possible types are listed
  in the union specifier below and the chorresponding classes
  (not in the same order or any particular order) are defined
  in the following enumeration.
*/

typedef struct Table {
  char *name;
  int class;
  union {
    void *item;
    int icode;
  } value;
} Table;

/*
 * Define macros used to access the table value field.
 */

#define TABITEM(tab)  ((tab)->value.item)
#define TABDESC(tab)  ((Descriptor *)(tab)->value.item)
#define TABINDX(tab)  ((Indexes *)(tab)->value.item)
#define TABFUNC(tab)  ((Functype *)(tab)->value.item)
#define TABOPER(tab)  ((Optype *)(tab)->value.item)
#define TABEXPR(tab)  ((Exprtype *)(tab)->value.item)
#define TABDOPAR(tab) ((Do_pars *)(tab)->value.item)
#define TABICODE(tab) ((tab)->value.icode)
#define TABSTR(tab)   ((char *)(tab)->value.item)
#define TABTAB(tab)   ((Table *)(tab)->value.item)

/*
  Define aliases for the various table types.
*/

enum {EMPTY,VAR,FUNC,OPER,CONST,BR_TRUE,BR_FALSE,BR_TO,BR_VIA,ABORT,END_LINK,
	IDO_PAR, DO_PAR, IDO_INI, DO_INI, START_EXPR, INDEX_EXPR, FN_RET,
	SUB_STRING, COMMAND, DECL, HASH,
	NUM_ARG, ARRAY_PTR, MODULE_SYM, HELP_SYM, ITOF, FTOI,
	ADD_OP, SUB_OP, MUL_OP, DIV_OP, POW_OP, GTE_OP, GT_OP,
	LT_OP, LTE_OP, EQ_OP, NE_OP, NO_OP, IADD_OP, ISUB_OP, IMUL_OP,
	IDIV_OP, IGTE_OP, IGT_OP, ILT_OP, ILTE_OP, IEQ_OP,
	INE_OP, SGTE_OP, SGT_OP, SLT_OP, SLTE_OP, SEQ_OP, SNE_OP, CAT_OP,
	NOT_OP, AND_OP, OR_OP, MINUS_OP, IMINUS_OP, REG_OP, NREG_OP};
/*
 * Declare the main table.
 */
extern Table **main_table;
extern int main_max;
extern int num_main;

/*
  Declare table handling functions.
*/

int up_shift(Table *stab[], int *ntab, int tab_size, int tab_pos);

char find_symbol(char *name, Table *stab[], int ntab, int *lolim, int *uplim);

int build_table(char *s_name[], Table *stab[], int nname, int class);

Table *install_new_symbol(Table *stab[], int max_sym, int *num_sym, char *sname, void *object, int class);

Table *table_alloc(int class, char *name);

Table *match_name(char *name);
int list_matches(int bot, int top, char *name);

#endif
