#include <stdio.h>
#include <stdlib.h>

#include "sphere.h"
#include "table.h"
#include "ops.h"
#include "logio.h"

/*
  This file contains the declaration of the operator symbol table.
  For each operator define its type in the following array.
*/

static Optype op_type[] = {
  {ADD_OP,   IADD_OP,   NO_OP,   NO_OP,    ADD      ,2, '*'},
  {SUB_OP,   ISUB_OP,   NO_OP,   NO_OP,    ADD      ,2, '*'},
  {MUL_OP,   IMUL_OP,   NO_OP,   NO_OP,    MULT     ,2, '*'},
  {DIV_OP,   IDIV_OP,   NO_OP,   NO_OP,    MULT     ,2, '*'},
  {POW_OP,   POW_OP,    NO_OP,   NO_OP,    POWER    ,2, '*'},
  {GTE_OP,   IGTE_OP,   SGTE_OP, NO_OP,    LOGIC    ,2, 'l'},
  {GT_OP,    IGT_OP,    SGT_OP,  NO_OP,    LOGIC    ,2, 'l'},
  {LT_OP,    ILT_OP,    SLT_OP,  NO_OP,    LOGIC    ,2, 'l'},
  {LTE_OP,   ILTE_OP,   SLTE_OP, NO_OP,    LOGIC    ,2, 'l'},
  {EQ_OP,    IEQ_OP,    SEQ_OP,  NO_OP,    LOGIC    ,2, 'l'},
  {NE_OP,    INE_OP,    SNE_OP,  NO_OP,    LOGIC    ,2, 'l'},
  {NO_OP,    NO_OP,     NO_OP,   NOT_OP,   UNARY    ,1, 'l'},
  {NO_OP,    NO_OP,     REG_OP,  NO_OP,    LOGIC    ,2, 'l'},  
  {NO_OP,    NO_OP,     NREG_OP, NO_OP,    LOGIC    ,2, 'l'},  
  {NO_OP,    NO_OP,     CAT_OP,  NO_OP,    CONCAT   ,2, 'c'},
  {NO_OP,    NO_OP,     NO_OP,   AND_OP,   AND      ,2, 'l'},
  {NO_OP,    NO_OP,     NO_OP,   OR_OP,    OR       ,2, 'l'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    OP_BR    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    CL_BR    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    OP_BR    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    CL_BR    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    COLON    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    COMMA    ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    EQUALS   ,2, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    FINISH   ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    OP_BRACE ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    CL_BRACE ,1, '*'},
  {NO_OP,    NO_OP,     NO_OP,   NO_OP,    CL_BR    ,1, '*'}
};

/*
  In the same order as the above array of types, define the array
  of operator symbols.
*/

static char *op_name[] = {
   "+",
   "-",
   "*",
   "/",
   "^",
   ">=",
   ">",
   "<",
   "<=",
   "==",
   "!=",
   "!",
   "~",
   "!~",
   "//",
   "&",
   "|",
   "(",
   ")",
   "[",
   "]",
   ":",
   ",",
   "=",
   "",
   "{",
   "}",
   "`"
};

/*
  Declare the operator table to have as many members as declared
  in f_name above.
*/

static Table *ops_table[sizeof(op_name)/sizeof(char *)];
static int num_ops = sizeof(op_name)/sizeof(char *);

/*
  Unary minus has to be treated separately since it has the same name as
  the binary minus operator.
*/

static Optype unmin_type = {MINUS_OP, IMINUS_OP, NO_OP, NO_OP,  UNARY, 1, '*'};
Table unminop = {"-", OPER, {&unmin_type}};


/*.......................................................................
  Build the operator symbol table.
*/
int build_ops(void)
{
  int i,ntab;
/*
  Loop for each operator.
*/
  for(i=0,ntab=0; i < num_ops; i++) {
    if(install_new_symbol(ops_table, num_ops, &ntab, op_name[i], &op_type[i], OPER) == NULL) {
      lprintf(stderr, "Error building operator symbol table\n");
      return -1;
    };
  };
  return no_error;
}

/*.......................................................................
  This routine attempts to match up to two characters in (char **s)
  with symbols in the operator symbol table. It first tries to match
  just the first character. If more than one match is found then the
  second character is used to resolve the ambiguity - if possible.
  The character pointer *s will be left pointing to the next character to
  be processed. On return (char *namebuf) will contain the name
  that was or wasn't matched and the table entry will be returned.
  If not found NULL will be returned.
*/
Table *find_ops(char **s, char *namebuf)
{
  static char retv;
  static int botbck,top,bot;
/*
  Copy the first character into namebuf and terminate it.
*/
  namebuf[0] = **s;
  namebuf[1] = '\0';
  (*s)++;
/*
  Search for it in the operator symbol-table.
*/
  retv = find_symbol(namebuf, ops_table, num_ops, &bot, &top);
/*
  If not found at all then signal it by returning NULL.
*/
  if(retv == 'n') return NULL;
/*
  If only one possible match exists then return it (exact match).
*/
  if(bot == top) return ops_table[top];
/*
  There is more than one table entry that could match the operator
  depending upon whether it is a single character or two character
  operator. Resolve the ambiguity by appending the next character to the
  operator in namebuf and searching for the two operator character
  sequence in ops_table.
*/
  botbck = bot;
  namebuf[1] = **s;
  namebuf[2] = '\0';
  retv = find_symbol(namebuf, &ops_table[bot], top-bot+1, &bot, &top);
/*
  If a match was made, report it as exact and increment the character pointer
  *s to point to the third character.
*/
  if(retv != 'n' && namebuf[1] != '\0') {
    (*s)++;
    return ops_table[botbck+bot];
  }
/*
  The second character didn't match any operator so return the single
  operator match.
*/
  else {
    bot = botbck;
    return ops_table[botbck];
  };
}
