#ifndef ops_h
#define ops_h

#ifndef table_h
#include "table.h"
#endif

/*
  This header contains the declaration of all operators available
  to the user. It also sets up operator precedences. The precedences
  are set up in the following enumeration - lowest precedence first.
*/

enum {CL_BRACE, OP_BRACE, FINISH, EQUALS, CL_BR, COMMA, COLON, OP_BR, OR, AND,
      LOGIC, ADD, MULT, POWER, UNARY, CONCAT, ARRAY, FN, NUM};

/*
  Each operator must be given a declarative structure. This contains a
  pointer to the actual function that handles the operation, an index
  to the relevant column and row of the precedence table, the
  number of arguments that it requires, and an indication of the data
  type(s) that it can handle ie numeric (n) or character (c).
*/

typedef struct {
  short int f_op,i_op,s_op,l_op; /* Oper code for each arg type */
  short int op_prec;             /* The operator precedence. */
  int narg;                      /* Its number of arguments (1 or 2) */
  char atyp;                     /* The result type (n,i,c,*)  */
} Optype;

extern struct Table unminop;

int build_ops(void);
struct Table *find_ops(char **s, char *namebuf);

#endif
