#ifndef run_h
#define run_h

extern int stack_ptr;   /* The compile stack pointer */
extern int run_ptr;     /* The run_stack pointer */
extern int expr_ptr;    /* The expr_stack pointer */
extern int in_run_mode; /* True if in execution phase */

#define MAXSTACK 10000

#ifndef table_h
struct Table;
#endif
extern struct Table *compile_stack[];

/*
  Where a function returns a single constant value in the middle of an
  array expression, the following structure is used to contain the
  return value. Thus the function only requires one invokation even
  though the rest of the expression is evaluated multiply for each
  element of the array expression. Prior to evaluation of the array
  expression around it, the run-time system makes (char skip) equal
  to 0. When execution proceeds the run-time system sees that skip is
  false and evaluates the function, and stores the return value in the
  value field. It then sets skip=1 such that on subsequent evaluation
  of other array elements it knows not to evaluate the function again, but
  to use the value in the value field. In order to be able to skip over
  the function argument expression, the value of skipby is the number
  of stack positions to skip to the next operand.
*/
typedef struct {
  char skip;
  short int skip_by;
  Descriptor type;
} Skipeval;

/*
  Provide a structure to be hooked onto the type field of a table structure
  on the stack, for use in initializing and maintaning user array index
  expressions and the variables they apply to. At run-time, before
  an expression is executed, the skip field is set to false (0).
  In this state, on encounterring this entry, the run-time system
  changes skip to true (1). start[3], end[3] and inc[3] specify which
  indexes for the array were specified by the user. Each element is
  either 0 for default, or equal to the number of the argument in which
  the index number resides. Thus inc[2]=1 specifies that the increment
  for the 3rd dimension was specified in the first user's argument.
  The address field points to the address of the last element processed
  and inc_add[3] the address increment for each axis.
    When this entry is encounterred again when skip is false (0) then
  the address increment is applied to the address field, using the
  declared type in the type field and skips the skip_by table entries
  holding the user expressions for the indexes.
*/
typedef struct {
  char nargs;
  char start[3];
  char end[3];
  char inc[3];
  char **ptr_to_elem_ptr;
  Descriptor *var;
} Indexes;
/*
  This structure is used to contain DO loop parameters and iteration
  state. It is initialized by the function DOinit(). The value of skipend
  is the number of table entries to skip to the next instruction after
  the end instruction of the loop. *value points to the location of the
  loop variable. At the start of each iteration *value=start+count*inc
  and when *value > end, the stack pointer is incremented by skipend.
  The multiplication in the *value assignment above is not as efficient
  as a simple increment by (inc) but obviates cumulative errors.
*/

typedef struct {
  short int skipend;   /* This is set at compile time. */
  int count;
  Equiv start;
  Equiv end;
  Equiv inc;
  void *value;
} Do_pars;

/*
  The following structure is used to carry information about expression
  types.
*/

typedef struct {
  char length;    /* The no. of stack entries used by the expression. */
  char type;      /* '*'=any, 'c'=char, 'f'=float, 'i'=int, 'l'=logic */
  char dim;       /* '*'=any,'0'=scaler,'1'=array,'2'=image,'3'=cube. */
  char access;    /* If 'v'=pass by value, 'r'=pass by reference */
} Exprtype;

void run_build(void);
int exe_control(short int start_ptr, short int end_ptr);

void compress_temp(short int ntab, char vtype, Equiv val);
void array_zap(short int ntab);
int stack_line(struct Table **fntst, char is_more, short int start_pos,
	       short int end_pos);
void found_op_err(struct Table *ttst);

size_t mem_size_of(char vtyp);

#endif
