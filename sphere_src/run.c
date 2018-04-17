#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "sphere.h"
#include "table.h"
#include "run.h"
#include "utils.h"

int stack_ptr=0;  /* The compile stack pointer */
int run_ptr=0;    /* The run_stack pointer */
int expr_ptr=0;   /* The expr_stack pointer */

/*
  The following structure array will hold a pointer into a user
  array and the amount by which to increment it along each dimension
  for each pass through the expression. vtyp holds the type of object
  pointed to by void *next.
*/
#define MAX_INDEXES 50
static struct {
        long addinc[3];
        char **ptr_to_elem_ptr;
} array_element[MAX_INDEXES];

/*
  Set up the run stack. This will contain a stack of intermediate
  scalar values during the execution of an arithmetic expression.
*/
#define MAXRUN 50
static Descriptor *run_stack[MAXRUN]; /* The run stack */
static Descriptor run_dsc[MAXRUN];    /* Empty descriptors for the run stack */
static Equiv temp[MAXRUN];  /* Re-useable value storage for the run stack. */

/*
  Set up the array stack. This holds the full array return values resulting
  from executing an expression.
*/
static Descriptor *expr_stack[MAXARG+1];

/*
  Keep a record of the number of indexes in use at a given time.
*/
static int num_indexes=0;

static int init_indices(char *name, Indexes *indval, long dims[3]);
static int exe_expr(long xyzmax[3]);

static int float_assign(Table *ttst, long addinc[], long ass_dims[],
			Descriptor *dtst);
static int int_assign(Table *ttst, long addinc[], long ass_dims[],
		      Descriptor *dtst);
static int logic_assign(Table *ttst, long addinc[], long ass_dims[],
			Descriptor *dtst);
static int char_assign(Table *ttst, long addinc[], long ass_dims[],
		       Descriptor *dtst);
static int pre_elemental_eval(Exprtype *expr_typ, short int start_ptr,
			      short int end_ptr, long xyzmax[]);
static void post_binop(char val_type, Equiv scalar_val);


/*.......................................................................
  This is the top level routine for controlling execution of instructions
  in the compile stack.
*/
int exe_control(short int start_ptr, short int end_ptr)
{
        int i;
	long dims[3];
	short sival,sidum;
	static Descriptor *tmp_dsc;
        Do_pars *dotst;
        Table *ttst;
	char vtyp;
/*
  Variables to record the 
*/
        long assign_addinc[3];
	long assign_dims[3];
/*
  Execute all commands up to the end_ptr stack position.
*/
        for(stack_ptr=start_ptr; stack_ptr <= end_ptr; stack_ptr++) {
          ttst = compile_stack[stack_ptr];
/*
  Determine the class of command received.
*/
          switch (ttst->class) {
/*
  New iteration of DO loop - increment value and check bounds.
*/
          case DO_PAR:
            dotst = TABDOPAR(ttst);
            *FLTPTR(dotst) = dotst->start.fval + dotst->count++ * dotst->inc.fval;
/*
  If the required end value of the DO loop has been exceeded, then
  exit the loop.
*/
            if((dotst->inc.fval > 0 && *FLTPTR(dotst) > dotst->end.fval) ||
	       (dotst->inc.fval < 0 && *FLTPTR(dotst) < dotst->end.fval) )
	      stack_ptr += dotst->skipend;
            break;
          case IDO_PAR:
            dotst = TABDOPAR(ttst);
            *INTPTR(dotst) = dotst->start.ival + dotst->count++ * dotst->inc.ival;
/*
  If the required end value of the DO loop has been exceeded, then
  exit the loop.
*/
            if((dotst->inc.ival > 0 && *INTPTR(dotst) > dotst->end.ival) ||
	       (dotst->inc.ival < 0 && *INTPTR(dotst) < dotst->end.ival) )
              stack_ptr += dotst->skipend;
            break;
/*
  A command function comes next.
*/
          case COMMAND:
/*
 * Step over the command description instruction.
 */
	    stack_ptr++;
/*
  Evaluate each argument expression of the function, stacking them on
  the array stack, and when all arguments have been evaluated, send them
  to the function.
*/
	    while(compile_stack[stack_ptr]->class == START_EXPR) {
	      for(i=0;i<3;i++) dims[i]=1;
	      if(exe_expr(dims) == -1)
		return -1;
	    };
/*
  The current compile stack entry is a NUM_ARG instruction, holding
  a record of the number of arguments to be stacked.
*/
	    sival = TABICODE(compile_stack[stack_ptr]);
/*
  Send the arguments to the function.
*/
	    if(TABFUNC(compile_stack[++stack_ptr])->fname(&expr_stack[expr_ptr-sival+1],(int) sival, NULL) == -1) {
	      lprintf(stderr,"Error occured in command: %s\n",
		      compile_stack[stack_ptr]->name);
	      return -1;
	    };
/*
  Zap the arguments from the array stack.
*/
	    array_zap(sival);
	    break;
/*
  A scalar variable follows for assignment - record it and
  record scalar dimensions and element increments.
*/
	  case VAR:
	    stack_ptr++;
            for(i=0;i<3;i++) dims[i]=1;
            if(exe_expr(dims) == -1)
              return -1;
            stack_ptr--;
/*
  Perform the assignment.
*/
	    switch (TABDESC(ttst)->atyp) {
	    case 'f':
	      *FLTPTR(TABDESC(ttst)) = *FLTPTR(expr_stack[expr_ptr]);
	      break;
	    case 'i':
	      *INTPTR(TABDESC(ttst)) = *INTPTR(expr_stack[expr_ptr]);
	      break;
	    case 'l':
	      *LOGPTR(TABDESC(ttst)) = *LOGPTR(expr_stack[expr_ptr]);
	      break;
	    case 'c':
	      if(string_copy(STRPTR(TABDESC(ttst)),
			     STRPTR(expr_stack[expr_ptr])) == -1)
		return -1;
	      break;
	    };
	    array_zap(1);
	    break;
/*
  An array variable expression follows - the left side of an
  assignment.
*/
          case ARRAY_PTR:
/*
  Ignore the ARRAY_PTR and get the stack_ptr position of the first
  assignment expression via the BR_TO instruction. Then get the INDEX_EXPR
  table entry - this contains all that we need to know.
*/
	    stack_ptr++;
	    sival = stack_ptr + TABICODE(compile_stack[stack_ptr]) + 1;
	    ttst = compile_stack[++stack_ptr];
/*
  Determine the type of the variable that is to be assigned.
*/
	    vtyp = TABINDX(ttst)->var->atyp;
/*
  Initialize the array indices.
*/
            for(i=0;i<3;i++) dims[i]=1;
	    if(init_indices(ttst->name, TABINDX(ttst), dims) == -1)
	      return -1;
/*
  Copy the indices offsets for the subsequent assignment and release
  the array_element[] slot.
*/
	    for(i=0;i<3;i++) {
	      assign_addinc[i] = array_element[num_indexes-1].addinc[i];
	      assign_dims[i]=dims[i];
	    };
	    num_indexes--;
/*
  Get the value(s) to be assigned - this may correspond to one
  array expression value or a number of scalar expression values to be
  assembled into an array prior to assignment.
*/
	    stack_ptr = sival;
	    sival = 0;
            while(stack_ptr < end_ptr && compile_stack[stack_ptr]->class == START_EXPR) {
	      sival++;
              for(i=0;i<3;i++) dims[i]=1;
              if(exe_expr(dims) == -1)
                return -1;
            };
            stack_ptr--;
/*
  One expression or many.
*/
	    if(sival != 1) {
/*
  A sequence of scalar expressions - assemble into an array before
  assignment.
*/
	      dims[0]=sival;
	      for(i=1;i<3;i++) dims[i]=1;
	      if( (tmp_dsc = descriptor_alloc(vtyp, '1', dims)) == NULL)
		return -1;
	      tmp_dsc->access = TEMP;
/*
  Copy the values into the new array.
*/
	      for(i=0;i < sival;i++) {
		switch (expr_stack[expr_ptr]->atyp) {
		case 'f':
		  FLTPTR(tmp_dsc)[i] = *FLTPTR(expr_stack[expr_ptr-sival+i+1]);
		  break;
		case 'i':
		  INTPTR(tmp_dsc)[i] = *INTPTR(expr_stack[expr_ptr-sival+i+1]);
		  break;
		case 'c':
		  if(string_copy(&STRPTR(tmp_dsc)[i],
		     STRPTR(expr_stack[expr_ptr-sival+i+1])) == -1) {
		    valof_free(tmp_dsc);
		    free(tmp_dsc);
		    return -1;
		  };
		  break;
		case 'l':
		  LOGPTR(tmp_dsc)[i] = *LOGPTR(expr_stack[expr_ptr-sival+i+1]);
		  break;
		};
	      };
/*
  Free the expression stack entries that held the individual values and
  copy the new concatenated version of the values onto the expression
  stack.
*/
	      array_zap(sival);
	      expr_stack[++expr_ptr] = tmp_dsc;
	    };
/*
  In specifying the variable to be assigned, the user had the option
  to specify array indexes - if this option was not taken then the
  variable should be re-declared with the dimensions of the assignment
  value.
*/
            if(TABINDX(ttst)->nargs == 0) {
              if(re_declare(TABINDX(ttst)->var, dims) == -1)
                return -1;
/*
  If the variable has been re-declared then *ttst->indx->ptr_to_elem_ptr
  will be the old address of the variable - reset it to point to the first
  element of the new array.
*/
	      *TABINDX(ttst)->ptr_to_elem_ptr = LOGPTR(TABINDX(ttst)->var);
/*
  Copy the new dimensions to assign_dims such that the required number
  of values are copied.
*/
              for(i=0;i<3;i++) assign_dims[i] = dims[i];
/*
  Change the addinc[] values to force an element step of one
  per assignment loop (in char_assign(), float_assign() etc..).
*/
              assign_addinc[0]=mem_size_of(vtyp);
              for(i=1;i<3;i++) assign_addinc[i] = 0;
            };
/*
  Check the dimensions of the assignment expression against those
  of the assignment variable expression.
*/
            for(i=0;i<3;i++) {
 	      if(dims[i] != 1 && assign_dims[i] != dims[i]) {
		array_zap(1);
		lprintf(stderr,"Illegal assignment due to differing array bounds.\n");
		return -1;
	      };
            };
/*
  Call the function responsible to assigning the current variable type.
*/
            switch (vtyp) {
            case 'f':
              if(float_assign(ttst,assign_addinc,assign_dims,expr_stack[expr_ptr]) == -1)
                return -1;
              break;
            case 'i':
              if(int_assign(ttst,assign_addinc,assign_dims,expr_stack[expr_ptr]) == -1)
                return -1;
              break;
            case 'l':
              if(logic_assign(ttst,assign_addinc,assign_dims,expr_stack[expr_ptr]) == -1)
                return -1;
              break;
            case 'c':
              if(char_assign(ttst,assign_addinc,assign_dims,expr_stack[expr_ptr]) == -1)
                return -1;
              break;
            };
/*
  Delete the assignment value.
*/
	    array_zap(1);
            break;
/*
  A variable expression follows - evaluate it.
*/
          case START_EXPR:
            for(i=0;i<3;i++) dims[i]=1;
            if(exe_expr(dims) == -1)
              return -1;
            stack_ptr--;
            continue;
          case BR_TRUE:
            if(*LOGPTR(expr_stack[expr_ptr]))
              stack_ptr += TABICODE(ttst);
            array_zap(1);
            break;
          case BR_FALSE:
            if(*LOGPTR(expr_stack[expr_ptr]) < 1)
              stack_ptr += TABICODE(ttst);
            array_zap(1);
            break;
          case BR_TO:
            stack_ptr += TABICODE(ttst);
            break;
          case BR_VIA:
            sival = TABICODE(ttst);
            sidum = TABICODE(compile_stack[stack_ptr+sival]);
            stack_ptr += sival + sidum;
            break;
/*
  Handle a declaratin of a user variable. The variable itself was created
  at compile time but we now wish to give it the number of elements per
  dimension specified by the user.
*/
          case DECL:
/*
  Get the user designations for the number of elements per dimension.
  Keep a record of the number of arguments recieved, in sival.
*/
            stack_ptr++;
            sival=0;
            while(compile_stack[stack_ptr]->class == START_EXPR) {
              for(i=0;i<3;i++) dims[i]=1;
              if(exe_expr(dims) == -1)
                return -1;
              sival++;
            };
/*
  Check the user designations of element numbers per dimension.
*/
            for(i=0;i<sival;i++) {
              if(*INTPTR(expr_stack[expr_ptr-sival+i+1]) < 1) {
                lprintf(stderr,"Illegal variable declaration - number of elements < 1");
                return -1;
              };
              dims[i]= *INTPTR(expr_stack[expr_ptr-sival+i+1]);
            };
            for(i=sival;i<3;i++) dims[i]=1;
/*
  Release the array stack entries.
*/
            array_zap(sival);
/*
  Re-declare the variable using the new dimensions.
*/
            if(re_declare(TABDESC(compile_stack[stack_ptr]), dims) == -1)
              return -1;
            break;
/*
  Initialize a DO loop, using the number of arguments specified in the
  type field of the current table entry. Operate on the following
  DO_PARS table entry.
*/
          case DO_INI:
            sival = TABICODE(ttst);
            dotst = TABDOPAR(compile_stack[stack_ptr+1]);
/*
  Get the values for the start and end DO variable limits.
*/
            dotst->start.fval = *FLTPTR(expr_stack[expr_ptr-sival+1]);
            dotst->end.fval   = *FLTPTR(expr_stack[expr_ptr-sival+2]);
/*
  If an increment was specified by the user then get it as well.
*/
	    dotst->inc.fval = (sival==2) ? 1.0 : *FLTPTR(expr_stack[expr_ptr]);
/*
  Test the sign of the increment against the sign of the specifed range.
*/
            if((dotst->inc.fval > 0 && dotst->end.fval < dotst->start.fval) ||
	       (dotst->inc.fval < 0 && dotst->end.fval > dotst->start.fval) ||
	       (dotst->inc.fval == 0) ) {
              lprintf(stderr,"Illegal DO step, %f for range %f -> %f \?.\n",
		      dotst->inc.fval, dotst->start.fval, dotst->end.fval);
              return -1;
            };
/*
  Initialize the iteration count to zero.
*/
            dotst->count = 0;
/*
  Zap the iteration parameters from the array stack.
*/
            array_zap(sival);
            break;
          case IDO_INI:
            sival = TABICODE(ttst);
            dotst = TABDOPAR(compile_stack[stack_ptr+1]);
/*
  Get the values for the start and end DO variable limits.
*/
            dotst->start.ival = *INTPTR(expr_stack[expr_ptr-sival+1]);
            dotst->end.ival   = *INTPTR(expr_stack[expr_ptr-sival+2]);
/*
  If an increment was specified by the user then get it as well.
*/
	    dotst->inc.ival = (sival==2) ? 1 : *INTPTR(expr_stack[expr_ptr]);
/*
  Test the sign of the increment against the sign of the specifed range.
*/
            if((dotst->inc.ival > 0 && dotst->end.ival < dotst->start.ival) ||
	       (dotst->inc.ival < 0 && dotst->end.ival > dotst->start.ival) ||
	       (dotst->inc.ival == 0) ) {
              lprintf(stderr,"Illegal DO step, %i for range %i -> %i \?.\n",
		      dotst->inc.ival, dotst->start.ival, dotst->end.ival);
              return -1;
            };
/*
  Initialize the iteration count to zero.
*/
            dotst->count = 0;
/*
  Zap the iteration parameters from the array stack.
*/
            array_zap(sival);
            break;
          case END_LINK: case EMPTY:
            break;
          case ABORT:
            return no_error;
          };
        };
        return no_error;
}


/*.......................................................................
  Re-declare the user variable (sent via the index_expr in indval)
  to have the dimensions given in dims[]. Normally this simply means
  copying dims[] into the descriptor version of adim[]. In the case where
  there are insufficient elements in the valof field then the value field
  is free()'d and re-allocated the required number of elements.
*/
int re_declare(Descriptor *dtst, long dims[3])
{
        static int num_new, num_now, i;
/*
  Calculate the total number of elements required and the current
  number in use.
*/
        num_new=num_now=1;
        for(i=0;i<3;i++) {
          num_new *= dims[i];
          num_now *= dtst->adim[i];
        };
/*
  If the variable type is string and the required number of elements is
  smaller than currently in use then first zap the current string values
  outside the required elements.
*/
        if(dtst->atyp == 'c' && num_new < num_now) {
          for(i=num_new; i<num_now;i++)
            char_free(&STRPTR(dtst)[i]);
        };
/*
  If more elements are required than can be accomodated by the current
  block of memory pointed to by the valof field, then free the current
  version (If not NO_DEL) and allocate a new block.
*/
        if(num_new > dtst->num_el) {
          if(dtst->access == NO_DEL) {
            lprintf(stderr,"Unable to allocate more memory for assignment to a\n");
            lprintf(stderr,"variable that has a system equivalent.\n");
            return -1;
          };
          if(dtst->num_el != 0)
            valof_free(dtst);
/*
  Allocate memory with the new variable dimensions and check for allocation
  errors. If one has occured then allocate a couple of elements
  for the variable, initialized to 0 and set up the dimensions accordingly.
*/
          if( (VOIDPTR(dtst) = valof_alloc(num_new, dtst->atyp)) == NULL) {
            lprintf(stderr,"Memory allocation failed.\n");
            dtst->num_el=2;
            dtst->adim[0]=2;
            for(i=1;i<3;i++)
              dtst->adim[i]=1;
            VOIDPTR(dtst) = calloc(2,sizeof(Equiv));
            return -1;
          };
/*
  Set up the descriptor so that it registers the new number of elements.
*/
          dtst->num_el=num_new;
        };
/*
  Give the descriptor its new dimensions.
*/
        for(i=0;i<3;i++) dtst->adim[i]=dims[i];
        return no_error;
}

/*.......................................................................
  Execute stack commands that make up an arithmetic, logical or character
  string expression. Return -1 on error. The return value is depositted
  in the array stack at expr_stack[expr_ptr], as a Descriptor structure.
*/
static int exe_expr(long xyzmax[3])
{
        char *elem_ptr;
        Exprtype *expr_typ;
        Table *ttst;
        short int sival,start_ptr,end_ptr, start_array;
        long xyz[3],dims[3];
        char start_index;
	int i;
/*
  Non-recursively used work variables.
*/
	static char *ctmp, var_typ;
        static size_t slen,tmpa,tmpb;
	static int itmp;
	static double fnum_a,fnum_b;
	static double ip;
	static Equiv scalar_val;
        static Descriptor d_ret_val={' ',TEMP,'0',1,{1,1,1}, &scalar_val};
/*
  Before proceding check that there is room for the expression result
  on the array stack and allocate a slot for the return value.
*/
	if(++expr_ptr > MAXARG) {
	  lprintf(stderr,"Sorry argument stack full\n");
	  return -1;
	};
/*
  Keep a record of the array stack pointer since pre_elemental_eval
  stores the return values of non-elemental array return functions
  on the expression stack and these will have to be deleted before
  exiting this funciton. Also make the return value NULL for the moment
  so that array_zap() won't attempt to zap a non-existent return value.
*/
	start_array=expr_ptr;
	expr_stack[start_array]=NULL;
/*
  Get a pointer to the expression type and determine the start and end
  stack positions over which the expression resides.
*/
        expr_typ = TABEXPR(compile_stack[stack_ptr]);
        end_ptr = stack_ptr + expr_typ->length;
/*
  Check for optional expression dimension specifiers given
  by the user.
*/
	sival=0;
	stack_ptr++;
	while(compile_stack[stack_ptr]->class == START_EXPR) {
/*
  Evaluate the scalar expression that specifies the number
  of elements on dimension 'sival'.
*/
	  for(i=0;i<3;i++) dims[i]=1;
	  if(exe_expr(dims) == -1) {
	    lprintf(stderr, "Error occurred in a {} dimensional cast\n");
	    return -1;
	  };
/*
  Apply it to the iteration counter xyzmax[].
*/
	  xyzmax[sival] = (long) *INTPTR(expr_stack[expr_ptr]);
	  if(xyzmax[sival] < 1) {
	    lprintf(stderr, "Illegal dimension specifier value: {%d}\n", xyzmax[sival]);
	    return -1;
	  };
	  array_zap(1);
	  sival++;
	};
/*
  The expression starts at the current stack position.
*/
        start_ptr = stack_ptr;
        start_index = num_indexes;
/*
  Resolve index expressions and evaluate non-elemental fuctions
  before proceding with evaluation of the elemental expression.
*/
	if(pre_elemental_eval(expr_typ, start_ptr, end_ptr, xyzmax) == -1)
	  return -1;
/*
  The expression may consist entirely of an array-return function, a
  reference to an array variable or a full descriptor - pass by name.
  These are non-elemental expressions and are resolved in pre_elemental_eval()
  and the return value of the current expression is already installed
  on the expression stack for return.
*/
	if(expr_typ->access != 'v') {
          stack_ptr = end_ptr+1;
          num_indexes = start_index;
          return no_error;
	};
/*
  Pass by value - ie evaluate the expression and return the answer
  as a temporary value on the array stack.
  First allocate memory for the return value.
*/
	if( ( expr_stack[start_array] = descriptor_alloc(expr_typ->type, expr_typ->dim, xyzmax)) == NULL)
	  return -1;
	expr_stack[start_array]->access = TEMP;
/*
  Get a pointer to the first element of the return array.
*/
	elem_ptr = (char *) VOIDPTR(expr_stack[start_array]);
/*
  The expression will be repeated for each element (x,y,z) of the array
  expression.
*/
    for(xyz[2]=0;xyz[2] < xyzmax[2];xyz[2]++) {
      for(xyz[1]=0;xyz[1] < xyzmax[1];xyz[1]++) {
        for(xyz[0]=0;xyz[0] < xyzmax[0];xyz[0]++) {
/*
  Initialize the stack pointer.
*/
          stack_ptr=start_ptr;
/*
  Execute all commands up to the end_ptr stack position.
*/
          do {
            ttst = compile_stack[stack_ptr];
/*
  Check if there is room on the run stack.
*/
            if( run_ptr+1 >= MAXRUN) {
              lprintf(stderr,"Sorry - run stack full - no more room to execute in.");
              lprintf(stderr,"Try shortening arithmetic expressions before retrying.");
              return -1;
            };
/*
  Determine the class of command received.
*/
            switch (ttst->class) {
/*
  A scalar variable or constant has been found - stack it.
*/
            case VAR: case CONST: case ARRAY_PTR: case FN_RET:
              *run_stack[++run_ptr] = *TABDESC(ttst);
              run_stack[run_ptr]->access = REF;
              break;
/*
  The current element along one of the array expression axes
  should be stacked as a scalar value.
*/
	    case HASH:
              run_stack[++run_ptr]->atyp = 'i';
              run_stack[run_ptr]->access = TEMP;
              *INTPTR(run_stack[run_ptr]) = xyz[TABICODE(ttst)];
	      break;
/*
  Pass arguments to an operator function and stack the return value on the
  run stack.
*/
	    case ADD_OP:
	      scalar_val.fval = *FLTPTR(run_stack[run_ptr-1]) +
		*FLTPTR(run_stack[run_ptr]);
	      post_binop('f',scalar_val);
	      break;
	    case SUB_OP:
	      scalar_val.fval = *FLTPTR(run_stack[run_ptr-1]) - *FLTPTR(run_stack[run_ptr]);
	      post_binop('f',scalar_val);
	      break;
	    case MUL_OP:
	      scalar_val.fval = *FLTPTR(run_stack[run_ptr-1]) * *FLTPTR(run_stack[run_ptr]);
	      post_binop('f',scalar_val);
	      break;
	    case DIV_OP:
	      if( *FLTPTR(run_stack[run_ptr]) == 0.0 ) {
		lprintf(stderr,"Division by zero error.\n");
		return -1;
	      };
	      scalar_val.fval = *FLTPTR(run_stack[run_ptr-1]) / *FLTPTR(run_stack[run_ptr]);
	      post_binop('f',scalar_val);
	      break;
	    case POW_OP:
	      var_typ = run_stack[run_ptr]->atyp;
	      switch (var_typ) {
	      case 'f':
		fnum_a = (double) *FLTPTR(run_stack[run_ptr-1]);
		fnum_b = (double) *FLTPTR(run_stack[run_ptr]);
		break;
	      case 'i':
		fnum_a = (double) *INTPTR(run_stack[run_ptr-1]);
		fnum_b = (double) *INTPTR(run_stack[run_ptr]);
		break;
	      };
/*
  If the operand to be raised is zero and the exponent is -ve then
  signal an error.
*/
	      if(fnum_a == 0.0 && fnum_b < 0.0) {
		lprintf(stderr,"Error raising 0 to a -ve power.\n");
		return -1;
	      }
/*
  Also if the operand is -ve and the exponent fractional, signal an
  error.
*/
	      else if(fnum_a < 0.0 && var_typ != 'i' && modf(fnum_b, &ip) != 0.0 ) {
		lprintf(stderr,"Error raising -ve number to a non-integral power power.\n");
		return -1;
	      };
/*
  All OK - perform the exponentiation.
*/
	      fnum_a = pow(fnum_a,fnum_b);
	      switch (var_typ) {
	      case 'f':
		scalar_val.fval = fnum_a;
		break;
	      case 'i':
		scalar_val.ival = fnum_a;
		break;
	      };
	      post_binop(var_typ,scalar_val);
	      break;
	    case GTE_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) >= *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case GT_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) > *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case LT_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) < *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case LTE_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) <= *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case EQ_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) == *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case NE_OP:
	      scalar_val.lval = *FLTPTR(run_stack[run_ptr-1]) != *FLTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case MINUS_OP:
	      scalar_val.fval = - *FLTPTR(run_stack[run_ptr]);
	      run_stack[run_ptr]->access = TEMP;
	      VOIDPTR(run_stack[run_ptr]) = &temp[run_ptr];
	      run_stack[run_ptr]->atyp = 'f';
	      temp[run_ptr]=scalar_val;
	      break;
	    case IADD_OP:
	      scalar_val.ival = *INTPTR(run_stack[run_ptr-1]) + *INTPTR(run_stack[run_ptr]);
	      post_binop('i',scalar_val);
	      break;
	    case ISUB_OP:
	      scalar_val.ival = *INTPTR(run_stack[run_ptr-1]) - *INTPTR(run_stack[run_ptr]);
	      post_binop('i',scalar_val);
	      break;
	    case IMUL_OP:
	      scalar_val.ival = *INTPTR(run_stack[run_ptr-1]) * *INTPTR(run_stack[run_ptr]);
	      post_binop('i',scalar_val);
	      break;
	    case IDIV_OP:
	      if( *INTPTR(run_stack[run_ptr]) == 0 ) {
		lprintf(stderr,"Division by zero error.\n");
		return -1;
	      };
	      scalar_val.ival = *INTPTR(run_stack[run_ptr-1]) / *INTPTR(run_stack[run_ptr]);
	      post_binop('i',scalar_val);
	      break;
	    case IGTE_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) >= *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case IGT_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) > *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case ILT_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) < *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case ILTE_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) <= *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case IEQ_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) == *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case INE_OP:
	      scalar_val.lval = *INTPTR(run_stack[run_ptr-1]) != *INTPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case IMINUS_OP:
	      scalar_val.ival = - *INTPTR(run_stack[run_ptr]);
	      run_stack[run_ptr]->access = TEMP;
	      VOIDPTR(run_stack[run_ptr]) = &temp[run_ptr];
	      run_stack[run_ptr]->atyp = 'i';
	      temp[run_ptr]=scalar_val;
	      break;
	    case SGTE_OP:
	      scalar_val.lval =  (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) >= 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case SGT_OP:
	      scalar_val.lval = (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) > 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case SLT_OP:
	      scalar_val.lval = (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) < 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case SLTE_OP:
	      scalar_val.lval = (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) <= 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case SEQ_OP:
	      scalar_val.lval = (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) == 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case SNE_OP:
	      scalar_val.lval = (char) (strcmp(*STRPTR(run_stack[run_ptr-1]),*STRPTR(run_stack[run_ptr])) != 0);
	      compress_temp(2,'l',scalar_val);
	      break;
	    case CAT_OP:
/*
  Determine the lengths of the two strings and thus the length of the
  addition.
*/
	      tmpa = strlen( *STRPTR(run_stack[run_ptr-1]));
	      tmpb = strlen( *STRPTR(run_stack[run_ptr]));
	      slen = tmpa + tmpb;
/*
  Allocate sufficient memory for the resultant string.
*/
	      if( (ctmp=stralloc(slen)) == NULL) {
		lprintf(stderr,
		   "Error concatenating: \"%.20s%s\"//\"%.20s%s\"\n",
		   *STRPTR(run_stack[run_ptr-1]), (tmpa>20) ? "...":"",
		   *STRPTR(run_stack[run_ptr]),  (tmpb>20) ? "...":"");
		return -1;
	      };
/*
  Copy each string individually.
*/
	      strncpy(ctmp, *STRPTR(run_stack[run_ptr-1]), tmpa);
	      strncpy(ctmp+tmpa, *STRPTR(run_stack[run_ptr]), tmpb);
/*
  Terminate the string.
*/
	      *(ctmp+slen) = '\0';
/*
  Install the string in the return descriptor.
*/
	      scalar_val.cptr = ctmp;
	      compress_temp(2,'c',scalar_val);
	      break;
	    case NOT_OP:
	      scalar_val.lval = ! *LOGPTR(run_stack[run_ptr]);
	      run_stack[run_ptr]->access = TEMP;
	      VOIDPTR(run_stack[run_ptr]) = &temp[run_ptr];
	      run_stack[run_ptr]->atyp = 'l';
	      temp[run_ptr]=scalar_val;
	      break;
	    case AND_OP:
	      scalar_val.lval = *LOGPTR(run_stack[run_ptr-1]) && *LOGPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case OR_OP:
	      scalar_val.lval = *LOGPTR(run_stack[run_ptr-1]) || *LOGPTR(run_stack[run_ptr]);
	      post_binop('l',scalar_val);
	      break;
	    case REG_OP:
	      scalar_val.lval = (char) match(*STRPTR(run_stack[run_ptr]), *STRPTR(run_stack[run_ptr-1]), &itmp);
	      if(itmp == 1)
		return -1;
	      compress_temp(2,'l',scalar_val);
	      break;
	    case NREG_OP:
	      scalar_val.lval = (char) !match(*STRPTR(run_stack[run_ptr]), *STRPTR(run_stack[run_ptr-1]), &itmp);
	      if(itmp == 1)
		return -1;
	      compress_temp(2,'l',scalar_val);
	      break;
	    case ITOF:
	      sival = run_ptr - TABICODE(ttst);
	      scalar_val.fval = (float) *INTPTR(run_stack[sival]);
	      VOIDPTR(run_stack[sival]) = &temp[sival];
	      run_stack[sival]->access = TEMP;
	      run_stack[sival]->atyp = 'f';
	      *FLTPTR(run_stack[sival]) = scalar_val.fval;
	      break;
	    case FTOI:
	      sival = run_ptr - TABICODE(ttst);
	      scalar_val.ival = (int) *FLTPTR(run_stack[sival]);
	      VOIDPTR(run_stack[sival]) = &temp[sival];
	      run_stack[sival]->access = TEMP;
	      run_stack[sival]->atyp = 'i';
	      *INTPTR(run_stack[sival]) = scalar_val.ival;
	      break;
/*
  Pass arguments to a function and stack its return value on the
  run stack.
*/
            case NUM_ARG:
/*
  The code in the type.icode field signals the number of arguments
  to be passed to the function.
*/
              sival= TABICODE(ttst);
/*
  The following instruction is the FUNC entry.
*/
	      ttst = compile_stack[++stack_ptr];
/*
  Send the arguments to the function.
*/
	      d_ret_val.atyp = *TABFUNC(ttst)->type;
              if( TABFUNC(ttst)->fname(&run_stack[run_ptr-sival+1], (int) sival, &d_ret_val) == -1) {
		lprintf(stderr,"Error occurred in function: %s().\n",ttst->name);
                return -1;
	      };
/*
  Zap the arguments from the run stack and stack the return value.
*/
              compress_temp(sival, d_ret_val.atyp, *EQUIVPTR(&d_ret_val));
              break;
            case BR_TRUE:
              if( *LOGPTR(run_stack[run_ptr]))
                stack_ptr += TABICODE(ttst);
              break;
            case BR_FALSE:
              if( *LOGPTR(run_stack[run_ptr]) < 1)
                stack_ptr += TABICODE(ttst);
              break;
	    case BR_TO:
	      stack_ptr += TABICODE(ttst);
	      break;
/*
  A user sub_string operation has been encounterred.
*/
            case SUB_STRING:
/*
  Get the number in the type field that signals which substring bound(s)
  were specified.
*/
              sival = TABICODE(ttst);
/*
  If there was a upper bound specified then copy it temporarily into tmpb
  and release the run_stack entry. sival=1 means that the lower bound was
  specified, sival=2 that the upper bound was specified and 3 that both
  were specified. Each bound is removed from the run stack when used, by
  compress_temp() - NB scalar_val is used for convenience - it is not
  actually used.
*/
              tmpb=0;
              if(sival > 1) {
                tmpb = (size_t) *INTPTR(run_stack[run_ptr]);
		compress_temp(1,' ',scalar_val);
              };
/*
  If the lower bound was specified then copy it into tmpa and release the
  run_stack entry. Otherwise assign the default lower bound of 1.
*/
              if(sival == 1 || sival == 3) {
                tmpa = (size_t) *INTPTR(run_stack[run_ptr]);
		compress_temp(1,' ',scalar_val);
              }
              else
                tmpa = 1;
/*
  Check relative user designations.
*/
              if(tmpa < 1 || (sival == 3 && tmpa > tmpb) ) {
                lprintf(stderr,"Illegal sub-string indices [%d:%d]\n",tmpa,tmpb);
                return -1;
              };
/*
  Find out the actual length of the string and enforce the maximum
  value of the upper bound.
*/
              slen = strlen(*STRPTR(run_stack[run_ptr]));
              tmpb = (tmpb == 0 || tmpb > slen) ? slen : tmpb;
/*
  Work out the length of the sub-string. If both the lower and upper
  requested bounds lie above the last character then we will create a
  null string, so make the length 0. Otherwise enforce a maximum upper
  bound equal to the length of the string and get the number of characters
  up to that point with respect to the requested lower bound.
*/
              slen = (tmpa > slen && tmpb >= slen) ? 0 : tmpb-tmpa+1;
/*
  Allocate the memory required to hold the string. (+ a null terminator '\0').
*/
              if( (ctmp = stralloc(slen)) == NULL) {
		lprintf(stderr,"Error occured while evaluating substring:\n\t\"%.40s...\"[%d:%d]\n", *STRPTR(run_stack[run_ptr]),tmpa,tmpb);
                return -1;
	      };
/*
  Copy the sub-string into it and terminate with a '\0'.
*/
              if(slen > 0)
                strncpy(ctmp, *STRPTR(run_stack[run_ptr])+tmpa-1, slen);
              *(ctmp+slen)='\0';
/*
  Zap the original string from the run stack and replace it with the
  sub_string.
*/
	      scalar_val.cptr = ctmp;
              compress_temp(1,'c',scalar_val);
	      break;
            default:
              lprintf(stderr,"syserr in exe_expr: unrecognised class: %d\n",ttst->class);
              return -1;
            };
          } while(stack_ptr++ < end_ptr);
/*
  Copy the new value into the next element of the return array.
*/
          switch (expr_typ->type) {
          case 'f':
            *((float *)elem_ptr) = *FLTPTR(run_stack[run_ptr]);
	    elem_ptr += sizeof(float);
            break;
          case 'i':
            *((int *)elem_ptr) = *INTPTR(run_stack[run_ptr]);
	    elem_ptr += sizeof(int);
            break;
          case 'l':
            *((char *)elem_ptr) = *LOGPTR(run_stack[run_ptr]);
	    elem_ptr += sizeof(char);
            break;
          case 'c':
/*
  When the string is already a TEMP, don't bother copying it.
*/
            if(run_stack[run_ptr]->access == TEMP) {
              *((char **)elem_ptr) = *STRPTR(run_stack[run_ptr]);
	      elem_ptr += sizeof(char **);
            }
/*
  If the string is not to be deleted then copy it rather than copying its
  pointer.  This prevents the problem that otherwise arrises when the
  return array is used to assign other elements of the same variable. eg.
  source(2:9)=source(1:8).
*/
            else {
              slen = strlen(*STRPTR(run_stack[run_ptr]));
              if( (ctmp=stralloc(slen)) == NULL)
                return -1;
              strcpy(ctmp, *STRPTR(run_stack[run_ptr]));
              *((char **)elem_ptr) = ctmp;
	      elem_ptr += sizeof(char **);
            };
            break;
          };
/*
  Decrement the run stack pointer now that the current value has been
  copied.  We can't use compress_temp() here since it would delete any
  TEMP class character strings, and these are required to be kept
  unscathed on the array stack for now.
*/
          VOIDPTR(run_stack[run_ptr]) = (void *) &temp[run_ptr];
          expr_stack[start_array]->access = run_stack[run_ptr]->access;
          run_ptr--;
/*
  If there are any further elements of the array expression to be evaluated
  then repeat the expression.
    First however, increment the memory offsets in the array index
  expressions parsed at the start of this routine. This will advance the
  pointers in each of the respective user arrays to the next required
  element.
*/
          for(i=start_index; i<num_indexes; i++) {
            *array_element[i].ptr_to_elem_ptr += array_element[i].addinc[0];
          };
        };
        for(i=start_index; i<num_indexes; i++) {
          *array_element[i].ptr_to_elem_ptr += array_element[i].addinc[1];
        };
      };
      for(i=start_index; i<num_indexes; i++) {
        *array_element[i].ptr_to_elem_ptr += array_element[i].addinc[2];
      };
    };
/*
  Restore the indexes pointer to the value it had on enterring this function.
*/
	num_indexes = start_index;
/*
  Remove the values of non-elemental array function returns from the
  expression stack.
*/
	array_zap(expr_ptr-start_array);
	return no_error;
}

/*.......................................................................
  A Private function of exe_expr() used to release the two run stack
  entries that held the arguments and install the operator result
  (scalar_val) and its type (val_type), on the run stack. This
  function must not be used when the arguments are strings - use
  compress_temp() for that.
*/
static void post_binop(char val_type, Equiv scalar_val)
{
        VOIDPTR(run_stack[run_ptr]) = &temp[run_ptr];
	run_stack[run_ptr]->access = TEMP;
	run_stack[--run_ptr]->access = TEMP;
	VOIDPTR(run_stack[run_ptr]) = &temp[run_ptr];
	run_stack[run_ptr]->atyp = val_type;
	temp[run_ptr]=scalar_val;
	return;
}

/*.......................................................................
  Zap (ntab) values from the run stack and decrement the stack pointer
  correspondingly. The only memory release that needs to be performed
  is memory allocated to character strings. These will only be deleted
  if their access class is TEMP (temporary). Next the value pointer
  is reset to point at the appropriate element of temp[].
  When this has been done, the value of vtyp is examined and if not
  ' ' then the run ptr is incremented again and the next entry assigned
  with the value of val.
*/
void compress_temp(short int ntab, char vtyp, Equiv val)
{
        static int i;
	static int last;
	last = run_ptr-ntab;
        for(i=run_ptr; i>last; i--) {
          if(run_stack[i]->atyp == 'c' && run_stack[i]->access == TEMP)
            char_free(STRPTR(run_stack[i]));
          VOIDPTR(run_stack[i]) = &temp[i];
	  run_stack[i]->access = TEMP;
        };
        run_ptr=i;
/*
  Now assign val to the next free entry.
*/
	if(vtyp != ' ') {
	  run_stack[++run_ptr]->atyp = vtyp;
	  *EQUIVPTR(run_stack[run_ptr]) = val;
	};
        return;
}

/*.......................................................................
  Allocate and return a pointer to an array of nchar+1 characters.
  nchar is the number of characters required minus the '\0' terminator.
  Whereas an nchar+1 char array is all that is required, nchar+2
  chars will be allocated to allow the pointer array-stepping idiom.
  Signals an error and returns NULL if memory allocation failed.
*/
char *stralloc(size_t nchar)
{
        static char *cptr;
        if( (cptr = (char *) calloc(nchar+2,sizeof(char))) == NULL)
          lprintf(stderr,"Memory allocation failed.");
        return cptr;
}


/*.......................................................................
  This routine is called by exe_control() to assign the return values of
  a user expression to the variable specification supplied by the user.
*/
static int float_assign(Table *ttst, long addinc[], long ass_dims[],
			Descriptor *dtst)
{
        static int i,j,k,el;
        static long num_el;
        static float *ass_el;
	static char *ass_var;
/*
  Get the pointer to the first element of the variable to be assigned.
*/
	ass_var = *TABINDX(ttst)->ptr_to_elem_ptr;
	num_el=dtst->num_el;
/*
  Index each element of the assignment variable using the three axis
  object increments provided in addinc[].
*/
	el=0;
	for(i=0;i<ass_dims[2];i++) {
	  for(j=0;j<ass_dims[1];j++) {
	    for(k=0;k<ass_dims[0];k++) {
	      ass_el = (float *) ass_var;
	      *ass_el = FLTPTR(dtst)[el];
	      if( (++el) >= num_el) el=0;
	      ass_var += addinc[0];
	    };
	    ass_var += addinc[1];
	  };
	  ass_var += addinc[2];
	};
        return no_error;
}

/*.......................................................................
  This routine is called by exe_control() to assign the return values of
  a user expression to the variable specification supplied by the user.
*/
static int int_assign(Table *ttst, long addinc[], long ass_dims[],
		      Descriptor *dtst)
{
        static int i,j,k,el;
        static long num_el;
        static int *ass_el;
	static char *ass_var;
/*
  Get the pointer to the first element of the variable to be assigned.
*/
	ass_var = *TABINDX(ttst)->ptr_to_elem_ptr;
	num_el=dtst->num_el;
/*
  Index each element of the assignment variable using the three axis
  object increments provided in addinc[].
*/
	el=0;
	for(i=0;i<ass_dims[2];i++) {
	  for(j=0;j<ass_dims[1];j++) {
	    for(k=0;k<ass_dims[0];k++) {
	      ass_el = (int *) ass_var;
	      *ass_el = INTPTR(dtst)[el];
	      if( (++el) >= num_el) el=0;
	      ass_var += addinc[0];
	    };
	    ass_var += addinc[1];
	  };
	  ass_var += addinc[2];
	};
	return no_error;
}

/*.......................................................................
  This routine is called by exe_control() to assign the return values of
  a user expression to the variable specification supplied by the user.
*/
static int logic_assign(Table *ttst, long addinc[], long ass_dims[],
			Descriptor *dtst)
{
        static int i,j,k,el;
        static long num_el;
        static char *ass_var;
/*
  Get the pointer to the first element of the variable to be assigned.
*/
	ass_var = *TABINDX(ttst)->ptr_to_elem_ptr;
	num_el= dtst->num_el;
/*
  Index each element of the assignment variable using the three axis
  object increments provided in addinc[].
*/
	el=0;
	for(i=0;i<ass_dims[2];i++) {
	  for(j=0;j<ass_dims[1];j++) {
	    for(k=0;k<ass_dims[0];k++) {
	      *ass_var = LOGPTR(dtst)[el];
	      if( (++el) >= num_el) el=0;
	      ass_var += addinc[0];
	    };
	    ass_var += addinc[1];
	  };
	  ass_var += addinc[2];
	};
	return no_error;
}

/*.......................................................................
  This routine is called by exe_control() to assign the return values of
  a user expression to the variable specification supplied by the user.
*/
static int char_assign(Table *ttst, long addinc[], long ass_dims[],
		       Descriptor *dtst)
{
        static int i,j,k,el;
        static long num_el;
        static Descriptor *datmp;
	static char *ass_var;
        static char **from_var, **ass_el, **from_el;
/*
  Get the relevant descriptor of the assignment variable.
*/
	datmp = TABINDX(ttst)->var;
/*
  Get the pointer to the first element of the variable to be assigned.
*/
	ass_var = *TABINDX(ttst)->ptr_to_elem_ptr;
/*
  And the pointer to the assignment value.
*/
        from_var = STRPTR(dtst);
	num_el = dtst->num_el;
/*
  Index each element of the assignment variable using the three axis
  object increments provided in addinc[].
*/
        el=0;
        for(i=0;i<ass_dims[2];i++) {
          for(j=0;j<ass_dims[1];j++) {
            for(k=0;k<ass_dims[0];k++) {
/*
  Get the pointer to the next element to be assigned.
*/
              ass_el = (char **) ass_var;
              from_el = (from_var + el);
	      if(string_copy(ass_el, from_el) == -1) {
/*
  If during allocation of the latest string, there is a memory allocation
  failure, zap all of the strings except one, via re-declaration to a one
  element array. This prevents one getting into the situation where
  allocation of one small string after possibly thousands being assigned
  fails, leaving all the rest of the succesfull strings allocated and only
  a few bytes of memory left, equal to the length of the final string.
*/
		for(i=0;i<3;i++) ass_dims[i]=1;
		re_declare(datmp, ass_dims);
		return -1;
	      };
/*
  Get the pointer to the next element of the array to be assigned.
*/
	      if( (++el) >= num_el) el=0;
              ass_var += addinc[0];
            };
            ass_var += addinc[1];
          };
          ass_var += addinc[2];
        };
        return no_error;
}

/*.......................................................................
  Zap (ntab) entries from the array stack and decrement the stack pointer
  correspondingly. If the access type of a particular entry is not TEMP
  then simply decrement the stack pointer. Otherwise:
        In the case of character string arrays, delete the
  strings first. In all other cases, delete the array and the descriptor.
*/
void array_zap(short int ntab)
{
        short int i;
        int last=expr_ptr-ntab;
        for(i=expr_ptr; i>last; i--) {
	  if(expr_stack[i] == NULL)
	    continue;
/*
  Only delete if the variable has the temporary access class.
*/
          switch(expr_stack[i]->access) {
	  case TEMP:
	    valof_free(expr_stack[i]);
/*
  Delete the descriptor if the access type was pass by reference or TEMP
  but definately leave it alone if it was pass by name since the descriptor
  belongs to a user variable. (NB. the absence of a break statement at the
  end of the TEMP: case is intentional.)
*/
	  case REF:
	    free(expr_stack[i]);
	    break;
	  case FN_ARRAY_VAL:
	    valof_free(expr_stack[i]);
	    break;
	  };
	};
/*
  Decrement the array stack pointer to reflect the deletions.
*/
        expr_ptr=expr_ptr-ntab;
        return;
}

/*.......................................................................
  Initialize the dimensions of a user variable before the first use in
  an expression. This requires the evaluation of index expressions
  written by the user and comparison with the array variable bounds.
  The number of elements on each dimension is compared to those sent
  in dims[3] and if these differ and the sent value is not 1, an error
  is signalled and -1 returned. The increments necessary to advance to
  the next object along each dimension is stored in the global array,
  array_element[num_indexes++]->addinc[3] and a pointer to the first
  element of the user array is deposited in the 'next' field of the
  same array.
*/
static int init_indices(char *name, Indexes *indval, long dims[3])
{
        static long vdim, new_dims[3], start[3],end[3],inc[3];
        static long mem_offset,meminc,dim;
	static long obsize;
	long inds[11];
        int i,j;
/*
  Check that there is room for the new indexes entry.
*/
	if(num_indexes >= MAX_INDEXES) {
	  lprintf(stderr,"Array index store overflowed with the addition of variable: %s\n",name);
	  return -1;
        };
/*
  Take local copies of any user arguments and then zap them from the
  array stack.
*/
        stack_ptr++;
        for(j=0;j < indval->nargs; j++) {
          for(i=0;i<3;i++) new_dims[i]=1;
          if(exe_expr(new_dims) == -1)
            return -1;
          inds[j] = (long) *INTPTR(expr_stack[expr_ptr]);
          array_zap(1);
        };
/*
  Find out the index limits - from the user arguments and from the
  variable declaration. In the Indexes start end and inc fields, a
  0 means - use the delcared bound, a number denotes the number of
  the user argument in which the bound resides.
*/
        for(i=0;i<3;i++) {
          if(indval->start[i] == 0)
            start[i] = 1;
          else
            start[i] = inds[indval->start[i]-1];

          if(indval->end[i] == 0)
            end[i] = indval->var->adim[i];
          else
            end[i] = inds[indval->end[i]-1];

          if(indval->inc[i] == 0)
            inc[i] = 1;
          else
            inc[i] = inds[indval->inc[i]-1];
        };
/*
  Initialize the address increments to zero.
*/
        for(i=0;i<3;i++)
          array_element[num_indexes].addinc[i] = 0;
/*
  Check that the bounds are sensible.
*/
        meminc=0;
        mem_offset = 0;
        for(i=0;i<3;i++) {
/*
  Check that the lower and upper bounds are between 1 and the declared
  size of the current axis of the user array.
*/
          if(end[i] < 1 || start[i] < 1) {
            lprintf(stderr,"Index specified below 1 for variable: %s.\n",name);
            return -1;
          };
          if(end[i] > indval->var->adim[i] || start[i] > indval->var->adim[i]) {
            lprintf(stderr,"Illegal request for element %d from %d elements on axis %d of variable: %s.\n",end[i],indval->var->adim[i],i,name);
	    return -1;
          };
/*
  A zero increment step is clearly illegal.
*/
          if(inc[i] == 0) {
            lprintf(stderr,"Zero array element step in index expression of variable: %s.\n",name);
              return -1;
          };
/*
  Test the sign of the increment against the sign of the specifed range.
*/
          if( (inc[i] > 0 && end[i] < start[i]) ||
              (inc[i] < 0 && end[i] > start[i]) ) {
            lprintf(stderr,"Illegal array index step: %d for range %d -> %d for variable: %s.\n",inc[i],start[i],end[i],name);
            return -1;
          };
/*
  Turn the increment on a given axis into the corresponding increment in
  objects (ie floats etc...). meminc is the amount incremented over the
  previous dimensions. NB each increment is the amount to increment
  from the last element+1 of the previous axes to the first element on the
  new axis. mem_offset is the object offset to the first element in the
  array.
*/
          if(i==0) {
            vdim = 1;
            j=0;
          }
          else
            vdim *= indval->var->adim[i-1];
/*
  Accumulate the object offset to the first element of the array.
*/
          mem_offset += (start[i]-1) * vdim;
/*
  Calculate the object increments necessary to step to the next element
  on the current axis, after the increments on the preceding axes have
  been completed.
*/
          if(start[i] != end[i]) {
            array_element[num_indexes].addinc[j] = vdim * inc[i] - meminc;
            new_dims[j] = 1+(end[i]-start[i])/inc[i];
/*
  NB. One extra increment is performed along each axis in the loops of
  incremented in exe_expr() at the same time that the next axis increment
  is performed.
*/
            meminc = vdim*inc[i]*new_dims[j];
            j++;
          };
        };
/*
  If there are any remaining dimensions that haven't been indexed,
  arrange their increments to point back to the 1st required element of the
  array. This makes it possible to use say a 1D array in a 2D expression.
  In such a case, the 1D array is visited for as many times as there
  are elements in the second dimension of the expression.
*/
	new_dims[j]=1;
	array_element[num_indexes].addinc[j]= -meminc;
        for(i=j+1;i<3;i++) {
          new_dims[i]=1;
          array_element[num_indexes].addinc[i]=0;
        };
/*
  Check the array bounds returned against those sent.
*/
        for(i=0;i<3;i++) {
          if(dims[i] != 1 && new_dims[i] != 1 && new_dims[i] != dims[i]) {
	    dim = indval->var->dim - '0';
	    lprintf(stderr,"The inclusion of array: %s(",name);
	    for(j=0; j<dim; j++) {
	      lprintf(stderr, "%d:%d:%d",start[j],end[j],inc[j]);
	      if(j < dim-1) lprintf(stderr,", ");
	    };
	    lprintf(stderr, ")\n in an array expression of dimensions (%d,%d,%d) doesn't make sense\n",dims[0],dims[1],dims[2]);
            return -1;
          }
/*
  In case the current dimension is scalar, copy the new bounds to dims[].
*/
          else if(new_dims[i] != 1)
            dims[i] = new_dims[i];
        };
/*
  Find out the size of one element of the variable with respect to
  one character.
*/
	obsize = mem_size_of(indval->var->atyp);
/*
  Scale the object increments and offset by the object size.
*/
	mem_offset *= obsize;
	for(i=0; i<3; i++)
          array_element[num_indexes].addinc[i] *= obsize;
/*
  Fill in the variable type and object offset fields in the current
  array_element entry.
*/
        array_element[num_indexes].ptr_to_elem_ptr = indval->ptr_to_elem_ptr;
        *array_element[num_indexes].ptr_to_elem_ptr = (char *)VOIDPTR(indval->var)+mem_offset;
        num_indexes++;
        stack_ptr--;
        return no_error;
}

/*.......................................................................
  Resolve the index increments and start addresses of all user variables
  and have return-once functions evaluated.
  Each set of index increments and the start address are stored in the
  static array_element[] array - see run.h. The static num_indexes counter
  contains the number of array elements in this array that are currently
  in use.
*/
static int pre_elemental_eval(Exprtype *expr_typ, short int start_ptr,
			      short int end_ptr, long xyzmax[])
{
        Table *ttst;
	Descriptor *dtmp;
	Indexes *indval;
	int i, num_args;
	long dims[3];
/*
  Look at each table entry in the current expression.
*/
        for(stack_ptr=start_ptr; stack_ptr<=end_ptr; stack_ptr++) {
          ttst=compile_stack[stack_ptr];
          switch(ttst->class) {
/*
  If a scalar variable or constant is the sole member of the expression
  then stack it in the expr_stack entry reserved by exe_expr() for the
  return value of the expression.
*/
	  case VAR:
	    if(expr_typ->access == 'N') {
	      expr_stack[expr_ptr] = TABDESC(ttst);
	      return no_error;
	    };
	    break;
	  case CONST:
	    if(expr_typ->access == 'V') {
	      expr_stack[expr_ptr] = TABDESC(ttst);
	      return no_error;
	    };
	    break;
/*
  The descriptor for a non-elemental, scalar return function
  has been encounterred - evaluate the function arguments and
  call the function giving the current descriptor as the return
  descriptor.
*/
	  case FN_RET:
/*
  Keep a record of the FN_RET descriptor and increment
  stack_ptr over the BR_TO instruction.
*/
	    dtmp = TABDESC(ttst);
	    stack_ptr += 2;
/*
  Get any optional function arguments.
*/
	    num_args=0;
	    while(compile_stack[stack_ptr]->class == START_EXPR) {
	      for(i=0;i<3;i++) dims[i]=1;
	      if(exe_expr(dims) == -1)
		return -1;
	      num_args++;
	    };
/*
  Send the arguments to the function.
*/
	    if( TABFUNC(compile_stack[stack_ptr])->fname(&expr_stack[expr_ptr-num_args+1], num_args, dtmp) == -1) {
	      lprintf(stderr,"Error occurred in function: %s().\n",compile_stack[stack_ptr]->name);
	      return -1;
	    };
/*
  Zap the arguments from the array stack.
*/
	    array_zap(num_args);
/*
  If the function call is the sole member of the expression then
  install the result on the array stack and return.
*/
	    if(expr_typ->access == 'V') {
	      expr_stack[expr_ptr] = dtmp;
	      return no_error;
	    };
	    break;
/*
  The ARRAY_PTR class entry holds a descriptor that will point to
  successive elements of an array during elemental expression
  evaluation. The array is either the return value of a return-once
  function or a user array variable. The following entry is a BR_TO entry
  used during elemental evaluation to skip index expressions etc.. we
  shall ignore this.  If the array value is a function return then
  argument expressions may follow before the function descriptor. In
  this case, and for the user array variable case, the next entry may
  be an INDEX_EXPR specifier, preceded by optional argument
  expressions.
*/
	  case ARRAY_PTR:
/*
  Increment the stack_ptr over the BR_TO instruction.
*/
	    stack_ptr += 2;
/*
  Get any optional function arguments.
*/
	    num_args=0;
	    while(compile_stack[stack_ptr]->class == START_EXPR) {
	      for(i=0;i<3;i++) dims[i]=1;
	      if(exe_expr(dims) == -1)
		return -1;
	      num_args++;
	    };
/*
  If the current entry is a function descriptor, have the function
  evaluated, giving the descriptor in the following INDEX_EXPR
  table entry as the return descriptor.
*/
	    if(compile_stack[stack_ptr]->class == FUNC) {
/*
  Get the index expression structure from the next table entry.
*/
	      indval = TABINDX(compile_stack[stack_ptr+1]);
	      if( TABFUNC(compile_stack[stack_ptr])->fname(&expr_stack[expr_ptr-num_args+1],num_args, indval->var) == -1) {
		expr_ptr--;
		lprintf(stderr,"Error occurred in function: %s().\n",compile_stack[stack_ptr]->name);
		return -1;
	      };
/*
  Zap the arguments from the array stack.
*/
	      array_zap(num_args);
/*
  Parse any user index specifications
  and set up the next entry of array_element[] to describe the address
  start position and increments necessary to step through the array
  elements.
*/
	      if(init_indices(compile_stack[stack_ptr++]->name, indval, xyzmax) == -1)
		return -1;
/*
  If the function call is the sole member of the expression, install
  the result on the expression stack entry reserved for return from
  exe_expr() and return - otherwise install it in a new entry on the
  expression stack to be zapped at the end of the elemental evaluation
  in exe_expr().
*/
	      if(expr_typ->access == 'V') {
		expr_stack[expr_ptr] = indval->var;
		return no_error;
	      };
/*
  Check if there is room for the return array descriptor on the array
  stack.
*/
	      if(++expr_ptr > MAXARG) {
		lprintf(stderr,"Sorry argument stack full\n");
		expr_ptr--;
		return -1;
	      };
/*
  Place the return array descriptor on the array stack such that its value
  will be deleted as soon as the current expression has been evaluated.
*/
	      expr_stack[expr_ptr] = indval->var;
	    }
/*
  Not a function - must be an array variable - the Indexes structure
  is held in the current table entry.
*/
	    else {
	      indval = TABINDX(compile_stack[stack_ptr]);
/*
  Parse any user index specifications
  and set up the next entry of array_element[] to describe the address
  start position and increments necessary to step through the array
  elements.
*/
	      if(init_indices(compile_stack[stack_ptr]->name, indval,  xyzmax) == -1)
		return -1;
	    };
/*
  If the expression had pass by reference or pass by name access then
  the variable should now be stacked for return, in the entry reserved
  in exe_expr on the array stack.
*/
	    switch (expr_typ->access) {
	    case 'N':
	      expr_stack[expr_ptr] = indval->var;
	      return no_error;
	      break;
	    case 'r': case 'V':
	      expr_stack[expr_ptr] = dtmp = TABDESC(ttst);
	      dtmp->num_el=1;
	      for(i=0;i<3;i++) {
		dtmp->num_el *= xyzmax[i];
		dtmp->adim[i] = xyzmax[i];
	      };
	      dtmp->dim = expr_typ->dim;
	      return no_error;
	      break;
	    };
            break;
          };
        };
/*
  Restore the stack pointer back to its entry value.
*/
        stack_ptr=start_ptr;
        return no_error;
}

/*.......................................................................
  Assign memory to the value field pointers of the run_stack descriptors
  from temp[MAXRUN] (See stack.h).
*/
void run_build(void)
{
        int i;
        for(i=0;i<MAXRUN;i++) {
          run_stack[i] = &run_dsc[i];
          VOIDPTR(run_stack[i]) = &temp[i];
          run_stack[i]->num_el=1;
          run_stack[i]->adim[0]=1;
          run_stack[i]->adim[1]=1;
          run_stack[i]->adim[2]=1;
          run_stack[i]->dim='0';
	  run_stack[i]->access = TEMP;
        };
}

/*.......................................................................
  Allocate a descriptor and an (empty) value field of the requested
  dimension. (char vtype) is used to specify the variable type
  ie. c=character string, n=floating point number, l=logical flag (char),
  i=integer constant. int adim,bdim,cdim are the three dimensions of the
  variable. Return the pointer to the entry if succesfull or a NULL
  pointer otherwise.
*/
Descriptor *descriptor_alloc(char vtype, char dim, long adim[3])
{
        static int nvals,i;
        static Descriptor *dtst;
/*
  Allocate memory for the descriptor.
*/
        dtst = (Descriptor *) malloc(sizeof(Descriptor));
        if(dtst == NULL) {
          lprintf(stderr,"Memory allocation failed.\n");
          return NULL;
        };
/*
  Copy the variable type, dimensional type and actual dimensions.
*/
        dtst->atyp    = vtype;
        dtst->dim = dim;
        for(i=0;i<3;i++) dtst->adim[i] = adim[i];
/*
  Determine the total number of elements required.
*/
        nvals = 1;
        for(i=0;i<3;i++)
          nvals *= dtst->adim[i];
        dtst->num_el = nvals;
/*
  Now allocate memory for the value fields, with the given variable type.
*/
        if( (VOIDPTR(dtst)=valof_alloc(nvals, vtype)) == NULL) {
          free(dtst);
          return NULL;
        };
/*
  Return a pointer to the new descriptor.
*/
        return dtst;
}

/*.......................................................................
  Allocate memory for a new table entry. Return a pointer to the
  entry or NULL on error. The class of the table entry and its name
  may be specified. If char *name is not a NULL pointer then memory will
  be allocated for it in the name field of the structure and the string
  pointed to by *name will be copied into it.
*/
Table *table_alloc(int class, char *name)
{
        Table *ttst;
/*
  Allocate memory for a table structure.
*/
        ttst = (Table *) malloc(sizeof(Table));
        if(ttst == NULL) {
          lprintf(stderr,"Memory alocation failed.\n");
          return NULL;
        };
/*
  Copy the class sent.
*/
        ttst->class = class;
        TABITEM(ttst) = NULL;
/*
  If a name was sent then allocate sufficent memory for a copy of it,
  copy the name into it and hook it onto the name field of the
  new table structure.
*/
        if(name == NULL) 
          ttst->name = NULL;
        else {
/*                 
  Allocate memory for the name.
*/
          ttst->name = stralloc(strlen(name));
          if(ttst->name == NULL) {
            free(ttst);
            lprintf(stderr,"Memory allocation failed.\n");
            return NULL;
          };
/*
  Install the name.
*/
          strcpy(ttst->name,name);
        };
        return ttst;
}

/*.......................................................................
  Return a pointer to 'nvals' values of the variable type encoded in
  'valtyp'. This is meant for use in allocating memory for user
  variables.
*/
void *valof_alloc(int nvals, char vartyp)
{
  return valof_realloc(NULL, vartyp, 0, nvals);
}

/*.......................................................................
 * Allocate or reallocate the memory of a user variable.
 *
 * Input:
 *  value    void *    The pointer to the current memory of the variable.
 *  vartyp   char      The type of the variable.
 *  n1        int      The current number of values in the array, if
 *                     value != NULL. Otherwise this argument is ignored.
 *  n2        int      The newly required number of values.
 * Output:
 *  return   void *    The reallocated memory, or NULL if the specified
 *                     pointer couldn't be reallocated, in which case
 *                     value remains a valid pointer to the original
 *                     array.
 */
void *valof_realloc(void *value, char vartyp, int n1, int n2)
{
  int i;
/*
 * If there is no existing array, the original number of elements is
 * zero.
 */
  if(!value)
    n1 = 0;
/*
 * If the array already has the required size, simply return the original
 * pointer.
 */
  if(value && n1 == n2) {
    return value;
  } else {
/*
 * Determine the size of the one element of the specified variable type.
 */
    size_t size = mem_size_of(vartyp);
    if(!size) {
      lprintf(stderr,"syserr: Unrecognised storage type in mem_size_of\n");
      return NULL;
    };
/*
 * If the variable type is string and the required number of elements is
 * smaller than currently in use, then first zap the current string values
 * beyond the required elements.
 */
    if(vartyp == 'c' && n2 < n1) {
      for(i=n2; i<n1;i++)
	char_free((char **) value);
    };
/*
 * Allocate or reallocate the memory as needed.
 */
    value = value ? realloc(value, size * (n2+1)) : malloc(size * (n2+1));
    if(!value) {
      lprintf(stderr,"Memory allocation failed.\n");
      return NULL;
    };
  };
/*
 * Initialize the uninitialized tail of the returned array.
 */
  switch (vartyp) {
  case 'c':
/*
 * Initialize all the newly allocated string pointers to point to a
 * null string.
 */
    {
      char **sptr = (char **) value;
      for(i=n1; i<n2; i++)
	*sptr++ = null_string;
    };
    break;
/*                                              
 * Initialize all of the newly allocated floats to 0.0.
 */
  case 'f':
    {
      float *fptr = (float *) value;
      for(i=n1; i<n2; i++)
	*fptr++ = 0.0f;
    };
    break;
/*
 * Initialize all of the newly allocated logical values to false.
 */
  case 'l':
    {
      char *lptr = (char *) value;
      for(i=n1; i<n2; i++)
	*lptr++ = 0;
    };
    break;
/*                                              
 * Initialize all of the newly allocated ints to 0.
 */
  case 'i':
    {
      int *iptr = (int *) value;
      for(i=n1; i<n2; i++)
	*iptr++ = 0;
    };
    break;
  default:
    lprintf(stderr,"syserr: Unrecognised storage type in valof_realloc\n");
    break;
  };
/*
 * Return a pointer to the new storage.
 */
  return value;
}

/*.......................................................................
  When memory is allocated for string user variables, each string pointer
  is initialised to point to the variable 'nul_string'. It is vital that
  subsequent free operations of the character strings do not attempt
  to free this pointer. In order to prevent that happening this buffer
  function has been provided, which checks the pointer sent, against
  null_string before free'ing it. Furthermore after zapping the string
  at the pointer sent, null_string is hooked onto it such that the
  elements of a user string array always have a value - hence the
  requirement that the pointer to the string be sent (char **).
*/
void char_free(char **cptr)
{
        if(*cptr != null_string) {
          free(*cptr);
          *cptr=null_string;
        };
        return;
}


/*.......................................................................
  Allocate memory for a copy of the string pointed to by *val
  at *var, having first released the string pointed to by *var.
  If a memory allocation failure occurs then the global null_string will
  be hung from *var and -1 returned.
*/
int string_copy(char **var, char **val)
{
        if(*var != *val) {
	  char_free(var);
	  if( (*var=stralloc(strlen(*val))) == NULL) {
	    *var=null_string;
	    return -1;
	  };
	  strcpy(*var,*val);
	};
	return no_error;
}
