#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "sig.h"
#include "sphere.h"
#include "lex.h"
#include "table.h"
#include "ops.h"
#include "utils.h"
#include "run.h"

static void lex_test(void);
static void clear_compile(short int start_ptr, short int end_ptr);

int in_run_mode=0;

char *null_string="";

int startup(Module **modules, int nmodule, const char *bootvar)
{
/*
  Build the operator symbol table.
*/
        if(build_ops() == -1)
          return -1;
/*
  Initilialize the run stack.
*/
        run_build();
/*
  Set up the signal handlers - interrupt handlers etc..
*/
        sig_init();
/*
  Initialise the user file-allocation-table.
*/
	fat_init();
/*
  Set up the main symbol table and do anything that the application
  writer deems neccessary to initialize their applications.
*/
	if(module_init(modules, nmodule))
	  return -1;
/*
 * Initialize the command input buffers.
 */
	if(com_init(bootvar))
	  return -1;
/*
  Call the control routine - This routine will only return to exit
  the program.
*/
        lex_test();
        return no_error;
}


/*.......................................................................
  Delete a variable-table entry. This entails deallocating memory
  held in the value field of the descriptor, deallocating memory
  held by the descriptor itself and finally deallocating memory held
  by the name string. The name and descriptor pointer fields will then
  be assigned NULL pointers. If this function is called for an entry
  with NULL fields already in place then none of the above will be
  performed. The exception is for an entry with a NULL name field
  but an active descriptor pointer field - the case for work
  variables used to contain temporary values during expression parsing.
  In this case the descriptor and its value will be removed as above.
*/
void var_free(Table *stab)
{
/*
  If the table current entry has a name delete it.
*/
        if( stab->name != NULL) {
          free(stab->name);
          stab->name = NULL;
        };
/*
  Give the table entry the EMPTY class.
*/
        stab->class = EMPTY;
/*
  See if the descriptor field is active before attempting to
  delete its subjects.
*/
        if(TABITEM(stab) != NULL) {
/*
  Free the value field and its contents (ie. if the variable
  is of string type, each string will have to be deleted
  first.
*/
          valof_free(TABDESC(stab));
/*                
  Then deallocate the memory held by the descriptor structure.
*/
          free(TABITEM(stab));
/*
  Assign a NULL pointer to the released entry.
*/
          TABITEM(stab) = NULL;
        };
        return;
}
/*.......................................................................
  Delete the value held in valof field of the descriptor sent as the
  sole argument. If the descriptor describes a string array, then
  in addition, each string will also be deallocated.
*/
void valof_free(Descriptor *dtst)
{
        static int i;
/*
  If the variable is a character string variable we must deallocate
  the memory assigned to each string before deallocating the value field.
*/
        if(dtst->atyp == 'c') {
          for(i=0;i < dtst->num_el; i++)
            char_free(&STRPTR(dtst)[i]);
        };
/*
  Deallocate the memory of the value field.
*/
        free(VOIDPTR(dtst));
        return;
}


/*.......................................................................
  Convert string, s, to lower case. Return the length of the string,
  the length not including the '\0' null terminator.
*/
int lowstr(char s[])
{
        static int i;
        for(i=0; (s[i] = (char) tolower((int) s[i])) != '\0'; i++);
        return i;
}


/*.......................................................................
  Test the lexical analysis module.
*/
static void lex_test(void)
{
        int i;
        Table *fntst;
        Table *ttst;
	Equiv eq_dummy;
        char no_more=0,exit_flag=0;
	extern char debug;
        stack_ptr = 0;
        while (!exit_flag) {
/*
  Compile the next line.
*/
          if(stack_line(&fntst,no_more,-1,-1) == -1) {
/*
  On error, close all command files, returning control to stdin,
  and find the next newline, ignoring anything left on the
  current line.
*/
            flush_input();
          }
/*
  Check for END_BLOCK type commands.
*/
          else if(fntst->class == FUNC && TABFUNC(fntst)->sub_class == END_BLOCK) {
            lex_err(comline.last);
            lprintf(stderr,"Unmatched '%s' statement found.\n",fntst->name);
            flush_input();
          }
          else {
/*
  Dump the contents of the compile stack to the terminal if the user
  variable, debug is set to true.
*/
            i=stack_ptr;
            if(debug) {
              printf("This is the content of the compile stack.\n");
              for(stack_ptr=0;stack_ptr<=i;stack_ptr++) {
                ttst = compile_stack[stack_ptr];
                found_op_err(ttst);
              };
              printf("Debug dump complete\n\n");
              debug=0;
            };
/*
  Execute block.
*/
            in_run_mode=1;
            if(exe_control(0,i-1) == -1) 
              flush_input();
            in_run_mode=0;
          };  
/*
  Release the memory held in constants and special classes of
  table entry on the compile stack.
*/
          clear_compile(0,stack_ptr-1);
/*
  Clear the run stack and array stack, releasing memory where
  possible.
*/
          compress_temp(run_ptr, ' ', eq_dummy);
          array_zap(expr_ptr);
/*
  Clear the interrupt flags.
*/
          no_error=0;
        };
        return;
}

/*.......................................................................
  Empty the compile stack between the entry limits sent. This includes
  deallocation of memory held in constants, branch type instructions
  etc.. The compile stack pointer is left at start_pos.
*/
static void clear_compile(short int start_ptr, short int end_ptr)
{
        Table *ttst;
	Indexes *indval;
/*
  Step along the compile stack and deallocate memory where applicable.
*/
        for(stack_ptr=start_ptr;stack_ptr<=end_ptr;stack_ptr++) {
          ttst = compile_stack[stack_ptr];
          if(ttst == NULL) continue;
          switch (ttst->class) {
          case CONST:
            free_const(ttst);
            break;
          case FN_RET:
	    if(TABDESC(ttst)->access == STACK)
	      valof_free(TABDESC(ttst));
	    free(TABITEM(ttst));
            free(ttst);
            break;
          case DO_PAR:
            free(TABITEM(ttst));
            free(ttst);
            break;
          case START_EXPR:
            free(TABITEM(ttst));
            free(ttst);
            break;
          case INDEX_EXPR:
	    indval = TABINDX(ttst);
	    switch (indval->var->access) {
	    case FN_ARRAY_REF: case FN_ARRAY_VAL:
	      free(indval->var);
	    };
            free(indval);
            free(ttst);
            break;
          case BR_TRUE: case BR_FALSE: case BR_TO: case BR_VIA: case HASH:
          case END_LINK: case DO_INI: case SUB_STRING:
	  case ARRAY_PTR: case NUM_ARG: case FTOI: case ITOF:
	  case COMMAND: case DECL:
          case ADD_OP: case SUB_OP: case MUL_OP: case DIV_OP: case POW_OP:
          case GTE_OP: case GT_OP: case LT_OP: case LTE_OP: case EQ_OP:
          case NE_OP: case NO_OP: case IADD_OP: case ISUB_OP: case IMUL_OP:
	  case IDIV_OP: case IGTE_OP: case IGT_OP:
          case ILT_OP: case ILTE_OP: case IEQ_OP: case INE_OP: case SGTE_OP:
	  case SGT_OP: case SLT_OP: case SLTE_OP: case SEQ_OP:
          case SNE_OP: case CAT_OP: case NOT_OP: case AND_OP: case OR_OP:
          case MINUS_OP: case IMINUS_OP: case REG_OP: case NREG_OP:
	    free(ttst);
	    break;
          };
        };
/*
  Leave the stack pointer at the next newly freed stack position.
*/
        stack_ptr = start_ptr;
        return;
}

/*.......................................................................
  This routine generates a variable-table entry for a temporary 
  scalar constant such as a constant written in an expression by the
  user. The routine returns a pointer to the new table entry.
  The value to be stored is sent via a (void *) pointer. The type of the
  variable should be sent as 'c' for character, 'f' for float or 'l' for
  logical in (char vtype).
*/
Table *store_const(char vtype, void *value)
{
        static Table *ttst;
        static Descriptor *dtst;
        static long adim[3]={1,1,1};
/*
  Allocate space for a new table structure.
*/
        if( (ttst = table_alloc(CONST,NULL)) == NULL)
          return NULL;
/*
  Allocate and initialize the descriptor.
*/
        if( (dtst=descriptor_alloc(vtype,'0',adim)) == NULL) {
          free(ttst);
          return NULL;
        };
/*
  Install the descriptor in the new table entry and give it the STACK
  access type (ie. don't zap until the compile stack has been finished with).
*/
        TABITEM(ttst) = dtst;
        dtst->access = STACK;
/*
  Now handle the different variable type values.
*/
        switch (vtype) {
        case 'c':
/*
  Allocate enough memory for the string that was sent and for a
  terminating '\0'. Install the address of the result in the first
  element of the array of pointers allocated above. 
*/
          if( (*STRPTR(dtst) = stralloc(strlen((char *) value))) == NULL)
            return NULL;
/*
  Copy the character string to the pointer in the descriptor value field.
*/
          strcpy(*STRPTR(dtst), (char *) value);
          break;
/*
  Copy a float into the value field.
*/
        case 'f':
          *FLTPTR(dtst) = *((float *) value);
          break;
/*                                              
  Copy a short int into the value field.
*/
        case 'i':
          *INTPTR(dtst) = *((int *) value);  
          break;
/*
  Copy a logical (char) into the value field.
*/
        case 'l':
         *LOGPTR(dtst) = *((char *) value);
          break;
        default:
          lprintf(stderr,"syserr: unrecognised storage type in store_const\n");
          var_free(ttst);
          break;
        };
/*
  Return a pointer to the new table entry.
*/
        return ttst;
}

/*.......................................................................
  Delete the table structure and value of a temporary constant.
*/
void free_const(Table *stab)
{
/*
  Delete the contents.
*/
        var_free(stab);
/*
  Now delete the table structure.
*/
        free(stab);
}

/*.......................................................................
  A utility routine for determining if an integer number is a power
  of two.
*/
int is_pow_of_two(int inum) {
        static double expon;
	expon = log((double) inum)/log(2.0);
	return (fabs(expon-((int) expon)) < 0.001);
}

/*.......................................................................
  Return the size of one element of a given variable type, in units of
  characters.
*/
size_t mem_size_of(char vtyp)
{
	switch(vtyp) {
	case 'l':
	  return 1;
	  break;
	case 'i':
	  return sizeof(int);
	  break;
	case 'f':
	  return sizeof(float);
	  break;
	case 'c':
	  return sizeof(char *);
	  break;
	default:
	  return 0;
	};
}
