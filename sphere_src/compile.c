#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "sphere.h"
#include "run.h"
#include "table.h"
#include "ops.h"
#include "lex.h"
#include "help.h"
#include "logio.h"

Table *compile_stack[MAXSTACK];

static int new_declare(Table *fntst);
static int stack_assign(Table *ttst, Table **optst);
static int stack_array_indexes(Exprtype *expr_typ, Table *array_ptr,
			       Descriptor *vardsc);
static Table *getoperator(void);
static int stack_block(Table **fntst, short int start_pos, short int end_pos);
static int while_block(Table *fntst, short int end_pos);
static int repeat_block(Table *fntst, short int end_pos);
static int if_block(Table *fntst, short int link_pos, short int start_pos,
		    short int end_pos);
static int do_block(Table *fntst, short int end_pos);
static int stack_instruct(int class, int icode);
static int link_expr_typ(Exprtype *expr_typ, short int start_expr);
static int stack_operand(Exprtype *expr_typ);
static int get_expr(char **opline, Table **ret_tab, int prec,
		    Exprtype *expr_typ);
static int stack_var(Table *ttst, Exprtype *expr_typ);
static int sub_index(short int ind[], int j, char *term, short int *nargs,
		     char *dims);
static int string_index(Exprtype *expr_typ);
static int stack_function(Table *ttst, int stmt, Exprtype *expr_typ);
static int stack_command(Table *fntst);
static int function_args(Table *fntst, int stmt, Exprtype *expr_typ,int *nargs);
static int stack(Table *ttst, int sptr);
static int stack_expr(char **opline, Table **ret_tab, int prec, char atyp,
		      char dim, char access, Exprtype *expr_typ);
static int check_expression(Exprtype *expr_typ, char atyp, char dim,
			    char access, char is_elemental);
static int stack_args(char *name, short int nmin, short int nmax, char type[],
		      char dim[], char access[], char once, int *nargs,
		      Exprtype *rtntype, Table **optst);
static int stack_lit(char **tmplin, int once, Table **optst,
		     Exprtype *expr_typ);

/*.......................................................................
  Invoke the lexical analyser to return the next operator, operand, or
  function. If an operand or function is returned then signal an error.
  and return a null pointer to the calling routine.
  If an Operator is returned then in turn pass it back to the calling
  routine.
*/
static Table *getoperator(void)
{
        static Table *ttst;
/*
  Get the next keyword via the lexical analyser.
*/
        if( (ttst = lex_expr(' ')) == NULL)
          return NULL;
        if(ttst->class == OPER)
          return  ttst;
/*
  Any other class of table entry is illegal - signal the error.
*/
        lex_err(comline.last);
        found_op_err(ttst);
        lprintf(stderr,"Where an operator was expected\n");
/*
 * Release memory.
 */
	switch(ttst->class) {
	case CONST:
	  free_const(ttst);
	  break;
	case HASH:
	  free(ttst);
	  break;
	};
        return NULL;
}

/*.......................................................................
  Parse a command line. A command line may consist of an assignment
  expression or a command plus arguments. It may also consist of an
  indirection to a command file in the format: @filename. If an END
  or ELSE statement is detected control is returned to the calling
  routine and with a pointer to the respective function table entry in
  *retfunc. If is_more is passed as true (1) then the previous line has
  not been fully parsed yet and another command is thought to reside in
  the command line. This is used for the statement-if command where
  a second command follows the if() command.
*/
int stack_line(Table **fntst, char is_more, short int start_pos, short int end_pos)
{            
        Table *ttst,*ttmp;
        Table *optst;
        short int link_pos;
        static char *cptr;
/*
  Read the next command line.
*/
        if(!is_more && newline() == -1)
	  return -1;
/*
  Get the a variable for assignment or a command.
*/
        if( (ttst=lex_expr(' ')) == NULL)
          return -1;
        *fntst =  ttst;
/*
  Decide what to do, from the table entry returned.
*/
        switch (ttst->class) {
        case FUNC:
/*  
  Certain special commands require separate attention - these are the
  block structure commands such as IF, WHILE, etc..
*/
          switch (TABFUNC(*fntst)->sub_class) {
          case END_BLOCK:
            return no_error;
            break;
          case START_BLOCK:
/*
  Without foreknowledge of the end stack position of a block of
  commands within a loop, it is difficult to link stack commands
  to the exit of the loop. The "break" command is a good example.
  It requires a branch stack offset to exit the loop and this should be
  set up at compile time, ie. now.
    In order to allow this, without an extra linking pass, the BR_VIA
  command is used instead. This is linked with an END_LINK table entry,
  which is itself placed before the loop. When the end of the new loop
  has been sighted, the offset of the END_LINK command is computed
  and placed in the END_LINK Table.type field. At run-time, branches
  to exit the loop go via the END_LINK table entry. All links are made
  via stack position offsets, rather than using absolute stack positions.
  This allows later rellocation of the compiled stack.
    Start by creating the END_LINK table entry.
*/
            link_pos = stack_ptr;
            if(stack_instruct(END_LINK,0) == -1)
              return -1;
/*------------------------------------ Individual block commands -----
  Handle compilation of a while loop.
*/
            if(strcmp((*fntst)->name,"while") == 0) {
              if(while_block(*fntst,link_pos) == -1)
                return -1;
            }
/*
  Handle compilation of a repeat-until loop.
*/
            else if(strcmp((*fntst)->name,"repeat") == 0) {
              if(repeat_block(*fntst,link_pos) == -1)
                return -1;
            }
/*
  Handle compilation of a "if-elseif-else-end if" block.
*/
            else if(strcmp((*fntst)->name,"if") == 0) {
              if(if_block(*fntst, link_pos, start_pos, end_pos) == -1)
                return -1;
            }
/*
  Handle compilation of a do-loop block.
*/
            else if(strcmp((*fntst)->name,"do") == 0) {
              if(do_block(*fntst, link_pos) == -1)
                return -1;
            };
/*--------------------------------------------------------------------
  Link up the END_LINK entry with the offset to the current stack position.
*/
            TABICODE(compile_stack[link_pos]) = stack_ptr-link_pos-1;
            return no_error;
            break;
/*
  The break statement is equivalent to the C break statement - it
  needs to be linked to exit the current loop. Link it via the
  END_LINK table entry at stack position end_pos.
*/
          case BRK_BLOCK:
            if(end_pos == -1) {
              lex_err(comline.last);
              lprintf(stderr,"%s statement outside any loop\?\n",ttst->name);
              return -1;
            };
            if(stack_instruct(BR_VIA,end_pos-stack_ptr) == -1)
              return -1;
            break;
/*
  Similarly, the continue statement works like the C continue statement.
  Link it to the start of the loop.
*/
          case CONT_BLOCK:
            if(start_pos == -1) {
              lex_err(comline.last);
              lprintf(stderr,"%s statement outside any loop\?\n",ttst->name);
              return -1;
            };
            if(stack_instruct(BR_TO,start_pos-stack_ptr-1) == -1)
              return -1;
            break;
/*
  The ABORT stack-command aborts execution at run-time.
*/
          case STOP_EXE:
	    if(stack_instruct(ABORT,0) == -1)
	      return -1;
            break;
/*
  A variable declaration function.
*/
          case DECLARE:
            if(new_declare(ttst) == -1)
              return -1;
            break;
/*
  Handle the whatis commands.
*/
	  case WHATVAR: case HELP:
	    if(comline.nxtc == '\0')
	      cptr = NULL;
	    else {
	      if((ttmp=lex_expr('n')) == NULL)
		return -1;
	      cptr = ttmp->name;
/*
  Ignore further input from the current sub-line.
*/
	      while(*comline.next != '\0')
		comline.next++;
	      comline.nxtc = '\0';
	    };
/*
  Send the symbol to the appropriate function.
*/
	    switch (TABFUNC(*fntst)->sub_class) {
	    case WHATVAR:
	      whatisvar(cptr);
	      break;
	    case HELP:
	      help(cptr);
	      break;
	    };
	    break;
          default:
/*
 * Stack a normal command or function.
 */
	    if(stack_command(*fntst))
	      return -1;
	    return no_error;
	    break;
          };
/*
  The next character should be the end of line '\0' character.
*/
          if(comline.nxtc != '\0') {
            lex_err(comline.last);
            lprintf(stderr,"Unexpected characters at end of line.\n");
            return -1;
          };
          return no_error;
          break;
/*
  A variable name at the start of a line signals an assignment.
*/
        case VAR:
	  if(stack_assign(ttst,&optst) == -1)
	    return -1;
	  break;
/*
  Check for other symbols that shouldn't be at the start of a command line.
*/
	default:
	  lex_err(comline.last);
	  found_op_err(ttst);
	  lprintf(stderr,"Where a command or assignment expression was expected.\n");
	  switch(ttst->class) {
	  case CONST:
	    free_const(ttst);
	    break;
	  case HASH:
	    free(ttst);
	    break;
	  };
	  return -1;
	};
/*
  The terminating operator should be a '\0' line terminator.
*/
	if(TABOPER(optst)->op_prec != FINISH) {
	  lex_err(comline.last);                                    
	  lprintf(stderr,"Unexpected argument before end of line.\n");
	  return -1;
	};
	return no_error;
}

/*.......................................................................
  This function parses command lines until either an END <WHATEVER>
  statement or an ELSE(<WHATEVER>) function is found (or an error occurs).
  When one of these is encounterred the function returns to the calling
  routine where the block of commands is properly linked and the
  <WHATEVER> expression is checked and parsed. This routine is called
  whenever a block of commands is required, such as the body of a FOR,
  WHILE, DO-WHILE loop or that of an (IF ELSE(WHATEVER) ELSE END IF).
*/
static int stack_block(Table **fntst, short int start_pos, short int end_pos)
{
	char no_more=0;
/*
  stack_line() only changes fntst if a function is found. Since we
  decide on when to exit this routine by the value of *fntst, assign
  a NULL pointer to *fntst.
*/
	*fntst = NULL;
/*
  Increment the record of nested user blocks (This is used when prompting,
  in newline()).
*/
	comline.nest_block++;
/*
  Parse command lines.
*/
	do {
	  if(stack_line(fntst, no_more, start_pos, end_pos) == -1)
	    return -1;
	} while(*fntst==NULL || !((*fntst)->class == FUNC && TABFUNC(*fntst)->sub_class == END_BLOCK));
/*
  Now that the latest user command block is complete, decrement the nest
  count.
*/
	comline.nest_block--;
	return no_error;
}

/*.......................................................................
  This function handles the while statement.
*/
static int while_block(Table *fntst, short int end_pos)
{
	static int nargs;
	Table *ttst;
	short int start_pos, stack_posb;
/*
  Keep a record of the current stack position - this will be required
  for linking the bottom of the loop back to the start.
*/
	start_pos = stack_ptr;
/*
  Then stack the logical argument of the while function.
*/
	if(function_args(fntst, 0, NULL, &nargs) == -1)
	  return -1;
/*
  Create a table entry which will later contain the offset to the
  instruction after the end of the while loop. Keep a record of
  the position in which it is stacked in stack_posb. This is the
  instruction to exit the loop when the while condition becomes
  false.
*/
	stack_posb = stack_ptr;
	if(stack_instruct(BR_FALSE,0) == -1)
	  return -1;
/*
  Now get the block of commands that make up the body of the while loop.
*/
	if(stack_block(&fntst,start_pos,end_pos) == -1)
	  return -1;
/*
  Test the end block function that terminated the command block.
*/
	if(strcmp(fntst->name,"end") != 0) {
	  lex_err(comline.last);
	  lprintf(stderr,"Unmatched %s statement in while loop.\n",fntst->name);
	  return -1;
	};
/*
  Check the next token is the function name "while".
*/
	if((ttst = lex_expr(' ')) == NULL)
	  return -1;
	if(strcmp(ttst->name,"while") != 0 || comline.nxtc != '\0') {
	  lex_err(comline.last);
	  lprintf(stderr,"Was expecting END WHILE;.\n");
	  return -1;
	};
/*
  Create a branch stack command to be loaded with a pointer to the
  start of the loop. Make it point to the start of the loop.
*/
	if(stack_instruct(BR_TO,start_pos-stack_ptr-1) == -1)
	  return -1;
/*
  Link up the BR_FALSE stack commands to point to the next stack posn.
*/
	TABICODE(compile_stack[stack_posb]) = stack_ptr-stack_posb-1;
	return no_error;
}

/*.......................................................................
  This function handles compilation of the repeat-until() statement.
*/
static int repeat_block(Table *fntst, short int end_pos)
{
	static int nargs;
	short int start_pos;
/*
  Keep a record of the current stack position - this will be required
  for linking the bottom of the loop back to the start.
*/
	start_pos = stack_ptr;
/*
  Check that there is no argument to the repeat function.
*/
	if(function_args(fntst, 0, NULL, &nargs) == -1)
	  return -1;
/*
  Now get the block of commands that make up the body of the while loop.
*/
	if(stack_block(&fntst,start_pos,end_pos) == -1)
	  return -1;
/*
  Test the end block function that terminated the command block.
*/
	if(strcmp(fntst->name,"until") != 0) {
	  lex_err(comline.last);
	  lprintf(stderr,"Unmatched %s statement in repeat loop.\n",fntst->name);
	  return -1;
	};
/*
  Stack the logical argument.
*/
	if(function_args(fntst, 0, NULL, &nargs) == -1)
	  return -1;
/*
  Create a branch stack command to be loaded with a pointer to the
  start of the loop.  Make it point to the start of the loop.
*/
	if(stack_instruct(BR_FALSE,start_pos-stack_ptr-1) == -1)
	  return -1;
	return no_error;
}

/*.......................................................................
  Handle compilation of the if_elseif..._else_end_if command.
*/
static int if_block(Table *fntst, short int link_pos, short int start_pos,
		    short int end_pos)
{
	static int nargs;
	Table *ttst;
	char is_more=1;
	short int is_if,is_end,is_else,is_elseif,was_else;
	short new_link=0;
/*
  Now get the block of commands that make up the conditionally executed
  sub-program.
*/
	was_else = is_end = 0;
	is_if = 1;
	for (;;) {
/*
  Work out which type of block function delimits the latest command block.
*/
	  is_elseif = (strcmp(fntst->name,"elseif") == 0);
	  is_else = (strcmp(fntst->name,"else") == 0);
/*
  If there was a previous command block terminated with an else
  statement, then was_else will be true. If so, only an "end if" command
  can follow the current block. Also check that the block delimiter is
  one of "else", "elseif", or "end if". is_if is 1 whenever if_block()
  is invoked. Reset it to 0 after the test.
*/
	  if( (was_else && !is_end) || !(is_if || is_else || is_elseif || is_end)) {
	    lex_err(comline.last);
	    lprintf(stderr,"Unmatched %s statement within if block.\n",fntst->name);
	    return -1;
	  };
	  is_if = 0;
/*
  The was_else flag is used to signal the past presence of an else
  statement. Set it here. Once set it should be made to stay set.
*/
	  was_else = was_else || is_else;
/*
  Check and stack arguments of the delimiting function.
*/
	  if(function_args(fntst, 0, NULL, &nargs) == -1)
	    return -1;
/*
  Terminate the for(;;) loop if the current command block was terminated with
  an "end if".
*/
	  if(is_end) break;
/*
  In order to be able to skip the next command block if if,elseif or else
  argument is false - set up a branch if false stack command. The stack
  offset to be skipped will be set when the end of the next block has been
  located. Don't do this if the latest end-block statement is an
  else statement since it has no such argument.
*/
	  if(!is_else) {
	    new_link = stack_ptr;
	    if(stack_instruct(BR_FALSE,0) == -1)
	      return -1;
	  };
/*
  If the next character on the command line is not the end of line
  character '\0', then assume that the user wants a statement-if,
  -elseif, or -else command. ie. if() command; as opposed to
  if();commands;end if. Get the command following the if command.
*/
	  if(comline.nxtc != '\0') {
	    if( stack_line(&fntst, is_more, start_pos, end_pos) == -1)
	      return -1;
/*
  Check to make sure that the last command found in stack_line was
  not an END_BLOCK type command. This prevents users from starting a new
  block inside the statement if command - a syntatically bad technique.
*/
	    if(fntst->class == FUNC && TABFUNC(fntst)->sub_class == END_BLOCK) {
	      lex_err(comline.last);
	      lprintf(stderr,"Unexpected '%s' statement detected in statement-if\n",fntst->name);
	      return -1;
	    };
/*
  Link up the branch if false stack entry to the current position.
*/
	    if(!was_else)
	      TABICODE(compile_stack[new_link]) = stack_ptr-new_link-1;
	    return no_error;
	  };
/*
  Get the next command block.
*/
	  if(stack_block(&fntst,start_pos,end_pos) == -1)
	    return -1;
/*
  Before getting the next elseif etc.. command, place a branch to
  exit loop command here. At run time, should the current command block
  be executed then this command will be encounterred before the
  next elseif etc... We don't need to bother with this if the current
  command block was terminated with an "end if".
*/
	  if( !(is_end=(strcmp(fntst->name,"end") == 0)) ) {
	    if(stack_instruct(BR_VIA,link_pos-stack_ptr) == -1)
	      return -1;
	  };
/*
  Now link up the last branch-if-false stack command to point to the
  current stack position.
*/
	  if(!(was_else && is_end))
	    TABICODE(compile_stack[new_link]) = stack_ptr-new_link-1;
/*
  Repeat the for(;;) loop for the next else,elseif or "end if" block.
*/
	};
/*
  An end function has been found - check that it "if", to form "end if".
*/
	if((ttst = lex_expr(' ')) == NULL)
	  return -1;
	if(strcmp(ttst->name,"if") != 0 || comline.nxtc != '\0') {
	  lex_err(comline.last);
	  lprintf(stderr,"Was expecting END IF;.\n");
	  return -1;
	};
/*
  If the last block had a branch-if-false stack command, link it to point
  to the current stack position.
*/
	if(!was_else)
	  TABICODE(compile_stack[new_link]) = stack_ptr-new_link-1;
	return no_error;
}

/*.......................................................................
  This function handles the do-end do statement.
*/
static int do_block(Table *fntst, short int end_pos)
{
	Exprtype expr_typ;
	Table *vtst;
	Table *optst;
	Table *ttst;
	char *tmplin;
	char atyp;
	Do_pars *dptst;
	int i,start_pos,narg=0;
/*
  The next token should be the DO-variable name.
*/
	if( (vtst=  lex_expr(' ')) == NULL)
	  return -1;
/*
  Check both that it is a scalar numeric variable.
*/
	atyp = TABDESC(vtst)->atyp;
	if(vtst->class != VAR || (atyp != 'f' && atyp != 'i') ||
	   TABDESC(vtst)->dim != '0') {
	  lex_err(comline.last);
	  lprintf(stderr,"DO requires a scalar numeric variable.\n");
	  return -1;
	};
/*
  Make sure that the variable is not a read-only parameter.
*/
	if(TABDESC(vtst)->access == R_ONLY) {
	  lex_err(comline.last);
	  lprintf(stderr,"Variable \"%s\" is a read-only parameter.\n",vtst->name);
	  return -1;
	};
/*
  The next operator should be an equals sign.
*/
	if(comline.nxtc != '=' || (optst = getoperator()) == NULL ||
	   (TABOPER(optst)->op_prec != EQUALS)) {
	  lex_err(comline.last);
	  lprintf(stderr,"Missing '=' in DO variable assignment\n");
	  return -1;
	};
/*
  Attempt to get up to 3 scalar, numeric arguments, separated by
  commas.
*/
	for(i=0;i < 3; i++) {
/*
  Stack the next argument expression.
*/
	  if( stack_expr(&tmplin, &optst, OP_BR, atyp, '0', 'v', &expr_typ) == -1)
	    return -1;
/*
  Check the delimiter is a comma or semi-colon.
*/
	  if(TABOPER(optst)->op_prec != COMMA && TABOPER(optst)->op_prec != FINISH) {
	    lex_err(comline.last);
	    lprintf(stderr,"Unexpected \"%s\" operator delimiting a DO parameter.\n",optst->name);
	    return -1;
	  };
/*
  Increment the argument count.
*/
	  narg++;
/*
  Stop when the '\0' character is encounterred.
*/
	  if(TABOPER(optst)->op_prec == FINISH) break;
	};
/*
  Make sure that at least two DO parameters were parsed.
*/
	if(narg < 2) {
	  lex_err(comline.last);                                    
	  lprintf(stderr,"DO requires at least 2 iteration parameters.\n");
	  return -1;
	};
/*
  The terminating operator should be a '\0' line terminator.
*/
	if(TABOPER(optst)->op_prec != FINISH) {
	  lex_err(comline.last);                                    
	  lprintf(stderr,"Unexpected argument before end of line.\n");
	  return -1;
	};
/*
  All syntax checks complete - now stack a DO_INI table entry and store the
  number of numeric expressions (narg) that it should take off the stack.
*/
	switch (atyp) {
	case 'f':
	  if(stack_instruct(DO_INI,narg) == -1)
	    return -1;
	  break;
	case 'i':
	  if(stack_instruct(IDO_INI,narg) == -1)
	    return -1;
	  break;
	};
/*
  Now create a DO loop parameter maintainance block. This has a
  type field that points to a Do_pars structure. The Do_pars
  structure has fields containing the stack offset to the exit of the
  loop, a count of the number of iterations completed, a record of the
  DO-variable start value, end value and step value and a pointer to
  the DO variable. At run-time, When the (above) DO_INI stack command
  is encounterred, the iteration parameters are set using the
  numeric expressions preceding it on the stack. The exit link field
  will be set shortly, as soon as the end instruction of the loop is
  sighted.
*/
	if( (dptst = (Do_pars *) malloc(sizeof(Do_pars))) == NULL) {
	  lprintf(stderr,"Memory allocation failed.\n");
	  return -1;
	};
	dptst->skipend = 0;
	if((ttst = table_alloc(DO_PAR,NULL)) == NULL) {
	  free(dptst);
	  return -1;
	};
	TABITEM(ttst) = dptst;
	VOIDPTR(TABDOPAR(ttst)) = VOIDPTR(TABDESC(vtst));
	switch (atyp) {
	case 'f':
	  ttst->class = DO_PAR;
	  break;
	case 'i':
	  ttst->class = IDO_PAR;
	  break;
	};
/*
  Keep a record of the current stack pointer before stacking the
  DO_PAR table entry. It will be used to link the loop end to here
  and vice versa.
*/
	start_pos = stack_ptr;
	if(stack(ttst,stack_ptr++) == -1) 
	  return -1;
/*
  Now get the block of commands that make up the body of the DO loop.
*/
	if(stack_block(&fntst,start_pos,end_pos) == -1)
	  return -1;
/*
  Test the end block function that terminated the command block.
*/
	if(strcmp(fntst->name,"end") != 0) {
	  lex_err(comline.last);
	  lprintf(stderr,"Unmatched %s statement in DO loop.\n",fntst->name);
	  return -1;
	};
/*
  Check the next token is the function name "DO".
*/
	if((ttst = lex_expr(' ')) == NULL)
	  return -1;
	if(strcmp(ttst->name,"do") != 0 || comline.nxtc != '\0') {
	  lex_err(comline.last);
	  lprintf(stderr,"Was expecting END DO;.\n");
	  return -1;
	};
/*
  Create a branch stack command to be loaded with a pointer to the
  start of the loop.  Make it point to the start of the loop.
*/
	if(stack_instruct(BR_TO,start_pos-stack_ptr-1) == -1)
	  return -1;
/*
  Link up the DO_PAR stack command to point to the next stack posn.
*/
	TABDOPAR(compile_stack[start_pos])->skipend = stack_ptr-start_pos-1;
	return no_error;
}

/*.......................................................................
  Create a special stack instruction command and place it at
  the current stack pointer. This is a table entry with a null name
  field, a class sent by the calling function eg. BR_TRUE, BR_FALSE or
  BR_TO, and a type field pointing to a a specified number of
  short int's. The latter is used to store things like stack branch
  offsets and should assigned by the calling routine - the memory for
  the (int's) is allocated herein.
*/
static int stack_instruct(int class, int icode)
{
	static Table *ttst;
/*
  Allocate memory for the table entry.
*/
	if( (ttst=table_alloc(class,NULL)) == NULL)
	  return -1;
/*
  Assign any instruction parameter sent, to the type field.
*/
	TABICODE(ttst) = icode;
/*
  Stack it at the current position.
*/
	if(stack(ttst,stack_ptr++) == -1)
	  return -1;
	return no_error;
}

/*.......................................................................
  Create a special stack instruction command and place it at
  the current stack pointer. This is a table entry with a null name
  field, a class sent by the calling function eg. BR_TRUE, BR_FALSE or
  BR_TO, and a type field pointing to a a specified number of
  short int's. The latter is used to store things like stack branch
  offsets and should assigned by the calling routine - the memory for
  the (int's) is allocated herein.
*/
static int link_expr_typ(Exprtype *expr_typ, short int start_expr)
{
	static Table *ttst;
/*
  Get a pointer to the stack entry to be updated.
*/
	ttst = compile_stack[start_expr];
/*
  Allocate memory for the Exprtype entry.
*/
	if( (TABITEM(ttst) = malloc(sizeof(Exprtype))) == NULL) {
	  free(ttst);
	  lprintf(stderr,"Memory allocation failure.\n");
	  return -1;
	};
/*
  Assign the sent Exptype version to the type field.
*/
	*TABEXPR(ttst) = *expr_typ;
/*
  Also assign the link stack offset to its length field.
*/
	TABEXPR(ttst)->length = stack_ptr-start_expr-1;
	return no_error;
}


/*.......................................................................
  Given any table entry, report (on stderr) what the table entry
  contains, including the values of constants.
*/
void found_op_err(Table *ttst)
{
	float *fval;
	char **cval, *lval;
	extern char debug;
	int i,adim,*ival;
	short int sival,sidum;
	Indexes *indval;
	Exprtype *expr_typ;
/*
  In debug mode print out the stack pointer from which ttst was
  extracted.
*/
	if(debug)
	  lprintf(stderr,"%d ",stack_ptr);
	if(ttst == NULL) {
	  lprintf(stderr,"Found null table pointer entry\n");
	  return;
	};
	switch (ttst->class) {             
	case VAR:
	  lprintf(stderr,"Found variable: %s\n",ttst->name);
	  break;
	case FUNC:
	  lprintf(stderr,"Found Function: %s()\n",ttst->name);
	  break;
	case OPER:
	  lprintf(stderr,"Found operator: %s\n",ttst->name);
	  break;
	case MODULE_SYM:
	  lprintf(stderr,"Found module help topic: %s\n",ttst->name);
	  break;
	case HELP_SYM:
	  lprintf(stderr,"Found help topic: %s\n",ttst->name);
	  break;
	case CONST:
	  lprintf(stderr,"Found constant: ");
/*
  Get type dependant versions of the value pointer and get the 1st
  dimension of the constant.
*/
	  fval = FLTPTR(TABDESC(ttst));
	  ival = INTPTR(TABDESC(ttst));
	  lval = LOGPTR(TABDESC(ttst));
	  cval = STRPTR(TABDESC(ttst));
	  adim = TABDESC(ttst)->adim[0];
/*
  Output the values to the terminal.
*/           
	  for(i=0;i < adim;i++) {
	    switch (TABDESC(ttst)->atyp) {
	    case 'f':
	      lprintf(stderr,"%.3f ", *(fval++));
	       break;
	    case 'i':
	      lprintf(stderr,"%3.3d ", *(ival++));
	      break;
	    case 'l':
	      lprintf(stderr,"%s ", *(lval++) ? "TRUE":"FALSE" );
	      break;
	    case 'c':
	      lprintf(stderr,"\"%s\"",*(cval++));
	      break;
	    };
	  };
	  lprintf(stderr,"\n");
	  break;
	case BR_TRUE: case BR_FALSE: case BR_TO: case BR_VIA: case END_LINK:
	case DO_INI: case ITOF: case FTOI: case SUB_STRING: case HASH:
	case NUM_ARG:
	  sival = TABICODE(ttst);
	  switch (ttst->class) {
	  case BR_TRUE:
	    lprintf(stderr,"Skip %d entries if true.\n", sival);
	    break;
	  case SUB_STRING:
	    lprintf(stderr,"Substring with '%d' index specifier.\n", sival);
	    break;
	  case BR_FALSE:
	    lprintf(stderr,"Skip %d entries if false.\n", sival);
	    break;
	  case BR_TO:
	    lprintf(stderr,"Skip %d entries.\n", sival);
	    break;
	  case BR_VIA:
	    if(sival == 0) break;
	    sidum = TABICODE(compile_stack[stack_ptr+sival]);
	    lprintf(stderr,"(wrt end_link @ %d) skip %d entries.\n", stack_ptr+sival, sival+sidum);
	    break;
	  case END_LINK:
	    lprintf(stderr,"end_link entry with skip %d to exit.\n",sival);
	    break;
	  case DO_INI:
	    lprintf(stderr,"Initialize DO loop, from %d arguments.\n",sival);
	    break;
	  case HASH:
	    lprintf(stderr, "#%d array index pseudo variable.\n",sival);
	    break;
	  case NUM_ARG:
	    lprintf(stderr, "The following function has %d arguments.\n",sival);
	    break;
	  case ITOF:
	    lprintf(stderr, "int->float argument %d arguments back\n",sival);
	    break;
	  case FTOI:
	    lprintf(stderr, "float->int argument %d arguments back\n",sival);
	    break;
	  };
	  break;
	case DO_PAR:
	  sival = TABDOPAR(ttst)->skipend;
	  if(sival == 0) break;
	  lprintf(stderr,"DO_PAR entry, including branch offset of %d entries.\n",sival);
	  break;
	case COMMAND:
	  lprintf(stderr, "The arguments of a command start here\n");
	  break;
	case FN_RET:
	  lprintf(stderr,"Entry for once-only function return.\n");
	  break;
	case INDEX_EXPR:
	  indval = TABINDX(ttst);
	  lprintf(stderr,"Array variable, type='%c', User indices in (%d) args:",indval->var->atyp,indval->nargs);
	  for(i=0;i<3;i++)
	    lprintf(stderr,"%d %d %d  ",indval->start[i],indval->end[i],indval->inc[i]);
	  lprintf(stderr,"\n");
	  break;
	case ARRAY_PTR:
	  lprintf(stderr,"Elemental pointer into: %s\n",ttst->name);
	  break;
	case ABORT:
	  lprintf(stderr,"Abort execution now.\n");
	  break;
	case START_EXPR:
	  expr_typ = TABEXPR(ttst);
	  lprintf(stderr,"Expression: spans %d entries: type: '%c', access='%c'.\n",
	    (int) expr_typ->length, expr_typ->type, expr_typ->access);
	  break;
	default:
	  lprintf(stderr,"syserror in found_op_err()\n");
	  break;
	};
	return;
}

/*.......................................................................
  Invoke the lexical analyser to compile the next operand.
  If an operator is found next then signal an error, and return
  -1 to the calling routine.
     An operand may consist of an optional unary operator followed
  by a constant, variable name (with optional index specifier), a
  function reference followed by a bracketed sequence of arguments,
  or a parenthetised expression.
*/                                                     
static int stack_operand(Exprtype *expr_typ)
{
	static int op_code;
	Table *ttst;
	char *tmplin;
	extern Table unminop;
/*
  Get the next keyword via the lexical analyser.
*/
	if( (ttst = lex_expr(' ')) == NULL)
	  return -1;
/*
  The only allowed operator types are unary operators and open parentheses.
*/
	switch (ttst->class) {
	case OPER:
	  switch( TABOPER(ttst)->op_prec) {
/*
  When a UNARY operator is encounterred, stack the following operand
  onto the run stack before stacking the unary operator. 
*/
	  case UNARY: case ADD:
	    tmplin = comline.last;
	    if( stack_operand(expr_typ) == -1)
	      return -1;
/*
  Check for the special case of unary + and unary -. The lexical analyser
  can not distinguish these from the binary operators, + and -.
  The unary minus operator has its own table entry (unminop) separate
  from the operator symbol table.
*/
	    switch (*ttst->name) {
	    case '+':
	      return no_error;
	    case '-':
	      ttst = &unminop;
	    };
/*
  Check that the operator supports an operand of the type found
  above.
*/
	    switch (expr_typ->type) {
	    case 'f':
	      op_code = TABOPER(ttst)->f_op;
	      break;
	    case 'i':
	      op_code = TABOPER(ttst)->i_op;
	      break;
	    case 'c':
	      op_code = TABOPER(ttst)->s_op;
	      break;
	    case 'l':
	      op_code = TABOPER(ttst)->l_op;
	      break;
	    };
/*
  Check that the operator has an instruction to deal with the type of
  operand given.
*/
	    if(op_code == NO_OP) {
	      lex_err(tmplin);
	      lprintf(stderr,"Illegal %s operand given to the unary %s operator.\n",
	       type_string(expr_typ->type), ttst->name);
	      return -1;
	    };
/*
  Stack the instruction.
*/
	    if(stack_instruct(op_code,0) == -1)
	      return -1;
	    break;
/*
  Round open bracket - the operand is an expression.
*/
	  case OP_BR:
	    if(*ttst->name != '(') {
	      lex_err(comline.last);
	      found_op_err(ttst);
	      lprintf(stderr,"Where an operand was expected\n");
	      return -1;
	    };
	    if(get_expr(&tmplin, &ttst, OP_BR, expr_typ) == -1)
	      return -1;
/*
  Check that the next symbol is a closing round bracket and signal
  an error if not.
*/
	    if( *ttst->name != ')') {
	      lex_err(comline.last);
	      found_op_err(ttst);
	      lprintf(stderr,"Where a ')' was expected\n");
	      return -1;
	    };
	    break;
/*
  Illegal operator encounterred. 
*/
	  default:
	    lex_err(comline.last);
	    found_op_err(ttst);
	    lprintf(stderr,"Where an operand was expected\n");
	    return -1;
	  };
/*
  Each of the above operand expressions was not formed of a single
  variable name - thus the operand type has the value access class.
*/
	  expr_typ->access = 'v';
	  break;
/*
  An variable name or constant followed by an optional () array
  index expression.
*/
	case VAR: case CONST:
	  if(stack_var(  ttst, expr_typ) != 0)
	    return -1;
	  break;
/*
  A #N array index pseudo variable.
*/
	case HASH:
	  expr_typ->type = 'i';
	  expr_typ->dim = '0';
	  expr_typ->access = 'v';
	  stack(ttst, stack_ptr++);
	  break;
/*
  A function name followed by (arguments...)?
*/
	case FUNC:
/*                                                   
  Make sure that the function returns scalar only.
*/
	  if(*TABFUNC(ttst)->type == ' ' ||
	     TABFUNC(ttst)->sub_class != NORM) {
	    lex_err(comline.last);
	    lprintf(stderr,"Illegal placement of command '%s' in the middle of an expression.\n",ttst->name);
	    return -1;
	  };
/*
  Stack the function and its arguments.
*/
	  if(stack_function(ttst, 0, expr_typ) == -1)
	    return -1;
	  break;
	case MODULE_SYM: case HELP_SYM:
	  lex_err(comline.last);
	  lprintf(stderr,"Illegal use of help topic '%s' as an operand\n",ttst->name); 
	  return -1;
	default:
	  lprintf(stderr,"Syserror: Unknown class in stack_operand.\n");
	  return -1;
	  break;
	};
/*
  Handle the sub-string operator - this is effectively an
  operator with up to 4 arguments, and needs special treatment.
*/
	if(expr_typ->type == 'c' && comline.nxtc == '[' && string_index(expr_typ) == -1)
	  return -1;
	return no_error;
}     

/*.......................................................................
  Get an arithmetic expression of unknown type (float or character) and
  return its type. An arithmetic expression is terminated by a close
  square or round bracket, a colon, a comma, or end of line. The actual
  (single character) terminator is returned via (char *term). An error
  condition is signalled by a -1 return value.
*/
static int get_expr(char **opline, Table **ret_tab, int prec,
		    Exprtype *expr_typ)
{
	Table *optst;
	int new_prec;
	short int stack_pos;
	char *newopline;
	Exprtype atyp,btyp;
	static int op_code;
/*
  An expression must start with an operand - stack the operand and
  keep a record of its type in atyp.
*/
	if( stack_operand(&atyp) == -1)
	  return -1;
/*
  This must be followed by an operator - get it.
*/
	if((optst = getoperator()) == NULL)
	  return -1;
/*
  An open bracket is illegal at this point.
*/
	if(TABOPER(optst)->op_prec == OP_BR) {
	  lex_err(comline.last);
	  found_op_err(optst);
	  lprintf(stderr,"Open bracket adjacent to operand.\n");
	  return -1;
	};
/*
  Keep a record of (Comline line) for future error reporting about
  this operator.
*/
	*opline = comline.last;                               
/*
  Test the precedence of the current binary operator against that with
  which this routine was called. If the new precedence is equal or lower
  then return to the calling routine.
*/
	for(;;) {
	  if((new_prec = TABOPER(optst)->op_prec) <= prec) {
	    *ret_tab = optst;
	    *expr_typ = atyp;       
	    return no_error;
	  }
	  else {
/*
  Otherwise, the new precedence is higher than the previous one
  so recursively call this routine to find the next operator with lower
  precedence before stacking the current operator. First of all check if
  the new operator is a logical AND or OR.  These are unique as binary
  operators since half of the time it is possible to tell the result
  from just the first operand, and it is usually assumed that the second
  expression won't be evaluated when this is the case - eg. the statement
  if(a == 0 & b/a == 1) would result in a divide by zero error if
  the second operand was (unnecessarily) evaluated when a=0. Thus
  it is necessary to conditionally be able to skip the second
  expression. In order to do this a special branch table entry should be
  created at this point, its stack position stored in stack_pos and
  after the new expression has been stacked via get_expr(), 
  the type field of the table entry will be assigned with the number
  of stack positions up to the end of the new expression.
*/
	    stack_pos = stack_ptr;
	    switch (TABOPER(optst)->op_prec) {
	    case AND:
	      if(stack_instruct(BR_FALSE,0) == -1)
		return -1;
	      break;
	    case OR: 
	      if(stack_instruct(BR_TRUE,0) == -1)
		return -1;
	      break;
	    };
/*
  Make sure that the current operator is binary before getting the
  next operand.
*/
	    if(TABOPER(optst)->narg != 2) {
	      lex_err(comline.last);
	      lprintf(stderr,"Illegal placement of '%s' operator within an operand.\n",optst->name);
	      return -1;
	    };
/*
  Now get the next operand expression.
*/
	    if(get_expr(&newopline, ret_tab, new_prec, &btyp) == -1)
	      return -1;
/*
  If one operand is a float and the other an int, stack an instruction
  to convert the int to a float.
*/
	    if(atyp.type == 'f' && btyp.type == 'i') {
	      if(stack_instruct(ITOF,0) == -1)
		return -1;
	      btyp.type = 'f';
	    }
	    else if(atyp.type == 'i' && btyp.type == 'f') {
	      if(stack_instruct(ITOF,1) == -1)
		return -1;
	      atyp.type = 'f';
	    }
/*
  It is assumed that binary operators require two
  operands of equal type. (We have already checked that the
  first operand has the correct type, above.)
*/
	    else if(btyp.type != atyp.type) {
	      lex_err(*opline);
	      lprintf(stderr,"Incompatible operands given: (%s) %s (%s)\n",
	       type_string(atyp.type), optst->name, type_string(btyp.type));
	      return -1;
	    };
/*
  Check the operand type against that specified in optst.
  '*' is a wild card that allows for any type of operand.
*/
	    switch (atyp.type) {
	    case 'f':
	      op_code = TABOPER(optst)->f_op;
	      break;
	    case 'i':
	      op_code = TABOPER(optst)->i_op;
	      break;
	    case 'c':
	      op_code = TABOPER(optst)->s_op;
	      break;
	    case 'l':
	      op_code = TABOPER(optst)->l_op;
	      break;
	    };
/*
  Check whether the operator has an instruction to handle the type
  of argument given.
*/
	    if(op_code == NO_OP) {
	      lex_err(*opline);
	      lprintf(stderr,"Illegal %s operand given to the binary %s operator.\n",
	       type_string(atyp.type), optst->name);
	      return -1;
	    };
/*
  Stack the instruction.
*/
	    if(stack_instruct(op_code,0) == -1)
	      return -1;
/*
  If the types are ok then promote the return type of the operator in optst
  to be the current operand type. Also carry over the current dimensional
  type. This should be the highest dimensional type of the two operands.
  Also, since we have found a binary expression, as opposed
  to a single variable, the access class is by value, ie. 'v'.
*/
	    if(TABOPER(optst)->atyp != '*')
	      atyp.type = TABOPER(optst)->atyp;
	    atyp.dim = (atyp.dim > btyp.dim) ? atyp.dim : btyp.dim;
	    atyp.access = 'v';
/*
  If the operator that was just stacked happened to be a logical AND
  or OR then link up the branch table entry created above to point
  to the current stack position.
*/
	    switch (TABOPER(optst)->op_prec) {
	    case AND: case OR:
	      TABICODE(compile_stack[stack_pos]) = stack_ptr-stack_pos-1;
	    };
/*
  Copy the returned operator and its line pointers to the current values
  in this routine.
*/
	    *opline = newopline;
	    optst = *ret_tab;
	  };
	};
}

/*.......................................................................
  Stack a constant or a variable with any user provided index expressions.
*/
static int stack_var(Table *ttst, Exprtype *expr_typ)
{
        Descriptor *vardsc, *tmpdsc;
	Table *tabtst;
	int skip_from;
/*
  Get the return expression type.
*/
	vardsc = TABDESC(ttst);
	expr_typ->type   = vardsc->atyp;
	expr_typ->dim    = vardsc->dim;
	expr_typ->access = (TABDESC(ttst)->access != R_ONLY) ? 'r':'V';
/*
  If the table entry sent is that of a constant - make the access
  type by value and return - constants are scalar and can't take
  indexes.
*/
	if(ttst->class == CONST) {
	  if(stack(ttst,stack_ptr++) == -1)
	    return -1;
	  expr_typ->access = 'V';
	  return no_error;
	};
/*
  A scalar variable - simply stack the variable descriptor.
*/
	if(vardsc->dim == '0') {
	  if(stack(ttst,stack_ptr++) == -1)
	    return -1;
	  expr_typ->access = 'N';
	}
/*
  If the variable is an array then create and stack a new descriptor
  that, during expression evaluation, will point to the latest element
  of the array.
*/
	else {
/*
  Get a table structure to point to the new descriptor.
*/
	  if( (tabtst=table_alloc(ARRAY_PTR,NULL)) == NULL)
	    return -1;
/*
  Create and stack the new descriptor.
*/
	  tmpdsc = (Descriptor *) malloc(sizeof(Descriptor));
	  if(tmpdsc == NULL) {
	    free(tabtst);
	    lprintf(stderr,"Memory allocation failed.\n");
	    return -1;
	  };
	  *tmpdsc = *vardsc;
	  VOIDPTR(tmpdsc) = NULL;
	  tmpdsc->num_el = 0;
	  tmpdsc->access = STACK;
	  TABITEM(tabtst) = tmpdsc;
	  tabtst->name = ttst->name;
	  if(stack(tabtst,stack_ptr++) == -1) {
	    free(tmpdsc);
	    free(tabtst);
	    return -1;
	  };
/*
  Stack a BR_TO instruction and keep a record of its stack position.
  This will point to the instruction following any user indexes etc...
  in order that the expression evaluation function will
  skip them after the indexes have been evaluated once.
*/
	  skip_from = stack_ptr;
	  if(stack_instruct(BR_TO,0) == -1)
	    return -1;
/*
  Stack any user indexes and a record of the variable descriptor.
*/
	  if(stack_array_indexes(expr_typ, tabtst, vardsc) == -1)
	    return -1;
/*
  Link up the BR_TO instruction to point to the next instruction.
*/
	  TABICODE(compile_stack[skip_from]) = stack_ptr-skip_from-1;
	};
	return no_error;
}


/*.......................................................................
  Stack user specified array index expressions of array variables and
  array function return values.
*/
static int stack_array_indexes(Exprtype *expr_typ, Table *array_ptr, Descriptor *vardsc)
{
	int i,j;
	short int nargs=0,indexes[9]={0,0,0,0,0,0,0,0,0};
	char term;
	Table *tabtst;
	char dims,ndim, first_dim;
	Indexes *indval;
/*
  Create the structure in which the index specifications are
  to be stored.
*/
	if( (indval = (Indexes *) malloc(sizeof(Indexes))) == NULL) {
	  lprintf(stderr,"Memory allocation failed.\n");
	  return -1;
	};
/*
  Make the indval descriptor field point the descriptor of the variable.
*/
	indval->var = vardsc;
/*
  Make the ptr_to_elem_ptr point to the valof field of the element
  pointer descriptor.
*/
	indval->ptr_to_elem_ptr = (char **) &VOIDPTR(TABDESC(array_ptr));
/*
  Create a table entry to hold the Indexes structure allocated above.
*/
	if( (tabtst = table_alloc(INDEX_EXPR,NULL)) == NULL) {
	  free(indval);
	  return -1;
	};
	TABITEM(tabtst) = indval;
	tabtst->name = array_ptr->name;
/*
  Stack it.
*/
	if(stack(tabtst,stack_ptr++) == -1) {
	  free(indval);
	  free(tabtst);
	  return -1;
	};
/*
  Ascertain the intrinsic dimension of the variable.
*/
	first_dim = ndim = dims = vardsc->dim - '0';
/*
  See if the next character is a '('. If so skip the '(' operator and
  parse the index expression that should follow it.
*/
	switch (comline.nxtc) {
	case '(':
	  getoperator();
/*
  Loop for each dimension's index expression separately.
*/
	  for(i=0,j=0; j < ndim; i += 3,j++) {
	    if( sub_index(indexes,i,&term,&nargs,&dims) == -1)
	      return -1;
/*
  Record the first dimension in which the user specified an index.
*/
	    if(first_dim == ndim && nargs != 0) first_dim = j;
/*
  Stop early if the termination of the expression is encounterred.
*/
	    if( term == ')')
	      break;
	  };
/*
  Check that the loop terminated because a ')' was found.
*/
	  if( term != ')') {
	    lex_err(comline.last);
	    lprintf(stderr,"Found operator '%c'\n",term);
	    lprintf(stderr,"Where ')' was expected.\n");
	    return -1;
	  };
/*
  If the variable is effectively scalar then it can be passed by reference
  as well as value (unless it is a parameter). Also if the variable is an
  array it can be passed by reference if all the required elements are
  contiguous. The more rapidly changing indices are the first ones.
  Thus one may pass the array by reference if the first dimensions are
  unspecified and the rest are single elements - otherwise by value.
*/
	  if(dims != 0 && first_dim != dims)
	    expr_typ->access = 'v';
	  break;
	default:
/*
  No index expression follows the variable name so mark it with the pass
  by name access class. (Unless it is a parameter or constant).
*/
	  if(expr_typ->access == 'r')
	    expr_typ->access = 'N';
	  else
	    expr_typ->access = 'V';
	  break;
	};
/*
  Fill in the index user argument specifier entries in the INDEX_EXPR
  table entry.
*/
	for(i=0,j=0;j<3;j++,i +=3){
	  indval->start[j]=indexes[i];
	  indval->end[j]=indexes[i+1];
	  indval->inc[j]=indexes[i+2];
	};
/*
  Store the number of arguments recieved, in its nargs field.
*/
	indval->nargs = nargs;
/*
  Turn the dimensional type into a character representation.
*/
	expr_typ->dim = dims+'0';
	return no_error;                                               
}

/*.......................................................................
  This function is called by stack_var() to handle the index expression
  for a single dimension.
*/
static int sub_index(short int ind[], int j, char *term, short int *nargs,
		     char *dims)
{
	Exprtype expr_typ;
	int i;
	Table *ttst;
	char *tmplin;
/*
  Loop for optional numeric expressions for the initial index,
  final index and index increment.
*/
	for(i=0; i < 3; i++) {
	  switch (comline.nxtc) {
/*
  Is there a numeric expression for the index, or is there a terminator?
*/
	  case ':': case ',': case ')':
	    ttst = getoperator();
	    *term = *ttst->name;
	    switch (*term) {
	    case ',': case ')':
	      return no_error;
	    };
	    break;
/*
  There must be a numeric expression.
*/
	  default:
	    if( stack_expr(&tmplin, &ttst, OP_BR, 'i', '0', 'v', &expr_typ) == -1)
	      return -1;
/*
  Keep a record of the number of numeric expressions to be evaluated,
  and the terminating operator.
*/
	    (*nargs)++;
	    *term = *ttst->name;
/*
  Decide what to do by looking at the terminator.
*/
	    switch (*term) {
/*
  Comma's delimit the current index entry, while ')' terminates the
  whole index expression. The existence of a user specified
  number is signalled by placing a zero in the pertinent index entry
  ind[j+i]. If only the first number has been specified before a
  comma or ) then only one element of the current dimension is required
  signal this by setting the upper element bound ie. ind[j+1] to -1.
  Also decrement the number of invoked array dimensions (*dims).
*/
	    case ',': case ')':
	      ind[j+i] = *nargs;
	      if(i==0) {
		ind[j+1] = *nargs;
		(*dims)--;
	      };
	      return no_error;
	    case ':':
	      ind[j+i] = *nargs;
	      break;
/*
  Illegal operator encounterred.
*/
	    default:
	      lex_err(comline.last);
	      lprintf(stderr,"Found '%c' where an operand or one of , : ) was expected\n",*term);
	      return -1;
	    };
	  };
	};
/*
  The above loop should terminate before this point is reached - otherwise
  an extra colon has been encounterred.
*/
	 lex_err(comline.last);
	 lprintf(stderr,"Too many ':'s found \n");
	 return -1;
}

/*.......................................................................
  This routine is called by stack_operand() to parse a substring index
  expression of a character variable. stack_operand() calls this routine
  when it finds a '[' following a character variable. The sub-string
  expression should be of the form:
	    [expression_a_opt :_opt expression_b_opt] 
  Where _opt means optional - under the condition that at least one of
  expression_a and _b must be present. eg. string[1:] and string[:10]
  and string[1:nchan*10] are legal, but not string[].
    As far as the stack is concerned, the result will be that either
  one or two numeric expressions will be stacked followed by a special
  SUB_STRING class table entry. The SUB_STRING table entry "type" field
  points to one short int number that in turn describes the relationship
  between the expressions and the substring indexes. If this number is
  1 then the first argument alone has a user expression associated with
  it - if 2 then the second argument alone has one and if 3, both
  indices are described by user expressions. The arguments must match
  the dimensional shape of the string array being operated upon, hence
  the requirement that the expression type of the variable be passed
  in expr_typ.
*/
static int string_index(Exprtype *expr_typ)
{
	Table *optst;
	Exprtype index_type;
	char *tmplin;
	static char ops[]={':',']'};
	int i,inds=0;
/*
  Skip the '[' operator.
*/
	getoperator();
/*
  Get up to 2 index arguments.
*/
	for(i=0;i<2;i++) {
/*
  Check if the next index is specified or left as default.
  For the i+1'th argument the appropriate terminating operator
  is held in ops[i].
*/
	  if(comline.nxtc == ops[i]) {
/*
  Skip the operator since we know what it is.
*/
	    getoperator();
	  }
/*
  Get the index expression.
*/
	  else {
	    if( get_expr(&tmplin, &optst, OP_BR, &index_type) == -1)
	      return -1;
	    if(check_expression(&index_type, 'i', expr_typ->dim, 'v', 1) == -1)
	      return -1;
/*
  Check that the terminating operator is the correct one.
*/
	    if(*optst->name != ops[i]) {
	      lex_err(comline.last);
	      lprintf(stderr,"Was expecting '%c'.\n",ops[i]);
	      return -1;
	    };
/*
  Increment the inds variable that describes which bound has been specified
  (1=first, 2=second, 3=both) by the appropriate amount.
*/
	    inds += i+1;
	  };
	};
/*
  Check that at least one bound was specified.
*/
	if(inds == 0) {
	  lex_err(comline.last);
	  lprintf(stderr,"Illegal null sub-string index encounterred.\n");
	  return -1;
	};
/*
  Stack the SUB_STRING table entry with 'inds' as its value.
*/
	if(stack_instruct(SUB_STRING,inds) == -1)
	  return -1;
/*
  A sub-string reference means that the variable can only be passed
  by value.
*/
	expr_typ->access = 'v';
	return no_error;
}

/*.......................................................................
 * Stack the arguments for a call to the function in Table *ttst. Then
 * stack a constant holding a record of the number of arguments parsed
 * and finally stack the function itself.
 *
 * Input:
 *  ttst        Table *  The symbol table entry for the function.
 *  stmt          int    If true, the function is to be treated as a
 *                       statement, where parentheses do not enclose
 *                       the argument list.
 * Input/output:
 *  expr_typ Exprtype *  Send the pointer to an expression-type descriptor.
 *                       On output the descriptor that this points to will
 *                       be assigned the return type of the function.
 * Output:
 *  return        int    0 - OK.
 *                      -1 - Error.
 */
static int stack_function(Table *ttst, int stmt, Exprtype *expr_typ)
{
	int nargs,i;
	Table *tabtst;
	Descriptor *tmpdsc, *retdsc;
	int skip_from;
/*
  Get the return expression type.
*/
	expr_typ->type   = *TABFUNC(ttst)->type;
	expr_typ->dim    = *TABFUNC(ttst)->dim;
	expr_typ->access = *TABFUNC(ttst)->access;
/*
 * Fill in the appropriate return type if the function is
 * marked as being callable either as a command or a function.
 */
	if(expr_typ->access == '?') {
	  if(stmt)
	    expr_typ->type = expr_typ->dim = expr_typ->access = ' ';
	  else 
	    expr_typ->access = 'v';
	};
/*
 * Initialize the argument counter.
 */
	nargs = 0;
/*
  If an "execute once" function has been encounterred that returns a value
  then create a descriptor to point into the return value array.
  If the fuction actually returns scalar then this will be the
  full descriptor for the return value. Otherwise it will
  just point into the return array an element at a time.
  After the return descriptor - a BR_TO instruction is stacked
  instructing the elemental pass to skip the function arguments
  etc..
*/
	if(TABFUNC(ttst)->once && expr_typ->type != ' ') {
/*
  Get a table structure to point to the return descriptor.
*/
	  if( (tabtst=table_alloc(ARRAY_PTR,NULL)) == NULL)
	    return -1;
/*
  Create and stack a null descriptor for the return value.
*/
	  tmpdsc = (Descriptor *) malloc(sizeof(Descriptor));
	  if(tmpdsc == NULL) {
	    free(tabtst);
	    lprintf(stderr,"Memory allocation failed.\n");
	    return -1;
	  };
/*
  The descriptor will only point at one element at a time - it is
  thus effectively scalar.
*/
	  tmpdsc->num_el = 1;
	  for(i=0; i<3; i++)
	    tmpdsc->adim[i]=1;
	  tmpdsc->dim = '0';
	  tmpdsc->access = NO_DEL;
	  TABITEM(tabtst) = tmpdsc;
	  tabtst->name = ttst->name;
	  if(stack(tabtst,stack_ptr++) == -1) {
	    free(tmpdsc);
	    free(tabtst);
	    return -1;
	  };
/*
  Stack a BR_TO instruction and keep a record of its stack position.
  This will point to the instruction following the function arguments,
  indexes etc... in order that the expression evaluation function will
  skip them after the function has been evaluated once.
*/
	  skip_from = stack_ptr;
	  if(stack_instruct(BR_TO,0) == -1)
	    return -1;
/*
  Stack the function arguments.
*/
	  if(function_args(ttst, stmt, expr_typ, &nargs) == -1)
	    return -1;
/*
  Give the return descriptor the storage type determined
  by function_args().
*/
	  tmpdsc->atyp   = expr_typ->type;
/*
  Stack the function descriptor.
*/
	  if(stack(ttst,stack_ptr++) == -1)
	    return -1;
/*
  If the function returns a scalar value then allocate memory for
  it now.
*/
	  if(expr_typ->dim == '0' && expr_typ->access != 'r') {
	    tabtst->class = FN_RET;
	    if( (VOIDPTR(tmpdsc)=valof_alloc(1,expr_typ->type)) == NULL)
	      return -1;
	    tmpdsc->access = STACK;
	  }
/*
  If the function returns an array then stack any user supplied indexes.
*/
	  else {
/*
  Allocate a descriptor for the return array.
*/
	    retdsc = (Descriptor *) malloc(sizeof(Descriptor));
	    if(retdsc == NULL) {
	      lprintf(stderr,"Memory allocation failed.\n");
	      return -1;
	    };
/*
  The FN_ARRAY_REF access class means that the value field should not
  be freed and that the descriptor should not be deleted until the
  compile stack has been finished. The FN_ARRAY_VAL class means the
  same except that the value should be deleted after the current expression
  has been evaluated.
*/
	    retdsc->num_el = 0;
	    for(i=0; i<3; i++) retdsc->adim[i]=1;
	    retdsc->dim = expr_typ->dim;
	    retdsc->access = (expr_typ->access == 'r') ? FN_ARRAY_REF:FN_ARRAY_VAL;
	    VOIDPTR(retdsc) = NULL;
	    retdsc->atyp = expr_typ->type;
	    if(stack_array_indexes(expr_typ, tabtst, retdsc) == -1) {
	      free(retdsc);
	      return -1;
	    };
	  };
/*
  Link up the BR_TO instruction to point to the next instruction.
*/
	  TABICODE(compile_stack[skip_from]) = stack_ptr-skip_from-1;
	}
/*
  Elemental function - simply stack the argument expressions
  the number of such arguments and the function descriptor.
*/
	else {
/*
  Stack the function arguments.
*/
	  if(function_args(ttst, stmt, expr_typ, &nargs) == -1)
	    return -1;
/*
  Stack a record of the number of arguments parsed.
*/
	  if(stack_instruct(NUM_ARG,nargs) == -1)
	    return -1;
/*
  Stack the function descriptor.
*/
	  if(stack(ttst,stack_ptr++) == -1)
	    return -1;
	};
	return no_error;
}


/*.......................................................................
 * Parse a function argument-list.
 *
 * Input:
 *  fntst       Table *  The symbol table entry of the function.
 *  stmt          int    If true the function must be compiled as a
 *                       statement. This means that parentheses enclosing
 *                       the argument list are not expected.
 * Input/output:
 *  expr_typ Exprtype *  Send the declaration of the return type of the
 *                       function. For certain function declarations the
 *                       return type is partially determined by the function
 *                       arguments. In such cases, *expr_typ will be modified.
 *  nargs         int *  On output the number of arguments stacked, will be
 *                       returned.
 * Output:
 *  return        int    0 - OK.
 */
static int function_args(Table *fntst, int stmt, Exprtype *expr_typ, int *nargs)
{
  Table *optst;  /* The operator that terminates the argument list */
  int no_args=0; /* True if the function has no arguments */
/*
 * Initialize the record of the number of arguments to zero.
 */
  *nargs = 0;
/*
 * Determine if the argument list is empty.
 */
  if(stmt) {
    no_args = comline.nxtc == '\0';
  } else {
/*
 * Consume the open parenthesis of the argument list, if there is one.
 */
    if(comline.nxtc == '(') {
      getoperator();
/*
 * If the open parenthesis is immediately followed by a matching close
 * parenthesis, then there are no arguments.
 */
      if(comline.nxtc == ')') {
	no_args = 1;
	getoperator();  /* Consume the redundant close parenthesis */
      };
    } else {
      no_args = 1;  /* No parentheses implies that there are no arguments */
    };
  };
/*
 * If there is no argument list, check that this is ok before returning.
 */
  if(no_args) {
    if(TABFUNC(fntst)->nmin != 0) {
      lex_err(comline.last);
      lprintf(stderr,"%s() requires at least %d argument(s).\n", fntst->name,
	      TABFUNC(fntst)->nmin);
      return -1;
    };
    return no_error;
  } else {
/*
 * Does this function take any arguments?
 */
    if(TABFUNC(fntst)->nmax == 0) {
      lex_err(comline.last);
      lprintf(stderr, "%s() expects no arguments.\n", (fntst)->name);
      return -1;
    };
/*
 * Parse and compile the arguments of the function.
 */
    if(stack_args(fntst->name, TABFUNC(fntst)->nmin,
		 TABFUNC(fntst)->nmax, TABFUNC(fntst)->type,
		 TABFUNC(fntst)->dim, TABFUNC(fntst)->access,
		 TABFUNC(fntst)->once, nargs, expr_typ, &optst) == -1)
      return -1;
/*
 * Check for the appropriate argument-list delimiter.
 */
    if(stmt) {
      if(*optst->name != '\0') {
	lex_err(comline.last);
	lprintf(stderr,"Unexpected characters at end of line.\n");
	return -1;
      };
    } else {
      if(*optst->name != ')') {
	lex_err(comline.last);
	lprintf(stderr,"Unmatched \')\' in the arguments of %s().\n",
		fntst->name);
	return -1;
      };
    };
  };
  return no_error;
}

/*.......................................................................
  This is a temporary substitute for the function that will stack the
  table entry sent, onto the compile stack. Presently the table entry
  will simply be reported to the user's terminal and constant entries
  will be deleted afterwards. (int sptr) is the pointer to the stack
  position to which the table entry should be written. The usual global
  stack pointer is int stack_ptr, thus to increment this when a new
  table entry is appended the call should be stack(ttst,stack_ptr++).
*/
static int stack(Table *ttst, int sptr)
{
/*
  Check that the compile stack can hold the new entry.
*/
	if(sptr >= MAXSTACK-1) {
	  lprintf(stderr,"Sorry - compile stack full - no more room to compile into.");
	  lprintf(stderr,"Try shortening the current block of commands before retrying.");
	  return -1;
	};
	compile_stack[sptr] = ttst;
	return no_error;
}

/*.......................................................................
  This function is a header routine to get_expr(). In addition to calling
  get_expr() to stack the expression, it stacks a START_EXPR table entry
  before it and links it up afterwards. It also parses the optional
  {} array expression dimension specifier.
*/
static int stack_expr(char **opline, Table **ret_tab, int prec, char atyp,
		      char dim, char access, Exprtype *expr_typ)
{
	short int expr_start;
	int nargs;
	static Table *ttmp;
/*
  Take a note of the stack position at which the START_EXPR table entry
  will be written.
*/
	expr_start = stack_ptr;
/*
  Allocate memory for an empty START_EXPR table entry and stack it.
*/
	if( (ttmp=table_alloc(START_EXPR,NULL)) == NULL)
	  return -1;
	if(stack(ttmp,stack_ptr++) == -1)
	  return -1;
/*
  Check for a '{' delimitting the start of a dimensional specifier.
*/
	if(comline.nxtc == '{') {
	  getoperator();
/*
  Stack up to three scalar arguments.
*/
	  if(stack_args("{N0,N1,N2}", 1,3,"iii","000","vvv",1,&nargs,NULL,&ttmp) == -1)
	    return -1;
/*
  Check that the arguments were terminated by a '}'.
*/
	  if(*ttmp->name != '}') {
	    lex_err(comline.last);
	    lprintf(stderr,"Unmatched terminator: %s\n", ttmp->name);
	    return -1;
	  };
	}
	else
	  nargs=0;
/*
  Evaluate the expression.
*/
	if( get_expr(opline, ret_tab, prec, expr_typ) == -1)
	  return -1;
/*
  If a expression dimensions specifier was parsed above and it had
  a higher dimensional type than the expression, then enforce the
  dimensional type of the specifier. This allows apparent scalar
  expressions to be turned into array expressions eg: '{10} 2'
  produces an array of 10 elements, each with the value 2.
*/
	if(expr_typ->dim-'0' < nargs) expr_typ->dim = nargs+'0';
/*
  If there was a dimensional cast then this turns the expression
  into a value.
*/
	if(nargs != 0) expr_typ->access = 'v';
/*
  Check the expression type against that required.
*/
	if(check_expression(expr_typ, atyp, dim, access, 0) == -1)
	  return -1;
/*
  Link the START_EXPR table entry to the last used stack position and
  place the expression type in its second int.
*/
	if(link_expr_typ(expr_typ,expr_start) == -1)
	  return -1;
	return no_error;
}


/*.......................................................................
  Check the return type of an expression (in expr_typ) 
  against the delarations in char atyp, *dim and access. Also append the
  numeric conversion operators to the stacked expressions where required.
  If the expression is an elemetal argument then int elem should be
  given the value 1, otherwise 0 and char *lastdim will be updated.
*/
static int check_expression(Exprtype *expr_typ, char atyp, char dim,
			    char access, char is_elemental)
{
/*
  Where necessary stack a conversion operator to turn a float expression
  into an int, and vice versa.
  Numeric conversions can only be applied when the required return type
  is a value.
*/
	if(access == 'v') {
	  switch (atyp) {
	  case 'f':
	    if(expr_typ->type == 'i') {
	      if(stack_instruct(ITOF,0) == -1) return -1;
	      expr_typ->type = atyp;
	      expr_typ->access = 'v';
	    };
	    break;
	  case 'i':
	    if(expr_typ->type == 'f') {
	      if(stack_instruct(FTOI,0) == -1) return -1;
	      expr_typ->type = atyp;
	      expr_typ->access = 'v';
	    };
	    break;
	  };
	};
/*
  Check that the storage class is now that requested.
*/
	switch (atyp) {
	case '*':
	  break;
/*
  'n' means numeric - both int and float are legal.
*/
	case 'n':
	  if(expr_typ->type == 'f' || expr_typ->type == 'i')
	    break;
	default:
	  if(expr_typ->type != atyp) {
	    lex_err(comline.last);
	    lprintf(stderr, "Illegal %s operand where a %s was expected\n",
		    type_string(expr_typ->type), type_string(atyp));
	    return -1;
	  };
	};
/*
  Check that the dimensional type is that required.
*/
	if(expr_typ->dim > dim && dim != '*' && !is_elemental) {
	  lex_err(comline.last);
	  lprintf(stderr,"Illegal %s expression where a %s was expected.\n",
		  dims_string(expr_typ->dim), dims_string(dim));
	  return -1;
	};
/*
  Check legal correspondence of access class with that required.
*/
	switch (access) {
/*
  Pass by value - anything else is legal.
*/
	case 'v':
	  break;
/*
  Pass by reference - only pass-by-value is illegal here.
*/
	case 'r':
	  switch(expr_typ->access) {
	  case 'r': case 'N':
	    break;
	  default:
	    lex_err(comline.last);
	    lprintf(stderr, "This argument ought to be a reference to a variable.\n");
	    return -1;
	  };
	  break;
/*
  Pass by name - only pass by name is legal.
*/
	case 'N':
	  if(expr_typ->access != 'N') {
	    lex_err(comline.last);
	    lprintf(stderr, "This argument ought to be the name of a variable.\n");
	    return -1;
	  };
	  break;
	};
        return no_error;
}

/*.......................................................................
  Parse and stack arguments from the command line. The name of the
  function, or command that requires the arguments is sent in
  'char *name'. The minimum and maximum number of arguments are specified
  in nmin and nmax. The argument types for the first nmin'th arguments
  are sent in 'char type[]', their required dimensions in 'char dim[]'
  and their access types in 'char access[]'. If 'char once' is true (1)
  then each argument will be stacked individually as a full independant
  array expression. Otherwise the arguments will be stacked as elemental
  arguments - ie at run time, the argument expressions will be evaluated
  together element by element. Each element of 'char type[], char dim[] and
  char access[]' corresponds to one argument. Element 1 (not element 0)
  corresponds to the 1st argument. Valid 'type's are '*' - any argument
  type, 'c' - character string argument type, 'f' - numeric argument
  type, and 'l' - logical argument type. Valid 'dim's are '0' - scalar,
  '1' - 1 dimensional array, '2', 2D array, '3' - 3D array and '*' -
  any dimensional type. Valid 'access' types are 'N' - pass by name
  (essentially this means that the variable descriptor will be sent to the
  user-function.) 'r' - pass by reference (send a temporary descriptor
  that points to some element or some (contiguous) elements of a user
  variable) , and 'v' - pass by value.
   The actual number of arguments is returned in 'short int *nargs'.
  If a parse error occurs, the function return value will be -1, otherwise
  it is 0.
    If the arguments are for an elemental function (once=0) then the
  dimensional type of the expression that the function reference resides
  in is sent via *lastdim. This should be '0' if the function reference
  is the first operand in the expression. Each argument will have its
  dimensional type checked against *lastdim and if *lastdim is '0' and
  a different dimensional type is encounterred *lastdim will be updated
  with the new type, for subsequent checking..Thus the elemental dimesional
  type is returned via *lastdim. The terminating operator is returned via
  Table **optst. The calling function is required to check
  the start and end terminators, and this routine should not be called
  when no arguments are given.
*/
static int stack_args(char *name, short int nmin, short int nmax, char type[],
		      char dim[], char access[], char once, int *nargs,
		      Exprtype *rtntype, Table **optst)
{
        short int i,j,num_spec;
        Exprtype argtype;
        char *tmplin;
        *nargs = 0;
/*
  Determine the number of argument-type sepcifiers sent.
*/
	num_spec = strlen(type)-1;
/*
  Parse up to the maximum legal number of arguments.
*/
        for(i=1;i <= MAXARG; i++) {
/*
  See if another argument is legal.
*/
	  if(*nargs > nmax) {
	    lex_err(comline.last);
	    lprintf(stderr, "You have given %d arguments to %s() - which accepts a maximum of %d.\n", *nargs, name, (int) nmax);
	    return -1;
	  };
/*
  The required argument type is held in type[j], access[j] and
  dim[j].
    For the n'th argument, if n is greater than the number of argument-type
  specifiers then use the last specifier in the strings. j will be used to
  index the required specifier in the strings.
*/
          j = (i < num_spec) ? i : num_spec;
/*
  If the function is of the once only variety, each argument is evaluated
  immediately, resulting in an array to be passed to the function. The
  run-time system will handle this by calling the expression
  executor to evaluate the new expression. It must be told how long the
  expression is, in order to know where to loop over for each element.
  Thus, we shall keep a record of the current stack position and store
  a START_EXPR table entry before each function argument, and when the
  argument expression has been parsed, the number of table entries
  involved will be written to the "type" field of the START_EXPR table entry.
  This is performed in stack_expr(). 
*/
	  if(type[j]=='C') {   /* Literal string argument */
	    if(stack_lit(&tmplin, once, optst, &argtype) == -1)
	      return -1;
          } else if(once) {
            if(stack_expr(&tmplin, optst, COMMA, type[j], dim[j], access[j],
			  &argtype) == -1)
              return -1;
          }
          else {
            if(get_expr(&tmplin, optst, COMMA, &argtype) == -1)
              return -1;
	    if(check_expression(&argtype, type[j], dim[j], access[j], 1) == -1)
	      return -1;
	    rtntype->dim = (argtype.dim > rtntype->dim) ? argtype.dim : rtntype->dim;
          };
/*
  Some functions have wild-card declarations for their return
  dimensional or/and storage types. In these cases the return
  type is that of the 1st argument - return these changes via
  rtntype.
*/
	  if(i==1 && rtntype != NULL) {
	    if(*dim == '*')
	      rtntype->dim = argtype.dim;
	    if(*type == '*')
	      rtntype->type = argtype.type;
	  };
/*
  Increment a record of the number of arguments read.
*/
	  (*nargs)++;
/*
  Stop when the last argument has been read.
*/
	  if(*(*optst)->name != ',')
            break;
        };
/*
  All the users arguments have been parsed. Check that at least the minimum
  number required have been parsed.
*/
        if(*nargs < nmin || *nargs > nmax) {
          lex_err(comline.last);
          lprintf(stderr, "You have given %d arguments to %s() - it requires ", *nargs, name);
          if(nmin != nmax)
            lprintf(stderr, "between %d and %d arguments\n", nmin, nmax);
          else
            lprintf(stderr, "%d argument(s)\n", nmin);
          return -1;
        };
        return no_error;
}
                                           
/*.......................................................................
  Handle user declaration of variables. The storage class identifiers
  are 'string', 'float', or 'logical'. These must be followed by one
  or more variable name identifiers in the following format.
   title, source(scalar_num_expr,...) where an identifier not followed
  by parentheses is taken as a scalar, otherwise it is taken as an array.
  The number of dimensions that the array has is determined from the
  number of arguments provided by the user. These arguments will be
  resolved at run time, into declarations of the number of elements
  required on each axis. At compile time the variables will be declared
  with just one element, but with a dimensional type of the kind required.
  Declaration, in this context means installing the variable into the
  user-variable symbol table. Re-declaration of a variable that already
  exists is only legal within its declared type - this includes the
  dimensional type. The function that will actually be responsible for
  declaring the number of elements per dimension at run time, is
  sent via 'Table *fntst'.
*/
static int new_declare(Table *fntst)
{
        static Table *ttst,*new_tab;
        static Descriptor *dtst;
        static char dim, atyp, match_typ;
        static int bot,top,tab_pos,nargs;
        static long adim[3]={1,1,1};
/*
  Determine the storage class required for the list of variables.
*/
	atyp = TABFUNC(fntst)->type[0];
/*
  Loop through the list of new variable names and dimensions given
  by the user.
*/
        for(;;) {
/*
  Stack a DECLARATION warning stack instruction.
*/
          if(stack_instruct(DECL,0) == -1)
            return -1;
/*
  Get the next variable-name identifier.
*/
          if( (ttst=lex_expr('n')) == NULL)
            return -1;
/*
  Find out where the new variable should be inserted in the main symbol table.
*/
          match_typ = find_symbol(ttst->name, main_table, num_main, &bot, &top);
/*
  Check that the symbol name doesn't conflict with a function name.
*/
          if(match_typ == 'e' && main_table[bot]->class != VAR) {
            lex_err(comline.last);
            lprintf(stderr, "Requested variable name '%s' clashes with a ",ttst->name);
	    switch(main_table[bot]->class) {
	    case FUNC:
	      lprintf(stderr, "function");
	      break;
	    case MODULE_SYM:
	      lprintf(stderr, "module help topic");
	      break;
	    case HELP_SYM:
	      lprintf(stderr, "help topic");
	      break;
	    };
	    lprintf(stderr, " of the same name.\n");
            free(ttst->name);
            return -1;
          };
/*
  Make it illegal to re-declare an existing variable to a different storage type.
*/
          if(match_typ == 'e') {
            if(TABDESC(main_table[bot])->atyp != atyp) {
              lex_err(comline.last);
              lprintf(stderr, "Illegal re-declaration of %s variable '%s' to %s type.\n",
               type_string(TABDESC(main_table[bot])->atyp), ttst->name, type_string(atyp));
              free(ttst->name);
              return -1;
            }
            else if(TABDESC(main_table[bot])->access == R_ONLY) {
              lex_err(comline.last);
              lprintf(stderr, "Illegal re-declaration of a read-only parameter.\n");
              free(ttst->name);
              return -1;
            };
          };
/*
  Parse between 0 and 3 numeric arguments enclosed by parentheses and
  separated by commas. At run time these will be
  used to declare the number of elements per dimension. Although the
  arguments actually follow the new variable name, at run time they
  will be sent to the appropriate function along with the descriptor
  of the variable (to be created next).
*/
          if(function_args(fntst, 0, NULL, &nargs) == -1) {
            free(ttst->name);
            return -1;
          };
/*
  The number of arguments given by the user may be used to determine the required
  dimensional type for it.
*/
          dim = (char) '0'+nargs;
/*
  If an exact match was found with an existing variable, check that it has the
  same dimensional type - otherwise re-declaration is illegal.
*/
          if(match_typ == 'e') {
            tab_pos=bot;
            if(TABDESC(main_table[tab_pos])->dim != dim) {
              lex_err(comline.last);
              lprintf(stderr, "Illegal re-declaration of the %s variable '%s' into a %s.\n",
               dims_string(TABDESC(main_table[tab_pos])->dim), ttst->name, dims_string(dim));
              free(ttst->name);
              return -1;
            };
            free(ttst->name);
          }
          else {
/*
  Create a table structure to be inserted in the symbol table.
*/
            if((new_tab=table_alloc(VAR, ttst->name)) == NULL) {
              free(ttst->name);
              return -1;
            };
/*
  Create a descriptor for the variable. Give it one element for now,
  but give it the appropriate dimensional type determined from the number
  of arguments received above.
*/
            if( (dtst=descriptor_alloc(atyp,dim,adim)) == NULL) {
              free(new_tab);
              free(ttst->name);
              return -1;
            };
/*
  Find the correct position for the new variable in the symbol table.
*/
            if(match_typ == 'n')
              tab_pos = top;
            else
              tab_pos = bot;
/*
  Shift the variables above the required position up to make way for the insertion.
*/
            if(up_shift(main_table, &num_main, main_max, tab_pos) == -1) {
              free(ttst->name);
              free(new_tab);
              valof_free(dtst);
              free(dtst);
              return -1;
            };
/*
  Install the variable name and descriptor.
*/
            TABITEM(new_tab) = dtst;
            dtst->access = RWD;
            main_table[tab_pos] = new_tab;
          };
/*
  Stack the variable on the compile stack.
*/
          if(stack(main_table[tab_pos], stack_ptr++) == -1)
            return -1;
          if(comline.nxtc != ',') break;
          getoperator();
        };
/*
  Make sure that the next character is a '\0' - ie that no more characters
  follow the declarations.
*/
        if(comline.nxtc != '\0') {
          lex_err(comline.next);
          lprintf(stderr, "Unexpected characters at end of line.\n");
          return -1;
        };
        return no_error;
}


/*.......................................................................
  Parse and compile an assignment.
*/
static int stack_assign(Table *ttst, Table **optst)
{
	static char *tmplin;
	static Exprtype expr_typa, expr_typb;
	static int num_args;
/*
  Stack the variable and any index expression.
*/
	if(stack_var(ttst,&expr_typa) != 0)
	  return -1;
/*
  The next operator should be an equals sign.
*/
	if(comline.nxtc != '=' || (*optst = getoperator()) == NULL ||
	   (TABOPER(*optst)->op_prec != EQUALS)) {
	  lex_err(comline.last);
	  lprintf(stderr,"Missing '=' in variable assignment\n");
	  return -1;
	};
/*
  Check that the variable to be assigned to, is not a read only parameter.
*/
	if(TABDESC(ttst)->access == R_ONLY) {
	  lex_err(comline.last);
	  lprintf(stderr,"Illegal assignment: %s is a read-only parameter.\n",ttst->name);
	  return -1;
	};
/*
  Get the expression to be assigned to the variable found above.
*/
	for(num_args=1 ;; num_args++) {
	  if(num_args > MAXARG) {
	    lex_err(comline.last);
	    lprintf(stderr, "No command or assignment may have more than %d arguments\n", MAXARG);
	    return -1;
	  };
/*
  Get the next expression.
*/
	  if(stack_expr(&tmplin,optst,COMMA,expr_typa.type, expr_typa.dim, 'v', &expr_typb) == -1)
	    return -1;
/*
  The assignment may be made either with a single expression with a
  dimensional type compatible with the assignment variable, or a list of
  scalar expressions seperated by commas.
*/
	  if(num_args == 1 && *(*optst)->name != ',') {
/*
  The case of a single expression.
  Don't repeat the loop for a second expression.
*/
	    break;
	  }
/*
  The case of a list of scalar expressions.
*/
	  else {
	    if(expr_typb.dim != '0') {
	      lex_err(comline.last);
	      lprintf(stderr, "Illegal non-scalar expression in assignment list.\n");
	      return -1;
	    };
	    if(expr_typa.dim == '0') {
	      lex_err(comline.last);
	      lprintf(stderr, "Illegal assignment list to scalar variable: %s\n",ttst->name);
	      return -1;
	    };
/*
  Stop parsing assignment expressions when one is found without a
  comma delimiter.
*/
	    if(*(*optst)->name != ',') {
	      expr_typb.dim = '1';
	      break;
	    };
	  };
	};
	return no_error;
}


/*.......................................................................
 * Stack a literal string as a once-only or elemental argument.
 */
static int stack_lit(char **tmplin, int once, Table **optst,
		     Exprtype *expr_typ)
{
  short int expr_start;
  static Table *ttmp;
/*
 * Take a note of the stack position at which the START_EXPR table entry
 * will be written.
 */
  expr_start = stack_ptr;
/*
 * If this is a imediate evaluation argument, allocate memory for an empty
 * START_EXPR table entry and stack it.
 */
  if(once) {
    if( (ttmp=table_alloc(START_EXPR,NULL)) == NULL)
      return -1;
    if(stack(ttmp,stack_ptr++) == -1)
      return -1;
  };
/*
 * Get the literal string or `` enclosed string expression.
 */
  if(comline.nxtc=='`') { /* `string-expression` in place of literal string */
    getoperator();        /* Skip the ` operator */
/*
 * Get the enclosed expression.
 */
    if(get_expr(tmplin, optst, CL_BR, expr_typ) == -1)
      return -1;
/*
 * Only allow string expressions between ``.
 */
    if(expr_typ->type != 'c') {
      lex_err(comline.last);
      lprintf(stderr, "Illegal %s operand where a %s was expected.\n",
	      type_string(expr_typ->type), type_string('c'));
      return -1;
    };
/*
 * Check for the closing ` character
 */
    if(*optst==NULL || (*optst)->name[0] != '`') {
      lex_err(comline.last);
      lprintf(stderr, "Missing \"`\"\n");
      return -1;
    };
  } else {
/*
 * Get the literal string.
 */
    if((ttmp = lex_expr('l')) == NULL)
      return -1;
/*
 * Stack the literal string constant.
 */
    if(stack_var(ttmp, expr_typ) != 0)
      return -1;
  };
/*
 * Link the START_EXPR table entry to the last used stack position and
 * place the expression type in its second int.
 */
  if(once) {
    if(link_expr_typ(expr_typ,expr_start) == -1)
      return -1;
  };
/*
 * Get the following operator character.
 */
  if((*optst=getoperator()) == NULL)
    return -1;
  return no_error;
}

/*.......................................................................
 * Stack a command or a function that can be used as a command.
 *
 * Input:
 *  fntst     Table *   The symbol-table entry for the function.
 * Output:
 *  return      int     0 - OK.
 */
static int stack_command(Table *fntst)
{
  Exprtype expr_typ; /* The descriptor of the function return type */
/*
 * Stack a COMMAND marker instruction.
 */
  if(stack_instruct(COMMAND,0) == -1)
    return -1;
/*
 * Stack the function if legal.
 */
  if(TABFUNC(fntst)->type[0] == ' ' || TABFUNC(fntst)->access[0] == '?') {
/*
 * Stack the function and its arguments.
 */
    if(stack_function(fntst, 1, &expr_typ) == -1)
      return -1;
  } else {
    lex_err(comline.last);
    lprintf(stderr, "Function %s() can not be called as a command.\n",
	    fntst->name);
    return -1;
  };
  return no_error;
}

