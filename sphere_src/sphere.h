#ifndef sphere_h
#define sphere_h

/*-----------------------------------------------------------*/
/* This header must be included by all user supplied modules */
/*-----------------------------------------------------------*/

/* -------  User variable types and method functions -------------------- */

/*
  Declare the variable symbol-table. The following is an array of
  pointers to up to MAXVAR possible user variables. This array
  will be filled in later, as the user declares the variables.
  The variable num_var will be incremented on each declaration to
  keep a record of the total number of user variables in existence.
*/

#define MAXVAR 100

/*
  Define the access classes - 
  0) RWD = declared at run time by user, read+write+delete.
  1) R_ONLY= parameters declared here - read access but no write or delete access.
  2) NO_DEL= variables declared here - read+write access but no delete access.
  3) STACK & TEMP are for internal use: STACK is for storage of constant
     values on the compile stack, TEMP is for temporaries on the run stack
*/
enum {RWD,R_ONLY,NO_DEL,TEMP,STACK,REF,FN_ARRAY_REF, FN_ARRAY_VAL, DESCR};


typedef union {          /* The pointer to each element of temp will be */
  float fval;      /* assigned to a run-stack value field.        */
  char lval;       /* They will be used to hold temporary values  */
  char *cptr;      /* during execution of the compile stack.      */
  int ival;
} Equiv;

/*
  Each variable and constant recognised by the handler is referred
  to by its descriptor. This consists of a storage type (atyp =
   '*'=any, 'i'=int, 'f'=float, 'c'=string, 'l'=logical), a dimensional
  type (dim =  '0'=scalar, '1'=1D_array, '2'=2D_array, '3'=3D_array),
  an access specifier (access =  R_ONLY, RWD, TEMP, NO_DEL,..),
  the total number of elements allocated to the variable, (num_el),
  the number of these currently in use on each dimension, (adim[3]),
  and a pointer to the first element of the variable (valof).
*/

typedef struct Descriptor {
  char atyp;          /* Type of variable (f,i, c or l) */
  char dim;           /* The dimensional type '0','1','2','3' */
  short access;       /* Access type - whether modifiable etc..*/
  long num_el;        /* The total number of elements. */
  long adim[3];       /* No. of elements per dimension. */
  void *value;
} Descriptor;         /*  element of the user-variable array */

/*.......................................................................
  Declare functions relevant to variable symbol-tables.
*/

#ifndef table_h
struct Table;
#endif

void var_free(struct Table *stab);
void valof_free(Descriptor *dtst);
struct Table *store_const(char vtype, void *value);
void free_const(struct Table *stab);
int re_declare(Descriptor *dtst, long dims[3]);
Descriptor *descriptor_alloc(char vtype, char dim, long adim[3]);
void *valof_alloc(int nvals, char vartyp);
void *valof_realloc(void *value, char vartyp, int n1, int n2);

/* Macros to access descriptor values. dsc should be a (Descriptor *) */

#define INTPTR(dsc) ((int *)(dsc)->value)
#define FLTPTR(dsc) ((float *)(dsc)->value)
#define STRPTR(dsc) ((char **)(dsc)->value)
#define LOGPTR(dsc) ((char *)(dsc)->value)
#define VOIDPTR(dsc) ((dsc)->value)
#define EQUIVPTR(dsc) ((Equiv *)(dsc)->value)

/* -----------  Run-time function types and method-functions ------------ */

/*
  Define the function sub-classes. The function sub-classes are NORM
  for normal user variables. The other types are for internal use
  as specialised stack commands.
*/

enum {BRK_BLOCK,CONT_BLOCK,START_BLOCK,END_BLOCK,STOP_EXE,WHATVAR,
      HELP,SYS,DECLARE,NORM};

#define MAXARG 40              /* The maximum number of arguments that a */
                               /* function may have. Be warned - This also */
                               /* determines the size of the array stack. */

/*
  A macro for declaring user functions.
*/
#define Template(a) int (a)(Descriptor *invals[],int npar,Descriptor *outvals)

/*
  A macro to count the number of strings in an array of strings.
*/
#define COUNT(a) sizeof(a)/sizeof(char *)

/*
  Each user function requires a structure. This contains a pointer
  to the actual function (Or a buffering function) - (*fname).
     When invoked, the function will be passed a pointer to an array of 
  descriptors containing the translated arguments that the
  user typed, an integer to specify how many arguments that entailed,
  and a pointer to the descriptor for the return value of the function.
     For the sake of the translator the required number of arguments
  will be specified via two integers - if a variable number of
  arguments can be sent then these two integers will have different
  values specifying the minimum and maximum number of arguments treatable.
     The type of each argument will be specified by two character arrays.
  Each character in each array describes one argument - the first array
  (atyp) specifies the type and the second (aarr) whether the argument
  should be an array. The zero'th element will specify the function
  return type.
*/

typedef struct {
  Template(*fname);       /* The user function */
  short int sub_class;    /* The subclass of function */
  short int nmin;         /* Min no. of arguments required */
  short int nmax;         /* Max no. of arguments required */
  char *type;     /* '*'=any, 'c'=char, 'n'=float, 'i'=int, 'l'=logic */
  char *dim;      /* '0'=scaler,'1'=array,'2'=image,'3'=cube. */
  char *access;   /* If 'v'=pass by value, 'r'=pass by reference */
  char once;      /* 'T', if any element of *type is other than scalar */
  struct Table *help;     /* Pointer to module help entry */
} Functype;

/* Closedown functions */

typedef enum {DO_EXIT, DO_QUIT} Exitcode;
#define EXITFN(fn) void (fn)(Exitcode code)
int add_exit_fn(EXITFN(*fn));
int closedown(int status, Exitcode code);

/* ----------- Module definition and method functions ------------------- */

/*
 * Define the structure that holds the members of a particular module.
 */

typedef struct {
  char *name;		/* Name of module */
  char *help_dir;	/* Directory of help files on this module */
  char **h_name;        /* Array of extra help-topic names */
  int h_num;            /* The number of elements in help_topics[] */
  Descriptor *v_type;	/* Pointer to variable declaration structures */
  char **v_name;	/* Array of variable names */
  int v_num;		/* Number of variables in module */
  Functype *f_type;	/* Pointer to Function declaration structures */
  char **f_name;	/* Array of function names */
  int f_num;		/* Number of functions in module */
  int (*begin)(void);   /* Optional intialization function */
  EXITFN(*end);         /* Optional cleanup function */
} Module;

int module_init(Module **modules, int nmodule);
int startup(Module **modules, int nmodule, const char *bootvar);
char *stralloc(size_t nchar);

int com_open(const char *filestr);  /* Have a comand file executed */

/* Get a command line from the user or file */

int lexgets(char *buff, int nmax, FILE *stream, char *prompt);
int push_command(FILE *fp, const char *comstr, const char *filename,
		 const char *argstr);
int pause_output(void);

/*
 * If no pgplot device is currently open, prompt the user for a device name.
 * On failure return -1.
 */
int make_open(void);

#ifndef FOPEN_MAX
#define FOPEN_MAX 10
#endif

extern int no_error;    /* Flag set in sig.c on interrupt or on user error */
extern char *null_string;

#endif
