#ifndef lex_h
#define lex_h

enum {MAX_LINE=512};  /* The maximum size of an input line of text. */

/*
  The Args type will be used to hold a record of arguments for the
  pre-processor. This will be used to keep a record of arguments
  to command files and when parsing macro definitions.
*/

typedef struct {
  char *arg_ptr;
  size_t arg_len;
} Args;

/*
  This is the structure to be used to maintain command lines.
  The '*last' and '*next' fields point to the last and next
  characters to be parsed in the command line. The nxtc field
  holds the next character in the line - ie the character pointed
  to by '*next'.
*/
typedef struct {
  char nxtc;       /* Command line look-ahead character */
  char *last;      /* Pointer to start of token in command line */
  char *next;      /* Pointer to next unused char in command line */
  int nest_block;  /* Command nesting level */
} Comline;

extern Comline comline;

/*
  Specify the legal number of nested command files and macro substitutions,
  and the maximum legal number of pre-processor arguments that they can have.
  Note - these should be kept small as they are used to declare multiple input
  buffers etc.. which are quite large.
*/

void lex_err(const char *line);

#ifndef table_h
struct Table;
#endif

struct Table *lex_expr(char optyp);

void flush_input(void);
int com_init(const char *bootenv);
int newline(void);

#endif
