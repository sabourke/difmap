#ifndef help_h
#define help_h

#include "pager.h"

#ifndef table_h
struct Table;
#endif

char *type_string(char type);
char *dims_string(char dim);
char *access_string(char access);
void whatisvar(char *name);
void help(char *name);
void apropos(char *name);
void list_modules(void);
void help_function(struct Table *symbol);
void help_module(struct Table *symbol);
void help_topic(struct Table *symbol);
int index_module(char *name);

#endif
