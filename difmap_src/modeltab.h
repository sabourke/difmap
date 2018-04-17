#ifndef modeltab_h
#define modeltab_h

#include "obs.h"
#include "model.h"
#include "chlist.h"

typedef struct ModelTable ModelTable;

ModelTable *new_ModelTable(unsigned size, unsigned blkfact);
ModelTable *del_ModelTable(ModelTable *mtab);
Model *add_ModelEntry(ModelTable *mtab, Model *model, Chlist *cl, Stokes pol,
		      float east, float north);
Model *rem_ModelEntry(ModelTable *mtab, Chlist *cl, Stokes pol,
		      float east, float north);
int have_ModelEntry(ModelTable *mtab, Chlist *cl, Stokes pol, int non_empty);
int clear_ModelTable(ModelTable *mtab);
int write_ModelTable(ModelTable *mtab, char *filename);
int read_ModelTable(ModelTable *mtab, char *filename);
int num_ModelTable_entries(ModelTable *mtab);

#endif
