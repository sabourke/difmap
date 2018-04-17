#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#include "modeltab.h"
#include "freelist.h"
#include "logio.h"

/*
 * Define a Model node entry.
 */
typedef struct ModelNode ModelNode;
struct ModelNode {
  ModelNode *next;      /* The next model node in a bucket list */
  Chlist *cl;           /* The list of selected channels */
  Stokes pol;           /* The selected polarization */
  Model *model;         /* The saved model */
};

/*
 * Each model hash-table bucket contains a linked list of model nodes that
 * hash to the same bucket.
 */
typedef struct {
  ModelNode *head;   /* The head of the bucket model-node list */
} ModelBucket;

/*
 * A model table consists of 'size' hash buckets. Note that the
 * ModelTable typedef for this struct is contained in modeltab.h.
 */
struct ModelTable {
  FreeList *node_mem;    /* ModelNode free-list */
  int size;              /* The number of hash buckets */
  ModelBucket *bucket;   /* An array of 'size' hash buckets */
  int nentry;            /* The number of entries in the table */
};

static ModelNode *del_ModelNode(ModelTable *mtab, ModelNode *node);
static ModelNode *new_ModelNode(ModelTable *mtab, Model *model, Chlist *cl,
				Stokes pol);
static ModelNode *find_ModelNode(ModelTable *mtab, ModelBucket *bucket,
				 Chlist *cl, Stokes pol, ModelNode **prev);
static ModelBucket *find_ModelBucket(ModelTable *mtab, Chlist *cl, Stokes pol);
static int write_ModelEntry(Chlist *cl, Stokes pol, Model *model,
			    char *filename, FILE *fp);
static int bad_mtab_file(FILE *fp, Model *model, Chlist *cl);

/*.......................................................................
 * Create a new hash table of models indexed by selected channel ranges
 * and polarizations.
 *
 * Input:
 *  size        unsigned   The size of the model hash table. Best
 *                         performance will be achieved if this is a
 *                         prime number.
 *  blkfact     unsigned   Model nodes are allocated in blocks. Specify
 *                         the blocking factor.
 * Output:
 *  return    ModelTable * The new model table, or NULL on error.
 */
ModelTable *new_ModelTable(unsigned size, unsigned blkfact)
{
  ModelTable *mtab;         /* The table to be returned */
  int i;
/*
 * Check arguments.
 */
  if(size <= 0) {
    lprintf(stderr, "new_ModelTable: Illegal table size (%d).\n", size);
    return NULL;
  };
/*
 * Allocate the container.
 */
  mtab = malloc(sizeof(ModelTable));
  if(!mtab) {
    lprintf(stderr, "new_ModelTable: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize
 * the container at least up to the point at which it can safely
 * be passed to del_ModelTable().
 */
  mtab->node_mem = NULL;
  mtab->size = size;
  mtab->bucket = NULL;
  mtab->nentry = 0;
/*
 * Allocate a freelist of ModelNode nodes.
 */
  mtab->node_mem = new_FreeList("new_ModelTable", sizeof(ModelNode), blkfact);
  if(!mtab->node_mem)
    return del_ModelTable(mtab);
/*
 * Allocate the array of 'size' hash buckets.
 */
  mtab->bucket = (ModelBucket *) malloc(sizeof(ModelBucket) * size);
  if(!mtab->bucket) {
    lprintf(stderr, "new_ModelTable: Insufficient memory for %d buckets.\n",
	    size);
    return del_ModelTable(mtab);
  };
/*
 * Initialize the bucket array.
 */
  for(i=0; i<size; i++)
    mtab->bucket[i].head = NULL;
/*
 * The table is ready for use - albeit currently empty.
 */
  return mtab;
}

/*.......................................................................
 * Delete a model table.
 *
 * Input:
 *  mtab   ModelTable *  The model table to be deleted.
 * Output:
 *  return ModelTable *  The deleted model table (always NULL).
 */
ModelTable *del_ModelTable(ModelTable *mtab)
{
  if(mtab) {
/*
 * Clear and delete the bucket array.
 */
    if(mtab->bucket) {
      clear_ModelTable(mtab);
      free(mtab->bucket);
      mtab->bucket = NULL;
    };
/*
 * Delete the freelist of model nodes.
 */
    mtab->node_mem = del_FreeList("del_ModelTable", mtab->node_mem, 1);
/*
 * Delete the container
 */
    free(mtab);
  };
  return NULL;
}

/*.......................................................................
 * Record a specified model along with the channel-range/polarization
 * selection that it is associated with. In order to be able to save
 * the model in an unshifted form, the option is provided to specify
 * shifts to be removed.
 *
 * Input:
 *  mtab   ModelTable *  The model table to insert the symbol into.
 *  model       Model *  The model to be recorded. A local copy of the
 *                       model will be recorded.
 *  cl         Chlist *  The list of selected channels associated with
 *                       the model. A local copy of the list will be made.
 *  pol        Stokes    The selected polarization associated with the
 *                       model.
 *  east        float    The eastward shift of the model to be removed
 *                       from the saved version of the model.
 *  north       float    The northward shift of the model to be
 *                       removed from the saved version of the model.
 * Output:
 *  return      Model *  The pointer to the recorded copy of the model,
 *                       or NULL on error.
 */
Model *add_ModelEntry(ModelTable *mtab, Model *model, Chlist *cl, Stokes pol,
		      float east, float north)
{
  ModelBucket *bucket;  /* The hash-bucket associated with the name */
  ModelNode *node;      /* The new node */
/*
 * Check arguments.
 */
  if(!mtab || !model || !cl) {
    lprintf(stderr, "add_ModelEntry: NULL argument(s).\n");
    return NULL;
  };
/*
 * Get the hash bucket of the specified selection.
 */
  bucket = find_ModelBucket(mtab, cl, pol);
/*
 * If a node with the same name already exists, replace its existing
 * model with the new one.
 */
  node = find_ModelNode(mtab, bucket, cl, pol, NULL);
  if(node) {
/*
 * Attempt to allocate a copy of the model.
 */
    Model *newmod = cpy_Model(model);
    if(!newmod)
      return NULL;
/*
 * Remove the specified eastward and northward shifts.
 */
    shiftmod(newmod, -east, -north);
/*
 * Replace the previous model with the new one.
 */
    node->model = del_Model(node->model);
    node->model = newmod;
    return model;
  };
/*
 * Allocate a new node containing copies of the model and channel list.
 */
  node = new_ModelNode(mtab, model, cl, pol);
  if(!node)
    return NULL;
/*
 * Install the node at the head of the hash-bucket list.
 */
  node->next = bucket->head;
  bucket->head = node;
/*
 * Remove the specified eastward and northward shifts.
 */
  shiftmod(node->model, -east, -north);
/*
 * Record the successful addition of a new entry.
 */
  mtab->nentry++;
  return node->model;
}

/*.......................................................................
 * Remove a given model from the table. In order to be able to move
 * the model to the parents coordinate frame, the option is provided
 * to specify shifts to be applied during removal.
 *
 * Input:
 *  mtab   ModelTable *  The model table to find the symbol in.
 *  cl         Chlist *  The channel list associated with the model.
 *  pol         Stokes   The polarization associated with the model.
 *  east        float    The eastward shift of the model to be
 *                       added to the returned model.
 *  north       float    The northward shift of the model to be
 *                       added to the returned model.
 * Output:
 *  return  ModelNode *  The removed model. If a model associated with
 *                       the requested channel ranges and polarization
 *                       isn't found, NULL is returned without an error
 *                       message being emitted.
 */
Model *rem_ModelEntry(ModelTable *mtab, Chlist *cl, Stokes pol, float east,
		      float north)
{
  Model *model = NULL;    /* The removed model */
  if(mtab && cl) {
    ModelBucket *bucket = find_ModelBucket(mtab, cl, pol);
    ModelNode *prev;   /* The node preceding the located node */
    ModelNode *node = find_ModelNode(mtab, bucket, cl, pol, &prev);
/*
 * Node found?
 */
    if(node) {
/*
 * Remove the node from the bucket list.
 */
      if(prev) {
	prev->next = node->next;
      } else {
	bucket->head = node->next;
      };
/*
 * Extract the model before deleting the node.
 */
      model = node->model;
      node->model = NULL;
/*
 * Delete the node.
 */
      (void) del_ModelNode(mtab, node);
/*
 * Record the removal of an entry from the table.
 */
      mtab->nentry--;
/*
 * Add the specified shift to the model.
 */
      shiftmod(model, east, north);
    };
  };
  return model;
}

/*.......................................................................
 * See if a model exists for the given channel range list and polarization.
 *
 * Input:
 *  mtab   ModelTable *  The model table to find the symbol in.
 *  cl         Chlist *  The channel list associated with the model.
 *  pol        Stokes    The polarization associated with the model.
 *  non_empty     int    If true, don't say that there is a model unless
 *                       it contains at least 1 component.
 * Output:
 *  return        int    0 - No model found.
 *                       1 - The table contains a model for the specified
 *                           channel ranges and polarization.
 */
int have_ModelEntry(ModelTable *mtab, Chlist *cl, Stokes pol, int non_empty)
{
/*
 * Lookup the model if possible.
 */
  if(mtab && cl) {
    ModelBucket *bucket = find_ModelBucket(mtab, cl, pol);
    ModelNode *node = find_ModelNode(mtab, bucket, cl, pol, NULL);
    if(node && (!non_empty || node->model->ncmp > 0))
      return 1;
  };
/*
 * No model found.
 */
  return 0;
}

/*.......................................................................
 * Private function used to allocate a model-table node.
 * The caller is responsible for checking that the specified symbol
 * is unique and for installing the returned entry in the table.
 *
 * Input:
 *  mtab   ModelTable *  The table to allocate the node for.
 *  model       Model *  The model to be recorded. A copy will be recorded.
 *  cl         Chlist *  The list of selected channels associated with
 *                       the model. A copy of this list will be recorded.
 *  pol         Stokes   The selected polarization associated with the
 *                       model.
 * Output:
 *  return  ModelNode *  The new node, or NULL on error.
 */
static ModelNode *new_ModelNode(ModelTable *mtab, Model *model, Chlist *cl,
				Stokes pol)
{
  ModelNode *node;  /* The new node */
/*
 * Allocate the new node from the free list.
 */
  node = new_FreeListNode("new_ModelNode", mtab->node_mem);
  if(!node)
    return NULL;
/*
 * Before attempting any operation that might fail, initialize the
 * new node at least up to the point at which it can safely be passed
 * to del_ModelNode().
 */
  node->next = NULL;
  node->model = NULL;
  node->cl = NULL;
  node->pol = pol;
/*
 * Allocate a new copy of the model.
 */
  node->model = cpy_Model(model);
  if(!node->model)
    return del_ModelNode(mtab, node);
/*
 * Allocate a new copy of the channel list.
 */
  node->cl = cpy_Chlist(cl);
  if(!node->cl)
    return del_ModelNode(mtab, node);
  return node;
}

/*.......................................................................
 * Private function used to delete a model-table node.
 * The node must have been removed from its list before calling this
 * function.
 *
 * Input:
 *  mtab   ModelTable *  The table for which the node was originally
 *                      allocated.
 *  node    ModelNode *  The node to be deleted.
 * Output:
 *  return  ModelNode *  The deleted node (always NULL).
 */
static ModelNode *del_ModelNode(ModelTable *mtab, ModelNode *node)
{
  if(node) {
/*
 * Delete the model and the channel list.
 */
    node->model = del_Model(node->model);
    node->cl = del_Chlist(node->cl);
    node->pol = NO_POL;
    node->next = NULL;
/*
 * Return the node to the free-list.
 */
    node = del_FreeListNode("del_ModelNode", mtab->node_mem, node);
  };
  return NULL;
}

/*.......................................................................
 * Private function to locate the hash bucket associated with a given
 * name.
 *
 * This uses a model-function described in the dragon-book
 * ("Compilers - Principles, Techniques and Tools", by Aho, Sethi and
 *  Ullman; pub. Adison Wesley) page 435.
 *
 * Input:
 *  mtab    ModelTable *   The table to look up the string in.
 *  cl          Chlist *   The list of channel ranges to look up.
 *  pol         Stokes *   The polarization to look up.
 * Output:
 *  return ModelBucket *   The located hash-bucket.
 */
static ModelBucket *find_ModelBucket(ModelTable *mtab, Chlist *cl, Stokes pol)
{
  unsigned long h = 0L;
  int range;
  for(range=0; range<cl->nrange; range++) {
    Chans *chans = cl->range + range;
    unsigned long ca = chans->ca;
    unsigned long cb = chans->cb;
    h = 65599UL * (65599UL * h + ca) + cb;  /* 65599 is a prime close to 2^16 */
  };
  return mtab->bucket + (h % mtab->size);
}

/*.......................................................................
 * Search for a given channel/polarization selection in the entries of
 * a specified bucket.
 *
 * Input:
 *  mtab     ModelTable *  The model-table being searched.
 *  bucket  ModelBucket *  The bucket to search (use find_ModelBucket()).
 *  cl           Chlist *  The channel list to find.
 *  pol          Stokes    The polarization to find.
 * Output:
 *  prev      ModelNode ** If prev!=NULL then the pointer to the node
 *                         preceding the located node in the list will
 *                         be recorded in *prev. This will be NULL either
 *                         if the name is not found or the located node is
 *                         at the head of the list of entries.
 * return     ModelNode *  The located model-table node, or NULL if not
 *                         found.
 */
static ModelNode *find_ModelNode(ModelTable *mtab, ModelBucket *bucket,
				 Chlist *cl, Stokes pol, ModelNode **prev)
{
  ModelNode *last;  /* The previously searched node */
  ModelNode *node;  /* The node that is being searched */
/*
 * Search the list for a node containing the specified polarization and
 * channel list.
 */
  for(last=NULL, node=bucket->head; node; last=node,node=node->next) {
    if(pol == node->pol && eq_Chlist(cl, node->cl))
      break;
  };
  if(prev)
    *prev = node ? last : NULL;
  return node;
}

/*.......................................................................
 * Empty a model-table by deleting all of its entries.
 *
 * Input:
 *  mtab    ModelTable *  The model table to clear.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int clear_ModelTable(ModelTable *mtab)
{
  int i;
  if(!mtab) {
    lprintf(stderr, "clear_ModelTable: NULL argument.\n");
    return 1;
  };
/*
 * Clear the contents of the bucket array.
 */
  for(i=0; i<mtab->size; i++) {
    ModelBucket *bucket = mtab->bucket + i;
/*
 * Delete the list of active model nodes from the bucket.
 */
    ModelNode *node = bucket->head;
    while(node) {
      ModelNode *next = node->next;
      (void) del_ModelNode(mtab, node);
      node = next;
    };
/*
 * Mark the bucket as empty.
 */
    bucket->head = NULL;
  };
/*
 * Mark the table as empty.
 */
  mtab->nentry = 0;
  return 0;
}

/*.......................................................................
 * Write the contents of a table of models to a file.
 *
 * Input:
 *  mtab   ModelTable *  The table to be saved.
 *  filename     char *  The name of the file to write the contents to.
 * Output:
 *  return        int    0 - OK.
 *                       1 - Error.
 */
int write_ModelTable(ModelTable *mtab, char *filename)
{
  FILE *fp;   /* The file pointer of the output file */
  int i;
/*
 * Check the arguments.
 */
  if(!mtab || !filename) {
    lprintf(stderr, "write_ModelTable: NULL argument(s).\n");
    return 1;
  };
/*
 * Attempt to open the output file.
 */
  fp = fopen(filename, "w");
  if(!fp) {
    lprintf(stderr, "Unable to open %s (%s).\n", filename, strerror(errno));
    return 1;
  };
/*
 * Scan the table for entries, writing each one to the file.
 */
  for(i=0; i<mtab->size; i++) {
    ModelBucket *bucket = mtab->bucket + i;
    ModelNode *node;
    for(node=bucket->head; node; node=node->next) {
      if(write_ModelEntry(node->cl, node->pol, node->model, filename, fp)) {
	fclose(fp);
	return 1;
      };
/*
 * Separate each entry from its successor by an empty line.
 */
      if(node->next && lprintf(fp, "\n") < 0) {
	lprintf(stderr,
		"Error writing to %s (%s).\n", filename, strerror(errno));
	fclose(fp);
	return 1;
      };
    };
  };
/*
 * Close the file.
 */
  if(fclose(fp) == EOF) {
    lprintf(stderr, "Error closing %s (%s).\n", filename, strerror(errno));
    return 1;
  };
  return 0;
}

/*.......................................................................
 * Write the contents of a single modeltable entry to a file.
 *
 * Input:
 *  cl     Chlist *  The channel range of the entry.
 *  pol    Stokes    The polarization of the entry.
 *  model   Model *  The model of the entry.
 *  filename char *  The name of the file being written to.
 *  fp       FILE *  The stream to write to.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int write_ModelEntry(Chlist *cl, Stokes pol, Model *model,
			    char *filename, FILE *fp)
{
/*
 * Check the arguments.
 */
  if(!cl || !model || !filename || !fp) {
    lprintf(stderr, "write_ModelEntry: NULL argument(s).\n");
    return 1;
  };
  if(pol==NO_POL) {
    lprintf(stderr, "write_ModelEntry: Polarization not specified.\n");
    return 1;
  };
/*
 * Don't write empty models.
 */
  if(model->ncmp < 1)
    return 0;
/*
 * Write the select command associated with the model.
 */
  if(lprintf(fp, "select %s, ", Stokes_name(pol)) < 0 ||
     write_Chlist(cl, fp, NULL) || lprintf(fp, "\n") < 0) {
    lprintf(stderr, "Error writing to %s (%s).\n", filename, strerror(errno));
    return 1;
  };
/*
 * Write the model.
 */
  if(wmodel(model, 0.0, 0.0, 0, 0.0, fp))
    return 1;
  return 0;
}

/*.......................................................................
 * Restore a model table from a file previously written by
 * write_ModelTable().
 * 
 * Input:
 *  mtab    ModelTable *  The model table to hold the restored table.
 *                        This will be cleared of existing entries before
 *                        reading the file.
 *  filename      char *  The name of the file to restore the table from.
 * Output:
 *  return         int    0 - OK.
 *                        1 - Error.
 */
int read_ModelTable(ModelTable *mtab, char *filename)
{
  FILE *fp=NULL;       /* The file stream connected to the input file */
  Model *model = NULL; /* The model of the latest entry */
  Chlist *cl = NULL;   /* The channel list of the latest entry */
  Stokes pol = NO_POL; /* The polarization of the latest entry */
  int nline=0;         /* The current line number */
/*
 * Check the arguments.
 */
  if(!mtab || !filename) {
    lprintf(stderr, "read_ModelTable: NULL argument(s).\n");
    return 1;
  };
/*
 * Start with an empty table.
 */
  if(clear_ModelTable(mtab))
    return 1;
/*
 * Open the file for read.
 */
  fp = fopen(filename, "r");
  if(!fp) {
    lprintf(stderr, "Unable to open %s (%s).\n", filename, strerror(errno));
    return 1;
  };
/*
 * Each entry in the file starts with a select command line followed by
 * the associated model.
 */
  while(!feof(fp)) {
    int c;            /* The latest character read from the file */
/*
 * Skip white-space to the first character of the next line.
 */
    do {
/*
 * Record the start of the new line.
 */
      nline++;
/*
 * Skip leading spaces and tabs.
 */
      do {c = getc(fp);} while(c==' ' || c=='\t');
/*
 * Skip full-line comments.
 */
      if(c=='!')
	do {c = getc(fp);} while(c!='\n' && c!=EOF);
/*
 * Skip empty lines.
 */
    } while(c=='\n');
/*
 * End of file?
 */
    if(feof(fp))
      break;
/*
 * Putback the character just read.
 */
    if(ungetc(c, fp) == EOF) {
      lprintf(stderr, "Error putting back character while reading %s.\n");
      return bad_mtab_file(fp, model, cl);
    };
/*
 * Only lines containing select commands start with 's' (all model file lines
 * start with ! or a signed number).
 */
    if(c=='s') {
      char polname[3];  /* The name of the polarization */
/*
 * This completes the previous entry, so add the previous model to the table,
 * and discard its model and channel list.
 */
      if(model) {
	if(add_ModelEntry(mtab, model, cl, pol, 0.0f, 0.0f) == NULL)
	  return bad_mtab_file(fp, model, cl);
/*
 * Report what has been read.
 */
	lprintf(stdout, "Read %d model components for stokes %s, channels ",
		model->ncmp, Stokes_name(pol));
	write_Chlist(cl, stdout, NULL);
	lprintf(stdout, "\n");
/*
 * The model and channel list were copied, so clear the model so that
 * it can be used again, and delete the channel list.
 */
	model = clr_Model(model);
	cl = del_Chlist(cl);
      };
/*
 * Create a new model and channel lists.
 */
      if(!model)
	model = new_Model();
      if(!model)
	return bad_mtab_file(fp, model, cl);
/*
 * Skip the 'select' command name, and read the polarization
 * name.
 */
      if(fscanf(fp, "select %2[^,]", polname) != 1) {
	lprintf(stderr, "Syntax error on line %d of %s.\n", nline, filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * Decode the polarization.
 */
      pol = Stokes_id(polname);
      if(pol==NO_POL) {
	lprintf(stderr, "Unknown polarization specified on line %d of %s.\n",
		nline, filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * The next non-space character should be a comma.
 */
      do {c = getc(fp);} while(c==' ' || c=='\t');
      if(c != ',') {
	lprintf(stderr, "Missing channel ranges on line %d of %s.\n",
		nline, filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * Read the channel ranges and add them to the channel-range list.
 */
      cl = read_Chlist(fp, filename, nline);
      if(!cl)
	return bad_mtab_file(fp, model, cl);
/*
 * No channel ranges?
 */
      if(!cl->nrange) {
	lprintf(stderr, "Missing channel ranges on line %d of %s.\n",
		nline, filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * Nothing should follow the channel list on the selection line.
 */
      do {c = getc(fp);} while(c==' ' || c=='\t');
      if(c != '\n') {
	lprintf(stderr, "Corrupt select specification on line %d of %s.\n",
		nline, filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * Read a model component line.
 */
    } else {
/*
 * If no select line has been seen yet, no model container will have been
 * allocated.
 */
      if(!model) {
	lprintf(stderr,
	     "Missing select line implies that %s isn't a multi-model file.\n",
	     filename);
	return bad_mtab_file(fp, model, cl);
      };
/*
 * Attempt to read the new component.
 */
      if(read_Modcmp(model, 0.0, 0.0, 1, filename, fp, &nline) == CMP_ERROR)
	return bad_mtab_file(fp, model, cl);
    };
  };
/*
 * To get here we must have reached the end of the file. Add the latest
 * model to the table and cleanup.
 */
  if(model && add_ModelEntry(mtab, model, cl, pol, 0.0f, 0.0f) == NULL)
    return bad_mtab_file(fp, model, cl);
  fclose(fp);
  model = clr_Model(model);
  cl = del_Chlist(cl);
  return 0;
}

/*.......................................................................
 * This is the private error return function of read_ModelTable().
 *
 * Input:
 *  fp         FILE *  The input stream.
 *  model     Model *  The model being composed.
 *  cl       Chlist *  The channel list being composed.
 * Output:
 *  return      int    The error return status of read_ModelTable().
 */
static int bad_mtab_file(FILE *fp, Model *model, Chlist *cl)
{
  fclose(fp);
  model = del_Model(model);
  cl = del_Chlist(cl);
  return 1;
}

/*.......................................................................
 * Return the count of the number of entries in the table.
 *
 * Input:
 *  mtab   ModelTable *  The table to count entries in.
 * Output:
 *  return        int    The number of entries. This will be zero if
 *                       mtab==NULL.
 */
int num_ModelTable_entries(ModelTable *mtab)
{
  return mtab ? mtab->nentry : 0;
}

