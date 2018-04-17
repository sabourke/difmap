#ifndef obedit_h
#define obedit_h

/* Define a node deferred edit operation lists */

typedef struct Edint {
  struct Edint *next;   /* Pointer to next edit in list of edits */
  short cif;            /* IF index used with 'selif' */
  short index;          /* Index of station or baseline */
  unsigned doflag  : 1; /* If set flag selected visibilities, else restore */
  unsigned selbase : 1; /* If set, only edit baseline 'index' */
  unsigned selstat : 1; /* If set, only edit baselines of station 'index' */
  unsigned selchan : 1; /* If set, only edit the source spectral-line channels*/
  unsigned selif   : 1; /* If set, only edit the current IF */
} Edint;

enum {EDBLK=256};

typedef struct Obedit {
  int nused;              /* Count of edit nodes that are in use */
  Edint *free;            /* Pointer to head of free list */
  struct edblock {
    Edint block[EDBLK];   /* Array of nodes for free list */
    struct edblock *next; /* Pointer to next block of free-list nodes */
  } blocks;               /* List of free-list blocks */
} Obedit;

Obedit *new_Obedit(Observation *ob);
Obedit *del_Obedit(Observation *ob);
int clr_Obedit(Observation *ob);
int app_Obedit(Observation *ob, int cif);

#endif
