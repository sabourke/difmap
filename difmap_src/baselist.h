#ifndef baselist_h
#define baselist_h

/*
 * Declare a type to be used to record a single baseline specification,
 * along with an indication of whether the baseline specification is
 * intended to be used to include or exclude baselines.
 */
typedef struct Basesel Basesel;
struct Basesel {
  int include;     /* If true include these baselines - if false exclude them */
  Basespec bs;     /* A baseline specification */
  Basesel *next;   /* The next baseline selection in the list */
};

/*
 * A container for the list of baseline groups.
 */
typedef struct Basegrp Basegrp;
struct Basegrp {
  Basesel *bsel;   /* The list of baseline selections. */
  int nnode;       /* The number of baseline selections in the list. */
  Basesel *tail;   /* The end of the Basesel list */
  Basegrp *next;   /* Pointer to the next selection list if part of a list */
};

/*
 * Declare a type to contain a list of Basegrp structures.
 */
typedef struct {
  Basegrp *bgrp;   /* The head of the list */
  Basegrp *tail;   /* The tail of the list */
  int nsel;
} Bgrplist;

/* Bgrplist constructor and destructor */

Bgrplist *new_Bgrplist(void);
Bgrplist *del_Bgrplist(Bgrplist *bsl);

/* Basegrp constructor and destructor */

Basegrp *new_Basegrp(void);
Basegrp *del_Basegrp(Basegrp *bgrp, Bgrplist *bsl);
Basegrp *clr_Basegrp(Basegrp *bgrp);

/* Append a baseline selection to a Basegrp list */

Basesel *add_Basesel(Observation *ob, Basegrp *bgrp, Basespec *bs, int include);

/* Discard a baseline selection node */

Basesel *del_Basesel(Basegrp *bgrp, Basesel *bsel);

/* Functions to write and read baseline specifications to and from strings */

int read_Basegrp(Observation *ob, Basegrp *bgrp, char *string, char **endp);
int write_Basegrp(Observation *ob, Basegrp *bgrp, int n, char *s);

/* Append a baseline group to a Bgrplist list */

Basegrp *add_Basegrp(Observation *ob, Bgrplist *bsl, Basegrp *bgrp,
		     char *bgrp_str);

/* Find out whether a baseline is included in a baseline group */

int in_Basegrp(Observation *ob, int isub, int base, Basegrp *bgrp);

/* Count the number of baselines in a group */

int size_Basegrp(Observation *ob, Basegrp *bgrp, int isub);

/* Search for the next included baseline in a given ordinal direction */

int srch_Basegrp(Observation *ob, Basegrp *bgrp, int forward,
		 int *s_isub, int *s_base);

/*
 * Declare a container of per-sub-array baseline index lists.
 */
typedef struct {   /* A baseline set for a single sub-array */
  int nbase;       /* The number of baselines in the list. */
  int *baselines;  /* An array of nbase baseline indexes. */
} Bsublist;

typedef struct {
  int ntotal;      /* The total number of baselines in the list */
  int nsub;        /* The number of sub-array entries in sub[]. */
  Bsublist *bsub;  /* Array of nsub baseline usability arrays */
} Baselist;

Baselist *new_Baselist(Observation *ob, Basegrp *bgrp);
Baselist *del_Baselist(Baselist *blist);

#endif
