#ifndef pollist_h
#define pollist_h

/*
 * Declare the type of a polarization linked-list node.
 */
typedef struct Polnode Polnode;
struct Polnode {
  Stokes pol;
  Polnode *next;
};

/*
 * Declare a polarization list type.
 */
typedef struct {
  Polnode *head;   /* The head of the list of polarizations */
  Polnode *tail;   /* The tail of the list of polarizations */
  int npol;        /* The number of polarizations in the list */
} Pollist;

Pollist *new_Pollist(void);
Pollist *del_Pollist(Pollist *pl);
Pollist *clr_Pollist(Pollist *pl);
Polnode *add_Polnode(Observation *ob, Pollist *pl, Stokes pol);
Polnode *del_Polnode(Pollist *pl, Polnode *polnode);

#endif
