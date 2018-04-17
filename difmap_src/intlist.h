#ifndef intlist_h
#define intlist_h

typedef struct Intbin { /* Container of list of groups in an integration bin */
  double ut;            /* Time-stamp of bin wrt fob->date (seconds) */
  long first;           /* The first group in the integration bin. */
  long last;            /* The last group in the integration bin */
  int isub;             /* Sub-array index of the integration */
  struct Intbin *next;  /* Pointer to a later integration bin */
} Intbin;

typedef struct Intlist Intlist;

Intlist *new_Intlist(int nsub, double origin, double binwid);
Intlist *del_Intlist(Intlist *igrid);
int add_group(Intlist *igrid, double ut, long group, int isub);

Intbin *nxt_Intbin(Intlist *igrid);
long nxt_group(Intbin *ibin);
int ibin_count(Intlist *igrid, int isub);

#endif
