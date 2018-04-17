#ifndef chlist_h
#define chlist_h

typedef struct Chans Chans;
struct Chans {
  int ca;
  int cb;
};

typedef struct {
  Chans *range;  /* Array of disjoint channel ranges in ascending order */
  int nrange;    /* The number of used elements in range[] */
  int ca;        /* Index of lowest channel in range[]... */
  int cb;        /* Index of highest channel in range[]... */
} Chlist;

Chlist *new_Chlist(void);
Chlist *sub_Chlist(Chlist *cl, int coff, int nchan);
Chlist *del_Chlist(Chlist *cl);
Chlist *cpy_Chlist(Chlist *cl);
int eq_Chlist(Chlist *cl1, Chlist *cl2);
int add_crange(Chlist *cl, int ca, int cb);
int lim_Chlist(Chlist *cl, int nchan);
int write_Chlist(Chlist *cl, FILE *fp, char *filename);
Chlist *read_Chlist(FILE *fp, char *filename, int nline);

#endif
