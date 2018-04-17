#ifndef clphs_h
#define clphs_h

/* Enumerate closure-phase flag states */

typedef enum {
  FLAG_CDEL=1,  /* Closure phase deleted */
  FLAG_CBAD=2,  /* Closure phase flagged */
  FLAG_CCOR=4   /* CLosure phase selfcal-correction flagged */
} Clflag;

/* Container to construct closure phase in */

typedef struct {
  float wt;       /* Weight of closure phase 1/variance */
  float ophs;     /* Observed closure phase (radians) */
  float mphs;     /* Model closure phase (radians) */
  int bad;        /* Flag status encoded as a bitmask union of Clflag values */
} Clphs;

Clphs *get_clphs(Trispec *ts, Visibility *vis);

#endif
