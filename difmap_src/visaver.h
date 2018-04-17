#ifndef visaver_h
#define visaver_h

/* Visibility averager class */

typedef struct Visaver Visaver;

/* Construct a new visibility averager object */

Visaver *new_Visaver(Dpage *dp, double avtime, int scatter);

/* Initialize to start a new output averaged integration */

int av_newint(Visaver *av, Visibility *vis, int nbase, long irec);

/* Accumulate visibility data weighted running means */

int av_dp(Visaver *av, float re, float im, float wt, long ivis);

/* Accumulate visibility U,V,W weighted running means */

int av_uvwt(Visaver *av, float uu, float vv, float ww, float wt, float inttim,
	    int base);

/* Complete the averaging of an integration */

int av_endint(Visaver *av);

/* Delete a visibility averager instance */

Visaver *del_Visaver(Visaver *av);

/*
 * Define a class to record the statistics required to calculate
 * the scatter of a single output visibility.
 */
typedef struct {
  long nsum;       /* Number of points in sum */
  float sqr_mean;  /* Sum of squared imaginary + squared real parts */
} Scatsum;

/*
 * Class to record intermediate sums required involved in calculating
 * the average U, V, and W coordinates of a single baseline.
 */
typedef struct {
  float wtsum;    /* Sum of weights used in weighted averages */
} Basesum;

/* Define an average object to retain the averaging state */

struct Visaver {
  long nvis;         /* The number of visibilities per integration record */
  int nbmax;         /* The max number of baselines per integration */
  int nbase;         /* The number of baselines in the current integration */
  Visibility *vis;   /* The array of nbase visibility descriptors */
  Scatsum *scatsum;  /* Array of 'nvis' visibility scatter sums */
  Basesum *basesum;  /* Array of 'nbmax' baseline averaging sums */
  Dpage *dp;         /* The output uvdata.scr file */
};

#endif
