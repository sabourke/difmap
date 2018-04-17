#ifndef visplot_h
#define visplot_h

/* Observation visibility plotting functions */

  /* Ampitude and/or phase sampling versus radius or radial vector */

int uvradplt(Observation *ob, Telspec *ts, int docurs, char *opts,
	     int doproj, float phi, float uvmin, float uvmax,
	     float ampmin, float ampmax, float phsmin, float phsmax,
	     int *modified);

  /* Telescope amplitude and/or phase corrections versus time */

int corplot(Observation *ob, Telspec *ts, int cif, int docurs, int *modified);

  /* Baseline amplitude and/or phase versus time */

int vedit(Observation *ob, Basespec *bs, int cif, int nrow, int npage,
	  int docurs, char *opts, int doscan, int doamp, int dophs,
	  int doflag, int domod, int dobars, int showall, int *modified);

  /* Sampling of the UV plane */

int uvplot(Observation *ob, Telspec *ts, int docurs, char *opts,
	   float umax, float vmax, int *modified);

  /* Sampling of stations vs. time */

int timplt(Observation *ob, Subspec *ss, int cif, int docurs,
	   char *opts, int *modified);

  /* Closure-phase plotting */

int clsplot(Observation *ob, Trispec *ts, int cif, int nrow, int npage,
	    int docurs, char *opts, int doscan, int doflag, int domod,
	    int dobars, int *modified);

#endif
