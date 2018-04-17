#ifndef visstat_h
#define visstat_h

/*
 * The ob_vis_stats() function fills containers of the following form
 * with the statistics of a specified observable.
 */
typedef struct {
  int nvis;         /* The number of visibilities used */
  double mean;      /* The mean of the observable */
  double sigma;     /* The standard error on the mean */
  double scatter;   /* The RMS scatter about the mean */
  double minval;    /* The minimum value of the observable */
  double maxval;    /* The maximum value of the observable */
} VisStat;

/*
 * Enumerate the types of observables that ob_vis_stats() recognizes.
 */
typedef enum {
  VS_AMP,         /* Visibility amplitudes */
  VS_PHS,         /* Visibility phases (radians) */
  VS_REAL,        /* The real parts of the visibilities */
  VS_IMAG,        /* The imaginary parts of the visibilities */
  VS_UMAG,        /* The magnitude of the U coordinate of the visibilities */
  VS_VMAG,        /* The magnitude of the V coordinate of the visibilities */
  VS_UVRAD        /* The UV radii of the visibilities */
} VisStatQty;

int ob_vis_stats(Observation *ob, VisStatQty qty, float uvmin, float uvmax,
		 VisStat *results);

#endif
