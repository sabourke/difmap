/* Ensure that obs.h is included before this */

int slfcal(Observation *ob, int isub, int doall, float gauval, float gaurad,
	   float solint, int doamp, int dophs, int dofloat, int mintel,
	   int doflag, int doone, float maxamp, float maxphs,
	   float uvmin, float uvmax, int *flagged);
