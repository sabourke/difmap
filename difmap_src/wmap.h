#ifndef wmap_h
#define wmap_h

/* Function used to write a map or beam to a new FITS file */

int w_MapBeam(Observation *ob, MapBeam *mb, int domap, const char *fname);

#endif
