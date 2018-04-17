#ifndef utbin_h
#define utbin_h

/* Define a struct for the return value of bintime() */

typedef struct {
  double beg_ut;   /* The start UT of the bin */
  double mid_ut;   /* The central UT of the bin */
  double end_ut;   /* The end UT of the bin */
} UTbin;

/* A utility function for defining averaging bins. */

UTbin *bintime(double origin, double ut, double binwid);

#endif
