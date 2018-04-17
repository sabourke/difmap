#ifndef scans_h
#define scans_h

/* Count the number of scans in a sub-array */

int nscans(Subarray *sub, double tsep);

/* Find the end of the next scan that follows integration uta of a sub-array */

int endscan(Subarray *sub, double tsep, int uta);

/* Determine the sum of the durations of all scans in a sub-array */

int timescans(Subarray *sub, double tsep);

/* Set the interscan delimiting gap (seconds) */

int scangap(Observation *ob, double gap, int isub);

/* Set the default interscan gap (seconds) */

#define DEFGAP 3600.0

#endif
