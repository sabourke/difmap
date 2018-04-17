#ifndef telspec_h
#define telspec_h

/* Enumerate telescope iterator operations */

typedef enum {FIND_FIRST, FIND_NEXT, SKIP_SUB, SKIP_TA, SKIP_TB, SKIP_TC} Findop;

/* Sub-array specification container */

typedef struct {
  short nfix;      /* The number of items explicitly specified (0 or 1) */
  short isub;       /* Sub-array index [specified or alluded] */
} Subspec;

Subspec *read_Subspec(Observation *ob, char *s, char **endp, int d_sub);
int write_Subspec(Observation *ob, Subspec *ss, int nref, int fixref,
		  int n, char *s);
Subspec *find_sub(Observation *ob, int nfix, int isub, int forward,
		  int nref, int fixref, int report);
int next_sub(Observation *ob, Findop oper, int forward, int nref, int fixref,
	     int report, Subspec *ss);

/* Telescope specification container */

typedef struct Telspec {
  short nfix;      /* The number of items explicitly specified (0 - 2) */
  short isub;       /* Sub-array index [specified or alluded] */
  short ta;         /* Index of telescope [specified or alluded] */
} Telspec;

Telspec *read_Telspec(Observation *ob, char *s, char **endp, int d_sub);
int write_Telspec(Observation *ob, Telspec *ts, int nref, int fixref,
		  int n, char *s);
Telspec *find_tel(Observation *ob, int nfix, int isub, int ta, int forward,
		  int nref, int fixref, int report);
int next_tel(Observation *ob, Findop oper, int forward, int nref, int fixref,
	     int report, Telspec *ts);

/* Baseline specification container */

typedef struct Basespec {
  short nfix;      /* The number of items explitictly specified (0-3) */
  short isub;       /* Sub-array index [specified or alluded] */
  short ta;         /* Index of the first telescope [specified or alluded] */
  short tb;         /* Index of the other telescope of baseline 'base' wrt ta */
  short base;       /* Index of the first baseline selected */
} Basespec;

/* Decode a baseline specification from a string or from stdin */

Basespec *read_Basespec(Observation *ob, char *s, char **endp, int d_sub);

/* Write a baseline specification string suitable for use by read_Basespec() */

int write_Basespec(Observation *ob, Basespec *bs, int nref, int fixref,
		   int n, char *s);

/* Search forward or backwards for a baseline meeting given requirements */

Basespec *find_base(Observation *ob, int nfix, int isub, int ta, int tb,
		    int forward, int nref, int allref, int fixref, int report);

/* Baseline iterator */

int next_base(Observation *ob, Findop oper, int forward, int nref, int allref,
	      int fixref, int report, Basespec *bs);

/* Find the index of the baseline associated with specified telescope */

int loc_base(Subarray *sub, int tel_a, int tel_b);

/* Container for return of closure triangle info from find_cls() */

typedef struct Trispec {
  short nfix;  /* The number of items explitictly specified (0-4) */
  short isub;   /* The index of the sub-array containing the triangle */
  short ta;     /* Index of the 1st telescope */
  short tb;     /* Index of the 2nd telescope tb > ta */
  short tc;     /* Index of the 3rd telescope tc > tb */
  struct {      /* Triangle baselines ta-tb, tb-tc, tc-ta */
    short base; /* Index of baseline */
    short sign; /* Sign of phase contribution */
  } b[3];
} Trispec;

Trispec *read_Trispec(Observation *ob, char *s, char **endp, int d_sub);

int write_Trispec(Observation *ob, Trispec *ts, int nref, int fixref,
		  int n, char *s);

Trispec *find_tri(Observation *ob, int nfix, int isub, int ta, int tb, int tc,
		  int forward, int nref, int allref, int fixref, int report);
int next_tri(Observation *ob, Findop oper, int forward, int nref, int allref,
	     int fixref, int report, Trispec *ts);

#endif
