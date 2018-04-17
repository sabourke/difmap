#ifndef fitkey_h
#define fitkey_h

/* FITS keyword value pair descriptor */

typedef struct Fitkey {
  char name[9];     /* Keyword name (Use "" for wildcard when reading) */
  int extn;         /* Postfixed keyword extension number (0 = no extension) */
  int keyid;        /* Use to enumerate keyword ID. */
  Fittype type;     /* Type of value required/received */
  void *value;      /* Pointer to keyword value */
  char *comment;    /* Header comment */
} Fitkey;

/* Macros for extracting values from the Fitkey value field */

#define KEYINT(a)  *((int *)(a).value)
#define KEYBOOL(a)  *((char *)(a).value)
#define KEYSTR(a)  ((char *)(a).value)
#define KEYFLT(a)  *((float *)(a).value)
#define KEYDBL(a)  *((double *)(a).value)
#define KEYCMP(a)  ((float *)(a).value)

/* Union for holding keyword values */

typedef union {
  double dval;   /* Double precision */
  float fval;    /* Float */
  int ival;      /* Integer */
  float cval[2]; /* Complex */
  char sval[81]; /* String */
  char bval;     /* Boolean ('T' or 'F') */
} Keyval;

typedef enum {  /* Return enumeration from key-value reading functions */
  KEY_FOUND=0,  /* Keyword matches one of those given */
  KEY_EOH,      /* No more header lines in header */
  KEY_UNKNOWN,  /* No match with given keywords */
  KEY_BAD       /* Error while reading or decomposing a header line */
} Keystat;

typedef enum {  /* Search request enumeration */
  NO_SEEK=0,    /* No search - return the next key in the header */
  EOH_SEEK,     /* Search up to the end of header (EOH) if necessary */
  LOOP_SEEK     /* If not found before EOH re-search from start of header */
} Seektype;

Keystat get_key(Fits *fits, Hdu *hdu, char *match, Fittype type,
		Seektype doseek, Fitkey *key);
Keystat next_key(Fits *fits, Hdu *hdu, Fitkey *keys, int nkey, Seektype doseek,
	     Fitkey *key);
Keystat read_key(Fits *fits, Hdu *hdu, Fitkey *keys, int nkey, Fitkey *key);

int putkey(Fits *fits, Hdu *hdu, Fitkey *key);

int wintkey(Fits *fits, Hdu *hdu, char *name, int extn, int value, char *comment);
int wfltkey(Fits *fits, Hdu *hdu, char *name, int extn, double value, char *comment);
int wlogkey(Fits *fits, Hdu *hdu, char *name, int extn, char value, char *comment);
int wcmpkey(Fits *fits, Hdu *hdu, char *name, int extn, float *value, char *comment);
int wstrkey(Fits *fits, Hdu *hdu, char *name, int extn, char *value, char *comment);
int wcomkey(Fits *fits, Hdu *hdu, char *name, int extn, char *value, char *comment);
int wvoidkey(Fits *fits, Hdu *hdu, char *name, int extn, char *comment);

#define end_hline(hdu)  (hdu)->endline
#define what_hline(hdu) (hdu)->nextline
int new_hline(Hdu *hdu, int iline);

#endif
