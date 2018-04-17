#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#include "sysfits.h"
#include "recio.h"
#include "fits.h"
#include "utils.h"
#include "fitkey.h"

/* Set up the HDU virtual function table */

extern Hdutab phdufns;
extern Hdutab ahdufns;
extern Hdutab bhdufns;
extern Hdutab uhdufns;

static struct Hdufns {
  Hdutype type;         /* The enumerated HDU type */
  char *name;           /* The name to recognise an HDU by */
  Hdutab *fns;          /* HDU type specific functions for the given HDU */
} hdufns[]={
  {F_UNKNOWN, "UNKNOWN", &uhdufns},
  {F_PRIMARY, "PRIMARY", &phdufns},
  {F_IMAGE,   "IMAGE",   &phdufns},
  {F_TABLE,   "TABLE",   &ahdufns},
  {F_BINTAB,  "BINTABLE",&bhdufns},
  {F_BINTAB,  "A3DTABLE",&bhdufns},
};

static int size_hdu(Hdu *hdu, int headrec);
static int len_data(Hdu *hdu);
static int len_header(Hdu *hdu);
static Hdutab *lookhdu(Hdutype type);
static Hdutype lookext(const char *xtension);
static long grp_size(Hdu *hdu);
static int adderr(Fits *fits, Hdu *hdu);
static Hdu *duperr(Fits *fits, Hdu *hdu);

/*.......................................................................
 * Look up an entry in hdutab[].
 *
 * Input:
 *  type   Hdutype   The enumerated HDU derived type.
 * Output:
 *  return  Hdutab * The corresponding HDU table entry or NULL on error.
 */
static Hdutab *lookhdu(Hdutype type)
{
  int i;
  for(i=0; i<sizeof(hdufns)/sizeof(struct Hdufns); i++) {
    if(type==hdufns[i].type)
      return hdufns[i].fns;
  };
  fprintf(stderr, "lookhdu: System error: un-handled HDU type\n");
  return NULL;
} 

/*.......................................................................
 * Look up the type of an HDU extension by its XTENSION keyword value.
 * Trailing spaces in the input string are ignored.
 *
 * Input:
 *  xtension  char * The extension type name to be looked up.
 * Output:
 *  return Hdutype   The corresponding HDU type - F_UNKNOWN if not found.
 */
static Hdutype lookext(const char *xtension)
{
  int i;
  const char *aptr; /* Pointer into 'xtension' */
  const char *bptr; /* Pointer into string being compared to 'xtension'. */
  for(i=0; i<sizeof(hdufns)/sizeof(struct Hdufns); i++) {
    aptr = xtension;
    bptr = &hdufns[i].name[0];
    while(*aptr && *bptr && *aptr == *bptr && *aptr != ' ') {
      aptr++;
      bptr++;
    };
    if(*bptr=='\0' && (*aptr=='\0' || *aptr==' '))
      return hdufns[i].type;
  };
  return F_UNKNOWN;
}

/*.......................................................................
 * Allocate memory for a new HDU of the appropriate derived type.
 * The only initialization performed is to record the HDU type and
 * NULL all pointers of the derived type so that the descriptor can
 * thereafter be deleted with del_Hdu().
 *
 * Input:
 *  type  Hdutype   The derived specialization of the base class HDU
 *                  descriptor. Recognised types:
 *                   F_PRIMARY   -   Primary header.
 *                   F_IMAGE     -   IUE IMAGE extension.
 *                   F_TABLE     -   ASCII table extension.
 *                   F_BINTAB    -   Binary table extension.
 * Output:
 *  return    Hdu * Generic pointer to the allocated (unitialized)
 *                  HDU descriptor, or NULL on memory allocation error.
 */
Hdu *new_Hdu(Hdutype type)
{
  Hdu *hdu;       /* Base class pointer to new descriptor */
  Hdutab *htab;   /* The virtual function-table entry for the HDU type */
/*
 * Look-up the HDU entry in hdutab.
 */
  htab = lookhdu(type);
  if(htab==NULL)
    return NULL;
/*
 * Allocate the HDU descriptor.
 */
  hdu = (Hdu *) malloc((*htab->sizefn)());
  if(hdu==NULL) {
    fprintf(stderr, "Insufficient memory for new HDU\n");
    return hdu;
  };
/*
 * Initialize the base-class parts of the descriptor.
 */
  hdu->type = type;
  hdu->dims = NULL;
  hdu->groups = 0;
  hdu->pcount = 0;
  hdu->gcount = 1;
  hdu->nrec = 0;
  hdu->headrec = hdu->datarec = 0;
  hdu->wnxtline = hdu->nextline = 0;
  hdu->endline = -1;
  hdu->pad = 0;
  hdu->grpsize = 0;
  hdu->complete = 0;
  hdu->state = HDU_DESCR;
  hdu->extname = NULL;
  hdu->extver = 0;
  hdu->extlevel = 1;
  hdu->next = NULL;
/*
 * Clear the derived parts of the HDU.
 */
  return (*htab->newfn)(hdu);
}

/*.......................................................................
 * Initialize the base class part of an HDU descriptor returned from
 * new_Hdu(). This includes allocating the hdu->dims[] array.
 *
 * Input:
 *  hdu       Hdu *  The base class pointer returned by new_Hdu().
 *  bitpix Bitpix    The number of bits per pixel in the data section.
 *                    B_CHAR  - Unsigned IEEE character or integer.
 *                    B_16INT - 16-bit IEEE integer.
 *                    B_32INT - 32-bit IEEE integer.
 *                    B_FLOAT - 32-bit IEEE single precision float.
 *                    B_DBLE  - 64-bit IEEE double precision float.
 *  dims      int *  An 'naxis' element array containing the dimensions
 *                   of each axis. Send NULL if this has already been
 *                   allocated and assigned to hdu->dims. In the latter
 *                   case the assigned hdu->dims array must have been
 *                   allocated via malloc since del_Hdu() will free it.
 *  naxis     int    The number of elements in 'dims'.
 *  groups    int    If true then this HDU contains random groups.
 *  pcount    int    Random parameter count - send zero if not relevant.
 *  gcount    int    Group count - send 1 if not relevant.
 *  extname  char *  The name to give this HDU. Send NULL if not known
 *                   or hdu->extname has already been assigned.
 *  extver    int    The version number of this extension.
 *  extlevel  int    Heirachical level - send 1 if not relevant.
 *  headrec   int    The FITS record number of the start of the HDU.
 *  endline   int    The required index of the END line.
 *                   This sets the header length.
 * Output:
 *  return    Hdu *  Base class pointer to the initialized descriptor, or NULL
 *                   if an error occurred and the descriptor was deleted.
 */
Hdu *ini_Hdu(Hdu *hdu, Bitpix bitpix, int *dims, int naxis,
	     int groups, int pcount, int gcount, char *extname, int extver,
	     int extlevel, int headrec, int endline)
{
  int i;
  if(hdu==NULL) {
    fprintf(stderr, "ini_Hdu: NULL Hdu descriptor recieved\n");
    return hdu;
  };
/*
 * Check and assign BITPIX.
 */
  switch(bitpix) {
  case B_CHAR: case B_16INT: case B_32INT: case B_FLOAT: case B_DBLE:
    hdu->bitpix = bitpix;
    break;
  default:
    fprintf(stderr, "ini_Hdu: Bad BITPIX value recieved: %d\n", (int) bitpix);
    return del_Hdu(hdu);
    break;
  };
/*
 * Initialize hdu->naxis and hdu->dims.
 */
  if(naxis<0 || naxis>999) {
    fprintf(stderr, "ini_Hdu: Illegal value of NAXIS: %d\n", naxis);
    return del_Hdu(hdu);
  };
  if(dims==NULL && hdu->dims==NULL && naxis!=0) {
    fprintf(stderr, "ini_Hdu: No \'dims\' array either sent or assigned\n");
    return del_Hdu(hdu);
  };
  if(naxis<0) {
    fprintf(stderr, "ini_Hdu: Illegal negative value of naxis\n");
    return del_Hdu(hdu);
  } else {
    hdu->naxis = naxis;
  };
/*
 * Allocate a new dims array?
 */
  if(naxis>0 && dims) {
    hdu->dims = (int *) malloc(naxis * sizeof(int));
    if(hdu->dims==0) {
      fprintf(stderr, "ini_Hdu: Insufficient memory to record dimensions\n");
      return del_Hdu(hdu);
    };
    for(i=0; i<naxis; i++)
      hdu->dims[i] = dims[i];
  };
/*
 * Initialize other members.
 */
  hdu->pcount = pcount;
  hdu->gcount = gcount;
  hdu->groups = groups;
  hdu->wnxtline = 0;
  hdu->nextline = 0;
  hdu->endline = endline;
/*
 * Copy extname if given, or if this is the primary header assign
 * the name "PRIMARY".
 */
  if( (extname && !(hdu->extname = fitsstr(extname))) ||
      (!extname && hdu->type==F_PRIMARY && !(hdu->extname=fitsstr("PRIMARY"))))
    return del_Hdu(hdu);
  hdu->extver = extver;
  hdu->extlevel = extlevel<1 ? 1 : extlevel;
/*
 * Work out the sizes of the hdu and its header in FITS 2880 byte
 * records.
 */
  if(size_hdu(hdu, headrec))
    return del_Hdu(hdu);
/*
 * Return the initialized descriptor.
 */
  return hdu;
}

/*.......................................................................
 * Determine and record the various sizes and record offsets that are
 * required before reading and writing can commence.
 *
 * Input:
 *  hdu      Hdu *  The descriptor to be set up.
 *  headrec  int    The record number of the header. datarec will be
 *                  calculated relative to this.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
static int size_hdu(Hdu *hdu, int headrec)
{
  int hlen;  /* Length of header in FITS records */
/*
 * Work out the sizes of the hdu and its header in FITS 2880 byte
 * records.
 */
  hdu->grpsize = grp_size(hdu);
  hdu->headrec = headrec;
  hlen = len_header(hdu);
  hdu->datarec = hdu->headrec + hlen;
  hdu->nrec = hlen + len_data(hdu);
  return 0;
}

/*.......................................................................
 * Release memory used in a Header-Data-Unit descriptor.
 *
 * Input:
 *  hdu     Hdu *  Base-class pointer to a derived Hdu descriptor previously
 *                 allocated in new_Hdu().
 * Output:
 *  return  Hdu *  Allways NULL. Use like:  hdu = del_Hdu(hdu);
 */
Hdu *del_Hdu(Hdu *hdu)
{
  Hdutab *htab;   /* HDU virtual function table entry */
/*
 * NOOP if there is no descriptor.
 */
  if(hdu==NULL)
    return hdu;
/*
 * Handle the base-class parts of the descriptor.
 */
  if(hdu->dims)
    free(hdu->dims);
  if(hdu->extname)
    free(hdu->extname);
/*
 * Delete the derived parts of the descriptor.
 */
  htab = lookhdu(hdu->type);
  if(htab!=NULL)
    (*htab->delfn)(hdu);
/*
 * Delete the container.
 */
  free(hdu);
  return NULL;
}

/*.......................................................................
 * Determine and record the characteristics of the next HDU in the FITS
 * file. If 'headrec' is 0 then the HDU will be assumed to be a primary
 * HDU. Otherwise an extension HDU will be expected.
 *
 * Input:
 *  fits     Fits *  The descriptor of the FITS file.
 *  headrec   int    The FITS record number where the new HDU starts.
 * Output:
 *  return    Hdu *  A pointer to the new Hdu descriptor or NULL if no valid
 *                   HDU was found.
 */
Hdu *get_Hdu(Fits *fits, int headrec)
{
  Hdu base;      /* Temporary base-class HDU */
  Hdu *hdu=NULL; /* Pointer to an HDU descriptor. */
  Hdutab *htab;  /* hdutab[] entry for derived HDU type */
  Fitkey key;    /* FITS keyword descriptor */
  int primary;   /* If true a primary header is required */
  int saveline;  /* Temporary storage of line number in header */
  int i;
/*
 * Optional keys.
 */
  enum {HDU_EXTNAME, HDU_EXTVER, HDU_EXTLEVEL, HDU_GROUPS, HDU_PCOUNT,
	HDU_GCOUNT};
  Fitkey opkeys[]={
    {"GROUPS",   0, HDU_GROUPS,   DAT_LOG, NULL, NULL},
    {"PCOUNT",   0, HDU_PCOUNT,   DAT_INT, NULL, NULL},
    {"GCOUNT",   0, HDU_GCOUNT,   DAT_INT, NULL, NULL},
    {"EXTNAME",  0, HDU_EXTNAME,  DAT_STR, NULL, NULL},
    {"EXTVER",   0, HDU_EXTVER,   DAT_INT, NULL, NULL},
    {"EXTLEVEL", 0, HDU_EXTLEVEL, DAT_INT, NULL, NULL},
  };
/*
 * Must the new header be a primary header?
 */
  primary = (headrec==0);
/*
 * Initialize a temporary HDU until we have sufficient information to
 * create a new one.
 */
  base.headrec = headrec;
  base.nrec = 0;
  base.endline = -1;
  base.nextline = 0;
  base.state = HDU_INFILE;
/*
 * Attempt to read the first record of the new HDU.
 */
  if(fits_read(fits, headrec, 0)) {
    if(primary)
      fprintf(stderr, "get_Hdu: Error reading first record of FITS file: %s\n",
	      fits->name);
    return NULL;
  };
/*
 * Primary header must have SIMPLE keyword.
 */
  if(primary) {
    if(get_key(fits, &base, "SIMPLE", DAT_LOG, NO_SEEK, &key) ||
       KEYBOOL(key)!='T') {
      fprintf(stderr, "get_Hdu: Non-standard FITS file\n");
      return NULL;
    };
    base.type = F_PRIMARY;
  }
/*
 * Extension headers must start with the XTENSION keyword.
 * If not found - assume that there is no new HDU to be read.
 */
  else {
    if(get_key(fits, &base, "XTENSION", DAT_STR, NO_SEEK, &key)!=KEY_FOUND) {
      rec_rewind(fits->rec);  /* Clear error status */
      return NULL;
    };
/*
 * Determine what type of extension follows, from the string value associated
 * with the XTENSION keyword.
 */
    base.type = lookext(KEYSTR(key));
    if(base.type==F_UNKNOWN) {
      fprintf(stderr, "get_Hdu: Warning: Unrecognised XTENSION=%s\n",
	      KEYSTR(key));
    };
  };
/*
 * Allocate the new HDU descriptor.
 */
  hdu = new_Hdu(base.type);
  if(hdu==NULL)
    return hdu;
  hdu->headrec = base.headrec;
  hdu->nextline = base.nextline;
  hdu->state = HDU_INFILE;
/*
 * The next keyword should be BITPIX.
 */
  if(get_key(fits, hdu, "BITPIX", DAT_INT, NO_SEEK, &key)) {
    fprintf(stderr, "Mandatory BITPIX keyword not in header\n");
    return NULL;
  } else {
    hdu->bitpix = (Bitpix) KEYINT(key);
  };
/*
 * The next keyword should be NAXIS.
 */
  if(get_key(fits, hdu, "NAXIS", DAT_INT, NO_SEEK, &key) || key.extn!=0) {
    fprintf(stderr, "Mandatory NAXIS keyword not in header\n");
    return NULL;
  };
  hdu->naxis = KEYINT(key);
/*
 * Allocate an array to store the NAXISn dimensions.
 */
  if(hdu->naxis > 0) {
    hdu->dims = (int *) malloc(sizeof(int) * hdu->naxis);
    if(hdu->dims==NULL) {
      fprintf(stderr, "Insufficient memory to record axis dimensions\n");
      return del_Hdu(hdu);
    };
  };
/*
 * Read the NAXISn keywords into hdu->dims[n-1].
 */
  for(i=1; i<=hdu->naxis; i++) {
    if(get_key(fits, hdu, "NAXIS", DAT_INT, NO_SEEK, &key) || key.extn != i) {
      fprintf(stderr, "Missing NAXIS%d in FITS header\n", i);
      return del_Hdu(hdu);
    };
    hdu->dims[i-1] = KEYINT(key);
  };
/*
 * Get the rest of the keywords that are needed by ini_hdu().
 * NOTE. ini_Hdu() requires hdu->endline to be initialized correctly.
 * The EOH_SEEK below accomplishes this as a side effect.
 */
  saveline = hdu->nextline;
  while(next_key(fits, hdu, opkeys, sizeof(opkeys)/sizeof(Fitkey), EOH_SEEK, &key) == KEY_FOUND) {
    switch(key.keyid) {
    case HDU_GROUPS:
      hdu->groups = KEYBOOL(key)=='T';
      break;
    case HDU_PCOUNT:
      hdu->pcount = KEYINT(key);
      break;
    case HDU_GCOUNT:
      hdu->gcount = KEYINT(key);
      break;
    case HDU_EXTNAME:
      hdu->extname = fitsstr(KEYSTR(key));
      break;
    case HDU_EXTVER:
      hdu->extver = KEYINT(key);
      break;
    case HDU_EXTLEVEL:
      hdu->extlevel = KEYINT(key);
      break;
    default:
      break;
    };
  };
/*
 * Now we know enough to initialize the base-class part
 * of the HDU descriptor.
 */
  hdu = ini_Hdu(hdu, hdu->bitpix, NULL, hdu->naxis, hdu->groups,
		hdu->pcount, hdu->gcount, NULL, hdu->extver,
		hdu->extlevel, hdu->headrec, hdu->endline);
  if(hdu==NULL)
    return hdu;
/*
 * Now initialize the derived parts of the HDU descriptor.
 */
  hdu->nextline = saveline;
  htab = lookhdu(hdu->type);
  if(htab==NULL || (*htab->getfn)(fits, hdu))
    return del_Hdu(hdu);
/*
 * Make sure that the index that records the position beyond which
 * records are blank in the file, points beyond the newly read HDU.
 * This assumes that the new HDU contains data.
 */
  if(fits->nullrec < hdu->headrec + hdu->nrec)
    fits->nullrec = hdu->headrec + hdu->nrec;
/*
 * Signal that this HDU is now ready for use and that the HDU descriptor
 * is hereafter read-only.
 */
  hdu->complete = 1;
/*
 * Return the new descriptor.
 */
  return hdu;
}

/*.......................................................................
 * Determine the number of records in the data segment of an HDU, as
 * implied by the NAXIS, NAXISn, PCOUNT and GCOUNT keyword information
 * stored in an HDU descriptor.
 *
 * Input:
 *  hdu    Hdu *  The descriptor of the HDU whose length is required.
 * Output:
 *  return int    The number of 2880 byte records used by the data segment.
 */
static int len_data(Hdu *hdu)
{
  long nbytes;  /* Number of bytes used in the data section of the HDU */
/*
 * Determine the number of bytes used in the data segment of the HDU.
 */
  nbytes = hdu->grpsize * hdu->gcount;
/*
 * Return the number of records used.
 */
  return (int) ((nbytes+2879L)/2880L);
}

/*.......................................................................
 * Determine the number of records in the header segment of an HDU, as
 * implied by the number of lines recorded in hdu->endline.
 *
 * Input:
 *  hdu    Hdu *  The descriptor of the HDU whose length is required.
 * Output:
 *  return int    The number of 2880 byte records used by the header
 *                segment.
 */
static int len_header(Hdu *hdu)
{
  long nbytes;  /* Number of bytes used in the header section of the HDU */
/*
 * Determine the number of bytes used in the header segment of the HDU.
 */
  if(hdu->endline < 0) {
    fprintf(stderr, "len_header: hdu->endline not initialized\n");
    return 0;
  };
/*
 * How many bytes are used up to and including the END line.
 */
  nbytes = (1+hdu->endline) * 80L;  /* 80 chars per line */
/*
 * Return the number of records used.
 */
  return (int) ((nbytes+2879L)/2880L);
}

/*.......................................................................
 * Convert from Bitpix to Fittype enumerators.
 *
 * Input:
 *  hdu         Hdu *  The descriptor of the HDU whose data-type is to be
 *                     found.
 * Output:
 *  return  Fittype    The data type of the HDU.
 */
Fittype dat_type(Hdu *hdu)
{
  switch(hdu->bitpix) {
  case B_CHAR:
    return DAT_CHR;
  case B_16INT:
    return DAT_SHT;
  case B_32INT:
    return (CHAR_BIT*sizeof(int)>=32) ? DAT_INT : DAT_LNG;
  case B_FLOAT:
    return DAT_FLT;
  case B_DBLE:
    return DAT_DBL;
  };
  fprintf(stderr, "dat_type: Unknown BITPIX type found\n");
  return DAT_NON;
}

/*.......................................................................
 * Compute the number of bytes per group ie the number of bytes between
 * the start of one group in the Hdu to the start of the next. This is
 * the sum of the number of bytes in the random-group section and the
 * number in the image section. If the header is not a random group
 * header then this is the number of bytes in the image.
 *
 * Input:
 *  hdu      Hdu *  The Hdu to be characterized.
 * Output:
 *  long  return    The number of bytes per group.
 */
static long grp_size(Hdu *hdu)
{
  long nbytes=1L;  /* Number of bytes used in the group */
  int i;
/*
 * Determine the number of bytes used in one group of the HDU.
 */
  if(hdu->naxis==0 || hdu->dims == NULL) {
    nbytes = 0L;
  } else if(hdu->groups) {
    for(i=1; i<hdu->naxis; i++)
      nbytes *= hdu->dims[i];
  } else {
    for(i=0; i<hdu->naxis; i++)
      nbytes *= hdu->dims[i];
  };
/*
 * Return the number of bytes.
 */
  return (abs((int) hdu->bitpix)/8L) * (hdu->pcount + nbytes);
}

/*.......................................................................
 * Add an HDU to a FITS file.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file to which the HDU is
 *                  to be added.
 *  hdu      Hdu *  An initialized HDU descriptor for the new HDU.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
ADDFN(add_Hdu)
{
  Hdu *tmphdu=NULL; /* An HDU currently in the FITS file */
  Hdutab *htab;     /* HDU virtual function table entry */
  int primary;      /* If true then the HDU is a primary HDU */
  int i;
/*
 * Sanity check the arguments.
 */
  if(fits==NULL || hdu==NULL) {
    fprintf(stderr, "add_Hdu: NULL %s descriptor intercepted\n",
	    fits==NULL ? "FITS" : "HDU");
    return 1;
  };
  if(hdu->state!=HDU_DESCR) {
    fprintf(stderr,
    "add_Hdu: Attempt to append aliased HDU from this (%s) or another file\n",
	    fits->name);
    return 1;
  };
/*
 * Was this FITS file opened for writing?
 */
  if(fits->readonly) {
    fprintf(stderr, "add_Hdu: FITS file not opened for writing\n");
    return 1;
  };
/*
 * Look up the HDU type.
 */
  htab = lookhdu(hdu->type);
  if(htab==NULL)
    return 1;
/*
 * Must the new header be a primary header?
 */
  primary = (fits->hdu==NULL);
  if(primary && hdu->type != F_PRIMARY) {
    fprintf(stderr, "add_Hdu: PRIMARY HDU required but not given\n");
    return 1;
  };
/*
 * Unknown HDU types can not be appended.
 */
  if(hdu->type == F_UNKNOWN) {
    fprintf(stderr, "add_Hdu: HDUs of unknown-type can not be written to FITS files\n");
    return 1;
  };
/*
 * Determine and record the FITS record at which the new HDU starts.
 */
  if(primary)                  /* Get primary header from next record */
    hdu->headrec = 0;
  else {                       /* Locate last HDU read */
    for(tmphdu=fits->hdu; tmphdu->next != NULL; tmphdu = tmphdu->next);
/*
 * Was that HDU complete?
 */
    if(!fits->complete || tmphdu->state!=HDU_INFILE) {
      fprintf(stderr, "add_Hdu: Call end_data() on previous HDU before appending a new one.\n");
      return 1;
    };
/*
 * Determine the start record of the new HDU header.
 */
    hdu->headrec = (tmphdu->headrec+tmphdu->nrec);
  };
/*
 * No header lines have been written yet.
 */
  hdu->endline = 0;
  hdu->nrec = 0;
/*
 * Force read_fits() to return blank header records, by ensuring that
 * it thinks that the new HDU is being appended at the end of the file.
 * This is necessary both because an existing file might contain
 * non-empty records after the final HDU, and because if a previous
 * add_Hdu() call failed it might have left some garbage beyond the end of
 * the last valid HDU.
 */
  fits->nullrec = hdu->headrec;
/*
 * The FITS file now contains an incomplete HDU. Mark this in the FITS
 * descriptor. This flag will not get reset until end_data() is called to
 * complete the HDU. It is used by fits_read() to determine whether to
 * attempt reads at the current EOF.
 */
  fits->complete = 0;
/*
 * Header records must be padded with ASCII blanks (ASCII 32).
 */
  fits->pad = 32;
/*
 * Ensure that the header line pointers are intialized.
 */
  hdu->nextline = hdu->wnxtline = 0;
/*
 * The HDU will not be ready for data access until the size of the
 * header becomes known when end_header() is called.
 */
  hdu->state = HDU_HEADER;
/*
 * Write the mandatory generic HDU keywords into the header.
 */
  if(primary) {
    if(wlogkey(fits, hdu, "SIMPLE", 0, 'T', "Standard FITS file"))
      return adderr(fits, hdu);
  } else {
    if(wstrkey(fits, hdu, "XTENSION", 0, fits->aips ? htab->aips:htab->name,
	       "FITS extension type"))
      return adderr(fits, hdu);
  };
  if(wintkey(fits, hdu, "BITPIX", 0, hdu->bitpix, "FITS data type"))
    return adderr(fits, hdu);
  if(wintkey(fits, hdu, "NAXIS", 0, hdu->naxis, "Dimensionality of array"))
    return adderr(fits, hdu);
  for(i=0; i<hdu->naxis; i++) {
    if(wintkey(fits, hdu, "NAXIS", i+1, hdu->dims[i], NULL))
      return adderr(fits, hdu);
  };
/*
 * Write the derived header keywords.
 */
  if((*htab->addfn)(fits,hdu))
    return adderr(fits, hdu);
/*
 * Append the HDU descriptor to the list of HDUs in the FITS file.
 */
  if(fits->hdu==NULL)
    fits->hdu = hdu;
  else {
    for(tmphdu=fits->hdu; tmphdu->next != NULL; tmphdu = tmphdu->next);
    tmphdu->next = hdu;
  };
/*
 * Terminate the list.
 */
  hdu->next = NULL;
/*
 * Signal that this HDU is ready for use and that the descriptor is
 * now read-only to users.
 */
  hdu->complete = 1;
  return 0;
}

/*.......................................................................
 * Clean up after an error during appending a new HDU to a FITS file.
 * This includes over-writing the start line of its header.
 */
static int adderr(Fits *fits, Hdu *hdu)
{
  hdu->state = HDU_DESCR;
/*
 * Over-write the initial header line with an empty comment line ie a
 * line filled with spaces.
 */
  wcomkey(fits, hdu, "", 0, "", NULL);
/*
 * All HDUs in the FITS file are complete now that the new one has been
 * aborted.
 */
  fits->complete = 1;
  return 1;
}

/*.......................................................................
 * Write the END line to a header to close it and ready the data section
 * of the HDU by filling in its size and record parameters and filling it
 * with ASCII NULs.
 *
 * Input:
 *  fits      Fits *   The descriptor of the FITS file.
 *  hdu        Hdu *   The HDU being processed.
 * Output:
 *  return     int     0 - OK.
 *                     1 - Error.
 */
int end_header(Fits *fits, Hdu *hdu)
{
/*
 * Valid operation?
 */
  if(fits->readonly) {
    fprintf(stderr, "end_header: Readonly file.\n");
    return 1;
  };
/*
 * Has the HDU been appended to a file.
 */
  if(hdu->state==HDU_DESCR) {
    fprintf(stderr, "end_header: Call add_Hdu() before this function.\n");
    return 1;
  };
/*
 * Has the HDU header already been closed?
 */
  if(hdu->next!=NULL || hdu->state!=HDU_HEADER) {
    fprintf(stderr, "end_header: The HDU header has apparently already been ended.\n");
    return 1;
  };
/*
 * Write the END line.
 */
  if(wvoidkey(fits, hdu, "END", 0, NULL))
    return 1;
/*
 * Fill any incomplete records of the header.
 */
  if(fits_flush(fits))
    return 1;
/*
 * Now that the header length is known we can determine the total
 * size of the header and the position of both the header and data
 * sections of the HDU.
 */
  if(size_hdu(hdu, hdu->headrec))
    return 1;
/*
 * Ensure that the first and subsequent records of the new data-section
 * are appropriately padded, byt telling fits_read() that the data section
 * starts at the end of the file.
 */
  fits->nullrec = hdu->datarec;
/*
 * The HDU is now ready for data I/O.
 */
  hdu->state = HDU_DATA;
  fits->pad = hdu->pad;
  return 0;
}

/*.......................................................................
 * Locate the descriptor of an HDU in a FITS file.
 *
 * Input:
 *  fits     Fits *  The descriptor of the FITS file.
 *  type  Hdutype    The enumerated type of the HDU being sought.
 *  extname  char *  The extension name of the HDU being sought. Send NULL
 *                   if any extension of type 'type' will do.
 *  extver    int    The version number of the extension being sought. Use
 *                   0 if the highest is required or -1 if the version
 *                   number is irrelevant.
 *  prev      Hdu *  The HDU to search forwards from. Send NULL to
 *                   start from the start HDU in the file. The first
 *                   hdu checked is prev->next.
 * Output:
 *  return    Hdu *  The descriptor of the requested HDU, or NULL if
 *                   not found.
 */
Hdu *find_hdu(Fits *fits, Hdutype type, char *extname, int extver, Hdu *prev)
{
  Hdu *hdu;       /* Generic HDU descriptor */
  Hdu *last=NULL; /* The last matching HDU descriptor */
  int lastver=0;  /* Record of highest version number found. */
/*
 * Loop through the extension HDUs looking for extensions
 * of the requested type.
 */
  for(hdu=prev?prev->next:fits->hdu; hdu!=NULL; hdu=hdu->next) {
/*
 * Is the type of the latest HDU in the set of requested types?
 */
    if(type!=F_ANY && !(type & hdu->type))
      continue;
/*
 * If a required extension name was sent compare it to that of the new HDU.
 * Ignore trailing spaces in the comparison.
 */
    if(extname && !matchstr(hdu->extname, extname, 0))
      continue;
/*
 * If an exact extension version number was requested check for it here.
 */
    if(extver>0) {
/*
 * Even if an exact match is found, don't return yet because there may be
 * more than one exact match, and the last such instance is probably the
 * most up to date.
 */
      if(hdu->extver==extver)
	last = hdu;
/*
 * Highest version number requested?
 */
    } else if(extver==0 && lastver <= hdu->extver) {
      lastver = hdu->extver;
      last = hdu;
/*
 * Will any version number do?
 */
    } else {
      return hdu;
    };
  };
/*
 * Return the descriptor of the last matching HDU in the file.
 */
  return last;
}

/*.......................................................................
 * Append a copy of an established HDU to the same or another FITS file.
 *
 * Input:
 *  afits    Fits *  The FITS file containing the HDU to be copied.
 *  ahdu      Hdu *  The HDU to be copied.
 *  bfits    Fits *  The FITS file to append the copy to.
 * Output:
 *  return    Hdu *  The pointer to the new descriptor of the copied HDU,
 *                   or NULL on error.
 */
Hdu *dup_Hdu(Fits *afits, Hdu *ahdu, Fits *bfits)
{
  Fitkey key;  /* Key descriptor used for modifying the EXTVER header line */
  Hdu *tail=NULL; /* Tail of list of HDUs in bfits->hdu */
  Hdu *tmphdu; /* Used to find preceding extension of same type */
  Hdu *newhdu; /* The duplicate HDU descriptor */
  int primary; /* True if the new HDU will be a primary HDU */
  int headrec; /* The record of the start of the new header */
  int saveline;/* Temporary storage of header-line write pointer */
  int extver;  /* Extension version number of new HDU */
  int i;
/*
 * Sanity check the arguments.
 */
  if(afits==NULL || ahdu==NULL || bfits==NULL) {
    fprintf(stderr, "dup_Hdu: NULL argument intercepted\n");
    return NULL;
  };
/*
 * Were the two FITS files opened with modes compatible with the copy
 * operation.
 */
  if(bfits->readonly) {
    fprintf(stderr,
	    "dup_Hdu: Destination FITS file \"%s\" only open for reading\n",
	    bfits->name);
    return NULL;
  };
/*
 * Only fully established HDUs can be copied.
 */
  if(ahdu->state!=HDU_INFILE) {
    fprintf(stderr, "dup_Hdu: The HDU to be copied is incomplete.\n");
    return NULL;
  };
/*
 * Must the new HDU be a primary header?
 */
  primary = (bfits->hdu==NULL);
/*
 * If a primary HDU is required - report error if the hew HDU is not
 * a primary HDU.
 */
  if(primary && ahdu->type != F_PRIMARY) {
    fprintf(stderr,
    "dup_Hdu: Can't append a non-primary HDU to a file with no primary HDU.\n");
    return NULL;
  };
/*
 * Check if a primary HDU has been given where an extension HDU is required.
 */
  if(!primary && ahdu->type==F_PRIMARY) {
    fprintf(stderr,
	    "dup_Hdu: Can't duplicate a primary HDU to a non-primary HDU\n");
    return NULL;
  };
/*
 * Determine and record the FITS record at which the new HDU will be
 * appended.
 */
  if(primary)                  /* Primary header starts in first record */
    headrec = 0;
  else {                       /* Locate last HDU read */
    for(tail=bfits->hdu; tail->next != NULL; tail = tail->next);
    headrec = (tail->headrec+tail->nrec);
/*
 * Is the file ready to receive a new HDU?
 */
    if(!bfits->complete || tail->state!=HDU_INFILE) {
      fprintf(stderr, "dup_Hdu: Call end_data() on previous HDU before appending a new one.\n");
      return NULL;
    };
  };
/*
 * Create a new HDU descriptor of the existing HDU.
 */
  newhdu = get_Hdu(afits, ahdu->headrec);
  if(newhdu==NULL)
    return newhdu;
  newhdu->next = NULL;
/*
 * Initialize the position and size info in the new HDU descriptor.
 */
  if(size_hdu(newhdu, headrec))
    return del_Hdu(newhdu);
/*
 * Flush any outstanding data to the output FITS file.
 */
  if(fits_flush(bfits))
    return duperr(bfits, newhdu);
/*
 * Copy the header and data of the original HDU to the new file.
 */
  for(i=0; i<ahdu->nrec; i++) {
/*
 * Read the next record of the input HDU.
 */
    if(fits_read(afits, ahdu->headrec + i, 1))
      return duperr(bfits, newhdu);
/*
 * Copy the record into the output FITS I/O buffer.
 */
    memcpy(bfits->buff, afits->buff, FITSLEN);
/*
 * Write the output record.
 */
    bfits->recnum = newhdu->headrec + i;
    bfits->modified = 1;
    if(fits_flush(bfits))
      return duperr(bfits, newhdu);
  };
/*
 * Determine what the extension number should be for the new HDU.
 */
  if(newhdu->type!=F_PRIMARY) {
    tmphdu = find_hdu(bfits, F_ANY, newhdu->extname, 0, NULL);
    extver = (tmphdu) ? tmphdu->extver+1 : 1;
/*
 * Does that of the new HDU need to be changed?
 */
    if(extver!=newhdu->extver) {
/*
 * Change the descriptor to show the new extension version.
 */
      newhdu->extver = extver;
/*
 * If the EXTVER keyword is found in the new HDU over-write it with the
 * new version.
 */
      if(get_key(bfits, newhdu, "EXTVER", DAT_INT, LOOP_SEEK, &key)==0) {
	saveline = newhdu->wnxtline;
	newhdu->wnxtline = newhdu->nextline-1;
	KEYINT(key) = extver;
	if(putkey(bfits, newhdu, &key))
	  return duperr(bfits, newhdu);
	newhdu->wnxtline = saveline;
/*
 * Write a new EXTVER header line.
 */
      } else {
	if(wintkey(bfits, newhdu, "EXTVER", 0, extver,
		   "Extension version number"))
	  return duperr(bfits, newhdu);
      };
    };
  };
/*
 * Append the HDU descriptor to the linked list of HDUs in bfits->hdu.
 */
  if(tail)
    tail->next = newhdu;
  else
    bfits->hdu = newhdu;
  return newhdu;
}

/*.......................................................................
 * Clean up after an error during appending a duplicated HDU to a FITS file.
 * This includes over-writing the start line of its header and deleting
 * the HDU descriptor.
 */
static Hdu *duperr(Fits *fits, Hdu *hdu)
{
/*
 * Over-write the initial header line with an empty comment line ie a
 * line filled with spaces.
 */
  wcomkey(fits, hdu, "", 0, "", NULL);
  return del_Hdu(hdu);
}

/*.......................................................................
 * Create a deep copy of a given HDU descriptor.
 *
 * NB. The following parameters will NOT be copied and will have the
 *     values given below.
 *
 *  state = HDU_DESCR;
 *  headrec = datarec = 0L;
 *  wnxtline = nextline = endline = 0;
 *  next = NULL;
 *
 * Input:
 *  hdu    Hdu *   The HDU descriptor to be copied.
 * Output:
 *  return Hdu *   The new copy of 'hdu'.
 */
Hdu *copy_Hdu(Hdu *hdu)
{
  Hdutab *htab;     /* The virtual function-table entry for the HDU type */
  Hdu *newhdu;      /* The returned copied HDU descriptor */
/*
 * Get the virtual function table for this HDU type and copy the
 * derived members.
 */
  htab = lookhdu(hdu->type);
  if(htab==NULL)
    return NULL;
/*
 * Call the appropriate copy function for the type of HDU.
 */
  newhdu = (*htab->copyfn)(hdu);
/*
 * Reset parameters that must not be copied.
 */
  if(newhdu) {
    newhdu->state = HDU_DESCR;
    newhdu->headrec = newhdu->datarec = 0L;
    newhdu->wnxtline = newhdu->nextline = 0;
    newhdu->next = NULL;
  };
/*
 * Return the copy.
 */
  return newhdu;
}

/*.......................................................................
 * User function, called to finish the data section of an HDU.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file in which the HDU
 *                  resides.
 *  hdu      Hdu *  The HDU to be completed. This must be the last HDU
 *                  in the FITS file.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int end_data(Fits *fits, Hdu *hdu)
{
  Hdutab *htab;     /* The virtual function-table entry for the HDU type */
/*
 * Check arguments.
 */
  if(fits==NULL || hdu==NULL) {
    fprintf(stderr, "end_data: NULL %s descriptor intercepted.\n",
	    fits==NULL ? "Fits":"Hdu");
    return 1;
  };
/*
 * Is the HDU already complete?
 */
  if(hdu->state == HDU_INFILE)
    return 0;
/*
 * Is the HDU ready to be completed?
 */
  if(hdu->state != HDU_DATA) {
    fprintf(stderr, "end_data: HDU is not ready to be completed.\n");
    return 1;
  };
/*
 * Get the size information of the completed HDU.
 */
  if(size_hdu(hdu, hdu->headrec))
    return 1;
/*
 * Fill any incomplete records of the data section of the HDU.
 */
  if(fits_flush(fits) || fits_pad(fits, hdu->headrec + hdu->nrec))
    return 1;
/*
 * Get the virtual function table for the given HDU.
 */
  htab = lookhdu(hdu->type);
  if(htab==NULL)
    return 1;
/*
 * Call the appropriate method function to have the HDU completed.
 */
  if(htab->endfn(fits, hdu))
    return 1;
/*
 * Mark the newly completed HDU and FITS file as complete.
 */
  hdu->state = HDU_INFILE;
  fits->complete = 1;
  return 0;
}

/*.......................................................................
 * Write the header keywords that record the name, version number and
 * hierachical level of an extension HDU.
 *
 * Input:
 *  fits    Fits *  The descriptor of the FITS file to which the HDU is
 *                  being added.
 *  hdu      Hdu *  The initialized HDU descriptor of the new HDU.
 * Output:
 *  return   int    0 - OK.
 *                  1 - Error.
 */
int w_extkeys(Fits *fits, Hdu *hdu)
{
/*
 * Check the HDU status.
 */
  if(hdu->state != HDU_HEADER) {
    fprintf(stderr, "w_extkeys: Header not open for writing.\n");
    return 1;
  };
/*
 * Write extension details only where given and relevant.
 */
  if(hdu->type!=F_PRIMARY && hdu->extname) {
/*
 * Write the extension name.
 */
    if(wstrkey(fits, hdu, "EXTNAME", 0, hdu->extname, "Extension name"))
      return 1;
/*
 * If a version number for this extension has not been provided,
 * determine the appropriate number from other extensions of the same
 * type that already reside in the FITS file.
 */
    if(!hdu->extver) {
      Hdu *tmphdu = find_hdu(fits, F_ANY, hdu->extname, 0, NULL);
      hdu->extver = tmphdu==NULL ? 1 : tmphdu->extver+1;
    };
/*
 * Write the extension version and hierachical level numbers.
 */
    if(wintkey(fits, hdu, "EXTVER", 0,hdu->extver,"Extension version number") ||
       wintkey(fits, hdu, "EXTLEVEL", 0, hdu->extlevel, "Hierarchical level"))
      return 1;
  };
  return 0;
}
