#ifndef uvpage_h
#define uvpage_h

/* Declare a type used to store one model visibility */

typedef struct {
  float amp;      /* Amplitude of the visibility */
  float phs;      /* Phase of the visibility */
} Mvis;

/* Declare a type used to store the state of a UV model scratch file */

typedef struct {
  Recio *rio;    /* Descriptor of binary record I/O descriptor */
  int ntime;     /* The number of integrations in each model */
  int nbase;     /* The number baselines covered by each model */
  int nif;       /* The number of IFs for which models are stored */
  int ioerr;     /* True after an I/O error */
  Mvis *mvis;    /* Array of nbase model visibilities */
} UVpage;

/* Open a uvmodel.scr paging file */

UVpage *new_UVpage(int ntime, int nbase, int nif);

/* Close a uvmodel.scr paging file */

UVpage *del_UVpage(UVpage *uvp);

/* Read a baseline's worth of visibilities from a given IF */

int uvp_read(UVpage *uvp, int ut, int cif);

/* Write a baseline's worth of visibilities for a given IF */

int uvp_write(UVpage *uvp, int ut, int cif);

/* Test for uvp I/O errors */

int uvp_error(UVpage *uvp, const char *fname);

/* Clear the uvp->mvis buffer */

int uvp_clear(UVpage *uvp);

/* Flush pending I/O to the paging file */

int uvp_flush(UVpage *uvp);

#endif
