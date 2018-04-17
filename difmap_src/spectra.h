#ifndef spectra_h
#define spectra_h

typedef struct {      /* The spectrum of one IF */
  Cvis *chan;         /* The spectrum channels */
  int nchan;          /* The number of entries in spectrum[] */
} IFspec;

typedef struct {
  int nbmax;          /* The dimension of baselines[] */
  int nbase;          /* The number of baselines sampled in baselines[] */
  int *baselines;     /* An array of the baselines to be sampled */
} Specsub;

typedef struct Spectrum Spectrum;
struct Spectrum {
  int uta,utb;       /* Start and end ob->rec[] integration indexes */
  Obpol obpol;       /* The polarization descriptor for the spectrum */
  float uvmin;       /* The minimum UV radius to sample */
  float uvmax;       /* The maximum UV radius to sample */
  Specsub *ssub;     /* Array of nsub baseline list arrays */
  int nsub;          /* The number of sub-arrays in ssub[] */
  IFspec *ifs;       /* 'nif' IF visibility spectra */
  int nif;           /* The number of IFs in the parent observation */
  int dovector;      /* True for a vector averaged spectrum, false for scalar */
  Spectrum *next;    /* The next spectrum in a list of spectra */
};

typedef struct {
  Observation *ob;   /* The observation of the spectra */
  Spectrum *head;    /* Head of list of spectra */
  Spectrum *tail;    /* End of list of spectra */
} Spectra;

Spectra *new_Spectra(Observation *ob);
Spectra *del_Spectra(Spectra *spectra);
int get_Spectra(Spectra *spectra);

Spectrum *add_Spectrum(Spectra *spectra, int dovector, Stokes stokes,
		       int uta, int utb, float uvmin, float uvmax,
		       Basegrp *bgrp);
void clr_Spectrum(Spectrum *spec);

int spc_set_bgrp(Spectra *spectra, Spectrum *spec, Basegrp *bgrp);
int spc_set_pol(Spectra *spectra, Spectrum *spec, Stokes pol);
int spc_set_ut(Spectra *spectra, Spectrum *spec, int uta, int utb);
int spc_set_avmode(Spectra *spectra, Spectrum *spec, int dovector);
int spc_set_uvrange(Spectra *spectra, Spectrum *spec, float uvmin, float uvmax);

Spectrum *rem_Spectrum(Spectra *spectra, Spectrum *spec);
Spectrum *del_Spectrum(Spectrum *spec);

#endif
