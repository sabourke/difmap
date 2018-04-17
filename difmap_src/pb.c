#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "logio.h"
#include "obs.h"
#include "freelist.h"
#include "cksum.h"
#include "telspec.h"
#include "model.h"
#include "pb.h"

/*
 * The voltage beam of the antenna is the voltage response of a
 * single antenna as a function of angular radius from the center
 * of the beam. The response is assumed to be circularly symmetric.
 */
struct VoltageBeam {
  VoltageBeam *next; /* The next voltage beam in the list */
  AntennaBeams *ab;  /* The parent antenna beams container */
  unsigned nref;     /* The reference count of the object */
  float *samples;    /* The array of voltage beam samples */
  int nsample;       /* The number of samples in samples[] */
  float binwidth;    /* The width of each sample on the sky (radians) */
  float freq;        /* The reference frequency of the beam (Hz) */
  unsigned long sum; /* A checksum of the contents of the samples[] array. */
};

/*
 * Set the number of VoltageBeam objects that are added to the freelist
 * whenever it becomes exhausted.
 */
#define VB_BLOCKING 20

/*
 * Encapsulate an ensemble of beams. To save memory, aliases of beam
 * objects with the same attributes are returned instead of making
 * new copies.
 */
struct AntennaBeams {
  FreeList *vbmem;   /* A freelist of VoltageBeam structures */
  FreeList *pbmem;   /* A freelist of PrimaryBeam structures */
  VoltageBeam *vbs;  /* The head of the list of voltage beam objects */
  CheckSum *cksum;   /* An object used to compute checksums */
  int total_nref;    /* The total number of references to voltage beams */
};

static int cant_set_antenna_beam(Observation *ob, VoltageBeam *vb);

/*.......................................................................
 * Create a new AntennaBeams object.
 *
 * Output:
 *  return  AntennaBeams *  The new object, or NULL on error.
 */
AntennaBeams *new_AntennaBeams(void)
{
  AntennaBeams *ab;  /* The object to be returned */
/*
 * Allocate the container.
 */
  ab = malloc(sizeof(AntennaBeams));
  if(!ab) {
    lprintf(stderr, "new_AntennaBeams: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_AntennaBeams().
 */
  ab->vbmem = NULL;
  ab->pbmem = NULL;
  ab->vbs = NULL;
  ab->cksum = NULL;
  ab->total_nref = 0;
/*
 * Allocate a freelist of VoltageBeam objects.
 */
  ab->vbmem = new_FreeList("new_AntennaBeams", sizeof(VoltageBeam),
			   VB_BLOCKING);
  if(!ab->vbmem)
    return del_AntennaBeams(ab);
/*
 * Allocate the resources needed to compute checksums.
 */
  ab->cksum = new_CheckSum();
  if(!ab->cksum)
    return del_AntennaBeams(ab);
  return ab;
}

/*.......................................................................
 * Delete an AntennaBeams object.
 *
 * Input:
 *  ab     AntennaBeams *  The object to be deleted.
 * Output:
 *  return AntennaBeams *  The deleted object (always NULL).
 */
AntennaBeams *del_AntennaBeams(AntennaBeams *ab)
{
  if(ab) {
    ab->vbmem = del_FreeList("del_AntennaBeams", ab->vbmem, 1);
    ab->pbmem = del_FreeList("del_AntennaBeams", ab->pbmem, 1);
    ab->vbs = NULL;          /* Deleted by deleting ab->vbmem */
    ab->cksum = del_CheckSum(ab->cksum);
    free(ab);
  };
  return NULL;
}

/*.......................................................................
 * Create and populate a new voltage beam container.
 *
 * Input:
 *  ab     AntennaBeams *  The parent ensemble of beams.
 *  samples       float *  An array of 'nsample' samples of the voltage
 *                         beam, where element i refers to the voltage
 *                         beam at i*binwidth. A copy is made of this
 *                         array. Note that the voltage beam beyond
 *                         samples*binwidth is assumed to be zero.
 *  nsample         int    The number of samples in samples[]. This
 *                         must be >= 2, to support linear interpolation.
 *  binwidth      float    The radial width covered by each sample in
 *                         samples[].
 *  freq          float    The frequency to which the numbers in samples[]
 *                         and binwidth refer (Hz). Numbers for other
 *                         frequencies will be computed by assuming that
 *                         binwidth increases linearly with increasing
 *                         frequency.
 *  nref       unsigned    Normally 1, this is the amount by which to
 *                         increment the reference count of the
 *                         object.  The caller is then guaranteed
 *                         that the returned object won't be deleted
 *                         until it has been passed at least nref
 *                         times to del_VoltageBeam(). If you
 *                         specify 0 here, be sure to call
 *                         dup_VoltageBeam() each time that you
 *                         place the returned pointer somewhere
 *                         where it will subsequently be passed to
 *                         del_VoltageBeam().
 * Output:
 *  return  VoltageBeam *  The new object, or NULL on error.
 */
VoltageBeam *new_VoltageBeam(AntennaBeams *ab, float *samples, int nsample,
			     float binwidth, float freq, unsigned nref)
{
  VoltageBeam *vb;    /* The object to be returned */
  unsigned long sum;  /* The checksum of the samples[] array */
  int i;
/*
 * Check the arguments.
 */
  if(!ab || !samples || nsample < 2 || binwidth <= 0.0) {
    lprintf(stderr, "new_VoltageBeam: Invalid argument(s).\n");
    return NULL;
  };
/*
 * Compute the checksum of the samples[] array.
 */
  sum = checksum_of_object(ab->cksum, samples, sizeof(*samples) * nsample);
/*
 * If we have already recorded an identical beam increment its reference
 * count and return it.
 */
  for(vb=ab->vbs; vb; vb=vb->next) {
    if(vb->sum == sum && vb->nsample == nsample && vb->binwidth == binwidth &&
       vb->freq == freq) {
      vb->nref += nref;
      ab->total_nref += nref;
      return vb;
    };
  };
/*
 * This is a new beam, so allocate a new container.
 */
  vb = malloc(sizeof(VoltageBeam));
  if(!vb) {
    lprintf(stderr, "new_VoltageBeam: Insufficient memory.\n");
    return NULL;
  };
/*
 * Before attempting any operation that might fail, initialize the
 * container at least up to the point at which it can safely be passed
 * to del_VoltageBeam().
 */
  vb->next = NULL;
  vb->ab = ab;
  vb->nref = nref;
  vb->samples = NULL;
  vb->nsample = nsample;
  vb->binwidth = binwidth;
  vb->freq = freq;
  vb->sum = sum;
/*
 * Allocate memory to record a copy of the array of samples.
 */
  vb->samples = (float *) malloc(sizeof(*vb->samples) * nsample);
  if(!vb->samples) {
    lprintf(stderr, "new_VoltageBeam: Insufficient memory.\n");
    return del_VoltageBeam(vb);
  };
/*
 * Copy the array of samples.
 */
  for(i=0; i<nsample; i++)
    vb->samples[i] = samples[i];
/*
 * Prepend the new object to the list of voltage beams.
 */
  vb->next = ab->vbs;
  ab->vbs = vb;
/*
 * Update the record of the number of references in use.
 */
  ab->total_nref += nref;
  return vb;
}

/*.......................................................................
 * Delete a VoltageBeam object.
 *
 * Input:
 *  vb     VoltageBeam *  The object to be deleted.
 * Output:
 *  return VoltageBeam *  The deleted object (always NULL).
 */
VoltageBeam *del_VoltageBeam(VoltageBeam *vb)
{
  AntennaBeams *ab;  /* The parent container object of the beam */
/*
 * Already deleted?
 */
  if(!vb)
    return NULL;
/*
 * Get the parent container.
 */
  ab = vb->ab;
/*
 * Decrement the reference count of the beam.
 */
  if(vb->nref != 0) {
    vb->nref--;
    ab->total_nref--;
  };
/*
 * Delete the object, if not already deleted, once its reference count
 * drops to zero.
 */
  if(vb->nref == 0) {
/*
 * Locate the object in the parent list.
 */
    VoltageBeam *node, *prev;
    for(prev=NULL, node=ab->vbs; node && node!=vb; prev=node, node=node->next)
      ;
    if(!node) {
      lprintf(stderr, "del_VoltageBeam: Object missing from parent list.\n");
      return NULL;
    };
/*
 * Remove the node from the list.
 */
    if(prev)
      prev->next = node->next;
    else
      ab->vbs = node->next;
/*
 * Delete the contents of the node.
 */
    if(vb->samples)
      free(vb->samples);
/*
 * Return the node to the freelist.
 */
    vb = del_FreeListNode("del_VoltageBeam", ab->vbmem, vb);
  };
  return NULL;
}

/*.......................................................................
 * Return a read-only duplicate of a VoltageBeam object. In reality the
 * same object is returned, but with its reference count incremented by
 * one.
 *
 * Input:
 *  vb     VoltageBeam *  The object to be duplicated.
 * Output:
 *  return VoltageBeam *  The duplicated object, or NULL if ab is NULL.
 */
VoltageBeam *dup_VoltageBeam(VoltageBeam *vb)
{
  if(vb) {
    vb->nref++;
    vb->ab->total_nref++;
  };
  return vb;
}

/*.......................................................................
 * Interpolate the given voltage beam (or return 1.0 if vb is NULL).
 *
 * Input:
 *  vb     VoltageBeam *   The voltage beam container, or NULL if there
 *                         is no voltage beam.
 *  radius       float     The radius at which to interpolate the voltage
 *                         beam (radians).
 *  freq         float     The frequency at which the beam is
 *                         required (Hz).
 * Output:
 *  return       float     The requested value of the voltage beam.
 */
float voltage_beam(VoltageBeam *vb, float radius, float freq)
{
  float fbin;    /* The floating point number of bins within 'radius' */
  int ia,ib;     /* The indexes of vb->samples[] that bracket the */
                 /*  requested radius. */
/*
 * If there is no voltage beam, assume that it has infinite extent.
 */
  if(!vb)
    return 1.0;
/*
 * Work out the two sample array indexes that bracket the desired radius.
 */
  fbin = radius / vb->binwidth * (freq / vb->freq);
  ia = floor(fbin);
  ib = ceil(fbin);
/*
 * Accomodate radii beyond the limits of the sampled beam.
 */
  if(ia < 0)
    return vb->samples[0];
  else if(ib >= vb->nsample)
    return 0.0;
  else if(ia==ib)
    return vb->samples[ia];
  return vb->samples[ia] + (fbin - ia)/(ib-ia) *
    (vb->samples[ib] - vb->samples[ia]);
}

/*.......................................................................
 * Change the voltage beam of one or more antennas.
 *
 * Input:
 *  ob    Observation *   The parent observation of the antennas.
 *  spec         char *   A string of one or more antenna specifications
 *                        separated by spaces. These specify which antennas
 *                        the beam refers to.
 *  samples      float *  An array of 'nsample' samples of the voltage
 *                        beam, where element i refers to the voltage
 *                        beam at i*binwidth. A copy is made of this
 *                        array. Note that the voltage beam beyond
 *                        samples*binwidth is assumed to be zero.
 *  nsample        int    The number of samples in samples[]. If non-zero,
 *                        this must be >= 2, to support linear interpolation.
 *                        To remove an existing beam, send 0, and note that
 *                        in this case the samples, binwidth and freq
 *                        arguments are completely ignored.
 *  binwidth     float    The radial width covered by each sample in
 *                        samples[].
 *  freq         float    The frequency to which the numbers in samples[]
 *                        and binwidth refer (Hz). Numbers for other
 *                        frequencies will be computed by assuming that
 *                        binwidth increases linearly with increasing
 *                        frequency.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int set_antenna_beam(Observation *ob, char *spec, float *samples, int nsample,
		     float binwidth, float freq)
{
  char *endp;      /* A pointer to the next unprocessed character in *spec */
  VoltageBeam *vb; /* The new antenna beam */
  int pass;        /* The number of the pass through the antenna */
                   /*  specifications. */
/*
 * Check the arguments.
 */
  if(!ob || !spec || (nsample > 0 && !samples)) {
    lprintf(stderr, "set_antenna_beam: NULL argument(s).\n");
    return 1;
  };
/*
 * Has a new beam description been provided?
 */
  if(nsample > 0) {
/*
 * Check the beam characteristics.
 */
    if(freq <= 0.0) {
      lprintf(stderr, "Invalid antenna beam frequency (%g).\n", freq);
      return 1;
    };
    if(binwidth <= 0.0) {
      lprintf(stderr, "Invalid sampling binwidth (%g).\n", binwidth);
      return 1;
    };
/*
 * Allocate a new beam description object.
 */
    vb = new_VoltageBeam(ob->ab, samples, nsample, binwidth, freq, 1);
    if(!vb)
      return 1;
  } else {
    vb = NULL;
  };
/*
 * Parse the list of telescope specifications twice, once to check it
 * for validity before doing anything that we can't back out of, then
 * a second time to associate the new beam with each of the selected
 * antennas.
 */
  for(pass=1; pass<=2; pass++) {
    char *s = spec;
/*
 * Parse the list of telescope specifications until the end of the
 * string is reached.
 */
    while(*s != '\0') {
/*
 * Read the latest telescope specification.
 */
      Telspec *ts = read_Telspec(ob, s, &endp, 0);
      if(!ts)
	return cant_set_antenna_beam(ob, vb);
/*
 * On the second pass iterate through the selected antennas, assigning
 * the new beam to each of them in turn.
 */
      if(pass==2) {
	Findop oper = FIND_FIRST;
	for(oper=FIND_FIRST; next_tel(ob, oper, 1, ts->nfix, 0, 0, ts)==0;
	    oper=FIND_NEXT) {
	  Station *tel = ob->sub[ts->isub].tel + ts->ta;
	  tel->vb = del_VoltageBeam(tel->vb);
	  tel->vb = dup_VoltageBeam(vb);
	};
      };
/*
 * The next character should either be a space, indicating that another
 * specification may follow, or the end of the string.
 */
      if(*endp != ' ' && *endp != '\t' && *endp != '\0') {
	lprintf(stderr, "Garbled telescope specification (%s).\n", endp);
	return cant_set_antenna_beam(ob, vb);
      };
/*
 * Skip spaces up to the next specification.
 */
      for(s=endp; *s==' ' || *s=='\t'; s++)
	;
    };
  };
/*
 * Delete our temporary reference to vb. Note that this won't delete
 * any references that were assigned to the selected antennas.
 */
  vb = del_VoltageBeam(vb);
  return 0;
}

/*.......................................................................
 * This is a private error cleanup return function of set_antenna_beam().
 * It releases temporary resources used by its caller, then returns
 * the error return code of its caller.
 */
static int cant_set_antenna_beam(Observation *ob, VoltageBeam *vb)
{
/*
 * Delete the temporary reference that the caller had to vb. Note that
 * this won't delete references that have already been assigned to
 * any of the selected antennas before this function was called.
 */
  vb = del_VoltageBeam(vb);
  return 1;
}

/*.......................................................................
 * Return a count of the current number of references to voltage beams.
 *
 * Input:
 *  ab    AntennaBeams *  The ensemble of beams to describe.
 * Output:
 *  return         int    The number of references to antenna voltage
 *                        beams.
 */
int count_antenna_beams(AntennaBeams *ab)
{
  return ab ? ab->total_nref : 0;
}

/*.......................................................................
 * Calculate the primary beam scale factor at a given radius from the
 * pointing center, for a given baseline and frequency.
 *
 * Input:
 *  sub          Subarray *  The parent subarray of the baseline.
 *  base              int    The baseline who's primary beam is to be
 *                           sampled.
 *  freq           double    The frequency at which to evaluate the
 *                           primary beam.
 *  radius          float    The radial distance of the target location,
 *                           (radians).
 * Output:
 *  return          float    The primary beam response for the given
 *                           position.
 */
float pb_bl_factor(Subarray *sub, int base, double freq, float radius)
{
  VoltageBeam *v1; /* The voltage beam container of the first antenna of */
                   /*  the baseline. */
  VoltageBeam *v2; /* The voltage beam container of the second antenna of */
                   /*  the baseline. */
  Baseline *b;     /* The baseline object corresponding to 'base' */
/*
 * Check the arguments.
 */
  if(!sub) {
    lprintf(stderr, "pb_bl_factor: NULL argument(s).\n");
    return 0.0;
  };
  if(base < 0 || base >= sub->nbase) {
    lprintf(stderr, "pb_bl_factor: Baseline index out of range.\n");
    return 0.0;
  };
/*
 * Get the baseline description.
 */
  b = sub->base + base;
/*
 * Get the voltage beams of the two antennas.
 */
  v1 = sub->tel[b->tel_a].vb;
  v2 = sub->tel[b->tel_b].vb;
/*
 * If voltage beams haven't been provided for either antenna, return
 * unit response.
 */
  if(!v1 && !v2)
    return 1.0;
/*
 * Compute the primary beam response at the computed radial offset.
 */
  return voltage_beam(v1, radius, freq) * voltage_beam(v2, radius, freq);
}

/*.......................................................................
 * Calculate the primary beam scale factor at a given radial distance
 * from the pointing center, averaged over all baselines, IFs and
 * subarrays.
 *
 * Input:
 *  ob        Observation *  The parent observation.
 *  radius          float    The radial distance of the target location,
 *                           (radians).
 * Input/Output:
 *  factor          float *  The primary beam response for the given
 *                           position will be assigned to *factor.
 * Output:
 *  return            int    0 - OK.
 *                           1 - Error.
 */
int pb_scale_factor(Observation *ob, float radius, float *factor)
{
  int isub;          /* The index of a subarray */
  double mean=0.0;   /* The weighted mean primary beam factor */
  double wtsum=0.0;  /* The sum of baseline weights */
/*
 * Check the arguments.
 */
  if(!ob_ready(ob, OB_SELECT, "pb_scale_factor"))
    return 1;
  if(!factor) {
    lprintf(stderr, "pb_scale_factor: NULL argument(s).\n");
    return 1;
  };
/*
 * Make sure that the per-baseline sums of visibility weights are up
 * to date.
 */
  if(update_baseline_weights(ob, -1))
    return 1;
/*
 * Loop over the subarrays.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Subarray *sub = ob->sub + isub;
    int base;
/*
 * Loop over the baselines of the subarray.
 */
    for(base=0; base<sub->nbase; base++) {
      Baseline *b = sub->base + base;
      int cif;
/*
 * Sum over the all of the IFs of the observation.
 */
      for(cif=0; cif<ob->nif; cif++) {
/*
 * Get the sum of weights of visibilities of baseline 'base' and IF 'cif'.
 */
	float wt = b->bwt[cif].wtsum;
/*
 * Are there any contributing visibilities?
 */
	if(wt > 0.0) {
/*
 * Compute the primary beam contribution of the current baseline at
 * the current IF frequency.
 */
	  float pb = pb_bl_factor(sub, base, ob->ifs[cif].freq, radius);
/*
 * Add the contribution of visibilities of the current baseline and
 * IF to the overall primary beam factor. The following statement
 * implements a weighted running mean.
 */
	  mean += (pb - mean) * wt / (wtsum += wt);
	};
      };
    };
  };
/*
 * Return the primary beam scale factor.
 */
  *factor = mean;
  return 0;
}

/*.......................................................................
 * Correct the flux of a model delta component for the effects of the
 * primary beam. On input the model component should model the flux in
 * a map that hasn't been corrected for the primary beam. On output the
 * flux will be the estimated true flux. The amplitude of the Fourier
 * transform of a delta function is constant for all visibilities. Thus
 * the the apparent map flux is related to the true sky flux by the
 * following sum over visibilties.
 *
 *                                 Sum(wt(sub,base,if) * pb(sub,base,if,x,y))
 * map_flux(x,y) = sky_flux(x,y) * ------------------------------------------
 *                                            Sum(wt(sub,base,if))
 *
 * Note that this is only correct for delta components, but other components
 * aren't rejected, since in some cases it may be a good approximation.
 *
 * Input:
 *  ob  Observation *  The parent observation.
 * Input/Output:
 *  cmp      Modcmp *  The delta component to be corrected.
 * Output:
 *  return      int    0 - OK.
 *                     1 - Error.
 */
int pb_correct_delta_cmp(Observation *ob, Modcmp *cmp)
{
  float factor;  /* The primary beam suppression factor */
/*
 * Check the arguments.
 */
  if(!ob || !cmp) {
    lprintf(stderr, "pb_correct_delta_cmp: NULL argument(s).\n");
    return 1;
  };
/*
 * Compute the combined primary beam factor for the component,
 * averaged over all baselines and IFs.
 */
  if(pb_scale_factor(ob, calc_pointing_offset(ob, cmp->x, cmp->y), &factor))
    return 1;
/*
 * Remove this factor from the flux of the component. Note that
 * pb_scale_factor() returns an error if factor would be zero,
 * so we don't have to check for zero-divide here.
 */
  if(factor==0.0)
    cmp->flux = 0.0;
  else
    cmp->flux /= factor;
  return 0;
}

/*.......................................................................
 * Change the primary beam of all baselines.
 *
 * Input:
 *  ob    Observation *   The parent observation of the antennas.
 *  samples      float *  An array of 'nsample' samples of the voltage
 *                        beam, where element i refers to the voltage
 *                        beam at i*binwidth. A copy is made of this
 *                        array. Note that the voltage beam beyond
 *                        samples*binwidth is assumed to be zero.
 *  nsample        int    The number of samples in samples[]. If non-zero,
 *                        this must be >= 2, to support linear interpolation.
 *                        To remove an existing beam, send 0, and note that
 *                        in this case the samples, binwidth and freq
 *                        arguments are completely ignored.
 *  binwidth     float    The radial width covered by each sample in
 *                        samples[].
 *  freq         float    The frequency to which the numbers in samples[]
 *                        and binwidth refer (Hz). Numbers for other
 *                        frequencies will be computed by assuming that
 *                        binwidth increases linearly with increasing
 *                        frequency.
 * Output:
 *  return        int     0 - OK.
 *                        1 - Error.
 */
int set_primary_beam(Observation *ob, float *samples, int nsample,
		     float binwidth, float freq)
{
  VoltageBeam *vb; /* The new antenna beam */
  int isub, itel;  /* The indexes of a subarray and telescope */
  float *vsamples; /* A dynamically allocated copy of the samples[] array */
                   /*  square rooted. */
  int i;
/*
 * Check the arguments.
 */
  if(!ob || (nsample > 0 && !samples)) {
    lprintf(stderr, "set_primary_beam: NULL argument(s).\n");
    return 1;
  };
/*
 * Has a new beam description been provided?
 */
  if(nsample > 0) {
/*
 * Check the beam characteristics.
 */
    if(freq <= 0.0) {
      lprintf(stderr, "Invalid primary beam frequency (%g).\n", freq);
      return 1;
    };
    if(binwidth <= 0.0) {
      lprintf(stderr, "Invalid sampling binwidth (%g).\n", binwidth);
      return 1;
    };
/*
 * Allocate a temporary array with the same number of elements as the
 * samples[] array.
 */
    vsamples = (float *) malloc(sizeof(*vsamples) * nsample);
    if(!vsamples) {
      lprintf(stderr, "primary_beam: Insufficient memory.\n");
      return 1;
    };
/*
 * Set vsamples[*] = sqrt(samples[*]);
 */
    for(i=0; i<nsample; i++)
      vsamples[i] = sqrt(fabs(samples[i]));
/*
 * Allocate a new antenna beam, using the square root of the array of
 * primary beam samples.
 */
    vb = new_VoltageBeam(ob->ab, vsamples, nsample, binwidth, freq, 1);
    if(!vb) {
      free(vsamples);
      return 1;
    };
/*
 * The vsamples[] array is no longer needed.
 */
    free(vsamples);
  } else {
    vb = NULL;
  };
/*
 * Install the new beam on all antennas.
 */
  for(isub=0; isub<ob->nsub; isub++) {
    Subarray *sub = ob->sub + isub;
    Station *tel = sub->tel;
    for(itel=0; itel<sub->nstat; itel++,tel++) {
      tel->vb = del_VoltageBeam(tel->vb);
      tel->vb = dup_VoltageBeam(vb);
    };
  };
/*
 * Delete our temporary reference to vb. Note that this won't delete
 * any references that were assigned to the selected antennas.
 */
  vb = del_VoltageBeam(vb);
  return 0;
}

