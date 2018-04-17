#ifndef pb_h
#define pb_h

/* Encapsulate an ensemble of antenna primary beams */

typedef struct AntennaBeams AntennaBeams;
AntennaBeams *new_AntennaBeams(void);
AntennaBeams *del_AntennaBeams(AntennaBeams *ab);

/*
 * Record the voltage beam of a single antenna.
 * Beware that two calls to this function with the same arguments will
 * return the same object, with its reference count incremented.
 */
typedef struct VoltageBeam VoltageBeam;
VoltageBeam *new_VoltageBeam(AntennaBeams *ab, float *samples, int nsample,
			     float binwidth, float freq, unsigned nref);
VoltageBeam *dup_VoltageBeam(VoltageBeam *vb);
VoltageBeam *del_VoltageBeam(VoltageBeam *vb);

/*
 * Interpolate the given volate beam (or return 1.0 if pb is NULL).
 */
float voltage_beam(VoltageBeam *pb, float freq, float radius);

/*
 * Return a count of the current number of references to voltage beams.
 */
int count_antenna_beams(AntennaBeams *ab);

#endif
