#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>

#include "utils.h"

#ifndef RAND_MAX
#define RAND_MAX LONG_MAX
#endif

/*.......................................................................
  This routine uses the machines random number generator (as defined in
  stdlib.h) to return a float uniform deviate between -1 and 1. In order
  that the results be free of sequential correlations, independant
  of the quality of the generator provided on any particular machine, 
  a shuffling table is established, as advised in numerical recipes.
  If the iseed argument is non-zero the native generator will be
  re-seeded and the shuffling table re-initialised. The same thing happens
  on the first call to the routine, as signalled by the is_first flag.
*/
float frand(unsigned int iseed)
{
	static int k,is_first=1;
	static float tmp, rand_num, stab[97];
	static double dtmp;
	const double maxval = RAND_MAX;
/*
  If this is the first time that the random number generator has been
  called or if a new non-zero seed has been parsed then seed the
  generator and initialise the shuffling table.
*/
        if(is_first || iseed != 0) {
	  iseed = (iseed == 0) ? 1U : iseed;
	  srand(iseed);
          is_first=0;
/*
  Get the system's random number generator warmed up.
*/
          for(k=0; k<97; k++)
            rand();
/*
  Fill the shuffling table with deviates, ranging between 0 and 1.0.
*/
          for(k=0; k<97; k++) {
	    dtmp = ((double) rand())/maxval;
	    stab[k] = (float) dtmp;
	  };
/*
  Get one extra random number.
*/
	  dtmp = ((double) rand())/maxval;
	  tmp = (float) dtmp;
        };
/*
  Use the latest deviate as a pointer into the shuffling table.
*/
	k=(int) 96.0*tmp;
	tmp=stab[k];
/*
  Make the new deviate cover the range -1 to +1.
*/
	rand_num = 2.0 * stab[k] - 1.0;
/*
  Having now used the deviate at element k of the shuffling table,
  replace it with a new one.
*/
	dtmp = ((double) rand())/maxval;
	stab[k] = (float) dtmp;
/*
  Return the new deviate.
*/
	return rand_num;
}

/*.......................................................................
  This function returns a uniformly distributed random number between
  -num and num, the first argument of the function.
*/
float uniform_rand(float num)
{
        static float rand_num;
	rand_num = num * frand(0U);
	return rand_num;
}

/*.......................................................................
  This function returns a random number, from a gaussian distribution of
  standard deviation, num.
*/
float gauss_rand(float num)
{
	static float rrad;
	static float aval, bval, rand_num;
/*
  Acquire two uniform random numbers between -1 and 1.
*/
	do {
	  aval=frand(0U);
	  bval=frand(0U);
/*
  The Box-Muller transformation to convert uniform to gaussian
  deviates requires that the two deviates be converted to a radius
  squared. The value of the radius must be less than one and not equal
  to zero.
*/
	  rrad = aval*aval + bval*bval;
	} while(rrad >= 1 || rrad == 0.0);
/*
  Apply the Box-Muller transformation to turn the random square radius
  found above into two Gaussian deviates.
*/
	rand_num = num * aval * sqrt(-2.0 * log((double) rrad)/rrad);
	return rand_num;
}

