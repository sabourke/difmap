#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

static void hartley(float array[], float work[], int num_el);

/*.......................................................................
  Evaluate the Fast Hartley transform of a two dimensional array,
  overwriting the data array with the transform. The reader is referred
  to the comments at the start of the 1-D fast_hartley_transform()
  function for the advantages and organisation of the FHT. The data
  array, data_array[] should actually be a 1-D array arranged as
  a contiguous 2-D array in FORTRAN sort order. The number of elements
  along each axis are sent as xnum and ynum. Both must be a power
  of two ie: xnum=2^I and ynum=2^J where I and J are integers.
*/
int two_dim_FHT(float data_array[], int xnum, int ynum, int forward)
{
        int i,j,y_element, inc_ad, inc_bc, half_xnum, half_ynum;
	float *work_array, *temp_array;
	float *ta,*tb,*tc,*td, temp;
/*
  Signal an error if the number of elements on either axis of the
  array to be transformed is not a power of two.
*/
	if(!is_pow_of_two(xnum) || !is_pow_of_two(ynum) ) {
	  lprintf(stderr, "Illegal array size (%d,%d) - not a powers of two - sent to\nthe Fast Hartley Transform function.\n", xnum, ynum);
	  return -1;
	};
/*
  Create a work array for the 1-D algorithm - equal to the length
  of the longest axis of the data array.
*/
	if((work_array = (float *) calloc(((xnum>ynum)?xnum:ynum)+1, sizeof(float))) == NULL) {
	  lprintf(stderr, "Insufficient memory available for work array of Hartley transform\n");
	  return -1;
	};
/*
  Since the second dimension of the data array is formed of
  non-contiguous elements, create another work array into which
  each y-array will be copied.
*/
	if((temp_array= (float *) calloc(ynum+1, sizeof(float))) == NULL) {
	  lprintf(stderr, "Insufficient memory available for work array of Hartley transform\n");
	  free(work_array);
	  return -1;
	};
/*
  Transform the 2nd dimension of the 2D data array first - if there
  is one.
*/
	if(ynum > 1) {
	  for(i=0; i<xnum; i++) {
/*
  Make a copy of the next y-array to be transformed.
*/
	    for(j=0,y_element=i; j<ynum; j++, y_element += xnum)
	      temp_array[j] = data_array[y_element];
/*
  Transform it.
*/
	    fast_hartley_transform(temp_array, work_array, ynum, forward);
/*
  Copy it back into the original array.
*/
	    for(j=0,y_element=i; j<ynum; j++, y_element += xnum)
	      data_array[y_element] = temp_array[j];
	  };
	};
/*
  The temporary array for extracting the y-arrays is no longer
  required.
*/
	free(temp_array);
/*
  We now have a 2-D array with transforms of each 1D-array along the
  y-axis. Now transform along the x-axis. NB. The Hartley transform is
  not really separable in x and y, however the result of transforming
  the x-axis transforms, along the y-axis, is a 2-D array that can
  easily be converted to the actual Hartley transform - see below.
*/
	for(i=0,j=0; i<ynum; i++,j += xnum)
	  fast_hartley_transform(&data_array[j], work_array, xnum, forward);
/*
  Transform completed - zap the work array.
*/
	free(work_array);
/*
  Convert the resulting (MxN) 2-D (T[u,v]) array to the 2-D Hartley transform
  as proscribed in Proc. IEEE, Vol 74, No. 9, September 1986, pp. 1282-1283.
  This gives:
   2.H(u,v) = T(u,v) + T(M-u,v) + T(u,N-v) - T(M-u,N-v)
  In order not to have to generate an intermediate (MxN) work array - the
  H(u,v) for each of (u,v), (M-u), (u,N-v) and (M-u,N-v) will be calculated
  simultaneously, via the intermediate result:
   E = 0.5 * ([T(u,v)+T(M-u,N-v)] - [T(M-u,v)+T(u,N-v)])
  And:
   H(u,v) -= E  and  H(M-u,v) += E  and  H(u,N-v) += E  and  H(M-u,N-v) -= E
  E is whenever one of u or v are 0 and whenever u = M/2 or v = N/2. Because
  of this and since four quadrants are handled at each iteration, the
  calculation need only be performed for u=1 to M/2-1 and v=1 to N/2-1.
  In our case i denotes u and j denotes v.
*/
/*
  Pre-calculate the often used half number of points on each axis
  and increments for stepping ta,td and tb,tc respectively, along the
  y-axis.
*/
	half_xnum = xnum/2;
	half_ynum = ynum/2;
	inc_ad = half_xnum+1;
	inc_bc = half_xnum+xnum-1;
/*
  Start by initializing pointers to the four T() components described
  above, to u=M/2-1, v=0 (ie. skip to the end of the v=0 iteration).
*/
	ta = data_array+half_xnum-1;
	tb = ta+2;
	tc = ta+ynum*xnum;
	td = tc+2;
/*
  Step along the y-axis.
*/
	for(j=1; j<half_ynum; j++) {
/*
  Increment the 4 pointers, to the next positions in the array along the
  y-axis.
*/
	  ta += inc_ad; tb += inc_bc; tc -= inc_bc; td -= inc_ad;
/*
  Deal with the x-axis points at the current y-axis position.
*/
	  for(i=1; i<half_xnum; i++) {
/*
  Evaluate E at the current coordinates and then increment the pointers
  to the next position on the x-axis.
*/
	    temp = ((*(++ta) + *(--td))  -  (*(--tb) + *(++tc)))/2.0;
/*
  Convert the current 4 points to their Hartley transform values.
*/
	    *ta -= temp;  *tb += temp;  *tc += temp;  *td -= temp;
	  };
	};
	return 0;
}

/*.......................................................................
  Perform the forward or reverse Fast-Hartley-Transform (FHT). The
  difference between the two is simply an extra factor of 1/num_el in
  the reverse transform (Much simpler than in the FFT). The FHT is
  superior to the FFT for a number of reasons - the reader is referred
  to the SECOND edition of "The Fourier Transform and its Applications",
  by R. N. Bracewell, for further details or (The same article) in
  Proc. IEEE, vol 72, no. 8, August 1984, pp. 1010-1018.
    The main advantage is that the Hartley transform is real, and thus
  requires no complex number manipulation and the turgid packing of
  their values into a single array, as required in most implementations
  of the FFT. In addition the real and imaginary values can still be
  recovered as the even and odd components of the tranform, so nothing
  is lost and a complex array can be transformed by first producing
  the respective even-odd real array. For applications such as
  convolution, it is vastly superior, since again no complex
  multiplications are required.
    The (1-Dimensional) data array is sent in data_array[] together
  with a work array of the same length (work_array[]). The number
  of elements in this array are given via num_el (num_el=2^N).
  If a forward transform is required do_forward_transform should
  be given a value > 0.
     The resulting transform is returned in data_array[] and
  is organised as follows: If the original data array had time
  increments of, dt, per channel: elements i=0,1,2..num_el/2 correspond
  to frequencies i/(dt*num_el), and subsequent elements correspond to
  frequencies (i-num_el)/(dt*num_el). For example if num_el=8 and dt=1,
  the frequencies of the channels would be:
    Element i:   0     1      2      3      4      5      6      7
    Frequency:   0    1/8    2/8    3/8    4/8   -3/8   -2/8   -1/8
  The real value corresponding to frequency f is the even part of the
  transform, H ie: [H(f)+H(-f)]/2, the imaginary value is minus the
  odd part: -[H(f)-H(-f)]/2 and its modulus is: [H(f)^2+H(-f)^2]/2.
    In terms of array elements these correspond to:
  real[i]    = (data_array[i]+data_array[i-num_el])/2
  imag[i]    = (data_array[i-num_el]-data_array[i])/2
  modulus[i] = (data_array[i]*data_array[i]+data_array[i-num_el]*data_array[i-num_el])/2
  Except for i=0 which is wholly real: real[0]=data_array[0], imag[0]=0.
*/
int fast_hartley_transform(float data_array[], float work_array[], int num_el, int do_forward_transform)
{
	int i;
/*
  If the data array has only one element then its transform is
  the same as the data array so no operations are required.
*/
	if(num_el == 1) return 0;
/*
  Signal an error if the number of elements in the array to be
  transformed is not a power of two.
*/
	if(!is_pow_of_two(num_el)) {
	  lprintf(stderr, "Illegal array length (%d) - not a power of two - sent to\nthe Fast Hartley Transform function.\n", num_el);
	  return -1;
	};
/*
  Perform the transform.
*/
	hartley(data_array, work_array, num_el);
/*
  Divide the transform by the number of elements in the array.
  This is only required on the forward transform.
*/
	if(do_forward_transform)
	  for(i=0;i<num_el;i++)
	    data_array[i] = data_array[i]/num_el;
	return 0;
}

/*.......................................................................
  Form the Hartley transform of a data array. The resulting transorm
  is actually double the value required and should be divided by two
  by the calling routine. The procedure is as follows:
    Rearrange the data array into two separate arrays in which the first
  array is formed of the odd elements and the second of the even
  elements. Then evaluate the Hartley transform of these two sub-arrays
  and combine their transforms to form the Hartley transform of the
  full array. The process is recursive - the recursion is terminated
  when the data array is made up of just two elements. Hence the
  original calling routine should divide the resulting transform by
  2 to the power of the number levels of recursion used - equal to the
  number of elements in the original array.
*/
static void hartley(float array[], float work[], int num_el)
{
        static int i,j,half_num;
	static double two_pi_div_n;
/*
  When only two elements are left to be transformed - finish.
*/
	if(num_el > 4) {
/*
  Sort even and odd elements into the lower and upper halves of the
  work array.
*/
	  half_num = num_el/2;
	  for(i=0,j=half_num; i < half_num; i++,j++) {
	    work[i] = array[i+i];
	    work[j] = array[i+i+1];
	  };
/*
  Now transform these two new sub-arrays.
*/
	  hartley(&work[half_num], &array[half_num], num_el/2);
	  hartley(work, array, num_el/2);
/*
  Having transformed the two new arrays, we need to combine them into
  one array, thus forming the transform of the data array that was passed
  to this routine. If we denote the upper and lower sub-arrays as A(i) and
  B(i) respectively, the combination follows as:
       H(i) = A(i) + B(i)*cos(2.PI.j/num_el) + C(i)*sin(2.PI.j/num_el)
    With C(i) = B(i)*(num_el-i) (except that B(0) when i=0).
  and where j=0,1,2,...num_el-1 and i increments by one as j does, but
  restarts from zero when it reaches num_el/2.
*/
	  half_num = num_el/2;
	  two_pi_div_n = 6.2831853072/num_el;
	  for(i=0,j=0; j<num_el; i++,j++) {
	    if(j == half_num) i=0;
	    array[j] = work[i] + work[i+half_num]*cos(two_pi_div_n*j) +
	     work[ ((i==0) ? half_num : (num_el-i))]*sin(two_pi_div_n*j);
	  };
	}
/*
  When the sub-array sent to this function has four elements the
  above sorting of even/odd numbers and the combination of these
  sub-arrays is an easy relation in which the cos's and sin's
  are 1,-1 or 0.
*/
	else if(num_el == 4) {
/*
  First copy the four elements of the data array into a work array.
*/
	  work[0]=array[0];work[1]=array[1];work[2]=array[2];work[3]=array[3];
/*
  For the special case of four elements the Hartley transform is as follows.
*/
	  array[0] = work[0] + work[1] + work[2] + work[3];
	  array[1] = work[0] + work[1] - work[2] - work[3];
          array[2] = work[0] - work[1] + work[2] - work[3];
	  array[3] = work[0] - work[1] - work[2] + work[3];
	}
/*
  When the sub-array sent to this routine is just a pair of elements it is
  easy to calculate the transorm of the pair. If the elements are (A,B)
  the resulting two element transform is (A+B,A-B). NB. In the following,
  array[0] in the second expression is A+B not A, hence the extra -B term.
  NB. This point will only be reached if the original transform array
  held just two elements since it is by-passed by the above if a four-element
  sub-array is encounterred.
*/
	else {
	  array[0] = array[0]+array[1];
	  array[1] = array[0]-array[1]-array[1];
	};
	return;
}

