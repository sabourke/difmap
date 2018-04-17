/*
 * Apply the shift theorem to switch the origin of a complex 2D
 * fourier transform between the center of the image and the corners.
 */
void fft_shift(float *image, int adim, int bdim);

/*
 * The same as fft_shift(), but designed for the real part
 * of a conjugate-symmetric.
 */
void cnj_shift(float *image, int adim, int bdim);

/*
 * Perform a 2D real or complex fft.
 */
void newfft(float *image, int adim, int bdim, int isign, int isreal,
	    int rescale);
