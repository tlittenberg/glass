/*
 * Copyright 2023 Tyson B. Littenberg
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *  http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/**
@file glass_math.h
\brief Common math functions used throughout GLASS library
 */

#ifndef math_h
#define math_h

struct CubicSpline
{
    int N;       //!<Number of grid points to be interpolated
    int nmin;    //!<Stored lower index of last call to interpolation
    int nmax;    //!<Stored upper index of last call to interpolation
    double *x;   //!<Independent variable of function to be interpolated
    double *y;   //!<Dependent variable of function to be interpolated
    double *d2y; //!<Second derivitives of function to be interpolated
};

struct CubicSpline* alloc_cubic_spline(int N);

/**
 \brief Computes cubic spline coefficients for input data {x,y}
 
 @param[in,out] spline cubic spline structure
 @param[in] x independent variable of interpolant
 @param[in] y dependent variable of interpolant
*/
void initialize_cubic_spline(struct CubicSpline *spline, double *x, double *y);

void free_cubic_spline(struct CubicSpline *spline);

/**
\brief GLASS implementation of solving for cubic spline interpolation coefficients
 
 Uses the tridiagonal algorithm to compute the second derivatives d2y of the input data y=f(x).
 It implicitly assumes 0 2nd deriviative at endpoints.
 
 Returns interpolated value y = f(x).
 
 @param[in] spline->N number of spline points
 @param[in] spline->x vector of independent-variable grid points
 @param[in] spline->y vector of dependent-variable grid points
 @param[out] spline->d2y second derivatives of \f$ f(x)\f$

 */
void spline_coefficients(struct CubicSpline *spline);

/**
\brief GLASS implementation of cubic spline interpolation
 
 Returns interpolated value y = f(x).
 
 @param[in] spline->N number of spline points
 @param[in] spline->x vector of independent-variable grid points
 @param[in] spline->y vector of dependent-variable grid points
 @param[in] spline->d2y second derivatives of \f$ f(x)\f$
 @param[in] x value of independent-variable where interpolated value is neede
 @return interpolated value \f$ y = f(x)\f$

 */
double spline_interpolation(struct CubicSpline *spline, double x);
double spline_interpolation_deriv(struct CubicSpline *spline, double x);
double spline_interpolation_deriv2(struct CubicSpline *spline, double x);
double spline_integration(struct CubicSpline *spline, double xi, double xf);

/**
 \brief Analytic in-place inversion of noise covariance matrix
 
 Substitutes contents of noise->Cij[n] with inverse, and sets deteriminent noise->detC[n]
 */
void invert_noise_covariance_matrix(struct Noise *noise);

/**
\brief Compute chirp mass from component masses
  
 @param m1 mass of primary
 @param m2 mass of secondary
 @return \f$ \mathcal{M} \equiv (m_1 m_2)^{3/5} (m_1+m_2)^{-1/5} \f$
 */
double chirpmass(double m1, double m2);

/**
\brief Compute symmetric mass ratio of binary

 @param Mchirp chirp mass \f$\mathcal{M}\f$ ( \f$ {\rm M}_\odot \f$ )
 @param Mtotal total mass \f$ M \f$( \f$ {\rm M}_\odot \f$ )
 @return \f$ \eta = \left( \frac{\mathcal{M}}{M} \right)^{5/3} \f$
 */
double symmetric_mass_ratio(double Mchirp, double Mtotal);

/**
\brief Compute component masses of binary
 
 Uses symmetric mass ratio \f$ \eta \f$ to compute mass difference
 \f$ \delta m = \sqrt{1-4\eta} \f$ such that \f$ m_{1,2}=M(1 \pm \delta m)/2\f$

 @param[in] Mchirp chirp mass \f$\mathcal{M}\f$ ( \f$ {\rm M}_\odot \f$ )
 @param[in] Mtotal total mass \f$ M \f$( \f$ {\rm M}_\odot \f$ )
 @param[out] m1 primary mass \f$ m_1 \f$ ( \f$ {\rm M}_\odot \f$ )
 @param[out] m2 secondary mass \f$ m_2 \f$( \f$ {\rm M}_\odot \f$ )
 */
void component_masses(double Mchirp, double Mtotal, double *m1, double *m2);

/**
\brief Compute GW amplitude from intrinsic parameters
 
 @param Mc chirp mass: \f$\mathcal{M}\ [{\rm M}_\odot]\f$
 @param f0 initial GW frequency: \f$ f_0\ [{\rm Hz}]\f$
 @param D luminosity distance: \f$ D_L [{\rm pc}]\f$
 @return \f$ \mathcal{A} = 2 \frac{ \mathcal{M}^{5/3} (\pi f_0)^{2/3} }{D_L} \f$
 */
double amplitude(double Mc, double f0, double D);

/** 
\brief Low PN approximation of current time given the frequency

 @param Mchirp chirp mass \f$\mathcal{M}\f$ ( \f$ {\rm M}_\odot \f$ )
 @param Mtotal total mass \f$ M \f$( \f$ {\rm M}_\odot \f$ )
 @param tc coalesence time (s)
 @param f gravitational wave frequency (Hz)
 @return t current time (s)
*/
double post_newtonian_time(double Mchirp, double Mtotal, double tc, double f);

/** 
\brief Low PN approximation of gravitational wave frequency given time until coalesence

 @param Mchirp chirp mass \f$\mathcal{M}\f$ ( \f$ {\rm M}_\odot \f$ )
 @param tc coalesence time (s)
 @param t current time (s)
 @return f gravitational wave frequency (Hz)
*/
double post_newtonian_frequency(double Mchirp, double tc, double t);

/**
\brief Compute integer powers of x by brute force multiplication
 
 Faster than pow() for small integers
  
 @param x variable
 @param n exponent
 @return \f$x^n =  x\times x \times\ ...\times x\f$
 */
double ipow(double x, int n);

/**
\brief Tukey window time series data
   
 @param data[in,out] time series to be windowed
 @param alpha[in] size of window filter in [s]
 @param N[in] size of time series
 */
void tukey(double *data, double alpha, int N);
double tukey_scale(double alpha, int N);

void detrend(double *data, int N, int Navg);

/**
\brief Rearrange output ofRFT
 
 Real and Imaginary components from  RFT functions are not ordered the way the rest of the package expects.
  
 @param x[out] array to be filled with ordered fourier coefficients
 @param x_packed[in] input array with ill-formatted fourier coefficients
 @param N[in] size of arrays
 */
void unpack_fft_output(double *x, double *x_packed, int N);

/**
\brief Wrappers to FFT functions
 
 In-place forward and reverse mixed radix Fourier transforms
  
 @param data[in/out] array to be transformed
 @param N[in] size of arrays
 */
void glass_forward_complex_fft(double *data, int N);
void glass_inverse_complex_fft(double *data, int N);
void glass_forward_real_fft(double *data, int N);
void glass_inverse_real_fft(double *data, int N);

/**
\brief Compute power of complex amplitude in single element of data
 
 Assumes the `data` array has alternating real and imaginary terms.
 
 @param data complex amplitude array
 @param n desired sample
 @return \f$P =  d^*_n d_n \f$
 */
double power_spectrum(double *data, int n);

/**
\brief Fourier-domain noise weighted inner product
  
 @param a complex amplitude array
 @param b complex amplitude array
 @param invC inverse covariance matrix
 @param N number of frequency bins in sum
 @return \f$(a|b) =  4 \sum_n a^*_n b_n C^-1_n \f$
 */
double fourier_nwip(double *a, double *b, double *invC, int N);
double wavelet_nwip(double *a, double *b, double *invC, int *list, int N);

/**
\brief Our implementation of the recursive binary search algorithm
   
 @param array list of values to search over
 @param nmin starting index of search
 @param nmax stopping index of search
 @param x value to search for in array
 @return index in array of x
 */
int binary_search(double *array, int nmin, int nmax, double x);

/**
\brief Computes eigenvectors and eigenvalues of matrix
   
 @param[in] N size of matrix
 @param[in] matrix the input \f$N\times N\f$ matrix
 @param[out] evector 2D array of eigenvector components
 @param[out] evalue 1D array of eigenvalues, corresponding with columns of evector
 */
void matrix_eigenstuff(double **matrix, double **evectors, double *evalues, int N);
/**
\brief Matrix inversion with replacement
   
 @param[in] N size of matrix
 @param[in,out] matrix the input \f$N\times N\f$ matrix, replaced with inverse
 */
void invert_matrix(double **matrix, int N);

/**
 \brief Matrix inversion preserving original and returning L and det

 @param[in] matrix input NxN matrix
 @param[out] inverse inverse of input matrix
 @param[out] L LU decomposition of input matrix
 @param[out] det deteriment of input matrix
 @param[in] N size of input matrxi
 */
void decompose_matrix(double **matrix, double **inverse, double **L, double *det, int N);

/**
\brief Matrix multiplication
   
 @param[in] N size of matrix
 @param[in] A \f$N\times N\f$ matrix \f$A\f$
 @param[in] B \f$N\times N\f$ matrix \f$B\f$
 @param[out] AB the matrix product \f$AB\f$
 */
void matrix_multiply(double **A, double **B, double **AB, int N);

/**
\brief Wrapper to LAPACK Cholesky decomposition routine
    
 @param[in] N size of matrix
 @param[in] A \f$N\times N\f$ matrix \f$A\f$
 @param[out] L Cholesky decomposition of \f$A = LL^{-1}\f$
 */
void cholesky_decomp(double **A, double **L, int N);

/**
 \brief Wrapper to `GLASS` cubic spline interpolation routines.

 @param[in] N number of spline points
 @param[in] x vector of independent-variable spline points
 @param[in] y vector of dependent-variable spline points
 @param[in] Nint number of interpolated points
 @param[out] xint vector of interpolated independent-variable points
 @param[out] yint vector of interpolated dependent-variable points
 */
void CubicSplineGLASS(int N, double *x, double *y, int Nint, double *xint, double *yint);

/**
\brief GLASS implementation of DBSCAN clustering algorithm
 
 Density based clustering algorithm implemented here for 1D
 data with simple Euclidean distance measure
  
 @param[in] X set of data points to cluster
 @param[in] eps maximum distance between two samples to be considered neighbors
 @param[in] min minimum number of samples in a neighborhood to be considered a cluster
 @param[out] C cluster assignments mapping C[n] = M means X[n] is assigned to cluster M
 @param[out] K total number of clusters found
 @param[in] size size of X
 */
void dbscan(double *X, double eps, int min, int C[], int *K, int size);

/** 
\brief Transform periodic to linear phase

 @param[in] N size of phase array
 @param[in,out] phase input from [0,2pi] replaced with unwrapped version
*/
void unwrap_phase(int N, double *phase);


double simpson_integration_3(double f0, double f1, double f2, double h);
double simpson_integration_5(double f0, double f1, double f2, double f3, double f4, double h);

/**
\brief wrapper for qsort() specific to integer arrays
*/
void integer_sort(int *x, int N);

/**
\brief wrapper for qsort() specific to double arrays
*/
void double_sort(double *x, int N);

/**
\brief wrapper for qsort() to get sorted indicies of input array
*/
void index_sort(int *index, double *data, int N);


void list_union(int *A, int *B, int NA, int NB, int *AUB, int *NAUB);

/**
 \brief return pdf gaussian at x
 
 @param[in] x value for which you want \f$ p(x) \f$
 @param[in] mean mean of the Gaussian distribution \f$(\mu)\f$
 @param[in] sigma standard deviaiton of the Gaussian distrubtion \f$(\sigma)\f$
 @return \f$ p(x) = \frac{1}{\sqrt{2\pi}\sigma} exp^{-\frac{1}{2}\frac{(x-\mu)^2}{\sigma^2} } \f$
 */
double gaussian_pdf(double x, double mean, double sigma);


/**
 \brief return variance of data vector x
 
 @param[in] x array
 @param[in] N size of array
 @return variance of x
 */
double get_variance(double *x, int N);

/**
 \brief return mean of data vector x
 
 @param[in] x array
 @param[in] N size of array
 @return mean of x
 */
double get_mean(double *x, int N);

/**
 \brief return quantile q of sorted data vector

 @param[in] data sorted data array
 @param[in] N size of sorted data array
 @param[in] q desired quantile
 @return value of data array d at quantile q
 */
double get_quantile_from_sorted_data(double *data, int N, double q);

/**
 \brief get minimum and maximum value of data vector
 
 @param[in] data unsorted data array
 @param[in] N size of data array
 @param[out] min minimum value of data array
 @param[out] max maximum vallue of data array
 */
void get_min_max(double *data, int N, double *min, double *max);

/**
 \brief GLASS implementation of the normalized incomplete Beta function
 
 which is \f$ I_x(a,b) = B_x(a,b) / B(a,b) \f$ computed using the relation
 \f$ I_x(a,b,x) = (1/a) x^a {}_2F_1(a,1-b,a+1,x)/B(a,b) \f$
 
 @param[in] a
 @param[in] b
 @param[in] x
 @return \f$ I_x(a,b) \f$
 */
double incomplete_beta_function(double a, double b, double x);

void extract_amplitude_and_phase(int Ns, double *As, double *Dphi, double *M, double *Mf, double *phiR);

#endif /* math_h */

