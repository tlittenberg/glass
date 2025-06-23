/*
 * Copyright 2019 Tyson B. Littenberg & Neil J. Cornish
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
 @file glass_ucb_waveform.h
 \brief Ultra Compact Binary waveform generator.
 */

#ifndef ucb_waveform_h
#define ucb_waveform_h

#include <stdio.h>

/**
\brief Signal to noise ratio in frequency domain
  
 Calls nwip() to compute inner product of waveform with itself summed over all frequencies and data channels `I`. If 4-link data, use single TDI channel `I` = `X`. If 6-link, use `I` = {`A`, `E`}
 
 @param source structure containing source parameters and waveform \f$h\f$
 @param noise structure containing noise model \f$S_n\f$
 @return \f$\rho =  \sqrt{\sum_I (h_I|h_I)} \f$
 */
double snr(struct Source *source, struct Noise *noise);

/**
 \brief Signal to noise ratio in wavelet domain
 @see snr()
 */
double snr_wavelet(struct Source *source, struct Noise *noise);

/**
\brief Analytic approximation to SNR
  
 Not exactly what is in the paper. Expression has been calibrated against snr().
 
 @param A gravitational wave amplitude \f$\mathcal{A}\f$
 @param Sn noise power spectral density
 @param Sf correction based on sources proximity to transfer frequency
 @param sqT \f$\sqrt{T}\f$
 @return \f$ \rho \approx  \frac{1}{2} \mathcal{A} \sqrt{T} S_f / S_n \f$
 */
double analytic_snr(double A, double Sn, double Sf, double sqT);

/**
\brief Compute prior on SNR
  
 Shouldn't this be in ucb_prior.h?
 Peak of distribution `SNRPEAK` \f$\rho_*\f$ defined in Constants.h
 
 @param SNR signal to noise ratio \f$\rho\f$ of source
 @return \f$ p(\rho) = \frac{3\rho}{4 \rho_*^2 \left(1 + \frac{\rho}{4\rho_*}
 \right)^{5}} \f$
 */
double snr_prior(double SNR);

/**
\brief Compute match between waveforms
   
 @param a waveform \f$h_a\f$
 @param b waveform \f$h_b\f$
 @return \f$  M = \frac{(h_a|h_b)}{\sqrt{(h_a|h_a)+(h_b|h_b)}} \f$
 */
double waveform_match(struct Source *a, struct Source *b, struct Noise *noise);

/**
\brief Compute match between waveforms in wavelet domain
 @see waveform_match()
 */
double waveform_match_wavelet(struct Source *a, struct Source *b, struct Noise *noise);

/**
\brief Compute distance between waveforms
   
 @param a waveform \f$h_a\f$
 @param b waveform \f$h_b\f$
 @return \f$  D = (h_a-h_b | h_a-h_b) \f$
 */
double waveform_distance(struct Source *a, struct Source *b, struct Noise *noise);


/**
 \brief Compute GR-driven frequency derivative from intrinsic parameters
 
 @param Mc chirp mass: \f$\mathcal{M}\ [{\rm M}_\odot]\f$
 @param f0 initial GW frequency: \f$ f_0\ [{\rm Hz}]\f$
 @return \f$ \dot f = \frac{96}{5} \pi^{8/3} \mathcal{M}^{5/3} f_0^{11/3} \ [{\rm Hz}\ {\rm s}^{-1}] \f$
 */
double ucb_fdot(double Mc, double f0);

/**
 \brief Compute chirp mass from frequency parameters
 
 @param f0 initial GW frequency: \f$ f_0\ [{\rm Hz}]\f$
 @param dfdt frequency derivative: \f$ \dot f \ {\rm s}^{-1}]\f$
 @return \f$\mathcal{M} = \left(\frac{\dot f} {\frac{96}{5} \pi^{8/3} f_0^{11/3}}\right)^{3/5} \ [{\rm M}_\odot]\f$
 */
double ucb_chirpmass(double f0, double dfdt);

/**
 \brief Compute luminosity distance assuming GR-driven orbital evolution
 
 @param f0 initial GW frequency: \f$ f_0\ [{\rm Hz}]\f$
 @param dfdt frequency derivative: \f$ \dot f \ [{\rm s}^{-1}]\f$
 @param A gravitational wave amplitude: \f$ \mathcal{A} \f$
 @return \f$D_L = \frac{5}{48} \frac{\dot f}{\pi^2 f_0^3 \mathcal{A}} \ [{\rm pc}]\f$
 */
double ucb_distance(double f0, double dfdt, double A);

/**
 \brief Compute UCB phase at given time
 
 @param t time (s)
 @param params[] source parameters
 @param T total observation time (s)
 @return \f$ \Phi = \varphi_0 + 2\pi\left( f_0 t + \frac{1}{2}\dot{f}t^2 + \frac{1}{6} \ddot{f}t^3 \right) \f$
 */
double ucb_phase(double t, double *params, double T);

/**
 \brief Compute UCB amplitude at given time
 
 @param t time (s)
 @param params[] source parameters
 @param T total observation time (s)
 @return \f$ A = A_0\left( 1 + \frac{2}{3}\frac{\dot{f}}{f_0}t \right) \f$
 */
double ucb_amplitude(double t, double *params, double T);

/**
 \brief Compute UCB time-dependent phase and amplitude
 
 @param[in] params[] source parameters
 @param[in] N number of time samples
 @param[in] times time grid for phase and amplitude calculation
 @param[out] phase array of \f$ \Phi(t) \f$
 @param[out] phase amp of \f$ A(t) \f$
 @param[in] T total observation time (s)
 */
void ucb_barycenter_waveform(double *params, int N, double *times, double *phase, double *amp, double T);

/**
 \brief computes Fisher Information Matrix for UCB waveform parameters Source::params
 
 Computes matrix elements numerically using central differencing
 \f$\Gamma_{ij} = \frac{\partial h}{\partial \theta_i} \frac{\partial h}{\partial \theta_j} \f$
 and stores in Source::fisher_matrix.
 Matrix eigenvectors and eigenvalues are then computed using matrix_eigenstuff() and stored in Source::fisher_evectr and Source::fisher_evalue, respectively.
 
 @param[in] Source::params
 @param[out] Source::fisher_matrix
 @param[out] Source::fisher_evectr
 @param[out] Source::fisher_evalue
 */
void ucb_fisher(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise);

/**
 \brief computes Fisher Information Matrix for UCB waveform parameters Source::params in wavelet domain
 @see ucb_fisher()
 */
void ucb_fisher_wavelet(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise);

/**
 \brief aligns generated UCB waveform with Data array
 
 The waveforms generated by ucb_waveform() are symmetric about the carrier frequency \f$f_0\f$.  To use the template waveform when analyzing the data it must be aligned with the frequency segment. This function computes the bandwidth of the source with parameters Source::params, and computes the offset between the first frequency bin of the template, Source::qmin, and the first frequency bin of the data Data::qmin.
 
 @param[in] Source::params
 @param[out] Source::BW bandwidth in frequency bins
 @param[out] Source::qmin absolute initial bin of template
 @param[out] Source::qmax absolute final bin of template
 @param[out] Source::imin relative initial bin of template w.r.t. Data
 @param[out] Source::imax realtive final bind of template w.r.t. Data
 */
void ucb_alignment(struct Orbit *orbit, struct Data *data, struct Source *source);

/**
 \brief computes frequency width of template based on SNR, Doppler spreading, sinc spreading from finite sampling, and frequency evolution.
 
 @param L average armlength
 @param fstar transfer frequency \f$ c / (2\pi L)\ [{\rm Hz}]\f$
 @param f central frequency \f$ f_0\ [{\rm Hz}]\f$
 @param fdot frequency derivative \f$ \dot f\ [{\rm Hz}\ {\rm s}^{-1}]\f$
 @param costheta cosine ecliptic co-latitude \f$ \cos\theta \f$
 @param A gravitational wave amplitude \f$ \mathcal{A}\f$
 @param T observation time \f$ T_{\rm obs}\ [{\rm s}]\f$
 @param N data length \f$ N [{\rm bins}]\f$
 @return BW bandwidth [bins]
 */
int ucb_bandwidth(double L, double fstar, double f, double fdot, double costheta, double A, double T, int N);

/**
 \brief Galactic binary waveform generator using fast-slow decomposition first described in <a href="https://journals.aps.org/prd/abstract/10.1103/PhysRevD.76.083006">Cornish and Littenberg, PRD 76, 083006</a>.

 Computes the frequency domain TDI response to a circular, slowly evolving, binary with parameters \p params.  The detector geometry is defined in \p Orbit.  The format of the TDI data, either "phase", "frequency", or "sangria", is specified by format.  The TDI response is defined in LISA_tdi() (phase), LISA_tdi_FF() (frequency), or LISA_tdi_Sangria() (frequency circa LDC2.1).
 
 The TDI response is returned for the Michelson-like channels (XYZ), or the two orthogonal channels A and E.  The T channel is practically a noise-only channel at typical UCB frequencies and is therefore neglected.
 
 The array params[] is required to contain parameters: \f$ f_0T,\cos\theta,\phi,\log\mathcal{A},\cos\iota,\psi,\varphi_0\f$,
 where \f$\theta\f$ is the ecliptic co-latitude and \f$\phi\f$ is the ecliptic longitude.
 Extra parameters for non-zero frequency evolution are optional: \f$ \dot{f}T^2, \ddot{f}T^3 \f$.
 
 
 @param[in] orbit LISA ephemerides
 @param[in] format TDI format, "phase" or "frequency" or "sangria"
 @param[in] T observation time \f$ T_{\rm obs}\ [{\rm s}]\f$
 @param[in] t0 start time of observations \f$ t_0\ [{\rm s}]\f$
 @param[in] params[] source parameters
 @param[in] NParams number of source parameters (7, 8, or 9)
 @param[in] BW source bandwidth [bins]
 @param[in] NI number of interferometer channels (1 for X, 2 for A,E, 3 for X,Y,Z)
 @param[out] X,Y,Z Michelson channels
 @param[out] A,E noise orthogonal TDI channels
 */
void ucb_waveform(struct Orbit *orbit, char *format, double T, double t0, double params[], int NParams, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI);

/**
 \brief Wavelet domain ultra compact binary waveform generators as first described in <a href="https://https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.124038">Cornish, PRD 102, 124038</a>.

 Computes the wavelet domain TDI response to a circular, slowly evolving, binary with parameters params.  The detector geometry is defined in Orbit.  The format of the TDI data is hard-coded to match the conventions of LDC2.1.
 
 The TDI response is returned for the Michelson-like channels X,Y,Z.

 The array params[] is required to contain parameters: \f$ f_0T,\cos\theta,\phi,\log\mathcal{A},\cos\iota,\psi,\varphi_0\f$,
 where \f$\theta\f$ is the ecliptic co-latitude and \f$\phi\f$ is the ecliptic longitude.
 Extra parameters for non-zero frequency evolution are optional: \f$ \dot{f}T^2, \ddot{f}T^3 \f$.
 
 @param[in] orbit LISA ephemerides
 @param[in] wdm defines wavelet basis
 @param[in] Tobs observation time \f$ T_{\rm obs}\ [{\rm s}]\f$
 @param[in] t0 start time of observations \f$ t_0\ [{\rm s}]\f$
 @param[in] params[] source parameters
 @param[in] wavelet_list list of active wavelet pixels for A and E channel
 @param[out] X,Y,Z wavelet domain TDI channels
 */
void ucb_waveform_wavelet(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z);

/**
 * @brief Wavelet domain UCB waveform using lookup table transform.
 * @see ucb_waveform_wavelet()
 */
void ucb_waveform_wavelet_tab(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z);


#endif /* ucb_waveform_h */
