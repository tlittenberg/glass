/*
 * Copyright 2025 Tyson B. Littenberg & Neil Cornish
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
@file glass_mbh_waveform.h
\brief Massive Black Hole binary waveform generator
 */

#ifndef glass_mbh_waveform_h
#define glass_mbh_waveform_h

#include <stdio.h>

/**
 \brief Compute MBH time/frequency-dependent phase and amplitude
 
 For frequency domain waveforms \p freq[] is input to the function and
 \p time[] is populated with \f$ t(f)\f$, as are \p amp[] with \f$ A(f)\f$ and \p phase[] with \f$ \Phi(f)\f$.
 For time domain waveforms \p time[] is input to the function and \p freq[]
 is populated with \f$ f(t) \f$ as well as \p amp[] with \f$ A(t)\f$ and \p phase[] with \f$ \Phi(t)\f$.
 
 @param[in] params[] source parameters
 @param[in] N number of time samples
 @param[in,out] freq frequency grid for phase and amplitude calculation
 @param[in,out] time time grid for phase and amplitude calculation
 @param[out] phase array of \f$ \Phi(t,f) \f$
 @param[out] amp array of \f$ A(t,f) \f$
 @param[in] model[] string for selecting waveform family. Supported options are "IMRPhenomD" and "IMRPhenomT"
 */
void mbh_barycenter_waveform(double *params, int N, double *freq, double *time, double *phase, double *amp, char model[]);

/**
 \brief Comput final spin of MBH merger remnant
 
 @param MBH parameter vector
 @return \f$ \chi_f \f$ using IMRPhenomD approximant
 */
double mbh_final_spin(double *params);

/**
 \brief Compute ringdown frequency of MBH merger remnant
 
 @param MBH parameter vector
 @return \f$ \omega_f \f$ using IMRPhenomD approximant
 */
double mbh_ringdown_frequency(double *params);


/**
 \brief Massive black hole binary merger waveform generator using fast-slow decomposition first described in <a href="https://arxiv.org/abs/2506.08093">Cornish and Littenberg, arXiv:2506.08093</a>.

 Computes the frequency domain TDI response to a circular, slowly evolving, binary with parameters params.  The detector geometry is defined in Orbit.  The conventions for the TDI data are assumed to be fractional frequency. The TDI response is defined in LISA_spline_response(). The TDI response is returned for the Michelson-like channels (XYZ).
 
 The array \p params[] is required to contain parameters:
 \f$
 \log M_c,
 \log M_t,
 \chi_1,
 \chi_2,
 \phi_c,
 t_c,
 \log D_L,
 \theta,
 \phi,
 \psi,
 \cos\iota\f$
 
 where
 \f$ M_c \f$ and \f$ M_t \f$ are the chirp mass and total mass in \f$ M_\odot \f$,
 \f$ \chi_i \f$ are the dimensionless spins of the component BHs, assumed to be aligned with the orbit,
 \f$ \phi_c\f$ and \f$ t_c \f$ are the phase and time at coalesence in seconds,
 \f$ D_L \f$ is the luminosity distance in Gpc,
 \f$\theta\f$ is the ecliptic latitude and \f$\phi\f$ is the ecliptic longitude in radians,
 \f$ \psi \f$ is the polarization angle in radians and \f$ \iota \f$ is the inclination angle.
 
 
 @param[in] orbit LISA ephemerides
 @param[in] wdm Wavelet basis
 @param[in] Tobs observation time \f$ T_{\rm obs}\ [{\rm s}]\f$
 @param[in] t0 start time of observations \f$ t_0\ [{\rm s}]\f$
 @param[in] params vector of MBH parameters
 @param[out] wavelet_list list of active wavelet pixels for the signal
 @param[out] Nwavelet number of active wavelet pixels for the signal
 @param[out] X,Y,Z Michelson channels
 */
void mbh_fd_waveform(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z);

/**
 \brief Same as mbh_fd_waveform() but for time-domain waveform.
 
 @see mbh_fd_wavefor()
 m
 */
void mbh_td_waveform(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z);

#endif /* glass_mbh_waveform_h */
