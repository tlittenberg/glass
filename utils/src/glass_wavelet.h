/*
 * Copyright 2024 Neil J. Cornish & Tyson B. Littenberg
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

#ifndef glass_wavelet_h
#define glass_wavelet_h

#define WAVELET_FILTER_CONSTANT 4  // filter steepness in frequency
#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

struct Wavelets
{
    /** @name define time-frequency-fdot grid */
     ///@{
    int NF;    //!<total number of frequency layers
    int NT;    //!<total number of time layers
    double cadence; //!<data sample cadence
    double dt; //!<wavelet pixel duration
    double df; //!<wavelet pixel bandwidth
    
    int frequency_steps; //!<frequency steps for wavelet filters
    int fdot_steps;      //!<fdot steps
    double d_fdot;       //!<fractional fdot increment
    ///@}

    int oversample;   //!<oversampling factor for basis functions
    
    double Omega;
    double dOmega;
    double domega;
    double inv_root_dOmega;
    double B;
    double A;
    double BW;
    double deltaf;
    double *fdot;
    double *window;
    double norm;
    
    int N;
    double T;

    int kmin; //!<minimum wavelet index
    int kmax; //!<maximum wavelet index
    
    int *n_table; //!< number of terms in the lookup table at each frequency
    double **table; //!< lookup table of wavelet coefficients
};

void initialize_wavelet(struct Wavelets *wdm, double T);
void wavelet_window_frequency(struct Wavelets *wdm, double *window, int Nlayers);
void wavelet_transform(struct Wavelets *wdm, double *data);
void wavelet_transform_inverse(struct Wavelets *wdm, double *data);
void wavelet_to_fourier_transform(struct Wavelets *wdm, double *data);
void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k);
void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k);
void wavelet_transform_by_layers(struct Wavelets *wdm, int jmin, int Nlayers, double *window, double *data);
void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax);
void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax);

#endif /* glass_wavelet_h */
