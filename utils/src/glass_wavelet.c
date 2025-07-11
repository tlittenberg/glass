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


#include "glass_utils.h"

struct TimeFrequencyTrack * malloc_time_frequency_track(struct Wavelets *wdm)
{
    struct TimeFrequencyTrack *track = malloc(sizeof(struct TimeFrequencyTrack));
    track->min_layer = 1;
    track->max_layer = wdm->NF;
    track->segment_size  = int_vector(wdm->NF);
    track->segment_midpt = int_vector(wdm->NF);
    return track;
}

void free_time_frequency_track(struct TimeFrequencyTrack *track)
{
    free_int_vector(track->segment_size);
    free_int_vector(track->segment_midpt);
    free(track);
}

static void setup_wdm_basis(struct Wavelets *wdm, int NF)
{
    wdm->NF = NF;
    wdm->cadence = WAVELET_DURATION/wdm->NF;
    wdm->Omega = M_PI/wdm->cadence;
    wdm->dOmega = wdm->Omega/(double)wdm->NF;
    wdm->inv_root_dOmega = 1.0/sqrt(wdm->dOmega);
    wdm->B = wdm->Omega/(double)(2*wdm->NF);
    wdm->A = (wdm->dOmega-wdm->B)/2.0;
    wdm->BW = (wdm->A+wdm->B)/M_PI;
}

static double phitilde(struct Wavelets *wdm, double omega)
{
    double insDOM = wdm->inv_root_dOmega;
    double A = wdm->A;
    double B = wdm->B;
    
    double x, y, z;
    
    z = 0.0;
    
    if(fabs(omega) >= A && fabs(omega) < A+B)
    {
        x = (fabs(omega)-A)/B;
        y = incomplete_beta_function(WAVELET_FILTER_CONSTANT, WAVELET_FILTER_CONSTANT, x);
        z = insDOM*cos(y*M_PI/2.0);
    }
    
    if(fabs(omega) < A) z = insDOM;
    
    return(z);
    
}

static void wavelet(struct Wavelets *wdm, int m, double *wave)
{
    
    int N = wdm->N;
    double dom = wdm->domega;
    double DOM = wdm->dOmega;
    
    double omega;
    double x, y, z;
    
    double *DE = (double*)malloc(sizeof(double)*(2*N));
    
    // zero and postive frequencies
    for(int i=0; i<=N/2; i++)
    {
        omega = (double)(i)*dom;
        
        y = phitilde(wdm,omega+(double)(m)*DOM);
        z = phitilde(wdm,omega-(double)(m)*DOM);
        
        x = y+z;
        
        REAL(DE,i) = M_SQRT1_2*x;
        IMAG(DE,i) = 0.0;
        
    }
    
    // negative frequencies
    for(int i=1; i< N/2; i++)
    {
        omega = -(double)(i)*dom;
        
        y = phitilde(wdm,omega+(double)(m)*DOM);
        z = phitilde(wdm,omega-(double)(m)*DOM);
        
        x = y+z;
        
        REAL(DE,N-i) = M_SQRT1_2*x;
        IMAG(DE,N-i) = 0.0;
        
    }
    
    glass_inverse_complex_fft(DE,N);

    for(int i=0; i < N/2; i++)
    {
        wave[i] = REAL(DE,N/2+i)/wdm->norm;
        wave[i+N/2] = REAL(DE,i)/wdm->norm;
    }

    free(DE);
}

static void wavelet_window_time(struct Wavelets *wdm)
{
    double *DX = (double*)malloc(sizeof(double)*(2*wdm->N));
    
    //zero frequency
    REAL(DX,0) =  wdm->inv_root_dOmega;
    IMAG(DX,0) =  0.0;
    
    for(int i=1; i<= wdm->N/2; i++)
    {
        int j = wdm->N-i;
        double omega = (double)(i)*wdm->domega;
        
        // postive frequencies
        REAL(DX,i) = phitilde(wdm,omega);
        IMAG(DX,i) =  0.0;
        
        // negative frequencies
        REAL(DX,j) =  phitilde(wdm,-omega);
        IMAG(DX,j) =  0.0;
    }
        
    glass_inverse_complex_fft(DX, wdm->N);

    wdm->window = (double*)malloc(sizeof(double)* (wdm->N));
    for(int i=0; i < wdm->N/2; i++)
    {
        wdm->window[i] = REAL(DX,wdm->N/2+i);
        wdm->window[wdm->N/2+i] = REAL(DX,i);
    }
    
    wdm->norm = sqrt((double)wdm->N * wdm->cadence / wdm->domega);

    free(DX);
}

static void wavelet_lookup_table(struct Wavelets *wdm)
{    
    // it turns out that all the wavelet layers are the same modulo a
    // shift in the reference frequency. Just have to do a single layer
    // we pick one far from the boundaries to avoid edge effects
    
    double *wave = (double*)malloc(sizeof(double)*(wdm->N));
    int ref_layer = wdm->NF/2;
    wavelet(wdm, ref_layer, wave);
    
    // The odd wavelets coefficienst can be obtained from the even.
    // odd cosine = -even sine, odd sine = even cosine
    // each wavelet covers a frequency band of width DW
    // execept for the first and last wasvelets
    // there is some overlap. The wavelet pixels are of width
    // DOM/PI, except for the first and last which have width
    // half that
    
    double f0 = ref_layer*wdm->df;
    
    #pragma omp parallel for
    for(int j=0; j<wdm->fdot_steps; j++)  // loop over f-dot slices
    {

        int NT = wdm->n_table[j];
        
        
        for(int n=0; n<NT; n++)  // loop of frequency slices
        {
            double f = f0 + ((double)(n-NT/2)+0.5)*wdm->deltaf;
            
            double real_coeff = 0.0;
            double imag_coeff = 0.0;
            
            for(int i=0; i<wdm->N; i++)
            {
                double t = ((double)(i-wdm->N/2))*wdm->cadence;
                double phase = PI2*f*t + M_PI*wdm->fdot[j]*t*t;
                real_coeff += wave[i]*cos(phase)*wdm->cadence;
                imag_coeff += wave[i]*sin(phase)*wdm->cadence;
            }
            wdm->table[j][2*n]   = real_coeff;
            wdm->table[j][2*n+1] = imag_coeff;
        }
    }
    
    free(wave);
}
void initialize_wavelet(struct Wavelets *wdm, double T)
{
    fprintf(stdout,"\n======= Initialize Wavelet Basis =======\n");
    
    wdm->NT = (int)ceil(T/WAVELET_DURATION);
    wdm->NF = WAVELET_DURATION/LISA_CADENCE;
    wdm->df = WAVELET_BANDWIDTH;
    wdm->dt = WAVELET_DURATION;

    setup_wdm_basis(wdm, wdm->NF);

    wdm->frequency_steps = 400;
    wdm->fdot_steps = 50;
    wdm->d_fdot = 0.1;
    wdm->oversample = 16.0;

    wdm->N = wdm->oversample * 2 * wdm->NF;
    wdm->T = wdm->N*wdm->cadence;

    wdm->domega = PI2/wdm->T;
    
    wdm->deltaf = wdm->BW/(double)(wdm->frequency_steps);

    wdm->fdot = malloc(wdm->fdot_steps*sizeof(double));

    /* Only needed when using the lookup table waveform generator
    wdm->table   = malloc(wdm->fdot_steps*sizeof(double *));
    wdm->n_table = malloc(wdm->fdot_steps*sizeof(int));

    double fdot_step = wdm->df/wdm->T*wdm->d_fdot; // sets the f-dot increment
        
    for(int n=0; n<wdm->fdot_steps; n++)
    {
        wdm->fdot[n] = -fdot_step*wdm->fdot_steps/2 + n*fdot_step;
        size_t N = (int)((wdm->BW+fabs(wdm->fdot[n])*wdm->T)/wdm->deltaf);
        if(N%2 != 0) N++; // makes sure it is an even number
        wdm->n_table[n] = N;
        wdm->table[n] = double_vector(2*N);
    }
        
    //stores lookup table of wavelet basis functions
    wavelet_lookup_table(wdm);
    */
    
    //stores window function and normalization
    wavelet_window_time(wdm);

    //set defaults for min and maximum pixels
    wavelet_pixel_to_index(wdm,0,1,&wdm->kmin);         //first pixel of second layer
    wavelet_pixel_to_index(wdm,0,wdm->NF-1,&wdm->kmax); //first pixel of last layer

    fprintf(stdout,"  Number of time pixels:        %i\n", wdm->NT);
    fprintf(stdout,"  Duration of time pixels:      %g [hr]\n", wdm->dt/3600);
    fprintf(stdout,"  Number of frequency layers:   %i\n", wdm->NF);
    fprintf(stdout,"  Bandwidth of frequency layer: %g [uHz]\n", wdm->df*1e6);
    fprintf(stdout,"\n========================================\n");
}

void wavelet_index_to_pixel(struct Wavelets *wdm, int *i, int *j, int k)
{
    int NT = wdm->NT;
    
    //which time
    *i = k%NT; 
    
    //which frequency
    *j = (k - (*i))/NT; //scary integer math
}

void wavelet_pixel_to_index(struct Wavelets *wdm, int i, int j, int *k)
{
    int NT = wdm->NT;
    
    *k = i + j*NT;
}

void wavelet_transform(struct Wavelets *wdm, double *data)
{
    //array index for tf pixel
    int k;
    
    //total data size
    int ND = wdm->NT*wdm->NF;
    
    //windowed data packets
    double *wdata = double_vector(wdm->N);

    //wavelet wavepacket transform of the signal
    double **wave = double_matrix(wdm->NT,wdm->NF);
    
    //normalization factor
    double fac = M_SQRT2*sqrt(wdm->cadence)/wdm->norm;
    
    //normalization fudge factor
    fac *= sqrt(wdm->cadence)/2;
        
    //do the wavelet transform by convolving data w/ window and FFT
    for(int i=0; i<wdm->NT; i++)
    {
        
        for(int j=0; j<wdm->N; j++)
        {
            int n = i*wdm->NF - wdm->N/2 + j;
            if(n < 0)   n += ND;  // periodically wrap the data
            if(n >= ND) n -= ND;  // periodically wrap the data
            wdata[j] = data[n] * wdm->window[j];  // apply the window
        }
                
        glass_forward_real_fft(wdata, wdm->N);

        //unpack Fourier transform
        wave[i][0] = wdata[0];
        for(int j=1; j<wdm->NF; j++)
        {
            int n = j*wdm->oversample;
            if((i+j)%2 ==0)
                wave[i][j] = wdata[2*n];
            else
                wave[i][j] = -wdata[2*n+1];
        }
    }
    
    //replace data vector with wavelet transform mapped from pixel to index
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            //get index number k for tf pixel {i,j}
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //replace data array
            data[k] = wave[i][j]*fac;
        }
    }
    
    free_double_vector(wdata);
    free_double_matrix(wave,wdm->NT);
}

void wavelet_tansform_inverse_fourier(struct Wavelets *wdm, double *data)
{
    int k;
    int N = wdm->NT*wdm->NF;
    double *phit  = double_vector(wdm->NT/2+1);
    double *row   = double_vector(wdm->NT*2);
    double *work  = double_vector(N);
    double Tobs   = N*wdm->cadence;
    double sign;


    for(int i=0; i<=wdm->NT/2; i++)
    {
        phit[i] = phitilde(wdm,i*PI2/Tobs);
    }

    for(int j=1; j<wdm->NF-1; j++)
    {
        if(j%2==0) sign =  1.0;
        else       sign = -1.0;

        for(int i=0; i<wdm->NT; i++)
        {
            REAL(row,i) = 0.0;
            IMAG(row,i) = 0.0;
            
            wavelet_pixel_to_index(wdm,i,j,&k);
                        
            if((i+j)%2==0)
            {
                REAL(row,i) = data[k];
            }
            else
            {
                if(j%2==0) IMAG(row,i) = -data[k];
                else       IMAG(row,i) =  data[k];
            }
        }
        
        glass_forward_complex_fft(row,wdm->NT);

        int jj = j*(wdm->NT/2);
        
        // negative frequencies
        for(int i=wdm->NT/2-1; i>0; i--)
        {
            double x = sign*phit[i];
            int kk = jj-i;
            work[kk] += x*REAL(row,wdm->NT-i);
            work[N-kk] += x*IMAG(row,wdm->NT-i);
        }
                
        // positive frequencies
        for(int i=0; i<wdm->NT/2; i++)
        {
            double x = sign*phit[i];
            int kk = i+jj;
            work[kk] += x*REAL(row,i);
            work[N-kk] += x*IMAG(row,i);
        }
    }
    
    //unpack work vector into real and imaginary parts consistent w/ GLASS conventions
    unpack_fft_output(data,work,N);

    //normalize -- wtf?
    double fft_norm = 2.*sqrt(M_PI/Tobs);
    for(int n=0; n<N; n++) data[n] *= fft_norm;


    free_double_vector(phit);
    free_double_vector(work);
    free_double_vector(row);
}

static void fourier_to_wavelet_transform_of_layer(struct Wavelets *wdm, double *window, double *data, int N, int layer)
{
    int i,n;
    double *wdata = double_vector(2*N);
    
    // window data
    for(i=-N/2; i<N/2; i++)
    {
        n = i+N/2;
        
        REAL(wdata,n) = 0.0;
        IMAG(wdata,n) = 0.0;
        
        if(n > 0 && n < N)
        {
            REAL(wdata,n) = data[2*n]   * window[abs(i)];
            IMAG(wdata,n) = data[2*n+1] * window[abs(i)];
        }
    }//end loop over window
    
    //get windowed data into time domain
    //to get FFTed
    glass_inverse_complex_fft(wdata, N);
    
    //normalization to match NJC's breadboard code
    for(int n=0; n<2*N; n++) wdata[n]/=(double)N;

    //populate layer
    for(n=0; n<N; n++)
    {
        if(layer%2 == 0)
        {
            if((n+layer)%2==0) data[n] = REAL(wdata,n);
            else               data[n] = IMAG(wdata,n);
        }
        else
        {
            if((n+layer)%2==0) data[n] =  REAL(wdata,n);
            else               data[n] = -IMAG(wdata,n);
        }
    }//end loop over time slices
    
    free_double_vector(wdata);
}

void wavelet_transform_by_layers(struct Wavelets *wdm, int jmin, int Nlayers, double *window, double *data)
{
    int i,j,k,m,n;
    
    // size of input data
    int N = (Nlayers+1)*wdm->NT;
    
    // transformed data
    double *data_wdm = double_vector(N);
    
    //windowed data
    double *wdata = double_vector(2*wdm->NT);

    double norm = 1.0/sqrt(0.5*N);
        
    double alpha = 8.0/wdm->NT;
    
    // FFT incoming data (timeseries)
    tukey(data, alpha, N);
    glass_forward_real_fft(data, N);
    
    // loop over frequency layers
    for(j=1; j<Nlayers+1; j++)
    {
        m = jmin + j - 1;
        
        // window data
        for(i=-wdm->NT/2; i<wdm->NT/2; i++)
        {
            n = i+wdm->NT/2;
            
            REAL(wdata,n) = 0.0;
            IMAG(wdata,n) = 0.0;
            
            k = i + j*wdm->NT/2;
            
            if(k > 0 && k < N/2)
            {
                REAL(wdata,n) = data[2*k]   * window[abs(i)];
                IMAG(wdata,n) = data[2*k+1] * window[abs(i)];
            }
        }//end loop over window
        
        glass_inverse_complex_fft(wdata, wdm->NT);
        
        // index magic
        for(i=0; i<wdm->NT; i++)
        {
            k = i*Nlayers + j - 1;
            
            if(m%2 == 0)
            {
                if((i+m)%2==0) data_wdm[k] = norm*REAL(wdata,i);
                else           data_wdm[k] = norm*IMAG(wdata,i);
            }
            else
            {
                if((i+m)%2==0) data_wdm[k] =  norm*REAL(wdata,i);
                else           data_wdm[k] = -norm*IMAG(wdata,i);
            }
        }//end loop over time slices
    }
    
    //replace input data w/ WDM'ed data
    memcpy(data, data_wdm, N*sizeof(double));
    
    free_double_vector(data_wdm);
    free_double_vector(wdata);
}


void wavelet_transform_inverse_time(struct Wavelets *wdm, double *data)
{
    int N = wdm->NT*wdm->NF;
    
    //transform data from WDM to frequency domain
    wavelet_tansform_inverse_fourier(wdm,data);
    
    //transform to time domain data
    glass_inverse_real_fft(data,N);
}

void wavelet_transform_from_table(struct Wavelets *wdm, double *phase, double *freq, double *freqd, double *amp, int *jmin, int *jmax, double *wave, int *list, int *rlist, int Nmax)
{
    
    int n, k, jj, kk;
    double dx, dy;
    double f, fdot;
    double fmid,fsam;
    double cos_phase, sin_phase, y, z, yy, zz;

    double df = wdm->deltaf;

    // maximum frequency and frequency derivative
    double f_max    = wdm->df*(wdm->NF-1);
    double fdot_max = wdm->fdot[wdm->fdot_steps-1];
    double fdot_min = wdm->fdot[0];
    double d_fdot   = wdm->fdot[1]-wdm->fdot[0]; // f-dot increment
    
    for(int i=0; i<wdm->NT; i++)
    {
        f     = freq[i];
        fdot  = freqd[i];
        
        //skip this step if f or fdot violate bounds
        if(f>=f_max || fdot>=fdot_max || fdot<=fdot_min)  continue;
        
        cos_phase = amp[i]*cos(phase[i]);
        sin_phase = amp[i]*sin(phase[i]);
        
        n = (int)floor((fdot-fdot_min)/d_fdot);  // lower f-dot layer
        dy = (fdot-fdot_min)/d_fdot - n;         // where in the layer
                                    
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
        
            // central frequency
            fmid = j*wdm->df;
                
            kk = (int)floor( ( f - (fmid + 0.5*df) )/df );
            fsam = fmid + (kk + 0.5)*df;
            dx = (f - fsam)/df; // used for linear interpolation
                
            // interpolate over frequency
            y = 0.0;
            z = 0.0;
            yy = 0.0;
            zz = 0.0;

            jj = kk + wdm->n_table[n]/2;
            if(jj>=0 && jj< wdm->n_table[n]-1)
            {
                y = (1.0-dx)*wdm->table[n][2*jj]   + dx*wdm->table[n][2*(jj+1)];
                z = (1.0-dx)*wdm->table[n][2*jj+1] + dx*wdm->table[n][2*(jj+1)+1];
            }

            jj = kk + wdm->n_table[n+1]/2;
            if(jj >=0 && jj < wdm->n_table[n+1]-1)
            {
                yy = (1.0-dx)*wdm->table[n+1][2*jj]   + dx*wdm->table[n+1][2*(jj+1)];
                zz = (1.0-dx)*wdm->table[n+1][2*jj+1] + dx*wdm->table[n+1][2*(jj+1)+1];
            }
                
            // interpolate over fdot
            y = (1.0-dy)*y + dy*yy;
            z = (1.0-dy)*z + dy*zz;

            // make sure pixel is in range
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {  
                int n = rlist[k - wdm->kmin];
                if(n<Nmax)
                {
                    if((i+j)%2 == 0) wave[n] =  (cos_phase*y - sin_phase*z);
                    else             wave[n] = -(cos_phase*z + sin_phase*y);
                }
                else
                {
                    //fprintf(stderr,"Warning, wavelet_transform_from_table tried accessing array out of bounds\n");
                    //fflush(stderr);
                }
            }
            
        } //loop over frequency layers
        
    } //loop over time steps
}

void active_wavelet_list(struct Wavelets *wdm, double *freqX, double *freqY, double *freqZ, double *fdotX, double *fdotY, double *fdotZ, int *wavelet_list, int *reverse_list, int *Nwavelet, int *jmin, int *jmax)
{
    
    int n;
    int k;
    int N;
    double fmx, fdmx, fdmn, dfd, HBW;
    double fmax, fmin;
    double fdotmax, fdotmin;
    int Xflag, Yflag, Zflag;
    
    double df = wdm->deltaf;
    double DF = wdm->df;
    double *fd = wdm->fdot;

    // maximum frequency and frequency derivative
    fmx  = (double)(wdm->NF-1)*DF;
    fdmx = fd[wdm->fdot_steps-1];
    fdmn = fd[0];
    dfd  = fd[1]-fd[0]; // f-dot increment
    
    N = 0;
    for(int i=0; i<wdm->NT; i++)
    {
        // check to see if any of the channels are ok
        Xflag = Yflag = Zflag = 0;
        if(freqX[i] < fmx) Xflag = 1;
        if(freqY[i] < fmx) Yflag = 1;
        if(freqZ[i] < fmx) Zflag = 1;

        // shut off any channel that does not have valid fdots
        if(fdotX[i] < fdmn || fdotX[i] > fdmx) Xflag = 0;
        if(fdotY[i] < fdmn || fdotY[i] > fdmx) Yflag = 0;
        if(fdotZ[i] < fdmn || fdotZ[i] > fdmx) Zflag = 0;

        // skip if no channels have valid values
        if(!Xflag && !Yflag && !Zflag) continue;

        /*  find the largest and smallest frequencies and frequency derivatives
        but only using the valid channels */
        fmin = 1;
        fmax = 0;
        fdotmin =  1;
        fdotmax = -1;

        if(Xflag)
        {
            if(freqX[i]>fmax) fmax=freqX[i];
            if(freqX[i]<fmin) fmin=freqX[i];
            if(fdotX[i]>fdotmax) fdotmax=fdotX[i];
            if(fdotX[i]<fdotmin) fdotmin=fdotX[i];
        }
        
        if(Yflag)
        {
            if(freqY[i]>fmax) fmax=freqY[i];
            if(freqY[i]<fmin) fmin=freqY[i];
            if(fdotY[i]>fdotmax) fdotmax=fdotY[i];
            if(fdotY[i]<fdotmin) fdotmin=fdotY[i];
        }
        
        if(Zflag)
        {
            if(freqZ[i]>fmax) fmax=freqZ[i];
            if(freqZ[i]<fmin) fmin=freqZ[i];
            if(fdotZ[i]>fdotmax) fdotmax=fdotZ[i];
            if(fdotZ[i]<fdotmin) fdotmin=fdotZ[i];
        }
       
        //skip if max/min fdot go out of bounds 
        if(fdotmax >= fdmx || fdotmin <= fdmn) continue;

        // lowest f-dot layer
        n = (int)(floor(fdotmin-fdmn/dfd));
        int NL = wdm->n_table[n];

        // highest f-dot layer
        n = (int)(floor(fdotmax-fdmn/dfd));
        int NH = wdm->n_table[n];

        // find which has the largest number of samples
        if(NL > NH) NH = NL;

        // half bandwidth of layer        
        HBW = 0.5*(NH-1)*df;
        
        // lowest frequency layer
        jmin[i] = (int)ceil((fmin-HBW)/DF);
        
        // highest frequency layer
        jmax[i] = (int)floor((fmax+HBW)/DF);   

        // skip any out-of-bounds layers
        if(jmin[i] < 0) jmin[i] = 0;
        if(jmax[i] > wdm->NF-1) jmax[i] = wdm->NF-1;
        
        for(int j=jmin[i]; j<=jmax[i]; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //check that the pixel is in range
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                wavelet_list[N]=k-wdm->kmin;
                reverse_list[k-wdm->kmin]=N;
                N++;
            }
        }  
    }
    *Nwavelet = N;
}

void wavelet_window_frequency(struct Wavelets *wdm, double *window, int Nlayers)
{
    int i;
    int N;
    double T;
    double domega;
    double omega;
    double norm=0.0;
    
    N = (Nlayers+1);

    //mini wavelet structure for basis covering just N layers
    struct Wavelets *wdm_temp = malloc(sizeof(struct Wavelets));
    setup_wdm_basis(wdm_temp, N);
    
    T = wdm->dt*wdm->NT;
    
    domega = PI2/T;

    //wdm window function
    for(i=0; i<=wdm->NT/2; i++)
    {
        omega = i*domega;
        window[i] = phitilde(wdm_temp,omega);
    }
    
    //normalize
    for(i=-wdm->NT/2; i<= wdm->NT/2; i++) norm += window[abs(i)]*window[abs(i)];
    norm = sqrt(norm/wdm_temp->cadence);
    
    for(i=0; i<=wdm->NT/2; i++) window[i] /= norm;
    
    free(wdm_temp);
    
}

void wavelet_transform_segment(struct Wavelets *wdm, int N, int layer, double *data)
{
    //wdm filter
    double *window = double_vector(N/2+1);
    
    //normalization factor
    double norm = 1./(2.0*wdm->NF*wdm->cadence);
    
    //layer-dependent factors
    double domega = PI2/(N*WAVELET_DURATION);
    
    //get filter
    for(int i=0; i<=N/2; i++)
        window[i] = norm * phitilde(wdm, i*domega);
        
    //wdm transfrom from frequency domain
    fourier_to_wavelet_transform_of_layer(wdm, window, data, N, layer);
    
    free_double_vector(window);
}
