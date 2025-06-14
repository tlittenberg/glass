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

#include <glass_utils.h>
#include "glass_mbh_IMRPhenom.h"
#include "glass_mbh_waveform.h"

double mbh_final_spin(double *params)
{
    double Mchirp = exp(params[0]); // chirp mass
    double Mtotal = exp(params[1]); // total mass
    double chi1   = params[2];      // Spin1
    double chi2   = params[3];      // Spin2

    double eta = symmetric_mass_ratio(Mchirp, Mtotal);

    return mbh_IMRPhenomD_final_spin_wrapper(eta, chi1, chi2);
}

double mbh_ringdown_frequency(double *params)
{
    double Mchirp = exp(params[0]); // chirp mass
    double Mtotal = exp(params[1]); // total mass
    double chi1   = params[2];      // Spin1
    double chi2   = params[3];      // Spin2
    
    return mbh_IMRPhenomD_ringdown_frequency_wrapper(Mchirp, Mtotal, chi1, chi2);
}

void mbh_barycenter_waveform(double *params, int N, double *freq, double *time, double *phase, double *amp, char model[])
{
        
    double Mchirp = exp(params[0]); // chirp mass
    double Mtotal = exp(params[1]); // total mass
    double dL     = exp(params[6]); // distance
    double chi1   = params[2];      // Spin1
    double chi2   = params[3];      // Spin2
    double tc     = params[5];      // merger time
    double phic   = params[4];      // merger phase
    
    if(!strcmp("IMRPhenomD",model)) mbh_IMRPhenomD_wrapper(Mchirp, Mtotal, chi1, chi2, dL, tc, phic, freq, time, phase, amp, N);
    if(!strcmp("IMRPhenomT",model)) mbh_IMRPhenomT_wrapper(Mchirp, Mtotal, chi1, chi2, dL, tc, phic, freq, time, phase, amp, N);
}

static void unpack_mbh_tf_data(double *data, double *time, double *omega, int N)
{
    for(int n=0; n<N; n++)
    {
        time[n]  = data[n];
        omega[n] = data[N+n];
    }
}

static void pack_mbh_tf_data(double *data, double *time, double *omega, int N)
{
    for(int n=0; n<N; n++)
    {
        data[n]   = time[n];
        data[N+n] = omega[n];
    }
}

static double * mbh_time_frequency_grid(double *params, int *N, int *Nc)
{
    // PhenomT dynamically figures out how big it is by
    // startint at the peak and working backwards
    //
    // The IMRPhenomT model defined merger at t=0
    // we have to offset this by the merger time tc
    // in the physical time array
    // The spacing in time is designe to give a maximum
    // phase shift from the previous point of dPhase radians
    // This is computed by dividing dPhase by the angular
    // frequency omega. Note that IMRPhenomT uses units where
    // the total mass = 1, so we have to divide by Mtotal in seconds
    // If deltaT = dPhase/omega exceeds dTmax (usually set to about a
    // day) then it gets set to dTmax. Note that once we get more than
    // t = 100 Mtotal from merger we start making dPhase larger to take
    // bigger steps. We work back from merger and forward from merger.
    //
    // The dynamically sized time and omega(t) arrays are packed into
    // data[2*N] with data[n=0..N-1] = t[n] and
    // data[n=N..2*N-1] = omega[n-N];
        
    //set up work space
    int n;
    int Nwork = 0;
    int Nmax  = 100000;
    double *work  = double_vector(Nmax);
    double *time  = double_vector(Nmax);
    double *omega = double_vector(Nmax);

    
    // unpack params vector
    double Mchirp = exp(params[0])*TSUN;
    double Mtotal = exp(params[1])*TSUN;
    double chi1   = params[2];
    double chi2   = params[3];
    double tc     = params[5];
    double eta    = symmetric_mass_ratio(Mchirp,Mtotal);

    // IMRPhonemT internals
    struct IMRPhenomT *IMRPT = setup_IMRPT(Mchirp, Mtotal, chi1, chi2);

    
    // start from merger and work backwards
    double t=0.0;
    double dPhase = 0.5;    //TODO: what is this 0.5
    double dt_max = 2.0e5;  // maximum time step for TDI extraction TODO: what is this?
    double dt;
    do
    {
        //get current time and angular frequency
        omega[Nwork] = mbh_IMRPhenomT_angular_frequency_wrapper(t, eta, Mtotal, IMRPT);
        time[Nwork]  = t + tc;

        //increment time
        dt = dPhase/omega[Nwork];
        if(dt > dt_max) dt = dt_max;
        t -= dt;

        //increase phase steps as we get into inspiral
        if(t < -100*Mtotal) dPhase*=1.1;

        Nwork++;
    } while (t > -(tc + 2.0e5)); //TODO: what is this 2e5  should it be dt_max?
    
    //store index of coalesence time
    *Nc = Nwork-1;

    // time order the omega and time arrays
    for(n=0; n<Nwork; n++) work[n] = time[Nwork-1-n];
    for(n=0; n<Nwork; n++) time[n] = work[n];
    for(n=0; n<Nwork; n++) work[n]  = omega[Nwork-1-n];
    for(n=0; n<Nwork; n++) omega[n] = work[n];

    // now pick up at merger and work forwards until tend
    dPhase = 0.5;
    t = dPhase/omega[Nwork-1];
    double tend = 500.0+1000.0*Mtotal; //TODO: what is this?

    do
    {
        //get current time and angular frequency
        omega[Nwork] = mbh_IMRPhenomT_angular_frequency_wrapper(t, eta, Mtotal, IMRPT);
        time[Nwork]  = t + tc;
        
        //increment time
        dt = dPhase/omega[Nwork];
        if(dt > dt_max) dt = dt_max;
        t += dt;
        
        //increase phase steps as we get into ringdown
        if(t>10*Mtotal) dPhase *= 1.1;
        
        Nwork++;
    } while (t-dt < tend);

    // package results and return
    *N = Nwork;
    

    // double the size of data array to hold both time and omega
    double *data = double_vector(Nwork*2);
    pack_mbh_tf_data(data, time, omega, Nwork);
    
    free_double_vector(work);
    free_double_vector(time);
    free_double_vector(omega);
    free_IMRPhenomT(IMRPT);
    
    return data;
}

void mbh_td_waveform(struct Orbit *orbit, double Tobs, double t0, double *params, double *X, double *Y, double *Z)
{
    int n;
    
    // get time and angular frequency grid
    int Nspline;
    int Nmerger;
    double *data = mbh_time_frequency_grid(params, &Nspline, &Nmerger);
     
    // extract time and angular frequency from packed data array
    double *time_ssb  = double_vector(Nspline);
    double *freq_ssb  = double_vector(Nspline);
    double *amp_ssb   = double_vector(Nspline);
    double *phase_ssb = double_vector(Nspline);
    unpack_mbh_tf_data(data, time_ssb, freq_ssb, Nspline);
     
    mbh_barycenter_waveform(params, Nspline, freq_ssb, time_ssb, phase_ssb, amp_ssb, "IMRPhenomT");
        
    // phase shift to phic at merger
    double phic = params[4];
    double dphi = phic - phase_ssb[Nmerger];
    for(n=0; n<Nspline; n++) phase_ssb[n] += dphi;
    

    /*
     Get spline interpolant for frequency and amplitude on time grid
     */
    struct CubicSpline *amp_ssb_spline   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_ssb_spline = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_ssb_spline,time_ssb,amp_ssb);
    initialize_cubic_spline(phase_ssb_spline,time_ssb,phase_ssb);
    
    // trim the edge of the interpolation domain so that we don't run off the end
    n=0;
    do
    {
        n++;
    } while (time_ssb[Nspline-1] - time_ssb[Nspline-n] < 2.*AU/CLIGHT);
    Nspline -= n;
    
    /*
     Resample SSB phase to reference spacecraft
     */
    double *phase_sc = double_vector(Nspline);
    double *time_sc  = double_vector(Nspline);
    
    // shift times from Barycenter to S/C 1 to get reference phase
    double costh = sin(params[7]); // EclipticLatitude
    double phi   = params[8];      // EclipticLongitude
    LISA_detector_time(orbit, costh, phi, time_ssb, Nspline, time_sc);
    
    // get signal phase at S/C 1
    for(int i=0; i< Nspline; i++)
        phase_sc[i] = spline_interpolation(phase_ssb_spline, time_sc[i]);
    
    /*
     Get TDI reponse for signal's phase and amplitude on time grid
     set by the frequency evolution
     */
    struct TDI *tdi_phase = malloc(sizeof(struct TDI));
    struct TDI *tdi_amp = malloc(sizeof(struct TDI));
    alloc_tdi(tdi_phase,Nspline,3);
    alloc_tdi(tdi_amp,Nspline,3);

    //extract rest of extrinsic parameters from MBH parameter vector
    double cosi  = params[10];
    double psi   = params[9];

    LISA_spline_response(orbit, time_ssb, Nspline, costh, phi, cosi, psi, amp_ssb_spline, NULL, phase_ssb_spline, phase_sc, tdi_amp, tdi_phase);
    
    /*
    FILE *out = fopen("PhenomT_TDI.dat","w");
    for(n=0; n<Nspline; n++)
        fprintf(out,"%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", time_ssb[n], phase_sc[n], tdi_phase->X[n], tdi_phase->Y[n], tdi_phase->Z[n], tdi_amp->X[n], tdi_amp->Y[n], tdi_amp->Z[n]);

    fclose(out);
    */
    
    // restore reference phase
    for(n=0; n<Nspline; n++)
    {
        tdi_phase->X[n] += phase_sc[n];
        tdi_phase->Y[n] += phase_sc[n];
        tdi_phase->Z[n] += phase_sc[n];
    }
    
    /*
    Interpolate amplitude and phase for instrument response of each TDI channel onto wavelet grid
    
    
    //set up time grid
    int N = (int)(Tobs/LISA_CADENCE);
    double *time = double_vector(N);
    for(n=0; n<N; n++) time[n] = t0 + n*LISA_CADENCE;
    
    double Amp,Phase;
    struct TDI *wave = malloc(sizeof(struct TDI));
    alloc_tdi(wave,N,3);

    struct CubicSpline *amp_interpolant   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_interpolant = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->X);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->X);

    for(n=0; n<N; n++)
    {
        if(time[n] < time_ssb[Nspline-1])
        {
            Amp    = spline_interpolation(amp_interpolant, time[n]);
            Phase  = spline_interpolation(phase_interpolant, time[n]);
            
            wave->X[n] = Amp*cos(Phase);
        }
    }
    
    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->Y);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->Y);

    for(n=0; n<N; n++)
    {
        if(time[n] < time_ssb[Nspline-1])
        {
            Amp    = spline_interpolation(amp_interpolant, time[n]);
            Phase  = spline_interpolation(phase_interpolant, time[n]);
            
            wave->Y[n] = Amp*cos(Phase);
        }
    }

    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->Z);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->Z);

    for(n=0; n<N; n++)
    {
        if(time[n] < time_ssb[Nspline-1])
        {
            Amp    = spline_interpolation(amp_interpolant, time[n]);
            Phase  = spline_interpolation(phase_interpolant, time[n]);
            
            wave->Z[n] = Amp*cos(Phase);
        }
    }
     
     for(n=0; n<N; n++)
     {
         X[n] = wave->X[n];
         Y[n] = wave->Y[n];
         Z[n] = wave->Z[n];
     }

     */
    /*
    out = fopen("PhenomT_wave.dat","w");
    for(n=0; n<N; n++)
        fprintf(out,"%.15e %.15e %.15e %.15e\n", time[n], wave->X[n],wave->Y[n], wave->Z[n]);
    fclose(out);
     
     free_double_vector(time);
     
     free_tdi(wave);
     
     free_cubic_spline(amp_interpolant);
     free_cubic_spline(phase_interpolant);

     */
    
    // clean your room
    free_double_vector(data);
    
    free_double_vector(time_ssb);
    free_double_vector(freq_ssb);
    free_double_vector(amp_ssb);
    free_double_vector(phase_ssb);
    
    free_cubic_spline(amp_ssb_spline);
    free_cubic_spline(phase_ssb_spline);
    
    free_double_vector(phase_sc);
    free_double_vector(time_sc);
    
    free_tdi(tdi_phase);
    free_tdi(tdi_amp);
    
    
}

static void get_mbh_frequency(double Tend, double *fnew, double *tf, double *Amp, double *Phase, double t, double fguess, double *params)
{
    int N=3;
    double *freq = double_vector(3);
    double *time = double_vector(3);
    double *phase= double_vector(3);
    double *amp  = double_vector(3);

    // unpack parameter vector
    double tc     = params[5];
    double Mtotal = exp(params[1])*TSUN;
    
    // scale of size for finite differences
    double epsilon = 1e-6/Mtotal;
    
    // check that the finite difference steps are for valid frequencies
    if(fguess < 1.0/Tend) fguess = 1.0/Tend;
    if(fguess-epsilon < 0.0) epsilon = 0.5*fguess;

    // actual step sizes for finite differences
    double v = 2*PI2*epsilon;
    double u = PI2*epsilon*epsilon;
    
    // set up small frequency grid for finite differencing
    freq[0] = fguess-epsilon;
    freq[1] = fguess;
    freq[2] = fguess+epsilon;
    
    // get waveform on that small grid
    mbh_barycenter_waveform(params, N, freq, time, phase, amp, "IMRPhenomD");

    // finite differencing to get the things we need
    double tnew    = (phase[2]-phase[0])/v + tc;
    double dtdf    = (phase[2]+phase[0]-2.0*phase[1])/u;
    double delta_t = t - tnew;
    double delta_f = delta_t/dtdf;
    
    // pass back the things we needed
    *tf    = tnew;
    *fnew  = fguess+delta_f;
    *Amp   = amp[1];
    *Phase = phase[1];
    
    free_double_vector(freq);
    free_double_vector(time);
    free_double_vector(phase);
    free_double_vector(amp);
}

static void mbh_frequency_bandwidth(double *params, double tstart, double tstop, double *fstart, double *fstop)
{

    // Here we space the frequency array to give approximately equal spacing in time
    // The dynamic frequency spacing is capped to be between 1e-6 and 1e-4 Hz
    double Amp, Phase;
    double fnew, tf;
    int i;

    // unpack parameter array
    double Mchirp = exp(params[0]);
    double tc     = params[5];
        
    // ringodwn frequency
    double f_ringdown = mbh_ringdown_frequency(params);

    // Nyquist frequency
    double f_nyquist  = 1./(2*LISA_CADENCE);

    // guess at fmin (where the signal starts at t=tstart)
    double Tseg = tstop - tstart;
    double fmin = post_newtonian_frequency(Mchirp,tc,tstart);
    if(fmin < 1.0/Tseg) fmin = 1.0/Tseg;
    
    // Iteratively determine starting frequency
    i = 0;
    do
    {
        get_mbh_frequency(Tseg, &fnew, &tf, &Amp, &Phase, tstart, fmin, params);
        if(fnew < 1.0/Tseg) fnew = 1.0/Tseg;
        fmin = fnew;

        i++;
    }while(i < 10 && fabs(tf-tstart) > 1.0 && fmin == fmin);

    // nan catcher
    if(fmin != fmin) fmin = 1.0/Tseg;
    if(fmin < 0.0)   fmin = 1.0/Tseg;
    
    // guess at fmax (twice the ringdown frequency)
    double fmax = 2.0*f_ringdown;
    
    /*
     if merger time is after the end of the observation time
     iteratively determine max frequency
     */
    if(tc > tstop)
    {
        // guess at fmax (where the signal stops at t=tstop
        fmax = post_newtonian_frequency(Mchirp,tc,tstop);
        if(fmax < fmin) fmax = fmin+1.0/Tseg;
        
        // Iteratively determine stopping frequency
        i = 0;
        do
        {
            get_mbh_frequency(Tseg, &fnew, &tf, &Amp, &Phase, tstop, fmax, params);
            if(fnew < fmin) fnew = fmin+1.0/Tseg;
            fmax = fnew;

            i++;
        }while(i < 10 && fabs(tf-tstop) > 1.0);
    }
    
    // nan catcher
    if(fmax != fmax)     fmax = 2.0*f_ringdown; //something went wrong, set it to ringdown
    if(fmax > f_nyquist) fmax = f_nyquist; //merges out of band, set it to nyquist
    if(fmax < fmin)      fmax = 2.0*fmin; // something went wrong, make it bigger than fmin
    
    *fstart = fmin;
    *fstop  = fmax;
}

static double * mbh_frequency_grid(double Tobs, double *params, int *Ngrid)
{
    // This subroutine sets up the frequency sample array,
    // the time-frequency map and finds the PhenomD amplitude and phase
    // Nmax is the size of the holder arrays. N is the actual size.
    
    //allocate workspace
    int Nmin  = 4;
    int Nmax  = 100000;
    int Nwork = 0;
    double *work = double_vector(Nmax);
    
    // unpack params vector
    double Mchirp = exp(params[0]);
    double Mtotal = exp(params[1]);
    double tc     = params[5];
    
    // pad the start so we have values to interpolate allowing for time delays
    double deltaT = 1.e5; //TODO: what is this 1e5?
    double tstop  = tc + 1.0e4+deltaT; //TODO: is that 1e4 really doing anything?
    double tstart = -1.0e4-deltaT; //TODO: what even is the 1e4?
    double fmin, fmax;
    mbh_frequency_bandwidth(params, tstart, tstop, &fmin, &fmax);
    
    // this can happen when tc is really small and the masses are large
    if(fmax < fmin) fmin = 0.5*fmax;
    
    // set max and min frequency step sizes for grid
    double dfmin = 1.0/Tobs;
    double dfmax = fmax/100.0; //TODO: 100?
    
    // default grid scale TODO: DT?
    double fac = deltaT * pow(8.0*M_PI, 8.0/3.0) * 3.0/40.0*pow(Mchirp*TSUN,5.0/3.0);
    
    // get adaptive grid spacing based on chirp rate and time to merger
    double f,df,t;
    work[0] = fmin;
    Nwork = 1;
    do
    {
        // frequency
        f = work[Nwork-1];
        
        // coarse sampling for TDI phase
        df = fac*pow(f,11.0/3.0);
        
        /*
         decrease step size near merger needed for
         interpolating the PhenomD phase and time.
         */
        t = tc - post_newtonian_time(Mchirp, Mtotal, tc, f);
        if(t < 2.0e6) df/=10; //TODO: what is this 2e6?
        
        // keep df in range
        if(df < dfmin) df = dfmin;
        if(df > dfmax) df = dfmax;
        
        // increment frequency
        f += df;
        
        // store new grid point
        work[Nwork] = f;
        Nwork++;
        
    }while(f < fmax);
    
    // minimum size of arrays for interpolation
    if(Nwork < Nmin)
    {
        Nwork = Nmin;
        df = (fmax-fmin)/(Nmin-1);
        for(int i=0; i< Nwork; i++) work[i] = fmin + i*df;
    }
    
    
    //allocate output array
    *Ngrid = Nwork;
    double *freq_grid = double_vector(*Ngrid);
    memcpy(freq_grid, work, *Ngrid*sizeof(double));
    
    //clean up work space and return frequency grid
    free_double_vector(work);
    return freq_grid;
}

void mbh_fd_waveform(struct Orbit *orbit, double Tobs, double t0, double *params, double *X, double *Y, double *Z)
{
    
    // get frequency grid
    int Nspline; //number of grid points in frequency
    double *freq_grid = mbh_frequency_grid(Tobs, params, &Nspline);    
    
    //get time, phase, amplitude on the same grid
    double *time_ssb  = double_vector(Nspline);
    double *amp_ssb   = double_vector(Nspline);
    double *phase_ssb = double_vector(Nspline);

    mbh_barycenter_waveform(params, Nspline, freq_grid, time_ssb, phase_ssb, amp_ssb, "IMRPhenomD");
    
    /*
     Get spline interpolant for frequency and amplitude on time grid
     */
    struct CubicSpline *amp_ssb_spline   = alloc_cubic_spline(Nspline);
    struct CubicSpline *freq_ssb_spline  = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_ssb_spline,time_ssb,amp_ssb);
    initialize_cubic_spline(freq_ssb_spline,time_ssb,freq_grid);
    
    /*
     get reference phase
     */
    // reference is just exp(2 pi i f t)
    double *phase_ref = double_vector(Nspline);
    for(int i=0; i<Nspline; i++) phase_ref[i] = PI2 * freq_grid[i] * time_ssb[i];

    /*
     Get TDI reponse for signal's phase and amplitude on time grid
     set by the frequency evolution
     */
    struct TDI *tdi_phase = malloc(sizeof(struct TDI));
    struct TDI *tdi_amp = malloc(sizeof(struct TDI));
    alloc_tdi(tdi_phase,Nspline,3);
    alloc_tdi(tdi_amp,Nspline,3);

    //extract extrinsic parameters from MBH parameter vector
    double costh = sin(params[7]);
    double phi   = params[8];
    double cosi  = params[10];
    double psi   = params[9];

    LISA_spline_response(orbit, time_ssb, Nspline, costh, phi, cosi, psi, amp_ssb_spline, freq_ssb_spline, NULL, phase_ref, tdi_amp, tdi_phase);

    // shift phase back while rectifying sign conventions w/ IMRPhenomD
    for(int i=0; i<Nspline; i++)
    {
        tdi_phase->X[i] = phase_ssb[i] - tdi_phase->X[i];
        tdi_phase->Y[i] = phase_ssb[i] - tdi_phase->Y[i];
        tdi_phase->Z[i] = phase_ssb[i] - tdi_phase->Z[i];
    }
    
    /*
     Interpolate amplitude and phase for instrument response of each TDI channel onto frequency grid
     

    int N = (int)(Tobs/LISA_CADENCE);
    double tc = params[5];    // merger time
    double delta_t = Tobs + LISA_CADENCE - tc;
    
    double Amp,Phase;
    struct TDI *wave = malloc(sizeof(struct TDI));
    alloc_tdi(wave,N,3);
    
    struct CubicSpline *amp_interpolant   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_interpolant = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_interpolant,   freq_grid, tdi_amp->X);
    initialize_cubic_spline(phase_interpolant, freq_grid, tdi_phase->X);
    
    for(int i=0; i<N/2; i++)
    {
        double f = i/Tobs;
        wave->X[2*i] = wave->X[2*i+1] = 0.0;
        
        if(f>freq_grid[0] && f<freq_grid[Nspline-1])
        {
            Amp   = spline_interpolation(amp_interpolant,f);   //amplitude
            Phase = spline_interpolation(phase_interpolant,f); //slow part of phase
            Phase = PI2 * f * delta_t - Phase;
            
            wave->X[2*i]   = Amp * cos(Phase);
            wave->X[2*i+1] = Amp * sin(Phase);
        }
    }
    
    initialize_cubic_spline(amp_interpolant,   freq_grid, tdi_amp->Y);
    initialize_cubic_spline(phase_interpolant, freq_grid, tdi_phase->Y);
    
    for(int i=0; i<N/2; i++)
    {
        double f = i/Tobs;
        wave->Y[2*i] = wave->Y[2*i+1] = 0.0;
        
        if(f>freq_grid[0] && f<freq_grid[Nspline-1])
        {
            Amp   = spline_interpolation(amp_interpolant,f);   //amplitude
            Phase = spline_interpolation(phase_interpolant,f); //slow part of phase
            Phase = PI2 * f * delta_t - Phase;
            
            wave->Y[2*i]   = Amp * cos(Phase);
            wave->Y[2*i+1] = Amp * sin(Phase);
        }
    }
    
    initialize_cubic_spline(amp_interpolant,   freq_grid, tdi_amp->Z);
    initialize_cubic_spline(phase_interpolant, freq_grid, tdi_phase->Z);
    
    for(int i=0; i<N/2; i++)
    {
        double f = i/Tobs;
        wave->Z[2*i] = wave->Z[2*i+1] = 0.0;
        
        if(f>freq_grid[0] && f<freq_grid[Nspline-1])
        {
            Amp   = spline_interpolation(amp_interpolant,f);   //amplitude
            Phase = spline_interpolation(phase_interpolant,f); //slow part of phase
            Phase = PI2 * f * delta_t - Phase;
            
            wave->Z[2*i]   = Amp * cos(Phase);
            wave->Z[2*i+1] = Amp * sin(Phase);
        }
    }

    for(int i=0; i<N; i++)
    {
        X[i] = wave->X[i];
        Y[i] = wave->Y[i];
        Z[i] = wave->Z[i];
    }
     
     free_tdi(wave);
     
     free_cubic_spline(amp_interpolant);
     free_cubic_spline(phase_interpolant);

     */
    
    free_double_vector(freq_grid);
    free_double_vector(time_ssb);
    free_double_vector(phase_ssb);
    free_double_vector(amp_ssb);
    free_double_vector(phase_ref);
    
    free_cubic_spline(amp_ssb_spline);
    free_cubic_spline(freq_ssb_spline);
    
    free_tdi(tdi_phase);
    free_tdi(tdi_amp);
}
