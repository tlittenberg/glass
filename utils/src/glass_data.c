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

#include "glass_utils.h"
#include "gitversion.h"

#define FIXME 0

#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 60


void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    fprintf(stdout, "\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}


void print_version(FILE *fptr)
{
    fprintf(fptr, "\n");
    fprintf(fptr, "=============== GLASS Version: ==============\n\n");
    //fprintf(fptr, "  Git remote origin: %s\n", GIT_URL);
    //fprintf(fptr, "  Git version: %s\n", GIT_VER);
    fprintf(fptr, "  Git commit: %s\n", GITVERSION);
    //fprintf(fptr, "  Git commit author: %s\n",GIT_AUTHOR);
    //fprintf(fptr, "  Git commit date: %s\n", GIT_DATE);
    fprintf(fptr, "\n=============================================\n\n");
}

void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain)
{
    
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);

    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

}

void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    /* Load spacecraft ephemerides */
    switch(flags->orbit)
    {
        case 0:
            initialize_analytic_orbit(orbit);
            break;
        case 1:
            initialize_numeric_orbit(orbit);
            break;
        default:
            fprintf(stderr,"unsupported orbit type\n");
            exit(1);
            break;
    }
    
    /* set approximate f/fstar for segment */
    data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);
}

//TODO: Move file pointers to Model instead of Chain structs
void initialize_chain(struct Chain *chain, struct Flags *flags, unsigned int *seed, const char *mode)
{
    int ic;
    int NC = chain->NC;
    char filename[MAXSTRINGSIZE];

    chain->index = calloc(NC,sizeof(int));
    chain->acceptance = calloc(NC,sizeof(double));
    chain->temperature = calloc(NC,sizeof(double));
    chain->avgLogL     = calloc(NC,sizeof(double));
    chain->dimension   = calloc(NC,sizeof(int *));
    for(ic=0; ic<NC; ic++)
    {
        chain->index[ic]=ic;
        chain->acceptance[ic] = 1.0;
        chain->temperature[ic] = pow(1.2,(double)ic);
        chain->avgLogL[ic] = 0.0;
        chain->dimension[ic] = calloc(flags->DMAX,sizeof(int));
        for(int id=0; id<flags->DMAX; id++) chain->dimension[ic][id] = 0;
    }
    //set hottest chain to ~infinite temperature
    if(NC>1) chain->temperature[NC-1] = 1e12;
    chain->logLmax = 0.0;
    
    chain->r = malloc(NC*sizeof(unsigned int *));
    
    for(ic=0; ic<NC; ic++)
    {
        //set seed
        chain->r[ic] = *seed;
        
        //evolve seed
        rand_r_U_0_1(seed);
    }
    
    if(!flags->quiet)
    {
        sprintf(filename,"%s/log_likelihood_chain.dat",chain->chainDir);
        chain->likelihoodFile = fopen(filename,mode);
        
        sprintf(filename,"%s/temperature_chain.dat",chain->chainDir);
        chain->temperatureFile = fopen(filename,mode);
    }
    
    chain->chainFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/model_chain.dat.0",chain->chainDir);
    chain->chainFile[0] = fopen(filename,mode);
    
    chain->parameterFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/parameter_chain.dat.0",chain->chainDir);
    chain->parameterFile[0] = fopen(filename,mode);
    
    chain->dimensionFile = malloc(flags->DMAX*sizeof(FILE *));
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        chain->dimensionFile[i]=NULL;
    }
    
    chain->noiseFile = malloc(NC*sizeof(FILE *));
    sprintf(filename,"%s/noise_chain.dat.0",chain->chainDir);
    chain->noiseFile[0] = fopen(filename,mode);
    
    if(flags->confNoise)
    {
        chain->foregroundFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/foreground_chain.dat.0",chain->chainDir);
        chain->foregroundFile[0] = fopen(filename,mode);
    }
    
    if(flags->calibration)
    {
        chain->calibrationFile = malloc(NC*sizeof(FILE *));
        sprintf(filename,"%s/calibration_chain.dat.0",chain->chainDir);
        chain->calibrationFile[0] = fopen(filename,mode);
    }
    
    if(flags->verbose)
    {
        for(ic=1; ic<NC; ic++)
        {
            sprintf(filename,"%s/parameter_chain.dat.%i",chain->chainDir,ic);
            chain->parameterFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/model_chain.dat.%i",chain->chainDir,ic);
            chain->chainFile[ic] = fopen(filename,mode);
            
            sprintf(filename,"%s/noise_chain.dat.%i",chain->chainDir,ic);
            chain->noiseFile[ic] = fopen(filename,mode);
        }
    }
}

void alloc_data(struct Data *data, struct Flags *flags)
{
    int NMCMC = flags->NMCMC;
        
    data->logN = log((double)(data->N*data->Nchannel));
    
    data->tdi   = malloc(sizeof(struct TDI));
    data->raw   = malloc(sizeof(struct TDI));
    data->dft   = malloc(sizeof(struct TDI));
    data->dwt   = malloc(sizeof(struct TDI));
    data->noise = malloc(sizeof(struct Noise));
            
    alloc_tdi(data->tdi, data->N, data->Nchannel);
    alloc_tdi(data->raw, data->N, data->Nchannel);
    alloc_tdi(data->dft, data->N, data->Nchannel);
    alloc_tdi(data->dwt, data->N, data->Nchannel);
    if(!strcmp(data->basis,"fourier")) alloc_noise(data->noise, data->NFFT, data->Nlayer, data->Nchannel);
    if(!strcmp(data->basis,"wavelet")) alloc_noise(data->noise, data->N, data->Nlayer, data->Nchannel);
    
    //reconstructed signal model
    int i_re,i_im;
    data->h_rec = malloc(data->N*sizeof(double **));
    data->h_res = malloc(data->N*sizeof(double **));
    if(!strcmp(data->basis,"fourier"))
    {
        data->r_pow = malloc(data->NFFT*sizeof(double **));
        data->h_pow = malloc(data->NFFT*sizeof(double **));
        data->S_pow = malloc(data->NFFT*sizeof(double **));
    }
    if(!strcmp(data->basis,"wavelet"))
    {
        data->r_pow = malloc(data->N*sizeof(double **));
        data->h_pow = malloc(data->N*sizeof(double **));
        data->S_pow = malloc(data->N*sizeof(double **));
    }    
    
    //number of waveform samples to save
    data->Nwave=100;
    
    //downsampling rate of post-burn-in samples
    data->downsample = NMCMC/data->Nwave;
    
    if(!strcmp(data->basis,"fourier"))
    {
        for(int i=0; i<data->NFFT; i++)
        {
            i_re = i*2;
            i_im = i_re+1;
            
            data->S_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->h_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->r_pow[i]    = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i_re] = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i_im] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i_re] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i_im] = malloc(data->Nchannel*sizeof(double *));
            for(int n=0; n<data->Nchannel; n++)
            {
                data->S_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->h_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->r_pow[i][n]    = calloc(data->Nwave,sizeof(double));
                data->h_rec[i_re][n] = calloc(data->Nwave,sizeof(double));
                data->h_rec[i_im][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i_re][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i_im][n] = calloc(data->Nwave,sizeof(double));
            }
        }
    }

    if(!strcmp(data->basis,"wavelet"))
    {
        for(int i=0; i<data->N; i++)
        {
            
            data->S_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->r_pow[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_rec[i] = malloc(data->Nchannel*sizeof(double *));
            data->h_res[i] = malloc(data->Nchannel*sizeof(double *));
            for(int n=0; n<data->Nchannel; n++)
            {
                data->S_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->r_pow[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_rec[i][n] = calloc(data->Nwave,sizeof(double));
                data->h_res[i][n] = calloc(data->Nwave,sizeof(double));
            }
        }
    }
    
    //Spectrum proposal
    data->p = calloc(data->N,sizeof(double));

    // Setup wavelet basis
    if(!strcmp(data->basis,"wavelet"))
    {
        data->wdm = malloc(sizeof(struct Wavelets));
        initialize_wavelet(data->wdm, data->T);
    }
}

void alloc_noise(struct Noise *noise, int N, int Nlayer, int Nchannel)
{
    noise->N = N;
    noise->Nlayer = Nlayer;
    noise->Nchannel = Nchannel;
    
    noise->eta = calloc(Nchannel*Nlayer,sizeof(double));

    noise->f = calloc(N,sizeof(double));

    noise->C    = malloc(Nchannel*sizeof(double **));
    noise->invC = malloc(Nchannel*sizeof(double **));
    
    for(int i=0; i<Nchannel; i++)
    {
        for(int j=0; j<Nlayer; j++) noise->eta[i*Nlayer+j] = 1.0;
        noise->C[i]    = malloc(Nchannel*sizeof(double *));
        noise->invC[i] = malloc(Nchannel*sizeof(double *));
        
        for(int j=0; j<Nchannel; j++)
        {
            noise->C[i][j]    = calloc(N,sizeof(double));
            noise->invC[i][j] = calloc(N,sizeof(double));
        }
    }

    noise->detC     = calloc(N,sizeof(double));
    noise->transfer = calloc(N,sizeof(double));
    
    int n;
    for(n=0; n<N; n++)
    {
        for(int i=0; i<Nchannel; i++) noise->C[i][i][n] = 1.0;
        for(int i=0; i<Nchannel; i++)
        {
            for(int j=i+1; i<Nchannel; i++)
            {
                noise->C[i][j][n] = 0.0;
                noise->C[j][i][n] = 0.0;
            }
        }
        noise->transfer[n] = 1.0;
    }
}

void alloc_calibration(struct Calibration *calibration)
{
    calibration->dampA = 0.0;
    calibration->dampE = 0.0;
    calibration->dampX = 0.0;
    calibration->dphiA = 0.0;
    calibration->dphiE = 0.0;
    calibration->dphiX = 0.0;
    calibration->real_dphiA = 1.0;
    calibration->real_dphiE = 1.0;
    calibration->real_dphiX = 1.0;
    calibration->imag_dphiA = 0.0;
    calibration->imag_dphiE = 0.0;
    calibration->imag_dphiX = 0.0;
}

//TODO: Expand copy_data() to include everything, and then replace where needed (NoiseWrapper.c, ...)
void copy_data(struct Data *origin, struct Data *copy)
{
    memcpy(copy->format, origin->format, sizeof(origin->format));
    memcpy(copy->basis, origin->basis, sizeof(origin->basis));
    memcpy(copy->fileName, origin->fileName, sizeof(origin->fileName));
    copy->T=origin->T;
    copy->sqT=origin->sqT;
    copy->N=origin->N;
    copy->NFFT=origin->NFFT;
    copy->Nlayer=origin->Nlayer;
    copy->Nchannel=origin->Nchannel;
    copy->qpad=origin->qpad;
    copy->cseed=origin->cseed;
    copy->nseed=origin->nseed;
    copy->iseed=origin->iseed;
    copy->t0=origin->t0;
    //TODO: need copy_wavelet
}

void copy_noise(struct Noise *origin, struct Noise *copy)
{
    copy->N = origin->N;
    copy->Nlayer = origin->Nlayer;
    copy->Nchannel = origin->Nchannel;
    memcpy(copy->eta,origin->eta,origin->Nchannel*origin->Nlayer*sizeof(double));

    memcpy(copy->f, origin->f, origin->N*sizeof(double));

    copy_Cij(origin->C, copy->C, origin->Nchannel, origin->N);
    copy_Cij(origin->invC, copy->invC, origin->Nchannel, origin->N);

    memcpy(copy->detC, origin->detC, origin->N*sizeof(double));
    memcpy(copy->transfer, origin->transfer, origin->N*sizeof(double));
}

void copy_Cij(double ***origin, double ***copy, int M, int N)
{
    for(int i=0; i<M; i++)
        for(int j=0; j<M; j++)
            memcpy(copy[i][j], origin[i][j], N*sizeof(double));
}

void copy_calibration(struct Calibration *origin, struct Calibration *copy)
{
    copy=origin;
    /*
    copy->dampA   = origin->dampA;
    copy->dampE   = origin->dampE;
    copy->dampX   = origin->dampX;
    copy->dphiA = origin->dphiA;
    copy->dphiE = origin->dphiE;
    copy->dphiX = origin->dphiX;
    copy->real_dphiA = origin->real_dphiA;
    copy->real_dphiE = origin->real_dphiE;
    copy->real_dphiX = origin->real_dphiX;
    copy->imag_dphiA = origin->imag_dphiA;
    copy->imag_dphiE = origin->imag_dphiE;
    copy->imag_dphiX = origin->imag_dphiX;
     */
}

void free_noise(struct Noise *noise)
{
    free(noise->eta);
    free(noise->f);
    for(int i=0; i<noise->Nchannel; i++)
    {
        for(int j=0; j<noise->Nchannel; j++)
        {
            free(noise->C[i][j]);
            free(noise->invC[i][j]);
        }
        free(noise->C[i]);
        free(noise->invC[i]);
    }
    free(noise->C);
    free(noise->invC);
    free(noise->detC);
    free(noise->transfer);
    free(noise);
}

void free_chain(struct Chain *chain, struct Flags *flags)
{
    free(chain->index);
    free(chain->acceptance);
    free(chain->temperature);
    free(chain->avgLogL);
    free(chain->dimension);
    free(chain->r);
    
    if(!flags->quiet)
    {
        fclose(chain->likelihoodFile);
        fclose(chain->temperatureFile);
    }

    fclose(chain->chainFile[0]);

    fclose(chain->parameterFile[0]);
    
    for(int i=0; i<flags->DMAX; i++)
    {
        /* only create these files when needed */
        if(chain->dimensionFile[i]!=NULL) fclose(chain->dimensionFile[i]);
    }
    free(chain->dimensionFile);

    fclose(chain->noiseFile[0]);
    
    if(flags->calibration)
    {
        fclose(chain->calibrationFile[0]);
        free(chain->calibrationFile);
    }

    
    if(flags->verbose)
    {
        for(int ic=1; ic<chain->NC; ic++)
        {
            fclose(chain->chainFile[ic]);
            fclose(chain->parameterFile[ic]);
            fclose(chain->noiseFile[ic]);
        }
    }
    free(chain->chainFile);
    free(chain->parameterFile);
    free(chain->noiseFile);
    
    free(chain);
}

void free_calibration(struct Calibration *calibration)
{
    free(calibration);
}

void ReadHDF5(struct Data *data, struct TDI *tdi, struct TDI *tdi_dwt, struct Flags *flags)
{
    /* LDASOFT-formatted structure for TDI data */
    struct TDI *tdi_td = malloc(sizeof(struct TDI));
        
    if(!strcmp(data->format,"frequency"))  LISA_Read_HDF5_LDC_RADLER_TDI(tdi_td, data->fileName);
    if(!strcmp(data->format,"sangria")) LISA_Read_HDF5_LDC_TDI(tdi_td, data->fileName, "/obs/tdi");
    
    
    /* Select time segment of full data set */
    double start_time = data->t0;
    double stop_time = start_time + data->T;
    double dt = tdi_td->delta;
    double Tobs = stop_time - start_time;
    int N = (int)floor(Tobs/dt);

    /* work space for selecting and transforming time series */
    double *X = malloc(N*sizeof(double));
    double *Y = malloc(N*sizeof(double));
    double *Z = malloc(N*sizeof(double));
    
    double *Xtime = malloc(N*sizeof(double));
    double *Ytime = malloc(N*sizeof(double));
    double *Ztime = malloc(N*sizeof(double));

    /* Allocate data->tdi structure for Fourier transform output */
    alloc_tdi(tdi, N, N_TDI_CHANNELS);
    alloc_tdi(tdi_dwt, N, N_TDI_CHANNELS);
    tdi->delta = 1./Tobs;

    /* Select requested time segment */
    int n_start = (int)floor(start_time/dt); // first sample of time segment
    
    for(int n=0; n<N; n++)
    {
        int m = n_start+n;
        Xtime[n] = tdi_td->X[m];
        Ytime[n] = tdi_td->Y[m];
        Ztime[n] = tdi_td->Z[m];
    }
    
    /* lets get rid of those black holes */
    if(flags->no_mbh)
    {
        struct TDI *tdi_td_mbhb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_mbhb, data->fileName, "/sky/mbhb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_mbhb->X[m];
            Ytime[n] -= tdi_td_mbhb->Y[m];
            Ztime[n] -= tdi_td_mbhb->Z[m];
        }
        free_tdi(tdi_td_mbhb);
    }
    
    /* lets get rid of the galaxy */
    if(flags->no_ucb)
    {
        struct TDI *tdi_td_dgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_dgb, data->fileName, "/sky/dgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_dgb->X[m];
            Ytime[n] -= tdi_td_dgb->Y[m];
            Ztime[n] -= tdi_td_dgb->Z[m];
        }
        free_tdi(tdi_td_dgb);
    
        struct TDI *tdi_td_igb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_igb, data->fileName, "/sky/igb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_igb->X[m];
            Ytime[n] -= tdi_td_igb->Y[m];
            Ztime[n] -= tdi_td_igb->Z[m];
        }
        free_tdi(tdi_td_igb);
    }

    /* lets get rid of the verification binaries */
    if(flags->no_vgb)
    {
        struct TDI *tdi_td_vgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_vgb, data->fileName, "/sky/vgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xtime[n] -= tdi_td_vgb->X[m];
            Ytime[n] -= tdi_td_vgb->Y[m];
            Ztime[n] -= tdi_td_vgb->Z[m];
        }
        free_tdi(tdi_td_vgb);
    }

    /* Detrend data */
    detrend(Xtime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    detrend(Ytime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    detrend(Ztime, N, (int)(FILTER_LENGTH/LISA_CADENCE));
    
    /* Tukey window time-domain TDI channels tdi_td */
    double alpha = (2.0*FILTER_LENGTH/Tobs);
    
    tukey(Xtime, alpha, N);
    tukey(Ytime, alpha, N);
    tukey(Ztime, alpha, N);
    
    /* Fourier transform time-domain TDI channels */
    for(int n=0; n<N; n++)
    {
        X[n] = Xtime[n];
        Y[n] = Ytime[n];
        Z[n] = Ztime[n];
    }

    glass_forward_real_fft(X,N);
    glass_forward_real_fft(Y,N);
    glass_forward_real_fft(Z,N);
    
    /* Normalize FD data */
    double rft_norm = sqrt(Tobs)/(double)N;
    
    /* Account for losses from windowing
    double tukey_norm = tukey_scale(alpha, N);
    rft_norm /= tukey_norm;
    */
    
    for(int n=0; n<N; n++)
    {
        tdi->X[n] = X[n]*rft_norm;
        tdi->Y[n] = Y[n]*rft_norm;
        tdi->Z[n] = Z[n]*rft_norm;
    }
    
    /* lets get rid of the high frequency ucbs */
    if(flags->no_ucb_hi)
    {
        double *Xgal = calloc(N,sizeof(double));
        double *Ygal = calloc(N,sizeof(double));
        double *Zgal = calloc(N,sizeof(double));

        struct TDI *tdi_td_dgb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_dgb, data->fileName, "/sky/dgb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xgal[n] += tdi_td_dgb->X[m];
            Ygal[n] += tdi_td_dgb->Y[m];
            Zgal[n] += tdi_td_dgb->Z[m];
        }
        free_tdi(tdi_td_dgb);
        
        struct TDI *tdi_td_igb = malloc(sizeof(struct TDI));
        LISA_Read_HDF5_LDC_TDI(tdi_td_igb, data->fileName, "/sky/igb/tdi");
        for(int n=0; n<N; n++)
        {
            int m = n_start+n;
            Xgal[n] += tdi_td_igb->X[m];
            Ygal[n] += tdi_td_igb->Y[m];
            Zgal[n] += tdi_td_igb->Z[m];
        }
        free_tdi(tdi_td_igb);
        
        tukey(Xgal, alpha, N);
        tukey(Ygal, alpha, N);
        tukey(Zgal, alpha, N);

        glass_forward_real_fft(Xgal, N);
        glass_forward_real_fft(Ygal, N);
        glass_forward_real_fft(Zgal, N);

        
        /* Allocate data->tdi structure for Fourier transform output */
        struct TDI *tdi_gal = malloc(sizeof(struct TDI));
        alloc_tdi(tdi_gal, N, N_TDI_CHANNELS);
        tdi_gal->delta = 1./Tobs;

        for(int n=0; n<N; n++)
        {
            tdi_gal->X[n] = Xgal[n] * rft_norm;
            tdi_gal->Y[n] = Ygal[n] * rft_norm;
            tdi_gal->Z[n] = Zgal[n] * rft_norm;
        }

        /* remove hi-f binaries */
        for(int n=0; n<N/2; n++)
        {
            double f = (double)n/Tobs;
            if(f>0.00504)
            {
                tdi->X[2*n]   -= tdi_gal->X[2*n];
                tdi->X[2*n+1] -= tdi_gal->X[2*n+1];

                tdi->Y[2*n]   -= tdi_gal->Y[2*n];
                tdi->Y[2*n+1] -= tdi_gal->Y[2*n+1];

                tdi->Z[2*n]   -= tdi_gal->Z[2*n];
                tdi->Z[2*n+1] -= tdi_gal->Z[2*n+1];

            }
        }
        
        free_tdi(tdi_gal);
        free(Xgal);
        free(Ygal);
        free(Zgal);
    }

    /* populate AET channels because I can't let go */
    for(int n=0; n<N; n++) XYZ2AET(tdi->X[n], tdi->Y[n], tdi->Z[n], &tdi->A[n], &tdi->E[n], &tdi->T[n]);


    /* Wavelet transform time-domain TDI channels */
    if(!strcmp(data->basis,"wavelet"))
    {
        for(int n=0; n<N; n++)
        {
            X[n] = Xtime[n];
            Y[n] = Ytime[n];
            Z[n] = Ztime[n];
        }

        wavelet_transform(data->wdm, X);
        wavelet_transform(data->wdm, Y);
        wavelet_transform(data->wdm, Z);
        
        for(int n=0; n<N; n++)
        {
            tdi_dwt->X[n] = X[n];
            tdi_dwt->Y[n] = Y[n];
            tdi_dwt->Z[n] = Z[n];
        }
        
        /* populate AET channels because I can't let go */
        for(int n=0; n<N; n++) XYZ2AET(tdi_dwt->X[n], tdi_dwt->Y[n], tdi_dwt->Z[n], &tdi_dwt->A[n], &tdi_dwt->E[n], &tdi_dwt->T[n]);
        
    }
    
    /* Free memory */
    free_tdi(tdi_td);
    free(X);
    free(Y);
    free(Z);
    free(Xtime);
    free(Ytime);
    free(Ztime);

}

void ReadASCII(struct Data *data, struct TDI *tdi)
{
    double f;
    double junk;
    
    FILE *fptr = fopen(data->fileName,"r");
    
    //count number of samples
    int Nsamples = 0;
    while(!feof(fptr))
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&junk,&junk,&junk,&junk);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
        Nsamples++;
    }
    rewind(fptr);
    Nsamples--;
    
    //load full dataset into TDI structure
    alloc_tdi(tdi, 2*Nsamples, 3);
    
    for(int n=0; n<Nsamples; n++)
    {
        int check = fscanf(fptr,"%lg %lg %lg %lg %lg",&f,&tdi->A[2*n],&tdi->A[2*n+1],&tdi->E[2*n],&tdi->E[2*n+1]);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",data->fileName);
            exit(1);
        }
        
    }
    fclose(fptr);
}

void ReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet) fprintf(stdout,"\n==== ReadData ====\n");
    
    /* load full dataset */
    struct TDI *tdi_full_dft = malloc(sizeof(struct TDI));
    struct TDI *tdi_full_dwt = malloc(sizeof(struct TDI));
    
    if(flags->hdf5Data)
        ReadHDF5(data,tdi_full_dft,tdi_full_dwt,flags);
    else
        ReadASCII(data,tdi_full_dft);
    
    
    /* select frequency segment */
    data->fmax = data->fmin + data->NFFT/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;
    
    //store frequency segment in TDI structure
    for(int n=0; n<data->N; n++)
    {
        int m = data->qmin*2+n;
        data->dft->X[n] = tdi_full_dft->X[m];
        data->dft->Y[n] = tdi_full_dft->Y[m];
        data->dft->Z[n] = tdi_full_dft->Z[m];
        data->dft->A[n] = tdi_full_dft->A[m];
        data->dft->E[n] = tdi_full_dft->E[m];
        data->dft->T[n] = tdi_full_dft->T[m];
    }
    
    /* select wavelet layers */
    if(!strcmp(data->basis,"wavelet"))
    {
        //set lmin and lmax to represent frequency layer instead of bin
        data->lmin = (int)floor(data->fmin/WAVELET_BANDWIDTH);
        data->lmax = data->lmin + data->Nlayer;

        printf("  Minimum frequency layer=%i, maximum layer=%i\n",data->lmin,data->lmax-1);
        printf("  fmin=%lg, fmax=%lg\n",data->fmin,data->fmax);

        //reset wavelet basis max and min ranges
        wavelet_pixel_to_index(data->wdm,0,data->lmin,&data->wdm->kmin);
        wavelet_pixel_to_index(data->wdm,0,data->lmax,&data->wdm->kmax);

        //store frequency segment in TDI structure
        for(int n=0; n<data->N; n++)
        {
            int m = data->wdm->kmin+n;
            data->dwt->X[n] = tdi_full_dwt->X[m];
            data->dwt->Y[n] = tdi_full_dwt->Y[m];
            data->dwt->Z[n] = tdi_full_dwt->Z[m];
            data->dwt->A[n] = tdi_full_dwt->A[m];
            data->dwt->E[n] = tdi_full_dwt->E[m];
            data->dwt->T[n] = tdi_full_dwt->T[m];
        }
    }
    
    /* copy correct representation of data into the main tdi structure */
    for(int n=0; n<data->N; n++)
    {
        if(!strcmp(data->basis,"fourier"))
        {
            data->tdi->X[n] = data->dft->X[n];
            data->tdi->Y[n] = data->dft->Y[n];
            data->tdi->Z[n] = data->dft->Z[n];
            data->tdi->A[n] = data->dft->A[n];
            data->tdi->E[n] = data->dft->E[n];
            data->tdi->T[n] = data->dft->T[n];
        }
        else if(!strcmp(data->basis,"wavelet"))
        {
            data->tdi->X[n] = data->dwt->X[n];
            data->tdi->Y[n] = data->dwt->Y[n];
            data->tdi->Z[n] = data->dwt->Z[n];
            data->tdi->A[n] = data->dwt->A[n];
            data->tdi->E[n] = data->dwt->E[n];
            data->tdi->T[n] = data->dwt->T[n];
        }
    }
     
    //free memory
    free_tdi(tdi_full_dft);
    free_tdi(tdi_full_dwt);
}

void GetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    double Spm, Sop;
    
    //if you are simulating/fitting the noise
    if(!flags->psd)
    {
        for(int n=0; n<data->NFFT; n++)
        {
            double f = data->fmin + (double)(n)/data->T;
            data->noise->f[n] = f;
            data->noise->transfer[n] = noise_transfer_function(f/orbit->fstar);

            if(strcmp(data->format,"phase")==0)
            {
                data->noise->C[0][0][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[1][1][n] = AEnoise(orbit->L, orbit->fstar, f);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise(data->T,f);
                    data->noise->C[1][1][n] += GBnoise(data->T,f);
                }
            }
            else if(strcmp(data->format,"frequency")==0)
            {
                get_noise_levels("radler",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else if(strcmp(data->format,"sangria")==0)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][1][n] = 0.0;
                if(flags->confNoise)
                {
                    data->noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f);
                    data->noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f);
                }
            }
            else
            {
                fprintf(stderr,"Unsupported data format %s\n",data->format);
                exit(1);
            }

            //TODO: 3-channel model only has support for Sangria data conventions
            if(data->Nchannel==3)
            {
                //TODO: Need a sqrt(2) to match Sangria data/noise
                get_noise_levels("sangria",f,&Spm,&Sop);
                data->noise->C[0][0][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][1][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[2][2][n] = XYZnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                data->noise->C[0][1][n] = data->noise->C[1][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[0][2][n] = data->noise->C[2][0][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);
                data->noise->C[1][2][n] = data->noise->C[2][1][n] = XYZcross_FF(orbit->L, orbit->fstar, f, Spm, Sop);

                if(flags->confNoise)
                {
                    double GBnoise=GBnoise_FF(data->T, orbit->fstar, f)/1.5; //GBnoise_FF() is hard-coded for AE channels
                    data->noise->C[0][0][n] += GBnoise;
                    data->noise->C[1][1][n] += GBnoise;
                    data->noise->C[2][2][n] += GBnoise;
                    data->noise->C[0][1][n] += -0.5*GBnoise;
                    data->noise->C[0][2][n] += -0.5*GBnoise;
                    data->noise->C[1][2][n] += -0.5*GBnoise;
                    data->noise->C[1][0][n] += -0.5*GBnoise;
                    data->noise->C[2][0][n] += -0.5*GBnoise;
                    data->noise->C[2][1][n] += -0.5*GBnoise;
                }
                
                /*normalize*/
                data->noise->C[0][0][n] /= 4.;
                data->noise->C[0][1][n] /= 4.;
                data->noise->C[0][2][n] /= 4.;
                data->noise->C[1][0][n] /= 4.;
                data->noise->C[1][1][n] /= 4.;
                data->noise->C[1][2][n] /= 4.;
                data->noise->C[2][0][n] /= 4.;
                data->noise->C[2][1][n] /= 4.;
                data->noise->C[2][2][n] /= 4.;
            }
        }
        
        invert_noise_covariance_matrix(data->noise);

    }
    //use PSD from file
    else
    {
        
        //parse input PSD file
        FILE *psdFile = fopen(flags->psdFile,"r");
        int lines=0;
        double f_temp, SnA_temp, SnE_temp;
        while(!feof(psdFile))
        {
            fscanf(psdFile,"%lg %lg %lg",&f_temp,&SnA_temp,&SnE_temp);
            lines++;
        }
        rewind(psdFile);
        lines--;
        
        double *f   = malloc(lines*sizeof(double));
        double *SnA = malloc(lines*sizeof(double));
        double *SnE = malloc(lines*sizeof(double));
        
        for(int l=0; l<lines; l++) fscanf(psdFile,"%lg %lg %lg",&f[l],&SnA[l],&SnE[l]);
        
        //interpolate input psd onto segment grid
        double *fint = malloc(data->NFFT*sizeof(double));
        for(int n=0; n<data->NFFT; n++) fint[n] = data->fmin + (double)(n)/data->T;

        CubicSplineGLASS(lines, f, SnA, data->NFFT, fint, data->noise->C[0][0]);
        CubicSplineGLASS(lines, f, SnE, data->NFFT, fint, data->noise->C[1][1]);
        
        free(f);
        free(SnA);
        free(SnE);
        fclose(psdFile);
    }
}

void AddNoise(struct Data *data, struct TDI *tdi)
{
    
    printf("   ...adding Gaussian noise realization\n");
    
    //set RNG for noise
    unsigned int r = data->nseed;
    
    double n_re[data->Nchannel];
    double n_im[data->Nchannel];
    double u_re[data->Nchannel];
    double u_im[data->Nchannel];
    
    //get LU decomposition of covariance matrix
    double **L = malloc(data->Nchannel*sizeof(double*));
    double **C = malloc(data->Nchannel*sizeof(double*));
    for(int i=0; i<data->Nchannel; i++)
    {
        L[i] = malloc(data->Nchannel*sizeof(double));
        C[i] = malloc(data->Nchannel*sizeof(double));
    }
    
    
    
    for(int n=0; n<data->NFFT; n++)
    {
        for(int i=0; i<data->Nchannel; i++)
        {
            u_re[i] = rand_r_N_0_1(&r);
            u_im[i] = rand_r_N_0_1(&r);
            n_re[i] = n_im[i] = 0.0;
        }
 
        // make sure both diagonals of the covariance matrix are filled
        for(int i=0; i<data->Nchannel; i++)
            for(int j=i; j<data->Nchannel; j++)
                C[i][j] = C[j][i] = data->noise->C[i][j][n];

        cholesky_decomp(C, L, data->Nchannel);

        // n = Lu
        for(int i=0; i<data->Nchannel; i++)
        {
            for(int j=0; j<data->Nchannel; j++)
            {
                n_re[i] += L[i][j]*u_re[j]/sqrt(2.);
                n_im[i] += L[i][j]*u_im[j]/sqrt(2.);
            }
        }
        
        switch(data->Nchannel)
        {
            case 1:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                break;
            case 2:
                tdi->A[2*n]   += n_re[0];
                tdi->A[2*n+1] += n_im[0];
                tdi->E[2*n]   += n_re[1];
                tdi->E[2*n+1] += n_im[1];
                break;
            case 3:
                tdi->X[2*n]   += n_re[0];
                tdi->X[2*n+1] += n_im[0];
                tdi->Y[2*n]   += n_re[1];
                tdi->Y[2*n+1] += n_im[1];
                tdi->Z[2*n]   += n_re[2];
                tdi->Z[2*n+1] += n_im[2];
                break;
        }
    }

    for(int i=0; i<data->Nchannel; i++)
    {
        free(L[i]);
        free(C[i]);
    }
    free(L);
    free(C);
}

void AddNoiseWavelet(struct Data *data, struct TDI *tdi)
{
    
    printf("   ...adding Gaussian noise realization\n");
    
    //set RNG for noise
    unsigned int r = data->nseed;
    
    double n[data->Nchannel];
    double u[data->Nchannel];
    
    //get LU decomposition of covariance matrix
    double **L = malloc(data->Nchannel*sizeof(double*));
    double **C = malloc(data->Nchannel*sizeof(double*));
    for(int i=0; i<data->Nchannel; i++)
    {
        L[i] = malloc(data->Nchannel*sizeof(double));
        C[i] = malloc(data->Nchannel*sizeof(double));
    }
    
    
    int k;
    struct Wavelets *wdm = data->wdm;
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=data->lmin; j<data->lmax; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            k-=data->wdm->kmin;

            for(int a=0; a<data->Nchannel; a++)
            {
                u[a] = rand_r_N_0_1(&r);
                n[a] = 0.0;
            }
    
            // make sure both diagonals of the covariance matrix are filled
            for(int a=0; a<data->Nchannel; a++)
                for(int b=a; b<data->Nchannel; b++)
                    C[a][b] = C[b][a] = data->noise->C[a][b][k];


            cholesky_decomp(C, L, data->Nchannel);

            // n = Lu
            for(int a=0; a<data->Nchannel; a++)
            {
                for(int b=0; b<data->Nchannel; b++)
                {
                    n[a] += L[a][b]*u[b];
                }
            }
            
            switch(data->Nchannel)
            {
                case 1:
                    tdi->X[k] += n[0];
                    break;
                case 2:
                    tdi->A[k] += n[0];
                    tdi->E[k] += n[1];
                    break;
                case 3:
                    tdi->X[k] += n[0];
                    tdi->Y[k] += n[1];
                    tdi->Z[k] += n[2];
                    break;
            }
        }
    }

    for(int i=0; i<data->Nchannel; i++)
    {
        free(L[i]);
        free(C[i]);
    }
    free(L);
    free(C);
}

void SimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags)
{
    if(!flags->quiet) fprintf(stdout,"\n==== SimulateData ====\n");
    struct TDI *tdi = data->tdi;

    //get max and min samples
    data->fmax = data->fmin + data->NFFT/data->T;
    data->qmin = (int)(data->fmin*data->T);
    data->qmax = data->qmin+data->NFFT;

    //Get noise spectrum for data segment
    GetNoiseModel(data,orbit,flags);
    
    //Add Gaussian noise to injection
    if(flags->simNoise) AddNoise(data,tdi);
    
    //print various data products for plotting
    print_data(data, flags);

}

void print_data(struct Data *data, struct Flags *flags)
{
    int k;
    FILE *fptr;
    char filename[256];
    struct TDI *tdi = NULL;
    
    
    /* Power spectra */
    tdi = data->dft;
    sprintf(filename,"%s/power_data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->A[2*i]*tdi->A[2*i]+tdi->A[2*i+1]*tdi->A[2*i+1], tdi->E[2*i]*tdi->E[2*i]+tdi->E[2*i+1]*tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg\n", f, tdi->X[2*i]*tdi->X[2*i]+tdi->X[2*i+1]*tdi->X[2*i+1], tdi->Y[2*i]*tdi->Y[2*i]+tdi->Y[2*i+1]*tdi->Y[2*i+1], tdi->Z[2*i]*tdi->Z[2*i]+tdi->Z[2*i+1]*tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);

    /* Scaleogram */
    if(!strcmp("wavelet",data->basis))
    {
        tdi = data->dwt;
        sprintf(filename,"%s/scaleogram_data.dat",data->dataDir);
        fptr=fopen(filename,"w");

        for(int j=data->lmin; j<data->lmax; j++)
        {
            for(int i=0; i<data->wdm->NT; i++)
            {
                wavelet_pixel_to_index(data->wdm,i,j,&k);
                k-=data->wdm->kmin;
                fprintf(fptr,"%lg %lg %.14e %.14e %.14e\n", i*data->wdm->dt, j*data->wdm->df + WAVELET_BANDWIDTH/2,tdi->X[k]*tdi->X[k], tdi->Y[k]*tdi->Y[k], tdi->Z[k]*tdi->Z[k]);
            }
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }

    tdi = data->dft;
    sprintf(filename,"%s/dft_data.dat",data->dataDir);
    fptr=fopen(filename,"w");
    for(int i=0; i<data->NFFT; i++)
    {
        double f = (double)(i+data->qmin)/data->T;
        switch(data->Nchannel)
        {
            case 1:
                fprintf(fptr,"%.12g %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1]);
                break;
            case 2:
                fprintf(fptr,"%.12g %lg %lg %lg %lg\n", f, tdi->A[2*i],tdi->A[2*i+1], tdi->E[2*i],tdi->E[2*i+1]);
                break;
            case 3:
                fprintf(fptr,"%.12g %lg %lg %lg %lg %lg %lg\n", f, tdi->X[2*i],tdi->X[2*i+1], tdi->Y[2*i],tdi->Y[2*i+1], tdi->Z[2*i],tdi->Z[2*i+1]);
                break;
        }
    }
    fclose(fptr);

    if(!strcmp("wavelet",data->basis))
    {
        tdi = data->dwt;
        sprintf(filename,"%s/dwt_data.dat",data->dataDir);
        fptr=fopen(filename,"w");
        for(int j=data->lmin; j<data->lmax; j++)
        {
            for(int i=0; i<data->wdm->NT; i++)
            {
                wavelet_pixel_to_index(data->wdm,i,j,&k);
                k-=data->wdm->kmin;
                fprintf(fptr,"%lg %lg %.14e %.14e %.14e\n", i*data->wdm->dt, j*data->wdm->df,tdi->X[k], tdi->Y[k], tdi->Z[k]);   
            }
            fprintf(fptr,"\n");
        }
        fclose(fptr);
    }
}

void print_wavelet_fourier_spectra(struct Data *data, struct TDI *tdi, char filename[])
{
    int k;
    struct Wavelets *wdm = data->wdm;
    int N = wdm->NF*wdm->NT;
    double T = N*LISA_CADENCE;

    double **freqData = double_matrix(3,N);
    double **waveData = double_matrix(3,N);

    //get TDI data into context
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                waveData[0][k] = tdi->X[k-wdm->kmin];
                waveData[1][k] = tdi->Y[k-wdm->kmin];
                waveData[2][k] = tdi->Z[k-wdm->kmin];
            }
        }
    }

    //wavelet to frequency
    memcpy(freqData[0],waveData[0],sizeof(double)*N);
    memcpy(freqData[1],waveData[1],sizeof(double)*N);
    memcpy(freqData[2],waveData[2],sizeof(double)*N);
    wavelet_tansform_inverse_fourier(wdm, freqData[0]);
    wavelet_tansform_inverse_fourier(wdm, freqData[1]);
    wavelet_tansform_inverse_fourier(wdm, freqData[2]);

    FILE *fptr=fopen(filename,"w");
    for(int n=0; n<N/2; n++) 
    {
        double f = (double)n/T;
        if(f>data->fmin && f < data->fmax)
        {
            fprintf(fptr,"%.14e %.14e %.14e %.14e %.14e %.14e %.14e\n", n/T, freqData[0][2*n],freqData[0][2*n+1],freqData[1][2*n],freqData[1][2*n+1], freqData[2][2*n],freqData[2][2*n+1]);
        }
    }
    fclose(fptr);

    free_double_matrix(freqData,3);
    free_double_matrix(waveData,3);
}

void wavelet_layer_to_fourier_transform(struct Data *data)
{
    int k;
    struct Wavelets *wdm = data->wdm;
    int N = wdm->NF*wdm->NT;

    double **freqData = double_matrix(3,N);
    double **waveData = double_matrix(3,N);

    /* store backup DWT data */
    memcpy(data->dwt->X,data->tdi->X,sizeof(double)*data->N);
    memcpy(data->dwt->Y,data->tdi->Y,sizeof(double)*data->N);
    memcpy(data->dwt->Z,data->tdi->Z,sizeof(double)*data->N);
    
    /* populate full DWT data with active layers */
    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=0; j<wdm->NF; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                waveData[0][k] = data->dwt->X[k-wdm->kmin];
                waveData[1][k] = data->dwt->Y[k-wdm->kmin];
                waveData[2][k] = data->dwt->Z[k-wdm->kmin];
            }
        }
    }

    /* copy full DWT into full DFT array for in-place transform */
    memcpy(freqData[0],waveData[0],sizeof(double)*N);
    memcpy(freqData[1],waveData[1],sizeof(double)*N);
    memcpy(freqData[2],waveData[2],sizeof(double)*N);
    
    /* in place DWT->DFT transform */
    wavelet_tansform_inverse_fourier(wdm, freqData[0]);
    wavelet_tansform_inverse_fourier(wdm, freqData[1]);
    wavelet_tansform_inverse_fourier(wdm, freqData[2]);

    /* copy active frequency bins into DFT struct */
    for(int n=0; n<data->NFFT; n++)
    {
        int m = n + data->qmin;
        data->dft->X[2*n]   = freqData[0][2*m];
        data->dft->X[2*n+1] = freqData[0][2*m+1];
        data->dft->Y[2*n]   = freqData[1][2*m];
        data->dft->Y[2*n+1] = freqData[1][2*m+1];
        data->dft->Z[2*n]   = freqData[2][2*m];
        data->dft->Z[2*n+1] = freqData[2][2*m+1];
    }
    
    /* Get A&E channels just in case... */
    for(int n=0; n<data->N; n++) XYZ2AE(data->dft->X[n],data->dft->Y[n],data->dft->Z[n],&data->dft->A[n],&data->dft->E[n]);

    free_double_matrix(freqData,3);
    free_double_matrix(waveData,3);
}

void print_glass_usage()
{
    fprintf(stdout,"\n");
    fprintf(stdout,"=============== GLASS Usage: ============== \n");
    fprintf(stdout,"REQUIRED:\n");
    fprintf(stdout,"\n");
    fprintf(stdout,"OPTIONAL:\n");
    fprintf(stdout,"  -h | --help        : print help message and exit         \n");
    fprintf(stdout,"  -v | --verbose     : enable verbose output               \n");
    fprintf(stdout,"  -q | --quiet       : restrict output                     \n");
    fprintf(stdout,"  -d | --debug       : leaner settings for quick running   \n");
    fprintf(stdout,"\n");
    
    //LISA
    fprintf(stdout,"       =========== LISA =========== \n");
    fprintf(stdout,"       --orbit       : orbit ephemerides file (2.5 GM MLDC)\n");
    fprintf(stdout,"       --channels    : # of channels [1->X,2->AE,3->XYZ](3)\n");
    fprintf(stdout,"       --phase       : phase data (fractional frequency)   \n");
    fprintf(stdout,"       --sangria     : use LDC Sangria TDI conventions     \n");
    fprintf(stdout,"\n");
    
    //Data
    fprintf(stdout,"       =========== Data =========== \n");
    fprintf(stdout,"       --data        : strain data file (ASCII)            \n");
    fprintf(stdout,"       --h5-data     : strain data file (HDF5)             \n");
    fprintf(stdout,"       --h5-no-mbh   : remove mbhs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-ucb   : remove ucbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-ucb-hi: remove high f ucbs from HDF5 data   \n");
    fprintf(stdout,"       --h5-no-vgb   : remove vgbs from HDF5 data          \n");
    fprintf(stdout,"       --h5-no-noise : remove noise from HDF5 data (TODO)  \n");
    fprintf(stdout,"       --psd         : psd data file (ASCII)               \n");
    fprintf(stdout,"       --samples     : number of DFT frequency bins (512)  \n");
    fprintf(stdout,"       --layers      : number of DWT frequency layers (1)  \n");
    fprintf(stdout,"       --padding     : number of bins padded on segment (0)\n");
    fprintf(stdout,"       --start-time  : initial time of epoch  (0)          \n");
    fprintf(stdout,"       --fmin        : minimum frequency                   \n");
    fprintf(stdout,"       --fmax        : maximum frequency (overrides --samples)\n");
    fprintf(stdout,"       --duration    : duration of epoch (31457280)        \n");
    fprintf(stdout,"       --sim-noise   : data w/out noise realization        \n");
    fprintf(stdout,"       --conf-noise  : include model for confusion noise   \n");
    fprintf(stdout,"       --noiseseed   : seed for noise RNG                  \n");
    fprintf(stdout,"\n");
    
    //Chain
    fprintf(stdout,"       ========== Chains ========== \n");
    fprintf(stdout,"       --steps       : number of mcmc steps (10000)        \n");
    fprintf(stdout,"       --chainseed   : seed for MCMC RNG                   \n");
    fprintf(stdout,"       --chains      : number of parallel chains (20)      \n");
    fprintf(stdout,"       --no-burnin   : skip burn in steps                  \n");
    fprintf(stdout,"       --resume      : restart from checkpoint             \n");
    fprintf(stdout,"       --threads     : number of parallel threads (max)    \n");
    fprintf(stdout,"       --prior       : sample from prior                   \n");
    fprintf(stdout,"       --no-rj       : turn off RJMCMC                     \n");
    fprintf(stdout,"\n");
    
    //Misc.
    fprintf(stdout,"       =========== Misc =========== \n");
    fprintf(stdout,"       --rundir      : top level run directory ['./']\n");
    fprintf(stdout,"       --match-in1   : input paramaters for overlap [filename] \n");
    fprintf(stdout,"       --match-in2   : output match values [filename] \n");
    fprintf(stdout,"\n");

    /*
    fprintf(stdout,"USEFUL TOBS:\n");
    fprintf(stdout,"   1 wk: %.0f\n",62914560./2./52.);
    fprintf(stdout,"   1 mo: %.0f \n",62914560./2./12.);
    fprintf(stdout,"   2 mo: %.0f \n",62914560./2./6.);
    fprintf(stdout,"   1 yr: %.0f (default)\n",62914560./2.);
    fprintf(stdout,"   2 yr: %.0f \n",62914560.);
    fprintf(stdout,"   5 yr: %.0f \n",5.*62914560./2.);
    fprintf(stdout,"  10 yr: %.0f \n",10.*62914560./2.);
     */
    fprintf(stdout,"\n");
}

void parse_data_args(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, char basis[])
{
    //copy argv since getopt permutes order
    char **argv_copy=malloc((argc+1) * sizeof *argv_copy);
    copy_argv(argc,argv,argv_copy);
    opterr=0; //suppress warnings about unknown arguments
    
    //Set defaults
    flags->rj          = 1;
    flags->help        = 0;
    flags->calibration = 0;
    flags->verbose     = 0;
    flags->quiet       = 0;
    flags->simNoise    = 0;
    flags->confNoise   = 0;
    flags->stationary  = 0;
    flags->burnin      = 1;
    flags->debug       = 0;
    flags->strainData  = 0;
    flags->hdf5Data    = 0;
    flags->psd         = 0;
    flags->orbit       = 0;
    flags->prior       = 0;
    flags->resume      = 0;
    flags->NMCMC       = 1000;
    flags->NBURN       = 1000;
    flags->threads     = omp_get_max_threads();
    sprintf(flags->runDir,"./");
    chain->NC          = 12;//number of chains
    int set_fmax_flag  = 0; //flag watching for if fmax is set by CLI
    
    /* Simulated data building blocks */
    flags->no_mbh = 0;
    flags->no_ucb = 0;
    flags->no_vgb = 0;
    flags->no_noise = 0;
    
    /*
     default data format is 'phase'
     optional support for 'frequency' a la LDCs
     */
    sprintf(data->format,"sangria");
    sprintf(data->basis,"%s",basis);

    data->T        = 31457280; /* one "mldc years" at 15s sampling */
    data->t0       = 0.0; /* start time of data segment in seconds */
    data->sqT      = sqrt(data->T);
    data->NFFT     = 512;
    data->Nlayer   = 1;
    data->Nchannel = 3; //1=X, 2=AE, 3=XYZ
    data->qpad     = 0;
    data->fmin     = 1e-4; //Hz
    
    data->cseed = 150914;
    data->nseed = 151226;
    data->iseed = 151012;

    if(!strcmp(data->basis,"fourier")) data->N = data->NFFT*2;
    if(!strcmp(data->basis,"wavelet"))
    {
        data->T = floor(data->T/WAVELET_DURATION)*WAVELET_DURATION;
        data->sqT = sqrt(data->T);
        data->Nlayer = 1;
        data->N = (int)floor(data->T/WAVELET_DURATION)*data->Nlayer;
        data->NFFT = data->N/2;
    }


    //Specifying the expected options
    static struct option long_options[] =
    {
        /* These options set a flag. */
        {"samples",    required_argument, 0, 0},
        {"layers",     required_argument, 0, 0},
        {"padding",    required_argument, 0, 0},
        {"duration",   required_argument, 0, 0},
        {"start-time", required_argument, 0, 0},
        {"orbit",      required_argument, 0, 0},
        {"chains",     required_argument, 0, 0},
        {"chainseed",  required_argument, 0, 0},
        {"noiseseed",  required_argument, 0, 0},
        {"data",       required_argument, 0, 0},
        {"h5-data",    required_argument, 0, 0},
        {"psd",        required_argument, 0, 0},
        {"fmin",       required_argument, 0, 0},
        {"fmax",       required_argument, 0, 0},
        {"channels",   required_argument, 0, 0},
        {"steps",      required_argument, 0, 0},
        {"threads",    required_argument, 0, 0},
        {"rundir",     required_argument, 0, 0},
        
        /* These options don’t set a flag.
         We distinguish them by their indices. */
        {"help",        no_argument, 0,'h'},
        {"verbose",     no_argument, 0,'v'},
        {"quiet",       no_argument, 0,'q'},
        {"debug",       no_argument, 0,'d'},
        {"resume",      no_argument, 0, 0 },
        {"sim-noise",   no_argument, 0, 0 },
        {"conf-noise",  no_argument, 0, 0 },
        {"stationary",  no_argument, 0, 0 },
        {"phase",       no_argument, 0, 0 },
        {"sangria",     no_argument, 0, 0 },
        {"prior",       no_argument, 0, 0 },
        {"no-burnin",   no_argument, 0, 0 },
        {"no-rj",       no_argument, 0, 0 },
        {"calibration", no_argument, 0, 0 },
        {"h5-no-mbh",   no_argument, 0, 0 },
        {"h5-no-ucb",   no_argument, 0, 0 },
        {"h5-no-ucb-hi",no_argument, 0, 0 },
        {"h5-no-vgb",   no_argument, 0, 0 },
        {"h5-no-noise", no_argument, 0, 0 },
        {0, 0, 0, 0}
    };
    
    int opt=0;
    int long_index=0;
    
    //Loop through argv string and pluck out arguments
    while ((opt = getopt_long_only(argc, argv_copy,"apl:b:", long_options, &long_index )) != -1)
    {
        switch (opt)
        {
                
            case 0:
                if(strcmp("samples",     long_options[long_index].name) == 0) data->NFFT        = atoi(optarg);
                if(strcmp("layers",      long_options[long_index].name) == 0) data->Nlayer      = atoi(optarg);
                if(strcmp("padding",     long_options[long_index].name) == 0) data->qpad        = atoi(optarg);
                if(strcmp("start-time",  long_options[long_index].name) == 0) data->t0          = (double)atof(optarg);
                if(strcmp("chains",      long_options[long_index].name) == 0) chain->NC         = atoi(optarg);
                if(strcmp("chainseed",   long_options[long_index].name) == 0) data->cseed       = (unsigned int)atoi(optarg);
                if(strcmp("noiseseed",   long_options[long_index].name) == 0) data->nseed       = (unsigned int)atoi(optarg);
                if(strcmp("injseed",     long_options[long_index].name) == 0) data->iseed       = (unsigned int)atoi(optarg);
                if(strcmp("sim-noise",   long_options[long_index].name) == 0) flags->simNoise   = 1;
                if(strcmp("conf-noise",  long_options[long_index].name) == 0) flags->confNoise  = 1;
                if(strcmp("stationary",  long_options[long_index].name) == 0) flags->stationary = 1;
                if(strcmp("prior",       long_options[long_index].name) == 0) flags->prior      = 1;
                if(strcmp("no-burnin",   long_options[long_index].name) == 0) flags->burnin     = 0;
                if(strcmp("no-rj",       long_options[long_index].name) == 0) flags->rj         = 0;
                if(strcmp("calibration", long_options[long_index].name) == 0) flags->calibration= 1;
                if(strcmp("resume",      long_options[long_index].name) == 0) flags->resume     = 1;
                if(strcmp("h5-no-mbh",   long_options[long_index].name) == 0) flags->no_mbh     = 1;
                if(strcmp("h5-no-ucb",   long_options[long_index].name) == 0) flags->no_ucb     = 1;
                if(strcmp("h5-no-vgb",   long_options[long_index].name) == 0) flags->no_vgb     = 1;
                if(strcmp("h5-no-ucb-hi",long_options[long_index].name) == 0) flags->no_ucb_hi  = 1;
                if(strcmp("h5-no-noise", long_options[long_index].name) == 0) flags->no_noise   = 1;
                if(strcmp("threads",     long_options[long_index].name) == 0) flags->threads    = atoi(optarg);
                if(strcmp("rundir",      long_options[long_index].name) == 0) strcpy(flags->runDir,optarg);
                if(strcmp("phase",       long_options[long_index].name) == 0) sprintf(data->format,"phase");
                if(strcmp("sangria",     long_options[long_index].name) == 0) sprintf(data->format,"sangria");
                if(strcmp("fmin",        long_options[long_index].name) == 0) sscanf(optarg, "%lg", &data->fmin);
                if(strcmp("fmax",        long_options[long_index].name) == 0)
                {
                    set_fmax_flag = 1;
                    sscanf(optarg, "%lg", &data->fmax);
                }
                if(strcmp("duration",    long_options[long_index].name) == 0)
                {   
                    data->T   = (double)atof(optarg);
                    if(!strcmp(data->basis,"wavelet")) data->T = floor(data->T/WAVELET_DURATION)*WAVELET_DURATION;
                    data->sqT = sqrt(data->T);
                }
                if(strcmp("steps",       long_options[long_index].name) == 0)
                {
                    flags->NMCMC = atoi(optarg);
                    flags->NBURN = flags->NMCMC;
                }
                if(strcmp("data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("h5-data", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->hdf5Data = 1;
                    flags->strainData = 1;
                    sprintf(data->fileName,"%s",optarg);
                }
                if(strcmp("psd", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->psd = 1;
                    sprintf(flags->psdFile,"%s",optarg);
                }
                if(strcmp("orbit", long_options[long_index].name) == 0)
                {
                    checkfile(optarg);
                    flags->orbit = 1;
                    sprintf(orbit->OrbitFileName,"%s",optarg);
                }
                if(strcmp("channels",long_options[long_index].name) == 0)
                {
                    data->Nchannel = (int)atoi(optarg);
                    if(data->Nchannel<1 || data->Nchannel>3)
                    {
                        fprintf(stderr,"Requested umber of channels (%i) not supported\n",data->Nchannel);
                        fprintf(stderr,"Use --channels 1 for X (Michelson) data\n");
                        fprintf(stderr,"    --channels 2 for AE data\n");
                        fprintf(stderr,"    --channels 3 for XYZ data\n");
                        exit(1);
                    }
                }
                break;
            case 'd' : flags->debug = 1;
                break;
            case 'h' : flags->help = 1;
                break;
            case 'v' : flags->verbose = 1;
                break;
            case 'q' : flags->quiet = 1;
                break;
            default:
                break;
        }
    }
    if(flags->cheat || !flags->burnin) flags->NBURN = 0;
    
    if(flags->verbose && flags->quiet)
    {
        fprintf(stderr,"--verbose and --quiet flags are in conflict\n");
        exit(1);
    }
    
    //Chains should be a multiple of threads for best usage of cores
    if(chain->NC % flags->threads !=0){
        chain->NC += flags->threads - (chain->NC % flags->threads);
    }
    
    //override size of data if fmax was requested
    if(set_fmax_flag)
    {
        if(!strcmp(data->basis,"fourier")) data->NFFT = (int)floor((data->fmax - data->fmin)*data->T); //number of frequency bins
        if(!strcmp(data->basis,"wavelet")) data->Nlayer = (int)( (ceil(data->fmax/WAVELET_BANDWIDTH) - floor(data->fmin/WAVELET_BANDWIDTH)) ); //number of frequency layers
    }

    //pad data
    if(!strcmp(data->basis,"fourier"))
    {
        data->NFFT += 2*data->qpad;
        data->fmin -= data->qpad/data->T;
    }
    if(!strcmp(data->basis,"wavelet")) 
    {
        data->Nlayer += 2; // pad with layers above and below
        data->fmin -= WAVELET_BANDWIDTH;
    }

    
    
    //map fmin to nearest bin
    if(!strcmp(data->basis,"fourier"))
    {
        data->fmin = floor(data->fmin*data->T)/data->T;
        data->fmax = data->fmin + (double)data->NFFT/data->T;
    }
    if(!strcmp(data->basis,"wavelet"))
    {
        data->fmin = floor(data->fmin/WAVELET_BANDWIDTH)*WAVELET_BANDWIDTH;
        data->fmax = data->fmin + (double)(data->Nlayer)*WAVELET_BANDWIDTH;

    }

    //after all of that resize data
    if(!strcmp(data->basis,"fourier")) data->N = data->NFFT*2;
    if(!strcmp(data->basis,"wavelet")) 
    {
        data->N = data->Nlayer*(int)floor(data->T / WAVELET_DURATION);
        data->NFFT = data->N/2.;
    }

    //Print version control
//    sprintf(filename,"glass.log");
//    FILE *runlog = fopen(filename,"w");
//    print_version(runlog);
    
    //Report on set parameters
//    if(!flags->quiet) print_run_settings(argc, argv, data, orbit, flags, stdout);
//    print_run_settings(argc, argv, data, orbit, flags, runlog);
    
//    fclose(runlog);

    //reset opt counter
    optind = 0;

    //free placeholder for argvs
    for(int i=0; i<=argc; i++)free(argv_copy[i]);
    free(argv_copy);

}

void copy_argv(int argc, char **argv, char **new_argv)
{
    for(int i = 0; i < argc; ++i)
    {
        size_t length = strlen(argv[i])+1;
        new_argv[i] = malloc(length);
        memcpy(new_argv[i], argv[i], length);
    }
    new_argv[argc] = NULL;
}

int checkfile(char filename[])
{
    FILE *fptr = fopen(filename, "r");
    if(fptr)
    {
        fclose(fptr);
        return 1;
    }
    else
    {
        fprintf(stderr,"File %s does not exist\n",filename);
        fprintf(stderr,"\n");
        exit(1);
    }
}
