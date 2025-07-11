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

#include <glass_utils.h>

#include "glass_ucb_model.h"
#include "glass_ucb_waveform.h"


double analytic_snr(double A, double Sn, double Sf, double sqT)
{
    return A*sqT*Sf/sqrt(Sn); //not exactly what's in paper--calibrated against (h|h)
}

double snr(struct Source *source, struct Noise *noise)
{
    double snr2=0.0;
    switch(source->tdi->Nchannel)
    {
        case 1: //Michelson
            snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->invC[0][0],source->tdi->N/2);
            break;
        case 2: //A&E
            snr2 += fourier_nwip(source->tdi->A,source->tdi->A,noise->invC[0][0],source->tdi->N/2);
            snr2 += fourier_nwip(source->tdi->E,source->tdi->E,noise->invC[1][1],source->tdi->N/2);
            break;
        case 3: //XYZ
            snr2 += fourier_nwip(source->tdi->X,source->tdi->X,noise->invC[0][0],source->tdi->N/2);
            snr2 += fourier_nwip(source->tdi->Y,source->tdi->Y,noise->invC[1][1],source->tdi->N/2);
            snr2 += fourier_nwip(source->tdi->Z,source->tdi->Z,noise->invC[2][2],source->tdi->N/2);
            snr2 += fourier_nwip(source->tdi->X,source->tdi->Y,noise->invC[0][1],source->tdi->N/2)*2.;
            snr2 += fourier_nwip(source->tdi->X,source->tdi->Z,noise->invC[0][2],source->tdi->N/2)*2.;
            snr2 += fourier_nwip(source->tdi->Y,source->tdi->Z,noise->invC[1][2],source->tdi->N/2)*2.;
            break;
    }
    
    return(sqrt(snr2));
}

double snr_wavelet(struct Source *source, struct Noise *noise)
{
    double snr2 = 0.0;
    snr2 += wavelet_nwip(source->tdi->X, source->tdi->X, noise->invC[0][0], source->list, source->Nlist);
    snr2 += wavelet_nwip(source->tdi->Y, source->tdi->Y, noise->invC[1][1], source->list, source->Nlist);
    snr2 += wavelet_nwip(source->tdi->Z, source->tdi->Z, noise->invC[2][2], source->list, source->Nlist);
    snr2 += wavelet_nwip(source->tdi->X, source->tdi->Y, noise->invC[0][1], source->list, source->Nlist)*2;
    snr2 += wavelet_nwip(source->tdi->X, source->tdi->Z, noise->invC[0][2], source->list, source->Nlist)*2;
    snr2 += wavelet_nwip(source->tdi->Y, source->tdi->Z, noise->invC[1][2], source->list, source->Nlist)*2;
    return sqrt(snr2);
}

double waveform_match_wavelet(struct Source *a, struct Source *b, struct Noise *noise)
{
    double aa = pow(snr_wavelet(a,noise),2);
    double bb = pow(snr_wavelet(b,noise),2);

    struct TDI *a_full = malloc(sizeof(struct TDI));
    struct TDI *b_full = malloc(sizeof(struct TDI));
    alloc_tdi(a_full,a->tdi->N,a->tdi->Nchannel);
    alloc_tdi(b_full,b->tdi->N,b->tdi->Nchannel);

    for(int n=0; n<a->Nlist; n++)
    {
        if(a->list[n]>0)
        {
            int k = a->list[n];
            a_full->X[k] = a->tdi->X[k];
            a_full->Y[k] = a->tdi->Y[k];
            a_full->Z[k] = a->tdi->Z[k];
        }
    }
    for(int n=0; n<b->Nlist; n++)
    {
        if(b->list[n]>0)
        {
            int k = b->list[n];
            b_full->X[k] = b->tdi->X[k];
            b_full->Y[k] = b->tdi->Y[k];
            b_full->Z[k] = b->tdi->Z[k];
        }
    }

    int *list = int_vector(a->Nlist+b->Nlist);
    int N;
    list_union(a->list,b->list,a->Nlist,b->Nlist,list,&N);

    double ab = 0.0;
    ab += wavelet_nwip(a_full->X, b_full->X, noise->invC[0][0], list, N);
    ab += wavelet_nwip(a_full->Y, b_full->Y, noise->invC[1][1], list, N);
    ab += wavelet_nwip(a_full->Z, b_full->Z, noise->invC[2][2], list, N);
    ab += wavelet_nwip(a_full->X, b_full->Y, noise->invC[0][1], list, N);
    ab += wavelet_nwip(a_full->X, b_full->Z, noise->invC[0][2], list, N);
    ab += wavelet_nwip(a_full->Y, b_full->Z, noise->invC[1][2], list, N);
    ab += wavelet_nwip(a_full->Y, b_full->X, noise->invC[1][0], list, N);
    ab += wavelet_nwip(a_full->Z, b_full->X, noise->invC[2][0], list, N);
    ab += wavelet_nwip(a_full->Z, b_full->Y, noise->invC[2][1], list, N);

    double match = ab/sqrt(aa*bb);

    free(list);
    free_tdi(a_full);
    free_tdi(b_full);

    return match;
}


double snr_prior(double SNR)
{
    //SNRPEAK defined in glass_constants.h
    double dfac  = 1.+SNR/(4.*SNRPEAK);
    double dfac5 = ipow(dfac,5);
    return (3.*SNR)/(4.*SNRPEAK*SNRPEAK*dfac5);
}

double waveform_match(struct Source *a, struct Source *b, struct Noise *noise)
{
    int N = a->tdi->N;
    int NFFT = 2*N;
    double match=0;
    
    double *a_A = calloc(NFFT,sizeof(double));
    double *a_E = calloc(NFFT,sizeof(double));
    double *b_A = calloc(NFFT,sizeof(double));
    double *b_E = calloc(NFFT,sizeof(double));
    
    int qmin = a->qmin - a->imin;
    
    
    //Align waveforms into arrays for summing
    for(int i=0; i<a->BW; i++)
    {
        int j = i+a->qmin-qmin;
        
        if(j>-1 && j<N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            a_A[j_re] = a->tdi->A[i_re];
            a_A[j_im] = a->tdi->A[i_im];
            a_E[j_re] = a->tdi->E[i_re];
            a_E[j_im] = a->tdi->E[i_im];
        }//check that index is in range
    }//loop over waveform bins
    
    //Align waveforms into arrays for summing
    for(int i=0; i<b->BW; i++)
    {
        int j = i+b->qmin-qmin;
        
        if(j>-1 && j<N)
        {
            int i_re = 2*i;
            int i_im = i_re+1;
            int j_re = 2*j;
            int j_im = j_re+1;
            
            b_A[j_re] = b->tdi->A[i_re];
            b_A[j_im] = b->tdi->A[i_im];
            b_E[j_re] = b->tdi->E[i_re];
            b_E[j_im] = b->tdi->E[i_im];
        }//check that index is in range
    }//loop over waveform bins
    
    
    double aa = fourier_nwip(a_A,a_A,noise->invC[0][0],N) + fourier_nwip(a_E,a_E,noise->invC[1][1],N);
    double bb = fourier_nwip(b_A,b_A,noise->invC[0][0],N) + fourier_nwip(b_E,b_E,noise->invC[1][1],N);
    double ab = fourier_nwip(a_A,b_A,noise->invC[0][0],N) + fourier_nwip(a_E,b_E,noise->invC[1][1],N);
    
    match = ab/sqrt(aa*bb);
    
    free(a_A);
    free(a_E);
    free(b_A);
    free(b_E);
    
    return match;
}

double waveform_distance(struct Source *a, struct Source *b, struct Noise *noise)
{
  int N = a->tdi->N;
  int NFFT = 2*N;

  double *a_A = calloc(NFFT,sizeof(double));
  double *a_E = calloc(NFFT,sizeof(double));
  double *b_A = calloc(NFFT,sizeof(double));
  double *b_E = calloc(NFFT,sizeof(double));

  int qmin = a->qmin - a->imin;


  //Align waveforms into arrays for summing
  for(int i=0; i<a->BW; i++)
  {
    int j = i+a->qmin-qmin;
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      a_A[j_re] = a->tdi->A[i_re];
      a_A[j_im] = a->tdi->A[i_im];
      a_E[j_re] = a->tdi->E[i_re];
      a_E[j_im] = a->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  //Align waveforms into arrays for summing
  for(int i=0; i<b->BW; i++)
  {
    int j = i+b->qmin-qmin;
    
    if(j>-1 && j<N)
    {
      int i_re = 2*i;
      int i_im = i_re+1;
      int j_re = 2*j;
      int j_im = j_re+1;
      
      b_A[j_re] = b->tdi->A[i_re];
      b_A[j_im] = b->tdi->A[i_im];
      b_E[j_re] = b->tdi->E[i_re];
      b_E[j_im] = b->tdi->E[i_im];
    }//check that index is in range
  }//loop over waveform bins

  
    double aa = fourier_nwip(a_A,a_A,noise->invC[0][0],N) + fourier_nwip(a_E,a_E,noise->invC[1][1],N);
    double bb = fourier_nwip(b_A,b_A,noise->invC[0][0],N) + fourier_nwip(b_E,b_E,noise->invC[1][1],N);
    double ab = fourier_nwip(a_A,b_A,noise->invC[0][0],N) + fourier_nwip(a_E,b_E,noise->invC[1][1],N);

  double distance = (aa + bb - 2*ab)/4.;

  free(a_A);
  free(a_E);
  free(b_A);
  free(b_E);
  
  return distance;
}

double ucb_fdot(double Mc, double f0)
{
    double f = f0;
    double M = Mc*TSUN;
    double Q = 19.2;//96./5.
    
    return Q*pow(pow(M_PI,8)*pow(M,5)*pow(f,11),1./3.);  
}
double ucb_chirpmass(double f0, double dfdt)
{
    double f = f0;
    double fd = dfdt;
    double pi83 = 21.170591578193; //pow(pi,8./3.)

    return pow(fd/(96./5.)/pi83/pow(f,11./3.), 3./5.)/TSUN;
}

double ucb_distance(double f0, double dfdt, double A)
{
    double f    = f0;
    double fd = dfdt;
    double amp   = A;
    return ((5./48.)*(fd/(M_PI*M_PI*f*f*f*amp))*CLIGHT/PC); //seconds  !check notes on 02/28!
}

double ucb_phase(double t, double *params, double T)
{
    double f0    = params[0]/T;
    double phi0  = params[6];
    double fdot  = params[7]/T/T;
    double fddot = 0.0;
    
    /*
     * LDC phase parameter in key files is
     * -phi0
     */
    return -phi0 + PI2*( f0*t + 0.5*fdot*t*t + 1.0/6.0*fddot*t*t*t );
}

double ucb_amplitude(double t, double *params, double T)
{
    double f0    = params[0]/T;
    double A0    = exp(params[3]);
    double fdot  = params[7]/T/T;

    return A0 * ( 1.0 + 2.0/3.0*fdot/f0*t );
}

void ucb_barycenter_waveform(double *params, int N, double *times, double *phase, double *amp, double T)
{
    for(int n=0; n<N; n++)
    {
        phase[n] = ucb_phase(times[n],params, T);
        amp[n]   = ucb_amplitude(times[n],params, T);
    }
}

void ucb_fisher(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise)
{
    //TODO:  ucb_fisher should compute joint Fisher
    int i,j,n;
        
    double epsilon    = 1.0e-6;
    //double invepsilon2= 1./(2.*epsilon);
    double invepsilon2= 1./(epsilon);
    double invstep;
    
    // Plus and minus parameters:
    double *params_p = calloc(UCB_MODEL_NP,sizeof(double));
    //double *params_m = calloc(NP,sizeof(double));
    
    // Plus and minus templates for each detector:
    struct Source *wave_p = malloc(sizeof(struct Source));
    //struct Source *wave_m = malloc(sizeof(struct Source));
    alloc_source(wave_p, data->N, data->Nchannel);

    //alloc_source(wave_m, data->N, data->Nchannel, NP);
    
    // TDI variables to hold derivatives of h
    struct TDI **dhdx = malloc(UCB_MODEL_NP*sizeof(struct TDI *));
    for(n=0; n<UCB_MODEL_NP; n++)
    {
        dhdx[n] = malloc(sizeof(struct TDI));
        alloc_tdi(dhdx[n], data->N, data->Nchannel);
    }
    
    /* assumes all the parameters are log or angle */
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        //step size for derivatives
        invstep = invepsilon2;
        
        // copy parameters
        for(j=0; j<UCB_MODEL_NP; j++)
        {
            wave_p->params[j] = source->params[j];
            //wave_m->params[j] = source->params[j];
        }
        
        // perturb parameters
        wave_p->params[i] += epsilon;
        //wave_m->params[i] -= epsilon;
        
	    // catch when cosine parameters get pushed out of bounds
        if(i==1 || i==4)
        {
            if(wave_p->params[i] > 1.0) wave_p->params[i] = 1.0;
            //if(wave_m->params[i] <-1.0) wave_m->params[i] =-1.0;
        }

        // complete info in source structure
        map_array_to_params(wave_p, wave_p->params, data->T);
        //map_array_to_params(wave_m, wave_m->params, data->T);
        
        // clean up TDI arrays, just in case
        for(j=0; j<data->N; j++)
        {
            wave_p->tdi->X[j]=0.0;
            wave_p->tdi->Y[j]=0.0;
            wave_p->tdi->Z[j]=0.0;
            wave_p->tdi->A[j]=0.0;
            wave_p->tdi->E[j]=0.0;
//            wave_m->tdi->X[j]=0.0;
//            wave_m->tdi->Y[j]=0.0;
//            wave_m->tdi->Z[j]=0.0;
//            wave_m->tdi->A[j]=0.0;
//            wave_m->tdi->E[j]=0.0;
        }
        
        // align perturbed waveforms in data array
        ucb_alignment(orbit, data, wave_p);
        //ucb_alignment(orbit, data, wave_m);
        
        // compute perturbed waveforms
        ucb_waveform(orbit, data->format, data->T, data->t0, wave_p->params, UCB_MODEL_NP, wave_p->tdi->X, wave_p->tdi->Y, wave_p->tdi->Z, wave_p->tdi->A, wave_p->tdi->E, wave_p->BW, wave_p->tdi->Nchannel);
        //ucb_waveform(orbit, data->format, data->T, data->t0, wave_m->params, UCB_MODEL_NP, wave_m->tdi->X, wave_m->tdi->Y, wave_m->tdi->Z, wave_m->tdi->A, wave_m->tdi->E, wave_m->BW, wave_m->tdi->Nchannel);

        // central differencing derivatives of waveforms w.r.t. parameters
        switch(source->tdi->Nchannel)
        {
            case 1:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    dhdx[i]->X[n] = (wave_p->tdi->X[n] - source->tdi->X[n])*invstep;
                }
                break;
            case 2:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    dhdx[i]->A[n] = (wave_p->tdi->A[n] - source->tdi->A[n])*invstep;
                    dhdx[i]->E[n] = (wave_p->tdi->E[n] - source->tdi->E[n])*invstep;
                }
                break;
            case 3:
                for(n=0; n<wave_p->BW*2; n++)
                {
                    dhdx[i]->X[n] = (wave_p->tdi->X[n] - source->tdi->X[n])*invstep;
                    dhdx[i]->Y[n] = (wave_p->tdi->Y[n] - source->tdi->Y[n])*invstep;
                    dhdx[i]->Z[n] = (wave_p->tdi->Z[n] - source->tdi->Z[n])*invstep;
                }
                break;

        }
    }
    
    // Calculate fisher matrix
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        for(j=i; j<UCB_MODEL_NP; j++)
        {
            switch(source->tdi->Nchannel)
            {
                case 1:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->X, dhdx[j]->X, noise->invC[0][0], wave_p->BW);
                    break;
                case 2:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->A, dhdx[j]->A, noise->invC[0][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->E, dhdx[j]->E, noise->invC[1][1], wave_p->BW);
                    break;
                case 3:
                    source->fisher_matrix[i][j] = fourier_nwip(dhdx[i]->X, dhdx[j]->X, noise->invC[0][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Y, dhdx[j]->Y, noise->invC[1][1], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Z, dhdx[j]->Z, noise->invC[2][2], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->X, dhdx[j]->Y, noise->invC[0][1], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->X, dhdx[j]->Z, noise->invC[0][2], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Y, dhdx[j]->Z, noise->invC[1][2], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Y, dhdx[j]->X, noise->invC[1][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Z, dhdx[j]->X, noise->invC[2][0], wave_p->BW);
                    source->fisher_matrix[i][j] += fourier_nwip(dhdx[i]->Z, dhdx[j]->Y, noise->invC[2][1], wave_p->BW);
                    break;
            }
            if(source->fisher_matrix[i][j]!=source->fisher_matrix[i][j])
            {
                fprintf(stderr,"WARNING: nan matrix element (line %d of file %s)\n",__LINE__,__FILE__);
                fprintf(stderr, "fisher_matrix[%i][%i], Snf=[%g,%g]\n",i,j,noise->C[0][0][data->NFFT/2],noise->C[1][1][data->NFFT/2]);
                for(int k=0; k<UCB_MODEL_NP; k++)
                {
                    fprintf(stderr,"source->params[%i]=%g\n",k,source->params[k]);
                }
                source->fisher_matrix[i][j] = 10.0;
            }
            source->fisher_matrix[j][i] = source->fisher_matrix[i][j];
        }
    }
    
    // Calculate eigenvalues and eigenvectors of fisher matrix
    matrix_eigenstuff(source->fisher_matrix, source->fisher_evectr, source->fisher_evalue, UCB_MODEL_NP);
    
    free(params_p);
    //free(params_m);
    free_source(wave_p);
    //free_source(wave_m);
    
    for(n=0; n<UCB_MODEL_NP; n++) free_tdi(dhdx[n]);
    free(dhdx);
}

void ucb_fisher_wavelet(struct Orbit *orbit, struct Data *data, struct Source *source, struct Noise *noise)
{
    //TODO:  ucb_fisher should compute joint Fisher
    int i,j,n;
        
    double epsilon    = 1.0e-6;
    double invepsilon2= 1./(epsilon);
    double invstep;
    
    // Plus and minus parameters:
    double *params_p = calloc(UCB_MODEL_NP,sizeof(double));
    
    // Plus and minus templates for each detector:
    struct Source *wave_p = malloc(sizeof(struct Source));
    alloc_source(wave_p, data->N, data->Nchannel);
    
    // TDI variables to hold derivatives of h
    struct TDI **dhdx = malloc(UCB_MODEL_NP*sizeof(struct TDI *));
    for(n=0; n<UCB_MODEL_NP; n++)
    {
        dhdx[n] = malloc(sizeof(struct TDI));
        alloc_tdi(dhdx[n], data->N, data->Nchannel);
    }
    
    /* assumes all the parameters are log or angle */
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        //step size for derivatives
        invstep = invepsilon2;
        
        // copy parameters
        for(j=0; j<UCB_MODEL_NP; j++)
        {
            wave_p->params[j] = source->params[j];
        }
        
        // perturb parameters
        wave_p->params[i] += epsilon;
        
	    // catch when cosine parameters get pushed out of bounds
        if(i==1 || i==4)
        {
            if(wave_p->params[i] > 1.0) wave_p->params[i] = 1.0;
        }

        // complete info in source structure
        map_array_to_params(wave_p, wave_p->params, data->T);
        
        // clean up TDI arrays, just in case
        for(j=0; j<data->N; j++)
        {
            wave_p->tdi->X[j]=0.0;
            wave_p->tdi->Y[j]=0.0;
            wave_p->tdi->Z[j]=0.0;
            wave_p->tdi->A[j]=0.0;
            wave_p->tdi->E[j]=0.0;
        }
        
        
        // compute perturbed waveforms
        ucb_waveform_wavelet(orbit,data->wdm,data->T, data->t0, wave_p->params, wave_p->list, &wave_p->Nlist, wave_p->tdi->X, wave_p->tdi->Y, wave_p->tdi->Z);

        // central differencing derivatives of waveforms w.r.t. parameters
        for(n=0; n<wave_p->Nlist; n++)
        {
            int k = wave_p->list[n];
            dhdx[i]->X[k] = (wave_p->tdi->X[k] - source->tdi->X[k])*invstep;
            dhdx[i]->Y[k] = (wave_p->tdi->Y[k] - source->tdi->Y[k])*invstep;
            dhdx[i]->Z[k] = (wave_p->tdi->Z[k] - source->tdi->Z[k])*invstep;
        }

    }
    
    // Calculate fisher matrix
    for(i=0; i<UCB_MODEL_NP; i++)
    {
        for(j=i; j<UCB_MODEL_NP; j++)
        {
            source->fisher_matrix[i][j]  = wavelet_nwip(dhdx[i]->X, dhdx[j]->X, noise->invC[0][0], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Y, dhdx[j]->Y, noise->invC[1][1], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Z, dhdx[j]->Z, noise->invC[2][2], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->X, dhdx[j]->Y, noise->invC[0][1], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->X, dhdx[j]->Z, noise->invC[0][2], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Y, dhdx[j]->Z, noise->invC[1][2], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Y, dhdx[j]->X, noise->invC[1][0], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Z, dhdx[j]->X, noise->invC[2][0], wave_p->list, wave_p->Nlist);
            source->fisher_matrix[i][j] += wavelet_nwip(dhdx[i]->Z, dhdx[j]->Y, noise->invC[2][1], wave_p->list, wave_p->Nlist);

            if(source->fisher_matrix[i][j]!=source->fisher_matrix[i][j])
            {
                fprintf(stderr,"WARNING: nan matrix element (line %d of file %s)\n",__LINE__,__FILE__);
                fprintf(stderr, "fisher_matrix[%i][%i], Snf=[%g,%g]\n",i,j,noise->C[0][0][data->N/2],noise->C[1][1][data->N/2]);
                for(int k=0; k<UCB_MODEL_NP; k++)
                {
                    fprintf(stderr,"source->params[%i]=%g\n",k,source->params[k]);
                }
                source->fisher_matrix[i][j] = 10.0;
            }
            source->fisher_matrix[j][i] = source->fisher_matrix[i][j];
        }
    }
    
    // Calculate eigenvalues and eigenvectors of fisher matrix
    matrix_eigenstuff(source->fisher_matrix, source->fisher_evectr, source->fisher_evalue, UCB_MODEL_NP);
    
    free(params_p);
    free_source(wave_p);
    
    for(n=0; n<UCB_MODEL_NP; n++) free_tdi(dhdx[n]);
    free(dhdx);
}

int ucb_bandwidth(double L, double fstar, double f, double fdot, double costheta, double A, double T, int N)
{
    int Nmin = 16;
    int Nmax = (int)pow(2,(int)log2((double)(N/2)));
    
    double sqT=sqrt(T);
    
    double sf = sin(f/fstar); //sin(f/f*)
    double sn = AEnoise(L,fstar,f);
    
    //Doppler spreading
    double sintheta = sin(acos(costheta));
    double bw = 2*T*((4.+PI2*f*(AU/CLIGHT)*sintheta)/YEAR + fabs(fdot)*T);
    int DS = (int)pow(2,(int)log2(bw-1)+1);
    if(DS > Nmax) DS = Nmax;
    if(DS < Nmin) DS = Nmin;
    
    
    //Sinc spreading
    double SNRm = analytic_snr(A, sn, sf, sqT);
    
    int SS = (int)pow(2,(int)log2(SNRm-1)+1);
    
    if(SS > Nmax) SS = Nmax;
    if(SS < Nmin) SS = Nmin;
    
    return (DS > SS) ? DS : SS; //return largest spread as bandwidth
}

void ucb_alignment(struct Orbit *orbit, struct Data *data, struct Source *source)
{
    map_array_to_params(source, source->params, data->T);
    
    source->BW   = 2*ucb_bandwidth(orbit->L, orbit->fstar, source->f0, source->dfdt, source->costheta, source->amp, data->T, data->NFFT);
    source->qmin = (int)(source->f0*data->T) - source->BW/2;
    source->qmax = source->qmin+source->BW;
    source->imin = source->qmin - data->qmin;
    source->imax = source->imin + source->BW;  
}

void ucb_waveform(struct Orbit *orbit, char *format, double T, double t0, double *params, int NParams, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI)
{
    /*   Indicies   */
    int i,j,n;
    /*   Carrier frequency bin  */
    long q;
    /*   Bandwidth      */
    int BW2   = BW*2;
    double invBW2 = 1./(double)BW2;
    
    /*   Gravitational Wave location vector   */
    double k[4];
    /*   Polarization basis tensors   */
    double eplus[4][4], ecross[4][4];
    /*   Spacecraft position and separation vector   */
    double *x, *y, *z;
    /*   Dot products   */
    double kdotx[4]={0},kdotr[4][4];
    /*   Convenient quantities   */
    double dplus[4][4],dcross[4][4];
    /*   GW source parameters   */
    double phi, psi, amp, Aplus, Across, f0, dfdt, d2fdt2, phi0;
    double costh, cosi, cosps, sinps;
    /*   Time and distance variables   */
    double t, xi[4] = {0};
    /*   Gravitational wave frequency & ratio of f and transfer frequency f*  */
    double f[4] = {0},fonfs[4] = {0};
    /*   LISA response to slow terms (Real & Imaginary pieces)   */
    //Static quantities (Re and Im)
    double DPr, DPi, DCr, DCi;
    //Time varrying quantities (Re & Im) broken up into convenient segments
    double TR[4][4], TI[4][4];
    //Miscellaneous constants used to speed up calculations
    double df;
    /*   Fourier coefficients before FFT and after convolution  */
    //Time series of slowly evolving terms at each vertex
    double *data12, *data13, *data21, *data23, *data31, *data32;
    //Fourier coefficients of slowly evolving terms (numerical)
    double a12[BW2+3], a13[BW2+3], a21[BW2+3], a23[BW2+3], a31[BW2+3], a32[BW2+3];
    //Package cij's into proper form for TDI subroutines
    double ***d;
    
    /*   Allocating Arrays   */
    x = calloc(4,sizeof(double));
    y = calloc(4,sizeof(double));
    z = calloc(4,sizeof(double));
    
    data12 = calloc((BW2+1),sizeof(double));
    data21 = calloc((BW2+1),sizeof(double));
    data31 = calloc((BW2+1),sizeof(double));
    data13 = calloc((BW2+1),sizeof(double));
    data23 = calloc((BW2+1),sizeof(double));
    data32 = calloc((BW2+1),sizeof(double));
    
    d = malloc(sizeof(double**)*4);
    for(i=0; i<4; i++)
    {
        d[i] = malloc(sizeof(double*)*4);
        for(j=0; j<4; j++)
        {
            	d[i][j] = calloc((BW2+1),sizeof(double));
        }
    }
    
    /*   Gravitational Wave source parameters   */
    
    f0     = params[0]/T;
    costh  = params[1];
    phi    = params[2];
    amp    = exp(params[3]);
    cosi   = params[4];
    psi    = params[5];
    phi0   = params[6];
    dfdt   = 0.0;
    d2fdt2 = 0.0;
    if(NParams>7)
        dfdt   = params[7]/(T*T);
    if(NParams>8)
        d2fdt2 = params[8]/(T*T*T);
    
    //Calculate carrier frequency bin
    q = (long)(f0*T);
    
    //Calculate cos and sin of sky position, inclination, polarization
    cosps	= cos(2.*psi);
    sinps	= sin(2.*psi);
    
    //Calculate GW polarization amplitudes
    Aplus  =  amp*(1.+cosi*cosi);
    Across = -amp*(2.0*cosi);
    
    df = PI2*(((double)q)/T);
    
    //Calculate constant pieces of transfer functions
    DPr =  Aplus*cosps;
    DPi = -Across*sinps;
    DCr = -Aplus*sinps;
    DCi = -Across*cosps;
    
    LISA_polarization_tensor(costh, phi, eplus, ecross, k);
    
    /* Main loop over signal bandwidth */
    for(n=1; n<=BW; n++)
    {
        //First time sample must be at t=0 for phasing
        t = t0 + T*(double)(n-1)/(double)BW;
        
        //Calculate position of each spacecraft at time t
        (*orbit->orbit_function)(orbit, t, x, y, z);
        
        //Form LISA detector tensor et al based on spacecraft and source location
        LISA_detector_tensor(orbit->L,eplus,ecross,x,y,z,k,dplus,dcross,kdotr);
        
        //Calculating LISA Transfer function
        for(i=1; i<=3; i++)
        {
            //Dot product of propogation vector with location of spacecrat i
            kdotx[i] = (x[i]*k[1]+y[i]*k[2]+z[i]*k[3])/CLIGHT;
            
            //Wave arrival time at spacecraft i
            xi[i] = t - kdotx[i];
            
            //Zeroeth order approximation to frequency at spacecraft i
            f[i] = f0;
            
            //First order in frequency
            if(NParams>7) f[i] += dfdt*xi[i];
            
            //Second order in frequency
            if(NParams>8) f[i] += 0.5*d2fdt2*xi[i]*xi[i];
            
            //Ratio of true frequency to transfer frequency
            fonfs[i] = f[i]/orbit->fstar;

            //Argument of complex exponentials
            /*
             * LDC phase parameter in key files is
             * -phi0, hence the -phi0 in arg2
             */
            double arg2 = PI2*f0*xi[i] - phi0 - df*t;
            
            
            //First order frequency evolution
            if(NParams>7) arg2 += M_PI*dfdt*xi[i]*xi[i];
            
            //Second order frequency evolution
            if(NParams>8) arg2 += (M_PI/3.0)*d2fdt2*xi[i]*xi[i]*xi[i];
            
            //Evolution of amplitude
            double aevol = 1.0;
            
            //First order amplitude evolution
            if(NParams>7) aevol += 0.66666666666666666666*dfdt/f0*xi[i];
            
            //Second order amplitude evolution
            //if(NParams>8) aevol += const.*d2fdt2*xi[i]*xi[i]/f0;
            
            for(j=1; j<=3; j++)
            {
                if(i!=j)
                {
                    //Argument of transfer function
                    /*
                     * Set to match Radler LDC convention
                     *
                     https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc
                     */
                    double arg1 = 0.5*fonfs[i]*(1.0 + kdotr[i][j]);
                    
                    //Transfer function
                    double sinc = 0.25*sin(arg1)/arg1;
                    
                    //Real and imaginary pieces of time series (no complex exponential)
                    double tran1r = aevol*(dplus[i][j]*DPr + dcross[i][j]*DCr);
                    double tran1i = aevol*(dplus[i][j]*DPi + dcross[i][j]*DCi);
                    
                    /*
                     * Set to match Sangria LDC convention
                     * which defines the GW as e(-i Phi)
                    */
                    //Real and imaginary components of complex exponential
                    double tran2r = cos(arg1 - arg2);
                    double tran2i = sin(arg1 - arg2);

                    //Real & Imaginary part of the slowly evolving signal
                    TR[i][j] = sinc*( tran1r*tran2r + tran1i*tran2i);
                    TI[i][j] = sinc*(-tran1r*tran2i + tran1i*tran2r);
                }
            }
        }
        
        //Fill  time series data arrays with slowly evolving signal->
        //dataij corresponds to fractional arm length difference yij
        j = 2*n;
        i = j-1;
        data12[i] = TR[1][2];   data21[i] = TR[2][1];   data31[i] = TR[3][1];
        data12[j] = TI[1][2];   data21[j] = TI[2][1];   data31[j] = TI[3][1];
        data13[i] = TR[1][3];   data23[i] = TR[2][3];   data32[i] = TR[3][2];
        data13[j] = TI[1][3];   data23[j] = TI[2][3];   data32[j] = TI[3][2];
    }
    
    /*   Numerical Fourier transform of slowly evolving signal */
    glass_forward_complex_fft(data12+1, BW);
    glass_forward_complex_fft(data21+1, BW);
    glass_forward_complex_fft(data31+1, BW);
    glass_forward_complex_fft(data13+1, BW);
    glass_forward_complex_fft(data23+1, BW);
    glass_forward_complex_fft(data32+1, BW);
     
    //Unpack arrays from fft and normalize
    for(i=1; i<=BW; i++)
    {
        j = i + BW;
        a12[i] = data12[j]*invBW2;  a21[i] = data21[j]*invBW2;  a31[i] = data31[j]*invBW2;
        a12[j] = data12[i]*invBW2;  a21[j] = data21[i]*invBW2;  a31[j] = data31[i]*invBW2;
        a13[i] = data13[j]*invBW2;  a23[i] = data23[j]*invBW2;  a32[i] = data32[j]*invBW2;
        a13[j] = data13[i]*invBW2;  a23[j] = data23[i]*invBW2;  a32[j] = data32[i]*invBW2;
    }
    
    /*   Renormalize so that the resulting time series is real   */
    for(i=1; i<=BW2; i++)
    {
        d[1][2][i] = a12[i];  d[2][1][i] = a21[i];  d[3][1][i] = a31[i];
        d[1][3][i] = a13[i];  d[2][3][i] = a23[i];  d[3][2][i] = a32[i];
    }
    
    /*   Call subroutines for synthesizing different TDI data channels  */
    if(strcmp("phase",format) == 0)
        LISA_tdi(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
    else if(strcmp("frequency",format) == 0)
        LISA_tdi_FF(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
    else if(strcmp("sangria",format) == 0)
        LISA_tdi_Sangria(orbit->L, orbit->fstar, T, d, f0, q, X-1, Y-1, Z-1, A-1, E-1, BW, NI);
    else
    {
        fprintf(stderr,"Unsupported data format %s",format);
        exit(1);
    }

    /*   Free Arrays   */
    free(x);
    free(y);
    free(z);
    
    free(data12);
    free(data21);
    free(data31);
    free(data13);
    free(data23);
    free(data32);
    
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            free(d[i][j]);
        }
        free(d[i]);
    }
    free(d);
    
    
    return;
}

static void ucb_wavelet_layers(double Tobs, double *params, struct Wavelets *wdm, int *jstart, int *jwidth)
{
    double fmin, fmax;
    double dfmin, dfmax;
    int jmin, jmax;
    
    //get start and stop frequencey of signal
    double fstart = params[0];
    double fstop  = (params[0]+params[7]*Tobs);//+0.5*params[8]*Tobs*Tobs);
    
    //find min and max frequencies including maximum possible Doppler shift
    if(fstart < fstop) //inspiral
    {
        fmin = fstart*(1.0-VEARTH);
        fmax = fstop *(1.0+VEARTH);
    }
    else //outspiral
    {
        fmin = fstop *(1.0-VEARTH);
        fmax = fstart*(1.0+VEARTH);
    }
    
    //frequency layer of start frequency
    int j = (int)rint(fmin/wdm->df);
    
    //get distance from top and bottom of frequency layer
    dfmin = fmin - j*wdm->df;
    dfmax = fmax - j*wdm->df;
    
    
    //find min and max layer with signal power, accounting for power bleeding above or below carrier
    jmin = j;
    jmax = j;
    
    if(dfmin < 0.0 && fabs(dfmin) > wdm->A/PI2) jmin = j-1;     // signal will extend into layer below
    if(dfmax > 0.0 && fabs(dfmax) > wdm->A/PI2) jmax = j+1;     // signal will extend into layer above
        
    *jstart = jmin;
    *jwidth = jmax-jmin+1;
}

/* Heterodyne wavelet transform */
void ucb_waveform_wavelet(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z)
{
    int Nspline = orbit->Norb;
    double dt = Tobs/(double)(Nspline-1);
    
    // get amplitude and phase at Barycenter on the orbit interpolation grid (with margin)...
    double *time_ssb  = double_vector(Nspline);
    double *amp_ssb   = double_vector(Nspline);
    double *phase_ssb = double_vector(Nspline);
    
    // store time array for full data on orbit cadence
    for(int i=0; i< Nspline; i++) time_ssb[i] = t0 + i*dt;
    
    //get ucb waveform on orbit grid
    ucb_barycenter_waveform(params, Nspline, orbit->t, phase_ssb, amp_ssb, Tobs);
    
    // get frequency layers containing signal
    int min_layer; //bottom layer
    int Nlayers;   //number of layers
    ucb_wavelet_layers(Tobs, params, wdm, &min_layer, &Nlayers);
    
    /*
    Get spline interpolant for SSB phase and amplitude
    */
    struct CubicSpline *amp_ssb_spline   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_ssb_spline = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_ssb_spline,orbit->t,amp_ssb);
    initialize_cubic_spline(phase_ssb_spline,orbit->t,phase_ssb);
    
    /*
     Resample SSB phase to reference spacecraft
     */
    double *phase_sc = double_vector(Nspline);
    double *time_sc  = double_vector(Nspline);
    
    // shift reference times from Barycenter to S/C 1
    double costh = params[1];
    double phi   = params[2];
    LISA_spacecraft_to_barycenter_time(orbit, costh, phi, time_ssb, time_sc, Nspline, -1);

    // get signal phase at S/C 1
    for(int i=0; i< Nspline; i++)
        phase_sc[i] = spline_interpolation(phase_ssb_spline, time_sc[i]);

    free(time_sc);
    
    /*
     Downsample waveform (i.e. shift to lower frequency layer)
     */
    int N_ds     = wdm->NT*(Nlayers+1); //number of downsampled data points
    double dt_ds = wdm->dt/(double)(Nlayers+1); //downsampled data cadence
   
    double *phase_ds  = double_vector(N_ds);  //downsampled phase
    double *time_ds   = double_vector(N_ds);  //downsampled time
    double *phase_het = double_vector(N_ds);  //heterodyne phase
    
    double f0 = (min_layer-1)*wdm->df; //"carrier" frequency
    
    for(int i=0; i<N_ds; i++)
    {
        time_ds[i]   = t0 + i*dt_ds;
        phase_het[i] = PI2 * f0 * time_ds[i];
    }
    
    // shift reference times from Barycenter to spacecraft 0
    time_sc = double_vector(N_ds);
    LISA_spacecraft_to_barycenter_time(orbit, costh, phi, time_ds, time_sc, N_ds, -1);
    
    for(int i=0; i<N_ds; i++)
        phase_ds[i] = spline_interpolation(phase_ssb_spline, time_sc[i]);

    free(time_sc);

    
    /*
    Get TDI responses back in terms of phase and amplitude
    */
    struct TDI *tdi_phase = malloc(sizeof(struct TDI));
    struct TDI *tdi_amp = malloc(sizeof(struct TDI));
    alloc_tdi(tdi_phase,Nspline,3);
    alloc_tdi(tdi_amp,Nspline,3);

    //extract remaining extrinsic parameters from UCB parameter vector
    double cosi  = params[4];
    double psi   = params[5]; 

    LISA_spline_response(orbit, time_ssb, Nspline, costh, phi, cosi, psi, amp_ssb_spline, NULL, phase_ssb_spline, phase_sc, tdi_amp, tdi_phase);

    /*
    Interpolate amplitude and phase for instrument response of each TDI channel onto wavelet grid
    */
    double Amp,Phase;
    struct TDI *wave = malloc(sizeof(struct TDI));
    alloc_tdi(wave,N_ds,3);

    struct CubicSpline *amp_interpolant   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_interpolant = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->X);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->X);

    for(int i=0; i<N_ds; i++)
    {
        Amp     = spline_interpolation(amp_interpolant, time_ds[i]);   //amplitude
        Phase   = spline_interpolation(phase_interpolant, time_ds[i]); //slow part of phase
        Phase  += phase_ds[i];                                         //carrier phase
        Phase  -= phase_het[i];                                        //remove heterodyne phase
        
        wave->X[i] = Amp*cos(Phase);
    }
    
    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->Y);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->Y);

    for(int i=0; i<N_ds; i++)
    {
        Amp    = spline_interpolation(amp_interpolant, time_ds[i]);
        Phase  = spline_interpolation(phase_interpolant, time_ds[i]);
        Phase += phase_ds[i];
        Phase -= phase_het[i];
        
        wave->Y[i] = Amp*cos(Phase);
    }
    
    initialize_cubic_spline(amp_interpolant,   time_ssb, tdi_amp->Z);
    initialize_cubic_spline(phase_interpolant, time_ssb, tdi_phase->Z);

    for(int i=0; i<N_ds; i++)
    {
        Amp    = spline_interpolation(amp_interpolant, time_ds[i]);
        Phase  = spline_interpolation(phase_interpolant, time_ds[i]);
        Phase += phase_ds[i];
        Phase -= phase_het[i];
        
        wave->Z[i] = Amp*cos(Phase);
    }
    
    /*
     Compute wavelet coefficients for signal's TDI response
     */
 
    // get freqeuncy wavelet window function for downsampled data
    double *window = double_vector((wdm->NT/2+1));
    wavelet_window_frequency(wdm, window, Nlayers);

    // wavelet transform on heterodyned data using downsampled windows.
    wavelet_transform_by_layers(wdm, min_layer, Nlayers, window, wave->X);
    wavelet_transform_by_layers(wdm, min_layer, Nlayers, window, wave->Y);
    wavelet_transform_by_layers(wdm, min_layer, Nlayers, window, wave->Z);

    /*
     Properly re-index to undo the heterodyning
    */
    int N=0;
    int k;

    for(int i=0; i<wdm->NT; i++)
    {
        for(int j=min_layer; j<min_layer+Nlayers; j++)
        {
            wavelet_pixel_to_index(wdm,i,j,&k);
            
            //check that the pixel is in range
            if(k>=wdm->kmin && k<wdm->kmax)
            {
                wavelet_list[N]=k-wdm->kmin;
                N++;
            }
        }
    }
    *Nwavelet = N;
        
    //insert non-zero wavelet pixels into correct indicies
    for(int n=0; n<*Nwavelet; n++)
    {
        X[wavelet_list[n]] = wave->X[n];
        Y[wavelet_list[n]] = wave->Y[n];
        Z[wavelet_list[n]] = wave->Z[n];
    }
    

    free_double_vector(time_ssb);
    free_double_vector(amp_ssb);
    free_double_vector(phase_ssb);

    free_cubic_spline(amp_ssb_spline);
    free_cubic_spline(phase_ssb_spline);

    free_double_vector(phase_sc);

    free_double_vector(phase_ds);
    free_double_vector(time_ds);
    free_double_vector(phase_het);

    free_tdi(tdi_phase);
    free_tdi(tdi_amp);
    free_tdi(wave);

    free_cubic_spline(amp_interpolant);
    free_cubic_spline(phase_interpolant);
    
    free_double_vector(window);
}

/* Lookup table wavelet transform */
void ucb_waveform_wavelet_tab(struct Orbit *orbit, struct Wavelets *wdm, double Tobs, double t0, double *params, int *wavelet_list, int *Nwavelet, double *X, double *Y, double *Z)
{
    /*
    Get waveform at solar system barycenter (SSB)
    */
    // each waveform type will want its own time spacing. For galactic binaries uniform spacing is fine
    // This section will need a flag to tell it what waveform type we are computing. For MBHMs the parameters
    // of the signal will impact the time spacing
    
    int Nspline = orbit->Norb;
    double dt = Tobs/(double)(Nspline-1);
            
    // get amplitude and phase at Barycenter on the orbit interpolation grid (with margin)...
    double *t         = malloc(sizeof(double)*Nspline);
    double *amp_ssb   = malloc(sizeof(double)*Nspline);
    double *phase_ssb = malloc(sizeof(double)*Nspline);
    
    // convert parameters
    params[0] = params[0]/Tobs;
    params[3] = exp(params[3]);
    params[7] = params[7]/(Tobs*Tobs);
    //params[8] = params[8]/(Tobs*Tobs*Tobs);
    
    ucb_barycenter_waveform(params, Nspline, orbit->t, phase_ssb, amp_ssb, Tobs);
    
    // convert back
    params[0] = params[0]*Tobs;
    params[3] = log(params[3]);
    params[7] = params[7]*(Tobs*Tobs);
    //params[8] = params[8]*(Tobs*Tobs*Tobs);
    
    /*
    Get spline interpolant for SSB phase and amplitude
    */
    struct CubicSpline *amp_ssb_spline   = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_ssb_spline = alloc_cubic_spline(Nspline);
    
    initialize_cubic_spline(amp_ssb_spline,orbit->t,amp_ssb);
    initialize_cubic_spline(phase_ssb_spline,orbit->t,phase_ssb);
    
    /*
    Interpolate phase at SSB now on the data's time grid.
    */
    for(int i=0; i< Nspline; i++)
    {
        t[i] = t0+(double)(i)*dt;
        phase_ssb[i] = spline_interpolation(phase_ssb_spline, t[i]);
    }

    /*
    Interpolate SSB phase, frequency, and frequency derivative on wdm time grid
    */
    double *time_wavelet_grid  = malloc(sizeof(double)*wdm->NT);
    double *phase_wavelet_grid = malloc(sizeof(double)*wdm->NT);
    double *freq_wavelet_grid  = malloc(sizeof(double)*wdm->NT);
    double *fdot_wavelet_grid  = malloc(sizeof(double)*wdm->NT);
    
    for(int i=0; i<wdm->NT; i++)
    {
        time_wavelet_grid[i]  = ((double)(i))*wdm->dt;  // time center of the wavelet pixels
        phase_wavelet_grid[i] = spline_interpolation(phase_ssb_spline, time_wavelet_grid[i]);
        freq_wavelet_grid[i]  = spline_interpolation_deriv(phase_ssb_spline, time_wavelet_grid[i])/PI2;
        fdot_wavelet_grid[i]  = spline_interpolation_deriv2(phase_ssb_spline, time_wavelet_grid[i])/PI2;
    }

    /*
    Get TDI response for signal's SSB phase and amplitude on spline grid
    */
    struct TDI *tdi_phase = malloc(sizeof(struct TDI));
    struct TDI *tdi_amp = malloc(sizeof(struct TDI));
    alloc_tdi(tdi_phase,Nspline,3);
    alloc_tdi(tdi_amp,Nspline,3);
    
    double costh = params[1]; 
    double phi   = params[2]; 
    double cosi  = params[4]; 
    double psi   = params[5]; 

    LISA_spline_response(orbit, t, Nspline, costh, phi, cosi, psi, amp_ssb_spline, NULL, phase_ssb_spline, phase_ssb, tdi_amp, tdi_phase);

    /*
    Interpolate amplitude and phase for instrument response of each TDI channel onto wavelet grid
    */
    struct TDI *phase = malloc(sizeof(struct TDI));
    struct TDI *freq  = malloc(sizeof(struct TDI));
    struct TDI *fdot  = malloc(sizeof(struct TDI));
    struct TDI *amp   = malloc(sizeof(struct TDI));

    alloc_tdi(phase,wdm->NT,3);
    alloc_tdi(freq, wdm->NT,3);
    alloc_tdi(fdot, wdm->NT,3);
    alloc_tdi(amp,  wdm->NT,3);

    struct CubicSpline *amp_interpolant = alloc_cubic_spline(Nspline);
    struct CubicSpline *phase_interpolant = alloc_cubic_spline(Nspline);
       
    initialize_cubic_spline(amp_interpolant,t,tdi_amp->X);
    initialize_cubic_spline(phase_interpolant,t,tdi_phase->X);

    for(int i=0; i< wdm->NT; i++)
    {
        amp->X[i]   = spline_interpolation(amp_interpolant, time_wavelet_grid[i]);
        phase->X[i] = spline_interpolation(phase_interpolant, time_wavelet_grid[i]) + phase_wavelet_grid[i];
        freq->X[i]  = spline_interpolation_deriv(phase_interpolant, time_wavelet_grid[i])/PI2 + freq_wavelet_grid[i];
        fdot->X[i]  = spline_interpolation_deriv2(phase_interpolant, time_wavelet_grid[i])/PI2 + fdot_wavelet_grid[i];
    }
    
    initialize_cubic_spline(amp_interpolant,t,tdi_amp->Y);
    initialize_cubic_spline(phase_interpolant,t,tdi_phase->Y);

    for(int i=0; i< wdm->NT; i++)
    {
        amp->Y[i]   = spline_interpolation(amp_interpolant, time_wavelet_grid[i]);
        phase->Y[i] = spline_interpolation(phase_interpolant, time_wavelet_grid[i])+phase_wavelet_grid[i];
        freq->Y[i]  = spline_interpolation_deriv(phase_interpolant, time_wavelet_grid[i])/PI2 + freq_wavelet_grid[i];
        fdot->Y[i]  = spline_interpolation_deriv2(phase_interpolant, time_wavelet_grid[i])/PI2 + fdot_wavelet_grid[i];
    }
    
    initialize_cubic_spline(amp_interpolant,t,tdi_amp->Z);
    initialize_cubic_spline(phase_interpolant,t,tdi_phase->Z);

    for(int i=0; i< wdm->NT; i++)
    {
        amp->Z[i]   = spline_interpolation(amp_interpolant, time_wavelet_grid[i]);
        phase->Z[i] = spline_interpolation(phase_interpolant, time_wavelet_grid[i])+phase_wavelet_grid[i];
        freq->Z[i]  = spline_interpolation_deriv(phase_interpolant, time_wavelet_grid[i])/PI2 + freq_wavelet_grid[i];
        fdot->Z[i]  = spline_interpolation_deriv2(phase_interpolant, time_wavelet_grid[i])/PI2 + fdot_wavelet_grid[i];
    }

    /*
    Wavelet transform of interpolated TDI channels
    */
    
    //minimum and maximum frequency layers for each TDI response
    int *min_layer = malloc(sizeof(int)*(wdm->NT));
    int *max_layer = malloc(sizeof(int)*(wdm->NT));

    int *reverse_list = malloc(sizeof(int)*(wdm->NF*wdm->NT));

    //get list of non-zero wavelet amplitudes for this signal
    active_wavelet_list(wdm, freq->X, freq->Y, freq->Z, fdot->X, fdot->Y, fdot->Z, wavelet_list, reverse_list, Nwavelet, min_layer, max_layer);

    //finally compute wavelet coefficients for signal's TDI response
    double *Xtemp = double_vector(*Nwavelet);
    double *Ytemp = double_vector(*Nwavelet);
    double *Ztemp = double_vector(*Nwavelet);

    wavelet_transform_from_table(wdm, phase->X, freq->X, fdot->X, amp->X, min_layer, max_layer, Xtemp, wavelet_list, reverse_list, *Nwavelet);
    wavelet_transform_from_table(wdm, phase->Y, freq->Y, fdot->Y, amp->Y, min_layer, max_layer, Ytemp, wavelet_list, reverse_list, *Nwavelet);
    wavelet_transform_from_table(wdm, phase->Z, freq->Z, fdot->Z, amp->Z, min_layer, max_layer, Ztemp, wavelet_list, reverse_list, *Nwavelet);


    //insert non-zero wavelet pixels into correct indicies
    for(int n=0; n<*Nwavelet; n++)
    {
        X[wavelet_list[n]] = Xtemp[n];
        Y[wavelet_list[n]] = Ytemp[n];
        Z[wavelet_list[n]] = Ztemp[n];
    }

    free_double_vector(Xtemp);
    free_double_vector(Ytemp);
    free_double_vector(Ztemp);


    free(t);
    free(amp_ssb);
    free(phase_ssb);

    free_cubic_spline(amp_ssb_spline);
    free_cubic_spline(phase_ssb_spline);

    free(time_wavelet_grid);
    free(phase_wavelet_grid);
    free(freq_wavelet_grid);
    free(fdot_wavelet_grid);

    free_tdi(tdi_phase);
    free_tdi(tdi_amp);

    free_tdi(phase);
    free_tdi(freq);
    free_tdi(fdot);
    free_tdi(amp);

    free_cubic_spline(amp_interpolant);
    free_cubic_spline(phase_interpolant);

    free(min_layer);
    free(max_layer);

    free(reverse_list);

}
