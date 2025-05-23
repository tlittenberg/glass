/*
 * Copyright 2019 Tyson B. Littenberg, Kristen Lackeos & Neil J. Cornish
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
#include "glass_ucb_io.h"
#include "glass_ucb_catalog.h"
#include "glass_ucb_prior.h"
#include "glass_ucb_waveform.h"


void alloc_entry(struct Entry *entry, int IMAX)
{
    entry->Nchain = 0;
    entry->source = malloc(IMAX*sizeof(struct Source*));
    entry->match  = malloc(IMAX*sizeof(double));
    entry->distance  = malloc(IMAX*sizeof(double));
    entry->stepFlag = calloc(IMAX,sizeof(int));
    entry->gmm = malloc(sizeof(struct GMM));
}

void free_entry(struct Entry *entry, int IMAX)
{
    for(size_t n=0; n<IMAX; n++) free_source(entry->source[n]);
    free(entry->source);
    free(entry->match);
    free(entry->distance);
    free(entry->stepFlag);
    for(size_t n=0; n<entry->gmm->NMODE; n++) free_MVG(entry->gmm->modes[n]);
    free(entry->gmm->modes);
    free(entry->gmm);
}

void create_empty_source(struct Catalog *catalog, int NFFT, int Nchannel)
{
    int N = catalog->N;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,1);
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->Nchain], NFFT, Nchannel);

    entry->match[entry->Nchain] = 1.0;
    entry->distance[entry->Nchain] = 0.0;
    
    entry->Nchain++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void create_new_source(struct Catalog *catalog, struct Source *sample, struct Noise *noise, int i, int IMAX, int NFFT, int Nchannel)
{
    int N = catalog->N;
    
    //allocate memory for new entry in catalog
    catalog->entry[N] = malloc(sizeof(struct Entry));
    struct Entry *entry = catalog->entry[N];
    
    alloc_entry(entry,IMAX);
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    alloc_source(entry->source[entry->Nchain], NFFT, Nchannel);
    
    //add sample to the catalog as the new entry
    copy_source(sample, entry->source[entry->Nchain]);
    
    //store SNR of reference sample to set match criteria
    entry->SNR = snr(sample,noise);
    
    entry->match[entry->Nchain] = 1.0;
    entry->distance[entry->Nchain] = 0.0;
    entry->stepFlag[i] = 1;
    
    entry->Nchain++; //increment number of samples for entry
    catalog->N++;//increment number of entries for catalog
}

void append_sample_to_entry(struct Entry *entry, struct Source *sample, int IMAX, int NFFT, int Nchannel)
{
    //malloc source structure for this sample's entry
    entry->source[entry->Nchain] = malloc(sizeof(struct Source));
    /* This is taking up too much memory
    alloc_source(entry->source[entry->Nchain], NFFT, Nchannel, NP);
    
    //copy chain sample into catalog entry
    copy_source(sample, entry->source[entry->Nchain]);
    */

    /* leaner way of storing source info */
    
    //only copy source parameters
    entry->source[entry->Nchain]->params=calloc(UCB_MODEL_NP,sizeof(double));
    memcpy(entry->source[entry->Nchain]->params, sample->params, UCB_MODEL_NP*sizeof(double));
    
    //need Tobs which isn't stored in entries
    double T = sample->params[0]/sample->f0;
    
    //get physical parameters for waveform calculations later
    map_array_to_params(entry->source[entry->Nchain], entry->source[entry->Nchain]->params, T);

    //increment number of stored samples for this entry
    entry->Nchain++;
}

void get_correlation_matrix(struct Data *data, struct Catalog *catalog, int *detection_index, int detections, int IMAX, double **corr)
{
    struct Entry *entry = NULL;
    
    /*
     compute mean and variance for each parameter
     */
    double **mean = malloc(detections*sizeof(double *));
    double **var = malloc(detections*sizeof(double *));
    
    for(int d=0; d<detections; d++)
    {
        mean[d] = calloc(UCB_MODEL_NP,sizeof(double));
        var[d] = calloc(UCB_MODEL_NP,sizeof(double));
        
        entry = catalog->entry[detection_index[d]];

        for(int n=0; n<UCB_MODEL_NP; n++)
        {
            double x;
            
            for(int i=0; i<entry->Nchain; i++)
            {
                x = entry->source[i]->params[n];
                mean[d][n] += x;
            }
            mean[d][n] /= (double)entry->Nchain;

            /*
             computing variance after mean (instead of using 1-pass method)
             because rounding error was non-negligible for frequency parameter
             std << mu
             */
            for(int i=0; i<entry->Nchain; i++)
            {
                x = entry->source[i]->params[n];
                var[d][n] += (x - mean[d][n])*(x - mean[d][n]);
            }
            
            var[d][n] /=(double)(entry->Nchain);
        }
    }
    
    /*
     compute correlation matrix
     */
    int N = detections*UCB_MODEL_NP;
    struct Entry *n_entry=NULL;
    struct Entry *m_entry=NULL;
    
    for(int n=0; n<N; n++)
    {
        for(int m=0; m<N; m++)
        {
            //which source row?
            int nd = n/UCB_MODEL_NP;

            //which source column?
            int md = m/UCB_MODEL_NP;

            //which parameter row?
            int nx = n - nd*UCB_MODEL_NP;
            
            //which parameter column?
            int mx = m - md*UCB_MODEL_NP;
                        
            //which entries?
            n_entry = catalog->entry[detection_index[nd]];
            m_entry = catalog->entry[detection_index[md]];
            

            /*
             for each element of the correlation matrix,
             sum over non-null chain samples
             */
            int corr_count=0;
            int ncount=0;
            int mcount=0;
                                     
            for(int i=0; i<IMAX; i++)
            {
                //is this a valid pairing?
                if(n_entry->stepFlag[i]*m_entry->stepFlag[i])
                {
                    double X = n_entry->source[ncount]->params[nx];
                    double Y = m_entry->source[mcount]->params[mx];
                    corr[n][m] += (X - mean[nd][nx]) * (Y - mean[md][mx]);
                    corr_count++;
                    
                }
                
                //advance n-counter
                if(n_entry->stepFlag[i]) ncount++;
                
                //advance m-counter
                if(m_entry->stepFlag[i]) mcount++;
                
            }
            corr[n][m] /= (double)corr_count;
            corr[n][m] /= sqrt(var[nd][nx]*var[md][mx]);;
        }
    }
}

int gaussian_mixture_model_wrapper(double **ranges, struct Flags *flags, struct Entry *entry, char *outdir, size_t NMODE, size_t NTHIN, unsigned int *seed, double *BIC)
{
    if(flags->verbose)fprintf(stdout,"Event %s, NMODE=%i\n",entry->name,(int)NMODE);
    
    // number of samples
    size_t NMCMC = entry->Nchain;
    
    // number of EM iterations
    size_t NSTEP = 100;
    
    // thin chain
    NMCMC /= NTHIN;
    
    struct Sample **samples = malloc(NMCMC*sizeof(struct Sample*));
    for(size_t n=0; n<NMCMC; n++)
    {
        samples[n] = malloc(sizeof(struct Sample));
        samples[n]->x = double_vector(UCB_MODEL_NP);
        samples[n]->p = double_vector(NMODE);
        samples[n]->w = double_vector(NMODE);
    }
    
    // covariance matrices for different modes
    struct MVG **modes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t n=0; n<NMODE; n++)
    {
        modes[n] = malloc(sizeof(struct MVG));
        alloc_MVG(modes[n],UCB_MODEL_NP);
    }
    
    // Logistic mapping of samples onto R
    double y;
    double pmin,pmax;
    double *y_vec = double_vector(NMCMC);
    
    /* parse chain file */
    double **params = malloc(UCB_MODEL_NP*sizeof(double *));
    for(size_t n=0; n<UCB_MODEL_NP; n++) params[n] = double_vector(NMCMC);
    double value[UCB_MODEL_NP];
    for(size_t i=0; i<NMCMC; i++)
    {
        value[0] = entry->source[i*NTHIN]->f0;
        value[1] = entry->source[i*NTHIN]->costheta;
        value[2] = entry->source[i*NTHIN]->phi;
        value[3] = log(entry->source[i*NTHIN]->amp);
        value[4] = entry->source[i*NTHIN]->cosi;
        value[5] = entry->source[i*NTHIN]->psi;
        value[6] = entry->source[i*NTHIN]->phi0;
        if(UCB_MODEL_NP>7)
            value[7] = entry->source[i*NTHIN]->dfdt;
        if(UCB_MODEL_NP>8)
            value[8] = entry->source[i*NTHIN]->d2fdt2;
        
        for(size_t n=0; n<UCB_MODEL_NP; n++)
        {
            params[n][i] = value[n];
        }
    }
    
    /* Use priors to set min and max of each parameter*/
    for(size_t n=0; n<UCB_MODEL_NP; n++)
    {
        // copy max and min into each MVG structure
        for(size_t k=0; k<NMODE; k++)
        {
            modes[k]->minmax[n][0] = ranges[n][0];
            modes[k]->minmax[n][1] = ranges[n][1];
        }
    }
    
    
    /* map params to R with logit function */
    for(size_t n=0; n<UCB_MODEL_NP; n++)
    {
        pmin = modes[0]->minmax[n][0];;
        pmax = modes[0]->minmax[n][1];
        logit_mapping(params[n], y_vec, pmin, pmax, NMCMC);
        
        for(size_t i=0; i<NMCMC; i++)
        {
            y = y_vec[i];
            samples[i]->x[n] = y;
        }
    }
    
    /* The main Gaussian Mixture Model with Expectation Maximization function */
    double logL;
    if(GMM_with_EM(modes,samples,UCB_MODEL_NP,NMODE,NMCMC,NSTEP,seed,&logL,BIC)) return 1;
    
    /* Write GMM results to binary for pick up by other processes */
    char filename[BUFFER_SIZE];
    sprintf(filename,"%s/%s_gmm.bin",outdir,entry->name);
    FILE *fptr = fopen(filename,"wb");
    fwrite(&NMODE, sizeof NMODE, 1, fptr);
    for(size_t n=0; n<NMODE; n++) write_MVG(modes[n],fptr);
    fclose(fptr);
    
    /* print 1D PDFs and 2D contours of GMM model */
    if(flags->verbose) print_model(modes, samples, UCB_MODEL_NP, NMODE, NMCMC, logL, *BIC, NMODE);
    
    /* clean up */
    for(size_t n=0; n<NMCMC; n++)
    {
        free_double_vector(samples[n]->x);
        free_double_vector(samples[n]->p);
        free_double_vector(samples[n]->w);
        free(samples[n]);
    }
    free(samples);
    
    // covariance matrices for different modes
    for(size_t n=0; n<NMODE; n++) free_MVG(modes[n]);
    free(modes);
    
    return 0;
}

