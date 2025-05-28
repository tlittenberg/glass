/*
 * Copyright 2019 Tyson B. Littenberg & Kristen Lackeos
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


/***************************  REQUIRED LIBRARIES  ***************************/

/**
 @file ucb_match.c
 \brief App for computing matches between UCB catalogs
 */

#include <glass_utils.h>
#include <glass_ucb.h>

// CODE USAGE:
// ucb_match --match-in1 /path/to/input1.dat --match-in2 /path/to/input2.dat --frac-freq --fmin 0.001249 --samples 512 --duration 62914560

static void print_usage()
{
    print_glass_usage();
    print_ucb_usage();
    exit(0);
}

/* ============================  MAIN PROGRAM  ============================ */

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================= UCB MATCH =================\n");

    if(argc==1) print_usage();

    FILE *chain_file1;
    FILE *chain_file2;
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    struct Data  *data  = malloc(sizeof(struct Data));
    
    
    //   Parse command line and set defaults/flags
    parse_data_args(argc,argv,data,orbit,flags,chain,"fourier");
    parse_ucb_args(argc,argv,flags);
    if(flags->help)print_usage();
    alloc_data(data, flags);
    data->qmin = (int)(data->fmin*data->T);
    
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
            return(1);
            break;
    }
    
    chain_file1 = fopen(flags->matchInfile1,"r");
    chain_file2 = fopen(flags->matchInfile2,"r");
    
    if ( chain_file1 == NULL )
    {
        printf("match-in1 is null\n");
    }
    
    if ( chain_file2 == NULL )
    {
        printf("match-in2 is null\n");
    }
    
    
    //allocate memory for two sources and noise
    struct Source *src1 = malloc(sizeof(struct Source));
    alloc_source(src1, data->N, 2);
    
    struct Source *src2 = malloc(sizeof(struct Source));
    alloc_source(src2, data->N, 2);
    
    struct Noise *noise = malloc(sizeof(struct Noise));
    alloc_noise(noise, data->NFFT, 1, 2);
    
    
    for(int n=0; n<data->N; n++)
    {
        src1->tdi->A[n] = 0.0;
        src1->tdi->E[n] = 0.0;
        src1->tdi->X[n] = 0.0;
        src2->tdi->A[n] = 0.0;
        src2->tdi->E[n] = 0.0;
        src2->tdi->X[n] = 0.0;
    }
    
    //Get noise spectrum for data segment
    for(int n=0; n<data->NFFT; n++)
    {
        double Spm, Sop;
        double f = data->fmin + (double)(n)/data->T;
        get_noise_levels("radler", f, &Spm, &Sop);

        noise->f[n] = f;
        if(strcmp(data->format,"phase")==0)
        {
            noise->C[0][0][n] = AEnoise(orbit->L, orbit->fstar, f);
            noise->C[1][1][n] = AEnoise(orbit->L, orbit->fstar, f);
        }
        else if(strcmp(data->format,"frequency")==0 || strcmp(data->format,"sangria")==0)
        {
            noise->C[0][0][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop)/sqrt(2.);
            noise->C[1][1][n] = AEnoise_FF(orbit->L, orbit->fstar, f, Spm, Sop)/sqrt(2.);
            noise->C[0][0][n] += GBnoise_FF(data->T, orbit->fstar, f)/sqrt(2.);
            noise->C[1][1][n] += GBnoise_FF(data->T, orbit->fstar, f)/sqrt(2.);

        }
        else
        {
            fprintf(stderr,"Unsupported data format %s",data->format);
            exit(1);
        }
    }
    invert_noise_covariance_matrix(noise);

    
    double max_match;
    double match;
    
    while(!feof(chain_file2))
    {
        for(int n=0; n<data->N; n++)
        {
            src2->tdi->A[n] = 0.0;
            src2->tdi->E[n] = 0.0;
        }
        
        scan_source_params(data, src2, chain_file2);
        ucb_alignment(orbit, data, src2);
        ucb_waveform(orbit, data->format, data->T, data->t0, src2->params, UCB_MODEL_NP, src2->tdi->X, src2->tdi->Y, src2->tdi->Z, src2->tdi->A, src2->tdi->E, src2->BW, 2);

        max_match=-INFINITY;
        while(!feof(chain_file1))
        {

            scan_source_params(data, src1, chain_file1);
            
            if( fabs(src1->params[0] - src2->params[0]) < 20.)
            {
                for(int n=0; n<data->N; n++)
                {
                    src1->tdi->A[n] = 0.0;
                    src1->tdi->E[n] = 0.0;
                }
                
                //Book-keeping of injection time-frequency volume
                ucb_alignment(orbit, data, src1);
                ucb_waveform(orbit, data->format, data->T, data->t0, src1->params, UCB_MODEL_NP, src1->tdi->X, src1->tdi->Y,src1->tdi->Z, src1->tdi->A, src1->tdi->E, src1->BW, 2);
                                
                match = waveform_match(src1, src2, noise);
                if(match>max_match) max_match=match;
            }
        }
        printf("%lg %lg %lg\n",max_match,snr(src2,noise),src2->f0);
        rewind(chain_file1);
    }

    fclose(chain_file1);
    fclose(chain_file2);
    
    return 0;
}
