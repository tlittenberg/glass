/*
 * Copyright 2021 Tyson B. Littenberg
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
 @file noise_spline_mcmc.c
 \brief Main function for stand-alone Noise spline model sampler 
 */

#include <glass_utils.h>
#include <glass_noise.h>

static void print_usage()
{
    print_glass_usage();
    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"noise_spline_mcmc --sim-noise --conf-noise --duration 7864320 --fmin 1e-4 --fmax 8e-3\n");
    fprintf(stdout,"\n");
    exit(0);
}

int main(int argc, char *argv[])
{
    fprintf(stdout, "\n============= NOISE SPLINE MCMC =============\n");

    time_t start, stop;
    start = time(NULL);
    char filename[128];
    
    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    /* Allocate data structures */
    struct Data *data   = malloc(sizeof(struct Data));
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    struct Chain *chain = malloc(sizeof(struct Chain));
    
    parse_data_args(argc,argv,data,orbit,flags,chain,"fourier");
    if(flags->help)print_usage();
    
    /*
     * Get Data
     */
    
    /* Setup output directories for data and chain files */
    sprintf(data->dataDir,"%s/data",flags->runDir);
    sprintf(chain->chainDir,"%s/chains",flags->runDir);
    sprintf(chain->chkptDir,"%s/checkpoint",flags->runDir);
    
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize data structures */
    alloc_data(data, flags);
    
    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);
    
    /* Initialize chain structure and files */
    initialize_chain(chain, flags, &data->cseed, "a");
    
    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    else if (flags->simNoise)
        SimulateData(data, orbit, flags);
    
    /*
     * Initialize Spline Model
     */
    int Nspline = 32+1;
    struct SplineModel **model = malloc(chain->NC*sizeof(struct SplineModel *));
    for(int ic=0; ic<chain->NC; ic++)
    {
        model[ic] = malloc(sizeof(struct SplineModel));
        initialize_spline_model(orbit, data, model[ic], Nspline);
    }
    
    sprintf(filename,"%s/initial_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->spline, filename);
    
    sprintf(filename,"%s/interpolated_spline_points.dat",data->dataDir);
    print_noise_model(model[0]->psd, filename);
    
    
    //MCMC
    sprintf(filename,"%s/chain_file.dat",chain->chainDir);
    FILE *chainFile = fopen(filename,"w");

    int numThreads;
    int step = 0;
    int NC = chain->NC;
    
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; step<flags->NMCMC;)
        {
            
            #pragma omp barrier
            
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                struct SplineModel *model_ptr = model[chain->index[ic]];
                for(int mc=0; mc<10; mc++)
                {
                    
                    if(rand_r_U_0_1(&chain->r[ic])<0.9)
                        noise_spline_model_mcmc(orbit, data, model_ptr, chain, flags, ic);
                    else
                        noise_spline_model_rjmcmc(orbit, data, model_ptr, chain, flags, ic);
                }
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            
            #pragma omp barrier
            
            if(threadID==0)
            {
                spline_ptmcmc(model, chain, flags);
                
                if(step%(flags->NMCMC/10)==0)printf("noise_spline_mcmc at step %i\n",step);
                
                if(step%(flags->NMCMC/10)==0)
                {
                    print_spline_state(model[chain->index[0]], chainFile, step);
                    
                    sprintf(filename,"%s/current_interpolated_spline_points.dat",data->dataDir);
                    print_noise_model(model[chain->index[0]]->psd, filename);
                    
                    sprintf(filename,"%s/current_spline_points.dat",data->dataDir);
                    print_noise_model(model[chain->index[0]]->spline, filename);


                }
                
                if(step%data->downsample==0 && step/data->downsample < data->Nwave)
                {
                    for(int n=0; n<data->N; n++)
                        for(int i=0; i<data->Nchannel; i++)
                            data->S_pow[n][i][step/data->downsample] = model[chain->index[0]]->psd->C[i][i][n];
                }

                step++;
                
                
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end of MCMC loop
        
    }// End of parallelization
    
    fclose(chainFile);
    
    sprintf(filename,"%s/final_spline_points.dat",data->dataDir);
    print_noise_model(model[chain->index[0]]->spline, filename);
    
    sprintf(filename,"%s/final_interpolated_spline_points.dat",data->dataDir);
    print_noise_model(model[chain->index[0]]->psd, filename);
    
    print_noise_reconstruction(data, flags);

    sprintf(filename,"%s/whitened_data.dat",data->dataDir);
    print_whitened_data(data, model[chain->index[0]]->psd, filename);

    for(int ic=0; ic<chain->NC; ic++) free_spline_model(model[ic]);
    free(model);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds\n",(double)(stop-start));
    
    
    return 0;
}

