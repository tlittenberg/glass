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
 @file vgb_mcmc.c
 \brief Main function for targeted binary sampler app `vgb_mcmc`
 */

/*  REQUIRED LIBRARIES  */

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_ucb.h>

static void print_usage()
{
    print_glass_usage();
    print_ucb_usage();

    fprintf(stdout,"EXAMPLE:\n");
    fprintf(stdout,"vgb_mcmc --known-sources /path/to/full_list.txt --quiet \n");
    fprintf(stdout,"\n");

    exit(0);
}

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================== VGB MCMC =================\n");

    time_t start, stop;
    start = time(NULL);
    
    char filename[MAXSTRINGSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    /* Allocate data structures */
    struct Flags *flags = malloc(sizeof(struct Flags));
    struct Orbit *orbit = malloc(sizeof(struct Orbit));

    /* get vgbmcmc-specific arguments before allocating Data structure*/
    parse_vgb_args(argc,argv,flags);
    if(flags->NVB==0)
    {
        fprintf(stdout, "ERROR: Verification binary list required\n");
        print_usage();
    }
    else
    {
        fprintf(stdout, "FOUND: %i binaries in %s\n", flags->NVB, flags->vbFile);
    }

    struct Data *data;
    struct Chain *chain;
    struct Source *inj;
    struct Data  **data_vec = malloc(flags->NVB*sizeof(struct Data*));
    struct Chain **chain_vec = malloc(flags->NVB*sizeof(struct Chain*));
    struct Source **inj_vec = malloc(flags->NVB*sizeof(struct Source*));

    for(int n=0; n<flags->NVB; n++)
    {
        data_vec[n] = malloc(sizeof(struct Data));
        chain_vec[n] = malloc(sizeof(struct Chain));
        inj_vec[n] = malloc(sizeof(struct Source));
    }
        
    data=data_vec[0];
    chain=chain_vec[0];
    parse_ucb_args(argc,argv,flags);
    parse_data_args(argc,argv,data,orbit,flags,chain,"fourier");
    if(flags->help) print_usage();
    
    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /*
     * Set custom flags for verification binary analysis
     *
     * We are using the injection infrastructure to keep track
     * of the known binary parameters
     */
    flags->knownSource = 1;
    flags->snrPrior    = 0;
    flags->fixSky      = 1;
    flags->fixFreq     = 1;
    flags->cheat       = 1; //initializes chain at injection values
    flags->NINJ        = flags->NVB;
    
    /* parse verification binary files */
    FILE *vbFile = fopen(flags->vbFile,"r");
    
    //strip off header
    char header[MAXSTRINGSIZE];
    if(fgets(header, MAXSTRINGSIZE, vbFile)==NULL)
    {
        fprintf(stderr,"Error reading %s\n",flags->vbFile);
        exit(1);
    }

    /* Initialize LISA orbit model */
    initialize_orbit(data, orbit, flags);

    
    /* initialize all data structures */

    /* Setup output directories for data and chain files */
    mkdir(flags->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    for(int n=0; n<flags->NVB; n++)
    {
        chain=chain_vec[n];
        data=data_vec[n];
        inj=inj_vec[n];

        if(n>0)
        {
            copy_data(data_vec[0],data);
            chain->NC = chain_vec[0]->NC;    //number of chains
        }

        alloc_source(inj, data->N, data->Nchannel);

        data->nseed+=n;

        char subDir[MAXSTRINGSIZE];
        sprintf(subDir,"%s/seg%02d",flags->runDir,n);
        mkdir(subDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        sprintf(data->dataDir,"%s/data",subDir);
        sprintf(chain->chainDir,"%s/chains",subDir);
        sprintf(chain->chkptDir,"%s/checkpoint",subDir);
        
        mkdir(data->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(chain->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        mkdir(chain->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

        
        /* Initialize data structures */
        alloc_data(data, flags);
        
        /* Get source from verification binary file */
        GetVerificationBinary(data, flags, inj, vbFile);
        
        /* Read strain data */
        if(flags->hdf5Data)
            ReadData(data,orbit,flags);

        /* Simulate strain data */
        else
            UCBInjectVerificationSet(data, orbit, flags, inj);
        
        /* Get noise model */
        GetNoiseModel(data,orbit,flags);

        /* Add Gaussian noise realization */
        if(flags->simNoise) AddNoise(data,data->tdi);

        /* set approximate f/fstar for segment */
        data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);

        //print various data products for plotting
        print_data(data, flags);
        
        //save parameters to file
        sprintf(filename,"%s/data/injection_parameters.dat",subDir);
        FILE *paramFile=fopen(filename,"w");
        print_source_params(data, inj, paramFile);
        fprintf(paramFile,"\n");
        fclose(paramFile);
        
        /* Initialize parallel chain */
        if(flags->resume)
            initialize_chain(chain, flags, &data->cseed, "a");
        else
            initialize_chain(chain, flags, &data->cseed, "w");
        
    }

    /* Setup the rest of the model */
    struct Prior *prior = NULL;
    struct Proposal **proposal = NULL;
    struct Model **trial = NULL;
    struct Model **model = NULL;
    struct Prior **prior_vec = malloc(flags->NVB*sizeof(struct Prior *));
    struct Proposal ***proposal_vec = malloc(flags->NVB*sizeof(struct Proposal **));
    struct Model ***trial_vec = malloc(flags->NVB*sizeof(struct Model**));
    struct Model ***model_vec = malloc(flags->NVB*sizeof(struct Model**));

    for(int n=0; n<flags->NVB; n++)
    {
        /* Initialize priors */
        prior_vec[n] = malloc(sizeof(struct Prior));
        
        /* Initialize MCMC proposals */
        proposal_vec[n] = malloc(UCB_PROPOSAL_NPROP*sizeof(struct Proposal*));
        initialize_vb_proposal(orbit, data_vec[n], prior_vec[n], chain_vec[n], flags, proposal_vec[n], DMAX);
        
        /* Initialize data models */
        trial_vec[n] = malloc(sizeof(struct Model*)*NC);
        model_vec[n] = malloc(sizeof(struct Model*)*NC);
        initialize_ucb_state(data_vec[n], orbit, flags, chain_vec[n], proposal_vec[n], model_vec[n], trial_vec[n], inj_vec);
    }
    
    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        fprintf(stdout,"\n=============== Checkpointing ===============\n");
        for(int n=0; n<flags->NVB; n++)
        {
            //check for files needed to resume
            FILE *fptr = NULL;
            int file_error = 0;
            
            for(int ic=0; ic<chain_vec[n]->NC; ic++)
            {
                sprintf(filename,"%s/chain_state_%i.dat",chain_vec[n]->chkptDir,ic);
                
                if( (fptr = fopen(filename,"r")) == NULL )
                {
                    fprintf(stderr,"Warning: Could not checkpoint run state\n");
                    fprintf(stderr,"         Parameter file %s does not exist\n",filename);
                    file_error++;
                    break;
                }
            }
            
            //if all of the files exist resume run from checkpointed state
            if(!file_error)
            {
                fprintf(stdout,"   Checkpoint files found. Resuming chain\n");
                restore_chain_state(orbit, data_vec[n], model_vec[n], chain_vec[n], flags, &mcmc_start);
            }
        }
        fprintf(stdout,"============================================\n\n");
    }
            
    /* Write example gb_catalog bash script in run directory */
    print_ucb_catalog_script(flags, data_vec[0], orbit);
    
    //For saving the number of threads actually given
    int numThreads;
    int mcmc = mcmc_start;
    #pragma omp parallel num_threads(flags->threads)
    {
        int threadID;
        //Save individual thread number
        threadID = omp_get_thread_num();
        
        //Only one thread runs this section
        if(threadID==0)  numThreads = omp_get_num_threads();
        
        #pragma omp barrier
        
        /* The MCMC loop */
        for(; mcmc < flags->NMCMC;)
        {
            
            if(threadID==0)
            {
                flags->burnin   = (mcmc<0) ? 1 : 0;
                flags->maximize = (mcmc<-flags->NBURN/2) ? 1 : 0;
            }
            
            #pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                //loop over verification binary segments
                for(int n=0; n<flags->NVB; n++)
                {
                    
                    struct Model *model_ptr = model_vec[n][chain_vec[n]->index[ic]];
                    struct Model *trial_ptr = trial_vec[n][chain_vec[n]->index[ic]];
                    
                    for(int steps=0; steps < 100; steps++)
                    {
                        ucb_mcmc(orbit, data_vec[n], model_ptr, trial_ptr, chain_vec[n], flags, prior_vec[n], proposal_vec[n], ic);
                    }//loop over MCMC steps
                    
                    //update fisher matrix for each chain
                    if(mcmc%100==0)
                    {
                        for(int i=0; i<model_ptr->Nlive; i++)
                        {
                            ucb_fisher(orbit, data_vec[n], model_ptr->source[i], data_vec[n]->noise);
                        }
                    }
                }
                
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            #pragma omp barrier
            if(threadID==0){
                
                for(int n=0; n<flags->NVB; n++)
                {
                    model = model_vec[n];
                    trial = trial_vec[n];
                    data = data_vec[n];
                    prior = prior_vec[n];
                    proposal = proposal_vec[n];
                    chain = chain_vec[n];
                    
                    ptmcmc(model,chain,flags);
                    adapt_temperature_ladder(chain, mcmc+flags->NBURN);
                    
                    print_chain_files(data, model, chain, flags, mcmc);
                    
                    //track maximum log Likelihood
                    if(mcmc%100)
                    {
                        if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
                    }

                    //store reconstructed waveform
                    if(!flags->quiet) print_waveform_draw(data, model[chain->index[0]], flags);

                    //update run status
                    if(mcmc%data->downsample==0)
                    {

                        if(!flags->quiet)
                        {
                            print_chain_state(data, chain, model[chain->index[0]], flags, stdout, mcmc); //writing to file
                            fprintf(stdout,"Sources: %i\n",model[chain->index[0]]->Nlive);
                            print_acceptance_rates(proposal, UCB_PROPOSAL_NPROP, 0, stdout);
                        }

                        //save chain state to resume sampler
                        save_chain_state(data, model, chain, flags, mcmc);

                    }

                    //dump waveforms to file, update avgLogL for thermodynamic integration
                    if(mcmc>0 && mcmc%data->downsample==0)
                    {
                        save_waveforms(data, model[chain->index[0]], mcmc/data->downsample);

                        for(int ic=0; ic<NC; ic++)
                        {
                            chain->dimension[ic][model[chain->index[ic]]->Nlive]++;
                            chain->avgLogL[ic] += model[chain->index[ic]]->logL + model[chain->index[ic]]->logLnorm;
                        }
                    }
                }
            mcmc++;
            }
            //Can't continue MCMC until single thread is finished
#pragma omp barrier
            
        }// end MCMC loop
        
    }// End of parallelization
    
    //print aggregate run files/results
    for(int n=0; n<flags->NVB; n++)
        print_waveforms_reconstruction(data_vec[n],flags);    
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    sprintf(filename,"%s/vb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
    
    
    return 0;
}



