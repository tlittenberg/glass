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
 @file ucb_mcmc.c
 \brief Main function for UCB sampler app `ucb_mcmc` 
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
    fprintf(stdout,"ucb_mcmc --inj [path to]/ldasoft/ucb/etc/sources/precision/PrecisionSource_0.txt --cheat\n");
    fprintf(stdout,"\n");

    exit(0);
}

/**
 * This is the main function
 *
 */
int main(int argc, char *argv[])
{
    fprintf(stdout, "\n================== UCB MCMC =================\n");

    time_t start, stop;
    start = time(NULL);
    char filename[MAXSTRINGSIZE];

    /* check arguments */
    print_LISA_ASCII_art(stdout);
    print_version(stdout);
    if(argc==1) print_usage();
    
    /* Allocate data structures */
    struct Flags  *flags = malloc(sizeof(struct Flags));
    struct Orbit  *orbit = malloc(sizeof(struct Orbit));
    struct Chain  *chain = malloc(sizeof(struct Chain));
    struct Data   *data  = malloc(sizeof(struct Data));
    
    /* Parse command line and set defaults/flags */
    sprintf(data->basis,"fourier");
    parse_data_args(argc,argv,data,orbit,flags,chain,"fourier");
    parse_ucb_args(argc,argv,flags);

    if(flags->help) print_usage();

    int NC = chain->NC;
    int DMAX = flags->DMAX;
    int mcmc_start = -flags->NBURN;
    
    /* Setup output directories for chain and data structures */
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

    /* read data */
    if(flags->strainData)
        ReadData(data,orbit,flags);
    
    /* noise model */
    GetNoiseModel(data,orbit,flags);

    struct Source **inj=NULL;
    if(flags->NINJ>0)
    {
        /* storage for injection parameters */
        inj = malloc(DMAX*sizeof(struct Source*));
        for(int n=0; n<DMAX; n++) inj[n] = malloc(sizeof(struct Source));
        
        /* Inject gravitational wave signal */
        if(flags->knownSource)
            UCBInjectVerificationSource(data,orbit,flags,inj[0]);
        else
            UCBInjectSimulatedSource(data,orbit,flags,inj);
        
        /* set approximate f/fstar for segment */
        data->sine_f_on_fstar = sin((data->fmin + (data->fmax-data->fmin)/2.)/orbit->fstar);
    }

    /* Add Gaussian noise realization */
    if(flags->simNoise) AddNoise(data,data->tdi);
    
    /* Store DFT copy of simulated data */
    copy_tdi(data->tdi, data->dft);

    /* print various data products for plotting */
    print_data(data, flags);
    
    /* Load catalog cache file for proposals/priors */
    struct Catalog *catalog=malloc(sizeof(struct Catalog));
    if(flags->catalog)
        UCBLoadCatalogCache(data, flags, catalog);

    
    /* Initialize data-dependent proposal */
    setup_frequency_proposal(data, flags);
    
    /* Initialize parallel chain */
    if(flags->resume)
        initialize_chain(chain, flags, &data->cseed, "a");
    else
        initialize_chain(chain, flags, &data->cseed, "w");
    
    /* Initialize priors */
    struct Prior *prior = malloc(sizeof(struct Prior));
    if(flags->galaxyPrior) set_galaxy_prior(flags, prior);
    if(flags->update) set_gmm_prior(flags, data, prior, catalog);

    /* Initialize MCMC proposals */
    struct Proposal **proposal = malloc(UCB_PROPOSAL_NPROP*sizeof(struct Proposal*));
    initialize_proposal(orbit, data, prior, chain, flags, catalog, proposal, DMAX);
    
    /* Test noise model */
    //test_noise_model(orbit);
    
    /* Initialize data models */
    struct Model **trial = malloc(sizeof(struct Model*)*NC);
    struct Model **model = malloc(sizeof(struct Model*)*NC);
    initialize_ucb_state(data, orbit, flags, chain, proposal, model, trial, inj);

    /* Start analysis from saved chain state */
    if(flags->resume)
    {
        fprintf(stdout,"\n=============== Checkpointing ===============\n");
        
        //check for files needed to resume
        FILE *fptr = NULL;
        int file_error = 0;
        
        for(int ic=0; ic<chain->NC; ic++)
        {
            sprintf(filename,"%s/checkpoint/chain_state_%i.dat",flags->runDir,ic);
            
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
            restore_chain_state(orbit, data, model, chain, flags, &mcmc_start);
        }
        fprintf(stdout,"============================================\n\n");
    }
    
    /*test proposals
     FILE *test=fopen("proposal_test.dat","w");
     for(int i=0; i<100000; i++)
     {
     double logP = draw_from_gmm_prior(data, model[0], model[0]->source[0], proposal[7], model[0]->source[0]->params, chain->r[0]);
     print_source_params(data, model[0]->source[0], test);
     fprintf(test,"%lg\n",logP);
     }
     fclose(test);*/
    //exit(1);
    
    //test covariance proposal
    if(flags->updateCov) test_covariance_proposal(data, flags, model[0], prior, proposal[8], &chain->r[0]);
    
    
    /* Write example gb_catalog bash script in run directory */
    print_ucb_catalog_script(flags, data, orbit);
    
    /* allow nested parallelization in mcmc loop (for rebuilding fstat proposal) */
    omp_set_max_active_levels(2);

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
                flags->maximize = 0;//(mcmc<-flags->NBURN/2) ? 1 : 0;
            }
            
            
            #pragma omp barrier
            // (parallel) loop over chains
            for(int ic=threadID; ic<NC; ic+=numThreads)
            {
                
                //loop over frequency segments
                struct Model *model_ptr = model[chain->index[ic]];
                struct Model *trial_ptr = trial[chain->index[ic]];
                copy_model(model_ptr,trial_ptr);

                for(int steps=0; steps < 500; steps++)
                {
                    //reverse jump birth/death or split/merge moves
                    if(rand_r_U_0_1(&chain->r[ic])<0.1 && flags->rj)
                    {
                        ucb_rjmcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                    }
                    //fixed dimension parameter updates
                    else
                    {
                            ucb_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, prior, proposal, ic);
                    }  

                    if( (flags->strainData || flags->simNoise) && !flags->psd)
                        noise_model_mcmc(orbit, data, model_ptr, trial_ptr, chain, flags, ic);

                }//loop over MCMC steps

                //update fisher matrix for each chain
                for(int n=0; n<model_ptr->Nlive; n++)
                    ucb_fisher(orbit, data, model_ptr->source[n], data->noise);
                
            }// end (parallel) loop over chains
            
            //Next section is single threaded. Every thread must get here before continuing
            #pragma omp barrier
            if(threadID==0){
                ptmcmc(model,chain,flags);
                adapt_temperature_ladder(chain, mcmc+flags->NBURN);
                
                print_chain_files(data, model, chain, flags, mcmc);
                
                //track maximum log Likelihood
                if(mcmc%100)
                {
                    if(update_max_log_likelihood(model, chain, flags)) mcmc = -flags->NBURN;
                }
                
                //store reconstructed waveform
                if(!flags->quiet) 
                {
                    print_waveform_draw(data, model[chain->index[0]], flags);
                    print_psd_draw(data, model[chain->index[0]], flags);
                }
                
                //update run status
                if(mcmc%data->downsample==0)
                {
                    
                    if(!flags->quiet)
                    {
                        print_chain_state(data, chain, model[chain->index[0]], flags, stdout, mcmc); //writing to file
                        fprintf(stdout,"Sources: %i/%i\n",model[chain->index[0]]->Nlive,model[chain->index[0]]->Neff-1);
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
                
                if(mcmc>-flags->NBURN+flags->NBURN/10. && model[0]->Neff < model[0]->Nmax && flags->rj)
                {
                    for(int ic=0; ic<NC; ic++) model[ic]->Neff++;
                    mcmc = -flags->NBURN;
                    
                    rebuild_fstatistic_proposal(orbit, data, model[chain->index[0]], flags, proposal[1]);
                }
                
                mcmc++;
            }
            //Can't continue MCMC until single thread is finished
            #pragma omp barrier
            
        }// end MCMC loop
    }// End of parallelization
    
    //store final state of sampler
    save_chain_state(data, model, chain, flags, mcmc);

    //print aggregate run files/results
    print_waveforms_reconstruction(data,flags);
    print_evidence(chain,flags);
    
    sprintf(filename,"%s/data/waveform_strain.dat",flags->runDir);
    FILE *waveFile = fopen(filename,"w");
    print_waveform_strain(data,model[chain->index[0]],waveFile);
    fclose(waveFile);


    sprintf(filename,"%s/avg_log_likelihood.dat",flags->runDir);
    FILE *chainFile = fopen(filename,"w");
    for(int ic=0; ic<NC; ic++) fprintf(chainFile,"%lg %lg\n",1./chain->temperature[ic],chain->avgLogL[ic]/(double)(flags->NMCMC/data->downsample));
    fclose(chainFile);
    
    //print total run time
    stop = time(NULL);
    
    printf(" ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    sprintf(filename,"%s/ucb_mcmc.log",flags->runDir);
    FILE *runlog = fopen(filename,"a");
    fprintf(runlog," ELAPSED TIME = %g seconds on %i thread(s)\n",(double)(stop-start),numThreads);
    fclose(runlog);
        
    return 0;
}

