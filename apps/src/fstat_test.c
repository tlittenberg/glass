//
//  fstat_test.c
//  ucb_wavelet_mcmc
//
//  Created by Tyson Littenberg on 10/8/24.
//

#include <glass_utils.h>
#include <glass_noise.h>
#include <glass_ucb.h>


int main(int argc, char *argv[])
{
    /* Allocate data structures */
    struct Flags  *flags_dft = malloc(sizeof(struct Flags));
    struct Flags  *flags_dwt = malloc(sizeof(struct Flags));
    struct Orbit  *orbit_dft = malloc(sizeof(struct Orbit));
    struct Orbit  *orbit_dwt = malloc(sizeof(struct Orbit));
    struct Chain  *chain_dft = malloc(sizeof(struct Chain));
    struct Chain  *chain_dwt = malloc(sizeof(struct Chain));
    struct Data   *data_dft  = malloc(sizeof(struct Data));
    struct Data   *data_dwt  = malloc(sizeof(struct Data));
    
    /* Parse command line and set defaults/flags */
    sprintf(data_dft->basis,"fourier");
    sprintf(data_dwt->basis,"wavelet");

    parse_data_args(argc,argv,data_dft,orbit_dft,flags_dft,chain_dft,"fourier");
    parse_ucb_args(argc,argv,flags_dft);
    parse_data_args(argc,argv,data_dwt,orbit_dwt,flags_dwt,chain_dwt,"wavelet");
    parse_ucb_args(argc,argv,flags_dwt);
    
    int DMAX = flags_dft->DMAX;

    /* Setup output directories for chain and data structures */
    sprintf(data_dft->dataDir,"%s/data_dft",flags_dft->runDir);
    sprintf(data_dwt->dataDir,"%s/data_dwt",flags_dwt->runDir);
    sprintf(chain_dft->chainDir,"%s/chains_dft",flags_dft->runDir);
    sprintf(chain_dft->chkptDir,"%s/checkpoint_dft",flags_dft->runDir);
    sprintf(chain_dwt->chainDir,"%s/chains_dwt",flags_dft->runDir);
    sprintf(chain_dwt->chkptDir,"%s/checkpoint_dwt",flags_dft->runDir);

    mkdir(flags_dft->runDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data_dft->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(data_dwt->dataDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain_dft->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain_dwt->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain_dft->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(chain_dwt->chkptDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    
    /* Initialize data structures */
    alloc_data(data_dft, flags_dft);
    alloc_data(data_dwt, flags_dwt);
    
    mkdir(chain_dft->chainDir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    /* Initialize LISA orbit model */
    initialize_orbit(data_dft, orbit_dft, flags_dft);
    initialize_orbit(data_dwt, orbit_dwt, flags_dwt);
    initialize_interpolated_analytic_orbits(orbit_dwt, data_dwt->T, data_dwt->t0);
    
    struct Source **inj_dft=malloc(DMAX*sizeof(struct Source*));
    struct Source **inj_dwt=malloc(DMAX*sizeof(struct Source*));
    
    for(int n=0; n<DMAX; n++)
    {
        inj_dft[n] = malloc(sizeof(struct Source));
        inj_dwt[n] = malloc(sizeof(struct Source));
    }
        
    UCBInjectSimulatedSource(data_dft,orbit_dft,flags_dft,inj_dft);
    UCBInjectSimulatedSource(data_dwt,orbit_dwt,flags_dwt,inj_dwt);
    
    //matrix to hold maximized extrinsic parameters
    double *Fparams = calloc(4,sizeof(double));
    
    double logL_X_dft;
    double logL_AE_dwt;
    double logL_AE_dft;
    
    double f = inj_dft[0]->f0;
    double fdot = inj_dft[0]->dfdt;
    double theta = acos(inj_dft[0]->costheta);
    double phi = inj_dft[0]->phi;
    
    get_Fstat_logL(orbit_dft, data_dft, f, fdot, theta, phi, &logL_X_dft, &logL_AE_dft, Fparams);
    logL_AE_dwt = get_Fstat_logL_wavelet(orbit_dwt, data_dwt, f, fdot, theta, phi);
    
    printf("logL_AE_dft = %lg, logL_AE_dwt = %lg\n", logL_AE_dft, logL_AE_dwt);
    printf("proposal_dft = %lg, proposal_dwt = %lg\n", exp(sqrt(2*logL_AE_dft)), exp(sqrt(2*logL_AE_dwt)));
    
    return 1;
}
