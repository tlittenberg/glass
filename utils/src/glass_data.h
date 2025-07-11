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

/**
@file glass_data.h
\brief Definition of GLASS data structures.
*/

#ifndef glass_data_h
#define glass_data_h

#define N_TDI_CHANNELS 3
#define FILTER_LENGTH 5e3 //seconds
#define MAXSTRINGSIZE 1024 //!<maximum number of characters for `path+filename` strings


#define WAVELET_DURATION 7680.0
#define WAVELET_BANDWIDTH 6.51041666666667e-5
//#define WAVELET_DURATION 20480.0 //!<duration of wavelet pixels [s]
//#define WAVELET_BANDWIDTH 2.44140625e-05 //!<bandwidth of wavelet pixels [Hz]
//#define WAVELET_DURATION 40960.0 //!<duration of wavelet pixels [s]
//#define WAVELET_BANDWIDTH 1.220703125e-05 //!<bandwidth of wavelet pixels [Hz]
//#define WAVELET_DURATION 81920.0 //!<duration of wavelet pixels [s]
//#define WAVELET_BANDWIDTH 6.103515625e-06 //!<bandwidth of wavelet pixels [Hz]

/*!
 * \brief Analaysis segment and meta data about size of segment, location in full data stream, and LISA observation parameters.
 *
 * The Data structure stores information about the data being analyzed,
 * including size of data, where it is located in the full LISA band,
 * information about gaps, etc.
 *
 * The structure stores the Fourier-domain TDI data channels, and a
 * fiducial model for the instrument noise. It also has memory allocated
 * for storing the reconstructed waveforms, residuals, and noise models.
 */
struct Data
{
    /** @name Size of Data and Model */
    ///@{
    int N;        //!<number of data points bins
    int NFFT;     //!<number of frequency bins
    int Nchannel; //!<number of data channels
    int Nlayer;   //!<number of frequency layers
    double logN;  //!<log total number of data points \f$ \log ( 2 \times N \times N_{\rm channel} \times N_{\rm T} )\f$
    ///@}

    /** @name Random Number Generator Seeds */
     ///@{
    unsigned int cseed; //!<seed for MCMC set by `[--chainseed=INT; default=150914]`
    unsigned int nseed; //!<seed for noise realization set by `[--noiseseed=INT; default=151226]`
    unsigned int iseed; //!<seed for injection parameters set by `[--injseed=INT; default=151012]`
    ///@}

    /** @name Time of Observations */
     ///@{
    double T; //!<observation time
    double sqT; //!<\f$\sqrt{T}\f$
    double t0;   //!<start times of segments
    ///@}

    /** @name Wavelet Basis */
    struct Wavelets *wdm;

    /** @name Analysis Frequency Segment */
     ///@{
    int qmin; //!<minimum DFT frequency bin of segment
    int qmax; //!<maximum DFT frequency bin of segment
    int qpad; //!<number of DFT frequency bins padding ends of segment
    
    int lmin; //!<minimum DWT frequency layer
    int lmax; //!<maximum DWT frequency layer

    double fmin; //!<minimum frequency of segment
    double fmax; //!<maximum frequency of segment
    double sine_f_on_fstar; //!<\f$sin(f * 2\pi L/c)\f$

    ///@}

    /** @name TDI Data and Noise */
     ///@{
    struct TDI *tdi; //!<TDI data channels as seen by sampler
    struct TDI *raw; //!<TDI data channels unaltered from input
    struct TDI *dft; //!<Fourier transform of TDI data channels
    struct TDI *dwt; //!<Wavelet transform of TDI data channels
    struct Noise *noise; //!<Reference noise model

    /**
     \brief Convention for data format
     
     format = "phase" for phase difference (distance)
     format = "frequency" for fractional frequency (velocity) **Use for matching LDC Radler**
     format = "sangria" for fractional frequency w/ LDC Sangria-era TDI & phase conventions
     */
    char format[16];
    char dataDir[MAXSTRINGSIZE]; //!<Directory for storing data files

    /** 
     \brief Basis for analysis

     basis = "fourier" for frequency-domain analysis
     basis = "wavelet" for WDM wavelet-domain analysis
    */
    char basis[16];

    //Spectrum proposal
    double *p; //!<power spectral density of data
    double pmax; //!<maximum power spectrial density
    double SNR2; //!<estimated \f${\rm SNR}^2\f$ of data
    ///@}

    /** @name Model Reconstructions */
     ///@{
    int Nwave; //!<Number of samples for computing posterior reconstructions
    int downsample; //!<Downsample factor for getting the desired number of samples

    double ***h_rec; //!<Store waveform reconstruction samples \f$ 2N \times N_\rm{channel} \times NMCMC \f$
    double ***h_res; //!<Store data residual samples \f$ 2N \times N_\rm{channel} \times NMCMC \f$
    double ***r_pow; //!<Store residual power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    double ***h_pow; //!<Store waveform power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    double ***S_pow; //!<Store noise power samples \f$ N \times N_\rm{channel} \times NMCMC \f$
    char fileName[MAXSTRINGSIZE]; //!<place holder for filnames
    ///@}

    /** @name Signal Injections */
     ///@{
    //int NP; //!<number of parameters of injection
    //struct Source *inj; //!<injected source structure
    ///@}

    
    /** @name Already known sources */
    ///@{
    //int Ncache; //!<number of sources in the cache file
    //char **cache; //!<contents of cache file
    //struct Catalog *catalog; //!< data and metadata for known sources
    ///@}
};

/*!
 *\brief Run flags set or changed from default values when parsing command line for `gb_mcmc` **and** `gb_catalog`.
 *
 *Descriptions of each flag includes default settings and command line arguments to change/set flags.
 *Boolean flags are defined as integers with `0==FALSE` and `1==TRUE`.
 *
 *To see how the defaults are set and then adjusted according to the command line, see parse().
 */
struct Flags
{
    /** @name Run Flags  */
     ///@{
    int verbose;    //!<`[--verbose; default=FALSE]`: increases file output (all chains, less downsampling, etc.)
    int quiet;      //!<`[--quiet; default=FALSE]`: decreases file output (less diagnostic data, only what is needed for post processing/recovery)
    int NMCMC;      //!<`[--steps=INT; default=100000]`: number of MCMC samples
    int NBURN;      //!<number of burn in samples, equals Flags::NMCMC.
    int NINJ;       //!<`[--inj=FILENAME]`: number of injections = number of `--inj` instances in command line
    int NVB;        //!<number of known binaries for `vb_mcmc`
    int DMAX;       //!<`[--sources=INT; default=10]`: max number of sources
    int simNoise;   //!<`[--sim-noise; default=FALSE]`: simulate random noise realization and add to data
    int stationary; //!<`[--stationary; default=FALSE]`: use stationary noise model in wavelet domain
    int fixSky;     //!<`[--fix-sky; default=FALSE]`: hold sky location fixed to injection parameters.  Set to `TRUE` if Flags::knownSource=`TRUE`.
    int fixFdot;
    int fixFreq;    //!<`[--fix-freq; default=FALSE]`: hold GW frequency fixed to injection parameters
    int galaxyPrior;//!<`[--galaxy-prior; default=FALSE]`: use model of galaxy for sky location prior
    int snrPrior;   //!<`[--snr-prior; default=FALSE]`: use SNR prior for amplitude
    int emPrior;    //!<`[--em-prior=FILENAME]`: use input data file with EM-derived parameters for priors.
    int knownSource;//!<`[--known-source; default=FALSE]`: injection is known binary, will need polarization and phase to be internally generated. Sets Flags::fixSky = `TRUE`.
    int detached;   //!<`[--detached; default=FALSE]`: assume binary is detached, fdot prior becomes \f$U[\dot{f}(\mathcal{M}_c=0.15),\dot{f}(\mathcal{M}_c=1.00)]\f$
    int strainData; //!<`[--data=FILENAME; default=FALSE]`: read data from ASCII file instead of simulate internally.
    int hdf5Data;   //!<'[--hdf5Data=FILENAME; default=FALSE]`: read data from LDC HDF5 file (compatible w/ Sangria dataset).
    int orbit;      //!<`[--orbit=FILENAME; default=FALSE]`: use numerical spacecraft ephemerides supplied in `FILENAME`. `--orbit` argument sets flag to `TRUE`.
    int prior;      //!<`[--prior; default=FALSE]`: set log-likelihood to constant for testing detailed balance.
    int debug;      //!<`[--debug; default=FALSE]`: coarser settings for proposals and verbose output for debugging
    int cheat;      //!<start sampler at injection values
    int burnin;     //!<`[--no-burnin; default=TRUE]`: chain is in the burn in phase
    int maximize;   //!<maximize over extrinsic parameter during burn in phase.
    int update;     //!<`[--update=FILENAME; default=FALSE]`: use Gaussian Mixture Model approximation to previous posterior as current prior.
    int updateCov;  //!<`[--update-cov=FILENAME; default=FALSE]`: updating fit from covariance matrix files built from chain samples, passed as `FILENAME`, used in draw_from_cov().
    int match;      //!<[--match=FLOAT; default=0.8]`: match threshold for chain sample clustering in post processing.
    int rj;         //!<--no-rj; default=TRUE]`: flag for determining if trans dimensional MCMC moves (RJMCMC) are enabled.
    int calibration;//!<`[--calibration; default=FALSE]`: flag for determining if model is marginalizing over calibration  uncertainty.
    int confNoise;  //!<`[--conf-noise; default=FALSE]`: include model of confusion noise in \f$S_n(f)\f$, either for simulating noise or as starting value for parameterized noise model.
    int resume;     //!<`[--resume; default=FALSE]`: restart sampler from run state saved during checkpointing. Starts from scratch if no checkpointing files are found.
    int catalog;    //!<`[--catalog=FILENAME; default=FALSE]`: use list of previously detected sources supplied in `FILENAME` to clean bandwidth padding (`gb_mcmc`) or for building family tree (`gb_catalog`).
    int grid;       //!<`[--ucb-grid=FILENAME; default=FALSE]`: flag indicating if a gridfile was supplied
    int threads;    //!<number of openMP threads for parallel tempering
    int psd;        //!<`[--psd=FILENAME; default=FALSE]`: use PSD input as ASCII file from command line
    int help;       //!<`[--help]`: print command line usage and exit
    ///@}

    
    /** @name Input File Names
     */
     ///@{
    char runDir[MAXSTRINGSIZE];       //!<store `DIRECTORY` to serve as top level directory for output files.
    char vbFile[MAXSTRINGSIZE];       //!<store `FILENAME` of list of known binaries `vb_mcmc`
    char ucbGridFile[MAXSTRINGSIZE];  //!<`[--ucb-grid=FILENAME]` frequency grid for multiband UCB analysis
    char **injFile;                   //!<`[--inj=FILENAME]`: list of injection files. Can support up to `NINJ=10` separate injections.
    char noiseFile[MAXSTRINGSIZE];    //!<file containing reconstructed noise model for `gb_catalog` to compute SNRs against.
    char cdfFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input chain file from Flags::update.
    char gmmFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input gmm file from Flags::update.
    char covFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input covariance matrix file from Flags::updateCov.
    char matchInfile1[MAXSTRINGSIZE]; //!<input waveform \f$A\f$ for computing match \f$(h_A|h_B)\f$
    char matchInfile2[MAXSTRINGSIZE]; //!<input waveform \f$B\f$ for computing match \f$(h_A|h_B)\f$
    char pdfFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input priors for Flags:knownSource.
    char psdFile[MAXSTRINGSIZE];      //!<store `FILENAME` of input psd file from Flags::psd.
    char catalogFile[MAXSTRINGSIZE];  //!<store `FILENAME` containing previously identified detections from Flags::catalog for cleaning padding regions
     ///@}
     
    /** @name Remove components of simulated data  */
     ///@{
    int no_mbh;
    int no_ucb;
    int no_ucb_hi;
    int no_vgb;
    int no_noise;
     ///@}
};

/**
 \brief Structure containing settings and housekeeping data for each of the parallel chains.
 */
struct Chain
{
    /// Number of chains
    int NC;
    
    /// Array containing current order of chains in temperature ladder
    int *index;
    
    /// Size of model being sampled by chain.  Depricated?
    int **dimension;
    
    /// Array tracking acceptance rate of parallel chain exchanges.
    double *acceptance;
    
    /// Array of chain temperatures
    double *temperature;
    
    /// Array storing \f$\langle \log p(d|\vec\theta)^{1/T} \rangle\f$ for thermodynamic integration.
    double *avgLogL;
    
    /// Annealing temperature for all chains used during burnin **DEPRICATED**
    double annealing;
    
    /// Store the maximum value of \f$\log p(d|\vec\theta)\f$ encountered by the chain.  Used for determining when burn-in is complete.
    double logLmax;
    
    /// Random Number Generator (RNG) seed
    unsigned int *r;
    
    /** @name Chain File Pointers
     By default only the cold chain `M=0` is saved.  When Flags::verbose = `TRUE` files for each of the parallel chain are written.
     */
    ///@{
    
    /**
     \brief Noise parameter chain file: `chains/noise_chain.dat.M`
     
     Columns:  `step | logL | logL_norm | etaA | eta E`
     */
    FILE **noiseFile;
    FILE **foregroundFile;
    
    /**
     \brief Markov chain state (iterations, likelihoods, etc.): `chains/model_chain.dat.M`
     
     Columns: `step | N_live | logL | logL_norm | t_0`
     */
    FILE **chainFile;
    
    /**
     \brief Calibration parameter chain file: `chains/calibration_chain.dat.M`
     
     Columns:
     */
    FILE **calibrationFile;
    
    /**
     
     \brief Signal parameter chain files for discrete models, only for cold chains: `chains/dimension_chain.dat.D`
     
     Columns: `f | fdot | A | cos_colat | long | cos_inc | psi | phi`
     */
    FILE **dimensionFile;
    
    /**
     \brief Full signal parameter chain files: `chains/dimension_chain.dat.M`

     Columns: `f | fdot | A | cos_colat | long | cos_inc | psi | phi`
     */
    FILE **parameterFile;
    
    /**
     \brief Log-likelhood values for each parallel chain.
     
     Columns: `step | logL_0 | logL_1 | ... | logL_NC-1`
     */
    FILE *likelihoodFile;
    
    /**
     \brief Temperature of each parallel chain to monitor adaptive temperature spacing.
     
     Columns: `step | 1/T_0 | 1/T_1 | ... | 1/T_NC-1`
     */
    FILE *temperatureFile;
    ///@}
    
    char chainDir[MAXSTRINGSIZE]; //!<store chain directory.
    char chkptDir[MAXSTRINGSIZE]; //!<store checkpoint directory.

};



/**
\brief Structure containing parameters and meta data for noise model in narrow-band analysis segment.
*/
struct Noise
{
    /// Number of data samples fit by noise model
    int N;

    /// Minimum index of wavelet pixels  
    int kmin; 

    /// Number of TDI channels
    int Nchannel;
    
    /// Number of frequency layers (for DWT data)
    int Nlayer;
    
    ///@name Constant Noise Parameters
    ///Each \f$\eta\f$ is a multiplier to the assumed noise level \f$S_n\f$ stored in Data structure. One per channel (`X` for 4-link, `A`,`E` or `X`,`Y`,`Z` for 6-link)
    ///@{
    double *eta;
    ///@}
    
    ///@name Noise Model
    ///Composite noise model to use over the analysis window \f$\eta_I \times Sn_I\f$
    ///@{
    double *f;
    
    double ***C;    //!<Covariance matrix
    double ***invC; //!<Inverse covariance matrix>
    double *detC;   //!<Determinent of covariance matrix
    
    double *transfer;
    ///@}
};

/**
\brief Structure containing calibration parameters
 */
struct Calibration
{
    ///@name Amplitude parameters for each TDI channel
    ///@{
    double dampA;
    double dampE;
    double dampX;
    double dampY;
    double dampZ;
    ///@}
    
    ///@name Overall phase parameters for each TDI channel
    ///@{
    double dphiA;
    double dphiE;
    double dphiX;
    double dphiY;
    double dphiZ;
    ///@}

    ///@name Phase correction to Re and Im part of TDI channels
    ///@{
    double real_dphiA;
    double real_dphiE;
    double real_dphiX;
    double real_dphiY;
    double real_dphiZ;
    double imag_dphiA;
    double imag_dphiE;
    double imag_dphiX;
    double imag_dphiY;
    double imag_dphiZ;
    ///@}

};

/**
 \brief Print command line arguments for GLASS apps
 */
void print_glass_usage();

/**
 \brief Show progress bar
 */
void printProgress (double percentage);

/**
 \brief Print git hash of running version
 */
void print_version(FILE *fptr);

/**
  \brief Name and create output directories
 */
void setup_run_directories(struct Flags *flags, struct Data *data, struct Chain *chain);

/**
 \brief Loads spacecraft orbits, either analytically or from tabulated ephemerides.
 */
void initialize_orbit(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Allocates and initializes Chain structure and prepares output files.
 */
void initialize_chain(struct Chain *chain, struct Flags *flags, unsigned int *seed, const char *mode);

/** @name Allocate memory for structures */
///@{
void alloc_data(struct Data *data, struct Flags *flags);
void alloc_noise(struct Noise *noise, int N, int Nlayer, int Nchannel);
void alloc_calibration(struct Calibration *calibration);
///@}

/**
 \brief Shallow copy of Data structure
 */
void copy_data(struct Data *origin, struct Data *copy);

/** @name Deep copy structure contents */
///@{
void copy_noise(struct Noise *origin, struct Noise *copy);
void copy_Cij(double ***origin, double ***copy, int M, int N);
void copy_calibration(struct Calibration *origin, struct Calibration *copy);
///@}

/** @name Free memory for structures */
///@{
void free_noise(struct Noise *noise);
void free_chain(struct Chain *chain, struct Flags *flags);
void free_calibration(struct Calibration *calibration);
///@}

/**
 \brief Reads data from external source input using `--data` flag
 */
void ReadData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/**
 \brief Reads LDC-formatted HDF5 data using `--h5-data` flag
 */
void ReadHDF5(struct Data *data, struct TDI *tdi, struct TDI *tdi_dwt, struct Flags *flags);

/**
 \brief Reads ASCII data using `--data` flag
 */
void ReadASCII(struct Data *data, struct TDI *tdi);

/**
 \brief Get theoretical noise PSDs for TDI channels
 */
void GetNoiseModel(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/** @name Add simulated Gaussian noise realization to data */
///@{
void AddNoise(struct Data *data, struct TDI *tdi);
void AddNoiseWavelet(struct Data *data, struct TDI *tdi);
///@}

/**
 \brief Generate noise-only simulated data
 */
void SimulateData(struct Data *data, struct Orbit *orbit, struct Flags *flags);

/** @name Wrapper functions that call data print functions */
///@{
void print_data(struct Data *data, struct Flags *flags);
void print_wavelet_fourier_spectra(struct Data *data, struct TDI *tdi, char filename[]);
///@}

/**
 \brief convert DWT of TDI data to DFT
 */
void wavelet_layer_to_fourier_transform(struct Data *data);

/**
 \brief Parse command line
 */
void parse_data_args(int argc, char **argv, struct Data *data, struct Orbit *orbit, struct Flags *flags, struct Chain *chain, char basis[]);

/**
 \brief copy `argv` string because `parse()` does not preserve order
 */
void copy_argv(int argc, char **argv, char **new_argv);

/**
 \brief Check if file exists and handle missing files cleanly
 */
int checkfile(char filename[]);



#endif /* glass_data_h */
