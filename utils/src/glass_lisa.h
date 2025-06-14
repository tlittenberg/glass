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

/**
@file glass_lisa.h
\brief Codes defining LISA instrument model
 
 Here are functions common to all parts of the LISA analysis including
  - Constellation configuration and orbit model
  - Instrument noise models
  - Methods for constructing TDI response to incident GWs.
*/

#ifndef lisa_h
#define lisa_h

/// Mean arm length of constellation (m) for baseline LISA configuration
#define LARM 2.5e9

/// Sample cadence for LISA (s)
#define LISA_CADENCE 5

/// Number of orbit ephemerides samples per year
#define LISA_ORBIT_SAMPLES_PER_YEAR 200

/** @name Component Noise Levels For Phase Data */
///@{

/// Photon shot noise power
#define SPS 8.321000e-23

/// Acceleration noise power
#define SACC 9.000000e-30

/// Position noise power when using phase data
#define SLOC 2.89e-24
///@}


/**
 * \brief Ephemerides of individual spacecraft and metadata for using orbits in waveform modeling.
 *
 * If numerical orbit files are provided, they are interpolated to the sample rate of the data using `GSL` cubic spline functions.
 *
 * If not, the eccentric inclined analytic model is computed once at the data sampling rate and stored.
 */
struct Orbit
{
    /// Filename input from `--orbit` command line argument when using numerical orbits
    char OrbitFileName[1024];
    
    /// Size of orbit arrays
    int Norb;
    
    /** @name Constellation Parameters */
    ///@{
    double L;        //!< Average armlength of constellation
    double fstar;    //!< Transfer frequency \f$f_* = 1/(L/c)\f$.
    double ecc;      //!< Eccentricity of spacecraft orbits
    double R;        //!< Distance to constellation guiding center from Sun (1 AU)
    double lambda_0; //!< Initial phase of constellation w.r.t. ecliptic
    double kappa_0;  //!< Initial phae of constellation guiding center   
    ///@}
    
    double *t;  //!<time step relative to start of mission (seconds)

    /** @name Spacecraft Ephemerides
     Cartesian ecliptic coordinates, in meters, where the x-y plane is the ecliptic.
     */
    ///@{
    double **x; //!<x-coordinate at each time step
    double **y; //!<y-coordinate at each time step
    double **z; //!<z-coordinate at each time step
    ///@}
    
    /** @name Derivatives of Orbits for Cubic Spline Interpolation
     Stores the derivatives for the cubic spline for the location of each spacecraft (\f$dx,dy,dz\f$) and some internal `GSL` workspace `acc`.
     */
    ///@{
    struct CubicSpline **dx; //!<spline derivatives in x-coordinate
    struct CubicSpline **dy; //!<spline derivatives in y-coordinate
    struct CubicSpline **dz; //!<spline derivatives in z-coordinate
    ///@}
    
    /**
     \brief Function prototyp for retreiving spacecraft locations.
     \param[in] time (double)
     \param[in] orbit data (Orbit*)
     \param[out] eccliptic cartesian location of each spacecraft
     
     If an orbit file is input this points to interpolate_orbits() which uses `GSL` cubic splines to interpolate the ephemerides at the needed time steps
     
     Otherwise this points to analytic_orbits() which is passed an arbitrary time \f$t\f$ and returns the spacecraft location.
     */
    void (*orbit_function)(struct Orbit*,double,double*,double*,double*);
};

/**
 \brief Structure for Time Delay Interferometry data and metadata
 
 Contains time or frequency series of TDI data channels (Michelson-like or orthogonal), and metadata about the number of channels in use, the sampling rate, and the size of the datastream.
 */
struct TDI
{
    /**
     @name Michelson TDI Channels
     */
    ///@{
    double *X; //!<X channel
    double *Y; //!<Y channel
    double *Z; //!<Z channel
    ///@}
    
    /**
     @name Noise-orthogonal TDI Channels
     */
    ///@{
    double *A; //!<A channel
    double *E; //!<E channel
    double *T; //!<T channel
    ///@}
    
    /// Number of data channels in use. 1 for 4-link, 2 for 6-link.
    int Nchannel;
    
    /// Size of data. Time samples or frequency bins.
    int N;
    
    /// Data cadence. \f$\Delta t\f$ for time-domain, \f$  \Delta f\f$ for frequency-domain.
    double delta;
};

/**
\brief Prints our beautiful LISA ASCII logo
 
\verbatim
                               OOOOO
                              OOOOOOO
                            11111OOOOO
 OOOOO            11111111    O1OOOOO
OOOOOOO  1111111             11OOOO
OOOOOOOO                    11
OOOOO1111                 111
  OOOO 1111             111
          1111    OOOOOO11
             111OOOOOOOOOO
               OOOOOOOOOOOO
               OOOOOOOOOOOO
               OOOOOOOOOOOO
                OOOOOOOOOO
                  OOOOOO
\endverbatim

 */
void print_LISA_ASCII_art(FILE *fptr);

/**
 \brief Numerical interpolation of spacecraft ephemerides using cubic spline
 */
void interpolate_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);

/**
 \brief Analytic function for spacecraft ephemerides
 
 Currently implemented as first order Keplerian elliptic inclined orbits.
 */
void analytic_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z);

/**
 \brief store metadata and assign orbit function to Orbit
 */
void initialize_analytic_orbit(struct Orbit *orbit);

/**
 \brief store meta data and prepare ephemerides interpolation
 
 - parse spacecraft ephemeris file in Orbit::OrbitFileName
 - compute cubic spline derivatives at data points
 - estimate average armlengths
 - store metadata
 */
void initialize_numeric_orbit(struct Orbit *orbit);

/**
 \brief store meta data and interpolant of analytic orbit ephemerides
 
 - pre-compute spacecraft ephemerides on time grid from analytic orbit model
 - set up interpolants for spacecraft ephemerides

 @param[in] Orbit structure
 @param[in] Tobs total observation time \f$[{\rm s}]\f$
 @param[in] t0 initial time of orbit grid \f$[{\rm s}]\f$
 */
void initialize_interpolated_analytic_orbits(struct Orbit *orbit, double Tobs, double t0);

/**
 \brief allocate memory for Orbit
 */
void alloc_orbit(struct Orbit *orbit, int Norb);

/**
 \brief free memory allocated for Orbit
 */
void free_orbit(struct Orbit *orbit);


/**
 @name  LISA Time Delay Interferometer functions
 
 TDI 1.5 assuming equal arm detector
 
 @param[in] L average LISA armlength \f$ L\ [{\rm s}]\f$
 @param[in] fstar LISA transfer frequency \f$ L/2\pi c\ [{\rm Hz}]\f$
 @param[in] T observation time \f$[{\rm s}]\f$
 @param[in] d gravitational wave phase at each vertex, at each frequency
 @param[in] f0 carrier frequency of signal \f$f_0\ [{\rm Hz}]\f$
 @param[in] q frequency bin with carrier frequency
 @param[in] NI number of interferometer channels (1 for 4-link, 2 for 6-link)
 @param[in] BW bandwidth of signal
 @param[out] M Michelson-like TDI channel (4-link)
 @param[out] A TDI A channel
 @param[out] E TDI E channel
*/
///@{

/// LISA TDI (phase)
void LISA_tdi(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI);
/// LISA TDI (frequency)
void LISA_tdi_FF(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI);
/// LISA TDI (LDC Sangria data)
void LISA_tdi_Sangria(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI);
///@}
void LISA_polarization_tensor_njc(double costh, double phi, double eplus[4][4], double ecross[4][4], double k[4]);


/**
 \brief Generic LISA response using interpolated signal frequency or phase and amplitude
 
 @param[in] Orbit structure
 @param[in] tarray time array for response
 @param[in] N size of tarray
 @param[in] costheta cosine ecliptic co-latitude
 @param[in] psi ecliptic longitude
 @param[in] cosi cosine inclination
 @param[in] psi polarization angle
 @param[in] amp_spline signal amplitude interpolant
 @param[in] freq_spline signal frequency interpolant
 @param[in] phase_spline signal phase interpolant
 @param[in] phase_ref reference phase (i.e. fast part of waveform)
 @param[out] amp TDI structure of amplitude response
 @param[out] phase TDI structure of phase response
 */
void LISA_spline_response(struct Orbit *orbit, double *tarray, int N, double costh, double phi, double cosi, double psi, struct CubicSpline *amp_spline, struct CubicSpline *freq_spline, struct CubicSpline *phase_spline, double *phase_ref, struct TDI *tdi_amp, struct TDI *tdi_phase);


/** @name  LISA Noise Model for equal arm, TDI1.5 configuration */
///@{
/// Pre-configured proof mass and optical path length noises
void get_noise_levels(char model[], double f, double *Spm, double *Sop);
/// A and E noise (phase)
double AEnoise(double L, double fstar, double f);
/// A and E noise (frequency)
double AEnoise_FF(double L, double fstar, double f, double Spm, double Sop);
/// T noise (frequency)
double Tnoise_FF(double L, double fstar, double f, double Spm, double Sop);
/// Michelson-like X,Y,Z channel noise (phase)
double XYZnoise(double L, double fstar, double f);
/// Michelson-like X,Y,Z channel noise (frequency)
double XYZnoise_FF(double L, double fstar, double f, double Spm, double Sop);
/// Michelson-like X,Y,Z channel cross spectra (frequency)
double XYZcross_FF(double L, double fstar, double f, double Spm, double Sop);
/// Confusion noise estimate for A,E channels (phase)
double GBnoise(double T, double f);
/// Confusion noise estimate for A,E channel (frequency)
double GBnoise_FF(double T, double fstar, double f);
/// Noise transfer function \f$ sin^2(f/f_*)\f$
double noise_transfer_function(double x);
///@}

/**
 \brief Print full-spectrum noise model to compare against other models/documents.
 
 @params[in] Orbit containing necessary LISA metadata
 @returns FILE (ascii) "psd.dat" with columns `f SnA SnE SnT`
 */
void test_noise_model(struct Orbit *orbit);

/**@name Memory handling for TDI structure
 */
///@{
/// allocate memory and initializ TDI structure
void alloc_tdi(struct TDI *tdi, int N, int Nchannel);
/// deep copy contents from origin into copy
void copy_tdi(struct TDI *origin, struct TDI *copy);
/// deep copy of segment of size `N` starting at `index`
void copy_tdi_segment(struct TDI *origin, struct TDI *copy, int index, int N);
/// free contents and overall TDI structure
void free_tdi(struct TDI *tdi);
///@}

/**
 \brief HDF5 parser for LISA Data
 
 Hard-coded for <a href="https://lisa-ldc.lal.in2p3.fr/challenge2">LDC Sangria</a> data format.
 Performs a deep copy of TDI time series data in HDF5 to input TDI structure.
 */
void LISA_Read_HDF5_LDC_TDI(struct TDI *tdi, char *fileName, const char *dataName);

/**
 \brief HDF5 parser for LISA Data
 
 Hard-coded for <a href="https://lisa-ldc.lal.in2p3.fr/challenge1">LDC Radler</a> data format.
 Performs a deep copy of TDI time series data in HDF5 to input TDI structure.
 */

void LISA_Read_HDF5_LDC_RADLER_TDI(struct TDI *tdi, char *fileName);

/**
 \brief Compute LISA detector tensor etc. based on spacecraft and source locations
 
 See Eqs 40 and 41 of <a href="https://doi.org/10.1103/PhysRevD.76.083006">Cornish & Littenberg, Phys. Rev. D76, 083006</a>.
 
 @param[in] L arm length
 @param[in] eplus plus polarization tensor
 @param[in] ecross cross polarization tensor
 @param[in] x,y,z cartesian ecliptic coordinates of each spacecraft \f$\vec{r}_i=\vec{x}_i+\vec{y}_i+\vec{z}_i\f$
 @param[in] k unit vector to source sky location to barycenter \f$\hat{k}\f$
 @param[out] dplus plus polarization detector tensor before applying transfer function
 @param[out] dcross cross polarization detector tensor before applying transfer function
 @param[out] kdotr \f$\hat{k} \dot \vec{r}_i\f$
 */
void LISA_detector_tensor(double L, double eplus[4][4], double ecross[4][4], double x[4], double y[4], double z[4], double k[4], double dplus[4][4], double dcross[4][4], double kdotr[4][4]);


/**
 \brief Compute LISA polarization tensor etc. based on source locations
 
 See Eq 30 <a href="https://doi.org/10.1103/PhysRevD.76.083006">Cornish & Littenberg, Phys. Rev. D76, 083006</a>.
 
 @param[in] costheta cosine of ecliptic co-latitude
 @param[in] phi ecliptic longitude
 @param[out] eplus cross polarization tensor
 @param[out] ecross cross polarization tensor
 @param[out] k unit vector to source location \f$\hat{k}\f$
 */
void LISA_polarization_tensor(double costh, double phi, double eplus[4][4], double ecross[4][4], double k[4]);

void LISA_detector_time(struct Orbit *orbit, double costh, double phi, double *time, int N,  double *time_sc);

/**
 \brief Convert noise-orthogonal AET channels from Michelson-like XYZ channels
 
 */
void XYZ2AE(double X, double Y, double Z, double *A, double *E);
void XYZ2AET(double X, double Y, double Z, double *A, double *E, double *T);

#endif /* lisa_h */
