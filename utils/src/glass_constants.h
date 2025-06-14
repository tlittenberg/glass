/*
 * Copyright 2019 Tyson B. Littenberg
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
 @file glass_constants.h
 \brief Mathematical and physical constants
 */


#ifndef constants_h
#define constants_h

#define FIXME 0
#define SNRPEAK 10 ///< Peak of SNR prior

/* --------------  MATHEMATICAL CONSTANTS  -------------- */
/* Some square roots */
#define SQ3 1.73205080757    ///< \f$\sqrt{3}\f$
#define SQ8 2.82842712474619 ///< \f$\sqrt{8}\f$

/* Pi's and frinds
 use math.h where possible
 M_PI for PI,
 M_PI_2 for PI/2
 M_PI_4 for PI/4
 */
#define PI2   6.283185307179586      ///< \f$2\pi\f$
#define RTPI  1.77245385090552       ///< \f$\sqrt{\pi}\f$
#define RT2PI 2.5066282746310005024  ///< \f$\sqrt{2\pi}\f$

/* Convert between angles*/
#define RAD2DEG 0.01745329251 ///< 1 deg [rad]
#define DEG2RAD 57.2957795131 ///< 1 rad [deg]

/* ----------------  NATURAL CONSTANTS  ----------------- */

//! Gravitational Constant \f$G\f$ [m\f$^3\f$ kg\f$^{-1}\f$ s\f$^{-2}\f$]
#define GNEWTON 6.67408e-11

//! Speed of light \f$c\f$ [m/s]
#define CLIGHT 299792458.

//! Mass of the Sun \f$M_\odot\f$ [kg]
#define MSUN 1.98848e30

//! Mass of the Sun \f$M_\odot G c^{-3}\f$ [s]
#define TSUN  4.92569043916e-6

//! Radius of the Sun \f$_\odot c^{-1}\f$ [s]
#define RSUN 2.32060541029

//! Parsec [m]
#define PC 3.0856775815e16

//! Year [s]
#define YEAR 3.15581497632e7

//! Day [s]
#define DAY 86400.0

//! Astronomical unit [m]
#define AU 1.49597870700e11

//! Velocity of the Earth around the Sun (relative to CLIGHT)
#define VEARTH 0.00010103671

#endif /* constants_h */
