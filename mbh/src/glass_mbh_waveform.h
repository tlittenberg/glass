/*
 * Copyright 2025 Tyson B. Littenberg & Neil Cornish
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
@file glass_mbh_waveform.h
\brief Massive Black Hole binary waveform generator
 */

#ifndef glass_mbh_waveform_h
#define glass_mbh_waveform_h

#include <stdio.h>


void mbh_barycenter_waveform(double *params, int N, double *freq, double *time, double *phase, double *amp, char model[]);

double mbh_final_spin(double *params);
double mbh_ringdown_frequency(double *params);

double get_PN_time(double f, double Mtotal, double eta, double tc);
double get_PN_frequency(double Mchirp, double tc, double t);

void mbh_fd_waveform(struct Orbit *orbit, double Tobs, double t0, double *params, double *X, double *Y, double *Z);
void mbh_td_waveform(struct Orbit *orbit, double Tobs, double t0, double *params, double *X, double *Y, double *Z);

#endif /* glass_mbh_waveform_h */
