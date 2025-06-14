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

#include <stdio.h>
#include <stdlib.h>

struct IMRPhenomT{};

struct IMRPhenomT * setup_IMRPT(double Mchirp, double Mtotal, double chi1, double chi2)
{
    printf("ERROR: setup_IMRPT() undefined\n");
    abort();
    return NULL;
}

void free_IMRPhenomT(struct IMRPhenomT *IMRPT)
{
    printf("ERROR: free_IMRPhenomT() undefined\n");
    abort();
}


double mbh_IMRPhenomD_ringdown_frequency_wrapper(double Mchirp, double Mtotal, double chi1, double chi2)
{
    printf("ERROR: mbh_IMRPhenomD_ringdown_frequency_wrapper() undefined\n");
    abort();
    return 0;
}

double mbh_IMRPhenomD_final_spin_wrapper(double eta, double chi1, double chi2)
{
    printf("ERROR: mbh_IMRPhenomD_final_spin_wrapper() undefined\n");
    abort();
    return 0;
}

double mbh_IMRPhenomT_angular_frequency_wrapper(double t, double eta, double Mtotal, struct IMRPhenomT *IMRPT)
{
    printf("ERROR: mbh_IMRPhenomT_angular_frequency_wrapper() undefined\n");
    abort();
    return 0;
}

void mbh_IMRPhenomD_wrapper(double Mchirp, double Mtotal, double chi1, double chi2, double dL, double tc, double phic, double *f, double *t, double *Phase, double *Amp, int N)
{
    printf("ERROR: mbh_IMRPhenomD_wrapper() undefined\n");
    abort();
}

void mbh_IMRPhenomT_wrapper(double Mchirp, double Mtotal, double chi1, double chi2, double dL, double tc, double phic, double *omega, double *t, double *Phase, double *Amp, int N)
{
    printf("ERROR: mbh_IMRPhenomT_wrapper() undefined\n");
    abort();
}
