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

#include "glass_utils.h"

double rand_r_U_0_1(unsigned int *seed)
{
   return (double)rand_r(seed) / (double)RAND_MAX;
}

double rand_r_N_0_1(unsigned int *seed)
{
    double u,v,s;
    do {
        u = 2*rand_r_U_0_1(seed)-1;
        v = 2*rand_r_U_0_1(seed)-1;
        s = u * u + v * v;
    } while (s>=1 || s==0);
    double z0 = u * sqrt(-2*log(s)/s);
    
    return z0;
}

int *int_vector(int N)
{
    int *v = malloc( N * sizeof(int) );
    return v;
}

void free_int_vector(int *v)
{
    free(v);
}

int **int_matrix(int N, int M)
{
    int **m = malloc( N * sizeof(int *));
    for(int i=0; i<N; i++) m[i] = int_vector(M);
    return m;
}

void free_int_matrix(int **m, int N)
{
    for(int i=0; i<N; i++) free_int_vector(m[i]);
    free(m);
}

double *double_vector(int N)
{
    double *v = calloc( N , sizeof(double) );
    return v;
}

void free_double_vector(double *v)
{
    free(v);
}

double **double_matrix(int N, int M)
{
    double **m = malloc( N * sizeof(double *));
    for(int i=0; i<N; i++) m[i] = double_vector(M);
    return m;
}

void free_double_matrix(double **m, int N)
{
    for(int i=0; i<N; i++) free_double_vector(m[i]);
    free(m);
}

double ***double_tensor(int N, int M, int L)
{
    double ***t = malloc( N * sizeof(double **));
    for(int i=0; i<N; i++) t[i] = double_matrix(M,L);
    return t;
}

void free_double_tensor(double ***t, int N, int M)
{
    for(int i=0; i<N; i++) free_double_matrix(t[i],M);
    free(t);
}

void astropy_pix2ang_ring(int nside, long ipix, double *theta, double *phi)
{
    int64_t xy = healpixl_ring_to_xy(ipix,nside);
    healpixl_to_radec(xy,nside,0.5,0.5,phi,theta);
    *theta = M_PI/2.0 - *theta;
}
