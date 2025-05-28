/*
 *  Copyright (C) 2023 Tyson B. Littenberg (MSFC-S.16)
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
 *  MA  02111-1307  USA
 */

//gcc ucb_reconstruction.c -lm -lgsl -lgslcblas -lhdf5 -lomp -L/Users/tlittenb/ldasoft/master/lib -I/Users/tlittenb/ldasoft/master/include -lglass_utils -lglass_ucb -lglass_noise -o ucb_reconstruction

#include <glass_utils.h>
#include <glass_ucb.h>

/*
 * Dependencies:
 * - GSL
 * - GSL CBLAS
 * - OpenMP
 * - HDF5
 * - ldasoft
 */

int main(int argc, char *argv[])
{
    /* usage */
    if(argc!=3)
    {
        fprintf(stdout,"Usage: gb_sim /path/to/injection_file Tobs[yr]\n");
        return 0;
    }
    
    /* parse command line */
    char filename[128];
    sprintf(filename,"%s",argv[1]);
    double T = atof(argv[2])*31457280.0;
    
    
    /* injection file */
    FILE *injectionFile = fopen(argv[1],"r");

    /* injection parameters */
    double f0,dfdt,theta,phi,amp,iota,phi0,psi;

    /* count sources in file */
    int N=0;
    while(!feof(injectionFile))
    {
        int check = fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&theta,&phi,&amp,&iota,&psi,&phi0);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",argv[1]);
            exit(1);
        }
        N++;
    }
    rewind(injectionFile);
    N--;
    fprintf(stdout,"%i injections in input file %s\n",N,argv[1]);
    
    /* set up data array */
    double fmin = 0;    //[Hz]
    double fmax = 0.1; //[Hz]
    int NFREQ = (int)((fmax - fmin)*T);
    
    /* TDI arrays */
    double *A = calloc(2*NFREQ,sizeof(double));
    double *E = calloc(2*NFREQ,sizeof(double));
    double *X = calloc(2*NFREQ,sizeof(double));
    double *Y = calloc(2*NFREQ,sizeof(double));
    double *Z = calloc(2*NFREQ,sizeof(double));
    
    /* parameter vector */
    double *params = malloc(sizeof(double)*8);
    
    /* set up LISA orbit structure */
    struct Orbit *orbit = malloc(sizeof(struct Orbit));
    initialize_analytic_orbit(orbit);
    
    /* now loop over injections and generate waveforms */
    double fstart = 1, fstop = 0; //for keeping track of total bandwidth of population from injections
    for(int n=0; n<N; n++)
    {
        fscanf(injectionFile,"%lg %lg %lg %lg %lg %lg %lg %lg",&f0,&dfdt,&amp,&phi,&theta,&iota,&psi,&phi0);
        
        /* reparameterize */
        params[0] = f0*T;
        params[1] = theta;
        params[2] = phi;
        params[3] = log(amp);
        params[4] = iota;
        params[5] = psi;
        params[6] = phi0;
        params[7] = dfdt*T*T;
        
        /* generate waveform */
        int BW = 2*galactic_binary_bandwidth(orbit->L, orbit->fstar, f0, dfdt, theta, amp, T, NFREQ); //bandwidth of waveform
        double *Atemp = malloc(sizeof(double)*2*BW);
        double *Etemp = malloc(sizeof(double)*2*BW);
        double *Xtemp = malloc(sizeof(double)*2*BW);
        double *Ytemp = malloc(sizeof(double)*2*BW);
        double *Ztemp = malloc(sizeof(double)*2*BW);
        galactic_binary(orbit, "sangria", T, 0.0, params, 8, Xtemp,Ytemp,Ztemp, Atemp, Etemp, BW, 2);
        
        
        /* insert waveform into full array */
        int qmin = (int)(f0*T) - BW/2;
        for(int i=0; i<BW; i++)
        {
            int j = i + qmin;
            
            A[2*j] += Atemp[2*i];
            A[2*j+1] += Atemp[2*i+1];

            E[2*j] += Etemp[2*i];
            E[2*j+1] += Etemp[2*i+1];

            X[2*j] += Xtemp[2*i];
            X[2*j+1] += Xtemp[2*i+1];

 	    Y[2*j] += Ytemp[2*i];
            Y[2*j+1] += Ytemp[2*i+1];

 	    Z[2*j] += Ztemp[2*i];
            Z[2*j+1] += Ztemp[2*i+1];
        }
        
        //keep track of minimum and maximum frequency of joint waveform
        double df = (double)BW/2./T;
        if(f0 - df < fstart) fstart = f0 - df;
        if(f0 + df > fstop)   fstop = f0 + df;
    
        
        //clean up memory for next round
        free(Atemp);
        free(Etemp);
        free(Xtemp);
        
        printProgress((double)(n+1)/(double)N);
        
    }
    fprintf(stdout,"\n");
    
    fprintf(stdout,"fmin = %lg, fmax = %lg\n",fstart,fstop);
    
    /* print injections to file */
    FILE *fptr = fopen("reconstruction_file.dat","w");
    for(int i = 0; i< (int)(fstart*T); i++) fprintf(fptr,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", (double)i/T, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
    for(int i = (int)(fstart*T); i <= (int)(fstop*T); i++)
    {
        fprintf(fptr,"%.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g\n", (double)i/T, X[2*i], X[2*i+1], Y[2*i], Y[2*i+1], Z[2*i], Z[2*i+1], A[2*i], A[2*i+1], E[2*i], E[2*i+1]);
    }
    
    fclose(fptr);

    return 0;
}
