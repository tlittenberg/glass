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


#include "glass_utils.h"

void print_LISA_ASCII_art(FILE *fptr)
{
    fprintf(fptr,"                                          \n");
    fprintf(fptr,"                                OOOOO     \n");
    fprintf(fptr,"                               OOOOOOO    \n");
    fprintf(fptr,"                             11111OOOOO   \n");
    fprintf(fptr,"  OOOOO            11111111    O1OOOOO    \n");
    fprintf(fptr," OOOOOOO  1111111             11OOOO      \n");
    fprintf(fptr," OOOOOOOO                    11           \n");
    fprintf(fptr," OOOOO1111                 111            \n");
    fprintf(fptr,"   OOOO 1111             111              \n");
    fprintf(fptr,"           1111    OOOOOO11               \n");
    fprintf(fptr,"              111OOOOOOOOOO               \n");
    fprintf(fptr,"                OOOOOOOOOOOO              \n");
    fprintf(fptr,"                OOOOOOOOOOOO              \n");
    fprintf(fptr,"                OOOOOOOOOOOO              \n");
    fprintf(fptr,"                 OOOOOOOOOO               \n");
    fprintf(fptr,"                   OOOOOO                 \n");
}

void alloc_orbit(struct Orbit *orbit, int Norb)
{
    //allocate memory for local workspace
    orbit->Norb = Norb;
    
    //allocate memory for orbit structure
    orbit->t  = double_vector(orbit->Norb);
    orbit->x  = double_matrix(3,orbit->Norb);
    orbit->y  = double_matrix(3,orbit->Norb);
    orbit->z  = double_matrix(3,orbit->Norb);
    orbit->dx = malloc(sizeof(struct CubicSpline *)*3);
    orbit->dy = malloc(sizeof(struct CubicSpline *)*3);
    orbit->dz = malloc(sizeof(struct CubicSpline *)*3);
    
    for(int i=0; i<3; i++)
    {
        orbit->dx[i] = alloc_cubic_spline(orbit->Norb);
        orbit->dy[i] = alloc_cubic_spline(orbit->Norb);
        orbit->dz[i] = alloc_cubic_spline(orbit->Norb);
    }
}

void interpolate_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z)
{
    int i;
    
    for(i=0; i<3; i++)
    {
        x[i+1] = spline_interpolation(orbit->dx[i], t);
        y[i+1] = spline_interpolation(orbit->dy[i], t);
        z[i+1] = spline_interpolation(orbit->dz[i], t);
    }
}

/* ********************************************************************* */
/*        Rigid approximation position of each LISA spacecraft           */
/* ********************************************************************* */
void analytic_orbits(struct Orbit *orbit, double t, double *x, double *y, double *z)
{
    
    double alpha = PI2*t/YEAR; // + kappa_0
    
    /*
     double beta1 = 0.;              //lambda_0
     double beta2 = 2.0943951023932; //2.*pi/3. + lambda_0;
     double beta3 = 4.18879020478639;//4.*pi/3. + lambda_0;
     */
    
    double sa = sin(alpha);
    double ca = cos(alpha);
    
    double sa2  = sa*sa;
    double ca2  = ca*ca;
    double saca = sa*ca;
    double AUca = AU*ca;
    double AUsa = AU*sa;
    
    double sb,cb;
    
    sb = 0.0;//sin(beta1);
    cb = 1.0;//cos(beta1);
    x[1] = AUca + orbit->R*(saca*sb - (1. + sa2)*cb);
    y[1] = AUsa + orbit->R*(saca*cb - (1. + ca2)*sb);
    z[1] = -SQ3*orbit->R*(ca*cb + sa*sb);
    
    sb = 0.866025403784439;//sin(beta2);
    cb = -0.5;//cos(beta2);
    x[2] = AUca + orbit->R*(saca*sb - (1. + sa2)*cb);
    y[2] = AUsa + orbit->R*(saca*cb - (1. + ca2)*sb);
    z[2] = -SQ3*orbit->R*(ca*cb + sa*sb);
    
    sb = -0.866025403784438;//sin(beta3);
    cb = -0.5;//cos(beta3);
    x[3] = AUca + orbit->R*(saca*sb - (1. + sa2)*cb);
    y[3] = AUsa + orbit->R*(saca*cb - (1. + ca2)*sb);
    z[3] = -SQ3*orbit->R*(ca*cb + sa*sb);

    // next order corrections (no longer equal arm)
     
    /*
    x[i] += 0.125*ec*ec*AU*(-10.*ca - 5.*ca*cb*cb + 3.*ca*ca*ca*cb*cb - 9.*ca*cb*cb*sa*sa - 10.*cb*sa*sb + 18.*ca*ca*cb*sa*sb - 6.*cb*sa*sa*sa*sb + 5.*ca*sb*sb - 3.*ca*ca*ca*sb*sb + 9.*ca*sa*sa*sb*sb);
    y[i] += 0.125*ec*ec*AU*(-10.*sa + 5.*cb*cb*sa + 9.*ca*ca*cb*cb*sa - 3.*cb*cb*sa*sa*sa - 10.*ca*cb*sb - 6.*ca*ca*ca*cb*sb + 18.*ca*cb*sa*sa*sb - 5.*sa*sb*sb - 9.*ca*ca*sa*sb*sb + 3.*sa*sa*sa*sb*sb);
    z[i] += sq3*AU*ec*ec*(1. + (sa*cb - sb*ca)*(sa*cb - sb*ca));
    */
}
void initialize_analytic_orbit(struct Orbit *orbit)
{
    //store armlenght & transfer frequency in orbit structure.
    orbit->L        = LARM;
    orbit->fstar    = CLIGHT/(2.0*M_PI*LARM);
    orbit->ecc      = LARM/(2.0*SQ3*AU);
    orbit->R        = AU*orbit->ecc;
    orbit->kappa_0  = 0.0;
    orbit->lambda_0 = 0.0;
    orbit->orbit_function = &analytic_orbits;
    
}

void initialize_numeric_orbit(struct Orbit *orbit)
{
    fprintf(stdout,"==== Initialize LISA Orbit Structure ====\n\n");
    
    int n,i,check;
    double junk;
    
    FILE *infile = fopen(orbit->OrbitFileName,"r");
    
    //how big is the file
    n=0;
    while(!feof(infile))
    {
        /*columns of orbit file:
         t sc1x sc1y sc1z sc2x sc2y sc2z sc3x sc3y sc3z
         */
        n++;
        check = fscanf(infile,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk,&junk);
        if(!check)
        {
            fprintf(stderr,"Failure reading %s\n",orbit->OrbitFileName);
            exit(1);
        }
    }
    n--;
    rewind(infile);
    
    //allocate orbit structure
    alloc_orbit(orbit, n);

    //allocate memory for local workspace
    double *t   = double_vector(orbit->Norb);
    double **x  = double_matrix(3,orbit->Norb);
    double **y  = double_matrix(3,orbit->Norb);
    double **z  = double_matrix(3,orbit->Norb);
    
    //read in orbits
    for(n=0; n<orbit->Norb; n++)
    {
        //First time sample must be at t=0 for phasing
        check = fscanf(infile,"%lg",&t[n]);
        for(i=0; i<3; i++) check = fscanf(infile,"%lg %lg %lg",&x[i][n],&y[i][n],&z[i][n]);
        
        if(!check)
        {
            fprintf(stderr,"Failure reading %s\n",orbit->OrbitFileName);
            exit(1);
        }

        orbit->t[n] = t[n];
        
        //Repackage orbit positions into arrays for interpolation
        for(i=0; i<3; i++)
        {
            orbit->x[i][n] = x[i][n];
            orbit->y[i][n] = y[i][n];
            orbit->z[i][n] = z[i][n];
        }
    }
    fclose(infile);
    
    //calculate derivatives for cubic spline
    for(i=0; i<3; i++)
    {
        initialize_cubic_spline(orbit->dx[i], t, orbit->x[i]);
        initialize_cubic_spline(orbit->dy[i], t, orbit->y[i]);
        initialize_cubic_spline(orbit->dz[i], t, orbit->z[i]);
    }
    
    //calculate average arm length
    printf("Estimating average armlengths -- assumes evenly sampled orbits\n\n");
    double L12=0.0;
    double L23=0.0;
    double L31=0.0;
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    for(n=0; n<orbit->Norb; n++)
    {
        x1 = orbit->x[0][n];
        x2 = orbit->x[1][n];
        x3 = orbit->x[2][n];
        
        y1 = orbit->y[0][n];
        y2 = orbit->y[1][n];
        y3 = orbit->y[2][n];
        
        z1 = orbit->z[0][n];
        z2 = orbit->z[1][n];
        z3 = orbit->z[2][n];
        
        
        //L12
        L12 += sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
        
        //L23
        L23 += sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) );
        
        //L31
        L31 += sqrt( (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3) );
    }
    L12 /= (double)orbit->Norb;
    L31 /= (double)orbit->Norb;
    L23 /= (double)orbit->Norb;
    
    printf("Average arm lengths for the constellation:\n");
    printf("  L12 = %g\n",L12);
    printf("  L31 = %g\n",L31);
    printf("  L23 = %g\n",L23);
    printf("\n");
    
    //are the armlenghts consistent?
    double L = (L12+L31+L23)/3.;
    printf("Fractional deviation from average armlength for each side:\n");
    printf("  L12 = %g\n",fabs(L12-L)/L);
    printf("  L31 = %g\n",fabs(L31-L)/L);
    printf("  L23 = %g\n",fabs(L23-L)/L);
    printf("\n");
    
    //store armlenght & transfer frequency in orbit structure.
    orbit->L     = L;
    orbit->fstar = CLIGHT/(2.0*M_PI*L);
    orbit->ecc   = L/(2.0*SQ3*AU);
    orbit->R     = AU*orbit->ecc;
    orbit->orbit_function = &interpolate_orbits;
    
    //free local memory
    free_double_vector(t);
    free_double_matrix(x,3);
    free_double_matrix(y,3);
    free_double_matrix(z,3);
    fprintf(stdout,"=========================================\n\n");
    
}

void initialize_interpolated_analytic_orbits(struct Orbit *orbit, double Tobs, double t0)
{
    fprintf(stdout,"==== Initialize LISA Orbit Structure ====\n\n");
    
    //store armlength & transfer frequency in orbit structure.
    orbit->L     = LARM;
    orbit->fstar = CLIGHT/(2.0*M_PI*LARM);
    orbit->ecc   = LARM/(2.0*SQ3*AU);
    orbit->R     = AU*orbit->ecc;
    orbit->orbit_function = &analytic_orbits;
    
    /*
    Setup coarse time sampling of slowly-evolving orbits 
    */
    int N = (int)(LISA_ORBIT_SAMPLES_PER_YEAR*Tobs/YEAR); 
    if (N < 20) N = 20; //but at least 20 samples

    double buffer = (Tobs)/(double)(N-1); // buffer for time sampling to keep interpolation away from edges of data
    double cadence = (Tobs + 2.0*buffer)/(double)(N-1);  // coarse sample cadence
    
    //allocate orbit structure
    alloc_orbit(orbit, N);

    //allocate memory for local workspace
    double t  = 0.0;
    double *x = double_vector(4);
    double *y = double_vector(4);
    double *z = double_vector(4);

    
    //compute orbits
    for(int n=0; n<orbit->Norb; n++)
    {
        t = t0 + n*cadence - buffer;
        (orbit->orbit_function)(orbit,t,x,y,z);

        //Repackage orbit positions into arrays for interpolation
        orbit->t[n] = t;
        for(int i=0; i<3; i++)
        {
            orbit->x[i][n] = x[i+1];
            orbit->y[i][n] = y[i+1];
            orbit->z[i][n] = z[i+1];
        }
    }
    
    //calculate derivatives for cubic spline
    for(int i=0; i<3; i++)
    {
        initialize_cubic_spline(orbit->dx[i],orbit->t,orbit->x[i]);
        initialize_cubic_spline(orbit->dy[i],orbit->t,orbit->y[i]);
        initialize_cubic_spline(orbit->dz[i],orbit->t,orbit->z[i]);
    }
    
    //calculate average arm length
    printf("Estimating average armlengths -- assumes evenly sampled orbits\n\n");
    double L12=0.0;
    double L23=0.0;
    double L31=0.0;
    double x1,x2,x3,y1,y2,y3,z1,z2,z3;
    for(int n=0; n<orbit->Norb; n++)
    {
        x1 = orbit->x[0][n];
        x2 = orbit->x[1][n];
        x3 = orbit->x[2][n];
        
        y1 = orbit->y[0][n];
        y2 = orbit->y[1][n];
        y3 = orbit->y[2][n];
        
        z1 = orbit->z[0][n];
        z2 = orbit->z[1][n];
        z3 = orbit->z[2][n];
        
        
        //L12
        L12 += sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
        
        //L23
        L23 += sqrt( (x3-x2)*(x3-x2) + (y3-y2)*(y3-y2) + (z3-z2)*(z3-z2) );
        
        //L31
        L31 += sqrt( (x1-x3)*(x1-x3) + (y1-y3)*(y1-y3) + (z1-z3)*(z1-z3) );
    }
    L12 /= (double)orbit->Norb;
    L31 /= (double)orbit->Norb;
    L23 /= (double)orbit->Norb;
    
    printf("Average arm lengths for the constellation:\n");
    printf("  L12 = %g\n",L12);
    printf("  L31 = %g\n",L31);
    printf("  L23 = %g\n",L23);
    printf("\n");
    
    //are the armlenghts consistent?
    double L = (L12+L31+L23)/3.;
    printf("Fractional deviation from average armlength for each side:\n");
    printf("  L12 = %g\n",fabs(L12-L)/L);
    printf("  L31 = %g\n",fabs(L31-L)/L);
    printf("  L23 = %g\n",fabs(L23-L)/L);
    printf("\n");
    
    
    //free local memory
    free(x);
    free(y);
    free(z);
    fprintf(stdout,"=========================================\n\n");
    
}

void free_orbit(struct Orbit *orbit)
{
    free_double_vector(orbit->t);
    free_double_matrix(orbit->x,3);
    free_double_matrix(orbit->y,3);
    free_double_matrix(orbit->z,3);

    for(int i=0; i<3; i++)
    {
        free_cubic_spline(orbit->dx[i]);
        free_cubic_spline(orbit->dy[i]);
        free_cubic_spline(orbit->dz[i]);
    }

    free(orbit->dx);
    free(orbit->dy);
    free(orbit->dz);
    
    free(orbit);
}

static void recursive_phase_evolution(double dre, double dim, double *cosPhase, double *sinPhase)
{
   /* Update re and im for the next iteration. */
   double cosphi = *cosPhase;
   double sinphi = *sinPhase;

   double newRe = cosphi*dre - sinphi*dim;
   double newIm = sinphi*dre + cosphi*dim;

   *cosPhase = newRe;
   *sinPhase = newIm;
   
}

static void double_angle(double cosA, double sinA, double *cos2A, double *sin2A)
{
    *cos2A = 2.*cosA*cosA - 1.0;
    *sin2A = 2.*sinA*cosA;
}

static void triple_angle(double cosA, double sinA, double cos2A, double *cos3A, double *sin3A)
{
    *cos3A = 2.*cosA*cos2A - cosA;
    *sin3A = 2.*sinA*cos2A + sinA;
}

void XYZ2AE(double X, double Y, double Z, double *A, double *E)
{
    /*
     * Conventions in ldc/common/series/tdi.py
     * A = (Z-X)/sqrt(2)
     * E = (X-2Y+Z)/sqrt(6)
     */
    double invSQ6 = 0.408248290463863;
    
    *A = (Z-X)*M_SQRT1_2;
    *E = (X-2*Y+Z)*invSQ6;

}

void XYZ2AET(double X, double Y, double Z, double *A, double *E, double *T)
{
    /*
     * Conventions in ldc/common/series/tdi.py
     * T = (X+Y+Z)/sqrt(3)
     */
    double invSQ3 = 0.577350269189626;

    XYZ2AE(X,Y,Z,A,E);

    *T = (X+Y+Z)*invSQ3;

    
}

void LISA_tdi(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI)
{
    int i,j,k;
    int BWon2 = BW/2;
    double fonfs;
    double c3, s3, c2, s2, c1, s1;
    double f;
    double sqT=sqrt(T);
    double invfstar = 1./fstar;
    double dPhi = invfstar/T;
    double cosdPhi = cos(dPhi);
    double sindPhi = sin(dPhi);
    double prefactor;
    
    /* Initialize recursion */
    i = 0;
    f = ((double)(q + i-1 - BWon2))/T;
    fonfs = f*invfstar;
    
    c1 = cos(fonfs);
    s1 = sin(fonfs);

    double_angle(c1,s1,&c2,&s2);
    triple_angle(c1,s1,c2,&c3,&s3);

    
    for(i=1; i<=BW; i++)
    {
        k = 2*i;
        j = k-1;
        
        f = ((double)(q + i-1 - BWon2))/T;
        fonfs = f*invfstar;
        
        prefactor = sqT;

        /* make use of recursion relationships & identities to get rid of trig calls
        c3 = cos(3.*fonfs);  c2 = cos(fonfs2);  c1 = cos(fonfs);
        s3 = sin(3.*fonfs);  s2 = sin(fonfs2);  s1 = sin(fonfs);
         */

        recursive_phase_evolution(cosdPhi, sindPhi, &c1, &s1);
        double_angle(c1,s1,&c2,&s2);
        triple_angle(c1,s1,c2,&c3,&s3);

        X[j] =	prefactor*((d[1][2][j]-d[1][3][j])*c3 + (d[1][2][k]-d[1][3][k])*s3 +
        (d[2][1][j]-d[3][1][j])*c2 + (d[2][1][k]-d[3][1][k])*s2 +
        (d[1][3][j]-d[1][2][j])*c1 + (d[1][3][k]-d[1][2][k])*s1 +
        (d[3][1][j]-d[2][1][j]));
        
        X[k] =	prefactor*((d[1][2][k]-d[1][3][k])*c3 - (d[1][2][j]-d[1][3][j])*s3 +
        (d[2][1][k]-d[3][1][k])*c2 - (d[2][1][j]-d[3][1][j])*s2 +
        (d[1][3][k]-d[1][2][k])*c1 - (d[1][3][j]-d[1][2][j])*s1 +
        (d[3][1][k]-d[2][1][k]));
        
        /* Michelson channel for single IFO
        M[j] = prefactor*(X[j]*cLS - X[k]*sLS);
        M[k] =-prefactor*(X[j]*sLS + X[k]*cLS);
         */
        
        //save some CPU time when only X-channel is needed
        if(NI>1)
        {
            Y[j] =	prefactor*((d[2][3][j]-d[2][1][j])*c3 + (d[2][3][k]-d[2][1][k])*s3 +
            (d[3][2][j]-d[1][2][j])*c2 + (d[3][2][k]-d[1][2][k])*s2+
            (d[2][1][j]-d[2][3][j])*c1 + (d[2][1][k]-d[2][3][k])*s1+
            (d[1][2][j]-d[3][2][j]));
            
            Y[k] =	prefactor*((d[2][3][k]-d[2][1][k])*c3 - (d[2][3][j]-d[2][1][j])*s3+
            (d[3][2][k]-d[1][2][k])*c2 - (d[3][2][j]-d[1][2][j])*s2+
            (d[2][1][k]-d[2][3][k])*c1 - (d[2][1][j]-d[2][3][j])*s1+
            (d[1][2][k]-d[3][2][k]));
            
            Z[j] =	prefactor*((d[3][1][j]-d[3][2][j])*c3 + (d[3][1][k]-d[3][2][k])*s3+
            (d[1][3][j]-d[2][3][j])*c2 + (d[1][3][k]-d[2][3][k])*s2+
            (d[3][2][j]-d[3][1][j])*c1 + (d[3][2][k]-d[3][1][k])*s1+
            (d[2][3][j]-d[1][3][j]));
            
            Z[k] =	prefactor*((d[3][1][k]-d[3][2][k])*c3 - (d[3][1][j]-d[3][2][j])*s3+
            (d[1][3][k]-d[2][3][k])*c2 - (d[1][3][j]-d[2][3][j])*s2+
            (d[3][2][k]-d[3][1][k])*c1 - (d[3][2][j]-d[3][1][j])*s1+
            (d[2][3][k]-d[1][3][k]));
            
            /* LDC conventions for A & E channels */
            XYZ2AE(X[j],Y[j],Z[j],&A[j],&E[j]);
            XYZ2AE(X[k],Y[k],Z[k],&A[k],&E[k]);
            
            A[k] = -A[k];
            E[k] = -E[k];
        }
    }
}

void LISA_tdi_FF(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI)
{
    int i,j,k;
    int BWon2 = BW/2;
    double fonfs,fonfs2;
    double c3, s3, c2, s2, c1, s1;
    double f;
    double sqT=sqrt(T);
    double invfstar = 1./fstar;
    double dPhi = invfstar/T;
    double cosdPhi = cos(dPhi);
    double sindPhi = sin(dPhi);
    double prefactor;
    
    
    /* Initialize recursion */
    i = 0;
    f = ((double)(q + i-1 - BWon2))/T;
    fonfs = f*invfstar;
    
    c1 = cos(fonfs);
    s1 = sin(fonfs);

    double_angle(c1,s1,&c2,&s2);
    triple_angle(c1,s1,c2,&c3,&s3);

    for(i=1; i<=BW; i++)
    {
        k = 2*i;
        j = k-1;
        
        f = ((double)(q + i-1 - BWon2))/T;
        fonfs = f*invfstar;
        fonfs2= 2.*fonfs;
        
        /* prefactor = sqrt(T) * 4 * (f/fstar) * sin(f/fstar) */
        prefactor = sqT*fonfs2;

        /* make use of recursion relationships & identities to get rid of trig calls
        c3 = cos(3.*fonfs);  c2 = cos(fonfs2);  c1 = cos(fonfs);
        s3 = sin(3.*fonfs);  s2 = sin(fonfs2);  s1 = sin(fonfs);
         */

        recursive_phase_evolution(cosdPhi, sindPhi, &c1, &s1);
        double_angle(c1,s1,&c2,&s2);
        triple_angle(c1,s1,c2,&c3,&s3);
        
        
        X[j] =	prefactor*((d[1][2][j]-d[1][3][j])*c3 + (d[1][2][k]-d[1][3][k])*s3 +
        (d[2][1][j]-d[3][1][j])*c2 + (d[2][1][k]-d[3][1][k])*s2 +
        (d[1][3][j]-d[1][2][j])*c1 + (d[1][3][k]-d[1][2][k])*s1 +
        (d[3][1][j]-d[2][1][j]));
        
        X[k] =	prefactor*((d[1][2][k]-d[1][3][k])*c3 - (d[1][2][j]-d[1][3][j])*s3 +
        (d[2][1][k]-d[3][1][k])*c2 - (d[2][1][j]-d[3][1][j])*s2 +
        (d[1][3][k]-d[1][2][k])*c1 - (d[1][3][j]-d[1][2][j])*s1 +
        (d[3][1][k]-d[2][1][k]));
        
        /* Michelson channels for single IFO
        M[j] = prefactor*(X[j]*cSL - X[k]*sSL);
        M[k] = prefactor*(X[j]*sSL + X[k]*cSL);
         */
        
        //save some CPU time when only X-channel is needed
        if(NI>1)
        {
            Y[j] =	prefactor*((d[2][3][j]-d[2][1][j])*c3 + (d[2][3][k]-d[2][1][k])*s3 +
            (d[3][2][j]-d[1][2][j])*c2 + (d[3][2][k]-d[1][2][k])*s2+
            (d[2][1][j]-d[2][3][j])*c1 + (d[2][1][k]-d[2][3][k])*s1+
            (d[1][2][j]-d[3][2][j]));
            
            Y[k] =	prefactor*((d[2][3][k]-d[2][1][k])*c3 - (d[2][3][j]-d[2][1][j])*s3+
            (d[3][2][k]-d[1][2][k])*c2 - (d[3][2][j]-d[1][2][j])*s2+
            (d[2][1][k]-d[2][3][k])*c1 - (d[2][1][j]-d[2][3][j])*s1+
            (d[1][2][k]-d[3][2][k]));
            
            Z[j] =	prefactor*((d[3][1][j]-d[3][2][j])*c3 + (d[3][1][k]-d[3][2][k])*s3+
            (d[1][3][j]-d[2][3][j])*c2 + (d[1][3][k]-d[2][3][k])*s2+
            (d[3][2][j]-d[3][1][j])*c1 + (d[3][2][k]-d[3][1][k])*s1+
            (d[2][3][j]-d[1][3][j]));
            
            Z[k] =	prefactor*((d[3][1][k]-d[3][2][k])*c3 - (d[3][1][j]-d[3][2][j])*s3+
            (d[1][3][k]-d[2][3][k])*c2 - (d[1][3][j]-d[2][3][j])*s2+
            (d[3][2][k]-d[3][1][k])*c1 - (d[3][2][j]-d[3][1][j])*s1+
            (d[2][3][k]-d[1][3][k]));
            
            /* LDC conventions for A & E channels */
            XYZ2AE(X[j],Y[j],Z[j],&A[j],&E[j]);
            XYZ2AE(X[k],Y[k],Z[k],&A[k],&E[k]);

        }
    }
}

void LISA_tdi_Sangria(double L, double fstar, double T, double ***d, double f0, long q, double *X, double *Y, double *Z, double *A, double *E, int BW, int NI)
{
    int i,j,k;
    int BWon2 = BW/2;
    double prefactor,fonfs;
    double c2, s2, c1, s1;
    double f;
    double sqT=sqrt(T);
    double invfstar = 1./fstar;
    double norm = 4.0*invfstar*sqT;
    double dPhi = invfstar/T;
    double cosdPhi = cos(dPhi);
    double sindPhi = sin(dPhi);
        
    /* Initialize recursion */
    i = 0;
    f = ((double)(q + i-1 - BWon2))/T;
    fonfs = f*invfstar;
    
    c1 = cos(fonfs);
    s1 = sin(fonfs);

    double_angle(c1,s1,&c2,&s2);

    for(i=1; i<=BW; i++)
    {
        k = 2*i;
        j = k-1;

        f = ((double)(q + i-1 - BWon2))/T;

        /* make use of recursion relationships & identities to get rid of trig calls
        c3 = cos(3.*fonfs);  c2 = cos(2.*fonfs);  c1 = cos(fonfs);
        s3 = sin(3.*fonfs);  s2 = sin(2.*fonfs);  s1 = sin(fonfs);*/
        recursive_phase_evolution(cosdPhi, sindPhi, &c1, &s1);
        double_angle(c1,s1,&c2,&s2);
                
        /* prefactor = sqrt(T) * 4 * (f/fstar) * sin(f/fstar) */
        prefactor = f * norm * s1;
        
        
        /* LDC Conventions for X,Y,Z channels (circa Sangria) */
        X[j] = prefactor*((d[1][2][j]-d[1][3][j])*c2 + (d[1][2][k]-d[1][3][k])*s2 +
                   (d[2][1][j]-d[3][1][j])*c1 + (d[2][1][k]-d[3][1][k])*s1);
        
        X[k] = prefactor*((d[1][2][k]-d[1][3][k])*c2 - (d[1][2][j]-d[1][3][j])*s2 +
                   (d[2][1][k]-d[3][1][k])*c1 - (d[2][1][j]-d[3][1][j])*s1);
        
        Y[j] = prefactor*((d[2][3][j]-d[2][1][j])*c2 + (d[2][3][k]-d[2][1][k])*s2 +
                   (d[3][2][j]-d[1][2][j])*c1 + (d[3][2][k]-d[1][2][k])*s1);
        
        Y[k] = prefactor*((d[2][3][k]-d[2][1][k])*c2 - (d[2][3][j]-d[2][1][j])*s2+
                   (d[3][2][k]-d[1][2][k])*c1 - (d[3][2][j]-d[1][2][j])*s1);
        
        Z[j] = prefactor*((d[3][1][j]-d[3][2][j])*c2 + (d[3][1][k]-d[3][2][k])*s2+
                   (d[1][3][j]-d[2][3][j])*c1 + (d[1][3][k]-d[2][3][k])*s1);
        
        Z[k] = prefactor*((d[3][1][k]-d[3][2][k])*c2 - (d[3][1][j]-d[3][2][j])*s2+
                   (d[1][3][k]-d[2][3][k])*c1 - (d[1][3][j]-d[2][3][j])*s1);
        

        /* LDC conventions for A & E channels */
        XYZ2AE(X[j],Y[j],Z[j],&A[j],&E[j]);
        XYZ2AE(X[k],Y[k],Z[k],&A[k],&E[k]);
    }
}

void LISA_detector_time(struct Orbit *orbit, double costh, double phi, double *time, int N,  double *time_sc)
{
    // sky location of source
    double sinth, cosph, sinph;
    
    // location of spacecraft 0
    double r[3];
    
    //Calculate cos and sin of sky position
    sinth = sqrt(1.0-costh*costh);
    cosph = cos(phi);
    sinph = sin(phi);
    
    //source propogation vector
    double k[3] = {-sinth*cosph,-sinth*sinph,-costh};
    
    for(int n=0; n<N; n++)
    {
        
        r[0] = spline_interpolation(orbit->dx[0],time[n])/CLIGHT;
        r[1] = spline_interpolation(orbit->dy[0],time[n])/CLIGHT;
        r[2] = spline_interpolation(orbit->dz[0],time[n])/CLIGHT;
        
        // k dot r_0
        double kdotr = 0.0;
        for(int i=0; i<3; i++) kdotr += k[i]*r[i];
        
        // shift time reference to spacecraft 0
        time_sc[n] = time[n] - kdotr;
    }
}

static void hplus_and_hcross(double t, double phase, double amp,  double Aplus, double Across, double cos2psi, double sin2psi,  double *hp, double *hc, double *hpf, double *hcf)
{
    double cp = cos(phase);
    double sp = sin(phase);

    *hp  = amp * ( Aplus*cos2psi*cp  + Across*sin2psi*sp );
    *hc  = amp * ( Across*cos2psi*sp - Aplus*sin2psi*cp  );
    
    *hpf = amp * (-Aplus*cos2psi*sp  + Across*sin2psi*cp );
    *hcf = amp * ( Across*cos2psi*cp + Aplus*sin2psi*sp  );              
}

static void LISA_TDI_spline(double *M, double *Mf, int a, int b, int c, double* tarray, int n, struct CubicSpline *amp_spline, struct CubicSpline *freq_spline, struct CubicSpline *phase_spline, double Aplus, double Across, double cos2psi, double sin2psi, double *App, double *Apm, double *Acp, double *Acm, double *kr, double *Larm)
{
    double t, f, amp, phase;
    double hp, hc, hpf, hcf;

    M[n] = 0.0;
    Mf[n] = 0.0;
    
    if(freq_spline) /* mbh */
    {
        amp = spline_interpolation(amp_spline,  tarray[n]);
        f   = spline_interpolation(freq_spline, tarray[n]);
    }

    t = tarray[n] - kr[a]-2.0*Larm[c]-2.0*Larm[b];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] += hp*Apm[b]+hc*Acm[b];
    M[n] -= hp*App[c]+hc*Acp[c];
    Mf[n] += hpf*Apm[b]+hcf*Acm[b];
    Mf[n] -= hpf*App[c]+hcf*Acp[c];
    
    t = tarray[n] - kr[b]-Larm[c]-2.0*Larm[b];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] -= hp*Apm[c]+hc*Acm[c];
    M[n] += hp*App[c]+hc*Acp[c];
    Mf[n] -= hpf*Apm[c]+hcf*Acm[c];
    Mf[n] += hpf*App[c]+hcf*Acp[c];
    
    t = tarray[n] - kr[c]-Larm[b]-2.0*Larm[c];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] += hp*App[b]+hc*Acp[b];
    M[n] -= hp*Apm[b]+hc*Acm[b];
    Mf[n] += hpf*App[b]+hcf*Acp[b];
    Mf[n] -= hpf*Apm[b]+hcf*Acm[b];
    
    t = tarray[n] - kr[a]-2.0*Larm[b];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] -= hp*Apm[b]+hc*Acm[b];
    M[n] += hp*Apm[c]+hc*Acm[c];
    Mf[n] -= hpf*Apm[b]+hcf*Acm[b];
    Mf[n] += hpf*Apm[c]+hcf*Acm[c];
    
    t = tarray[n] - kr[a]-2.0*Larm[c];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] += hp*App[c]+hc*Acp[c];
    M[n] -= hp*App[b]+hc*Acp[b];
    Mf[n] += hpf*App[c]+hcf*Acp[c];
    Mf[n] -= hpf*App[b]+hcf*Acp[b];
    
    t = tarray[n] - kr[c]- Larm[b];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] -= hp*App[b]+hc*Acp[b];
    M[n] += hp*Apm[b]+hc*Acm[b];
    Mf[n] -= hpf*App[b]+hcf*Acp[b];
    Mf[n] += hpf*Apm[b]+hcf*Acm[b];
    
    t = tarray[n] - kr[b]- Larm[c];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] += hp*Apm[c]+hc*Acm[c];
    M[n] -= hp*App[c]+hc*Acp[c];
    Mf[n] += hpf*Apm[c]+hcf*Acm[c];
    Mf[n] -= hpf*App[c]+hcf*Acp[c];
    
    t = tarray[n] - kr[a];
    if(phase_spline) /* ucb */
    {
        amp   = spline_interpolation(amp_spline,   t);
        phase = spline_interpolation(phase_spline, t);
    }
    if(freq_spline) phase = PI2*f*t; /* mbh */
    hplus_and_hcross(t, phase, amp, Aplus, Across, cos2psi, sin2psi, &hp, &hc, &hpf, &hcf);
    M[n] -= hp*Apm[c]+hc*Acm[c];
    M[n] += hp*App[b]+hc*Acp[b];
    Mf[n] -= hpf*Apm[c]+hcf*Acm[c];
    Mf[n] += hpf*App[b]+hcf*Acp[b];    
}

void LISA_polarization_tensor_njc(double costh, double phi, double eplus[4][4], double ecross[4][4], double k[4])
{
    /*   Gravitational Wave basis vectors   */
    double u[4],v[4];

    /*   Sky location angles   */
    double sinth = sqrt(1.0 - costh*costh);
    double cosph = cos(phi);
    double sinph = sin(phi);

    /*   Tensor construction for building slowly evolving LISA response   */
    //Gravitational Wave source basis vectors
    u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
    k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
    
    //GW polarization basis tensors
    /*
     * NJC convention:
     */
    for(int i=1;i<=3;i++)
    {
        for(int j=1;j<=3;j++)
        {
            eplus[i][j]  = u[i]*u[j] - v[i]*v[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }

}

void LISA_spline_response(struct Orbit *orbit, double *tarray, int N, double costh, double phi, double cosi, double psi, struct CubicSpline *amp_spline, struct CubicSpline *freq_spline, struct CubicSpline *phase_spline, double *phase_ref, struct TDI *tdi_amp, struct TDI *tdi_phase)
{
    /* work space for TDI response */
    struct TDI *R = malloc(sizeof(struct TDI));
    struct TDI *Rf = malloc(sizeof(struct TDI)); //_f has a phase flip for building up the full response later?
    alloc_tdi(R,N,3);
    alloc_tdi(Rf,N,3);

    /* aliases for TDI response */
    double *X = R->X;
    double *Y = R->Y;
    double *Z = R->Z;
    double *Xf = Rf->X;
    double *Yf = Rf->Y;
    double *Zf = Rf->Z;

    /*   Indicies    */
    int i,j;

    double k[4]; 
    double kdotr[3], kdotn[3];
    double dplus[3], dcross[3];
    
    /*   Polarization basis tensors   */
    double eplus[4][4], ecross[4][4];
    
    double cos2psi, sin2psi;
    double t;
    double Aplus, Across;
    double *App, *Apm, *Acp, *Acm;
    
    /*   Allocating Arrays   */
            
    App = malloc(sizeof(double)*3); 
    Apm = malloc(sizeof(double)*3);
    Acp = malloc(sizeof(double)*3);
    Acm = malloc(sizeof(double)*3);
    
    cos2psi = cos(2.*psi);  
    sin2psi = sin(2.*psi);
    
    Aplus = 0.5*(1.+cosi*cosi);
    Across = -cosi;

    LISA_polarization_tensor(costh, phi, eplus, ecross, k);


    double x[3][3];
    double n[3][3];
    double L[3];

    for(int m=0; m<N; m++)
    {
        
        t = tarray[m];

        // position vectors for each spacecraft (converted to seconds)
        for(i=0; i<3; i++)
        {
            x[i][0] = spline_interpolation(orbit->dx[i],t)/CLIGHT;
            x[i][1] = spline_interpolation(orbit->dy[i],t)/CLIGHT;
            x[i][2] = spline_interpolation(orbit->dz[i],t)/CLIGHT;
        }

        // arm vectors
        for(i=0; i<3; i++)
        {
            n[0][i] = x[1][i] - x[2][i];
            n[1][i] = x[2][i] - x[0][i];
            n[2][i] = x[0][i] - x[1][i];
        }

        //arm lengths
        for(i=0; i<3; i++) 
        {
            L[i]=0.0;
            for(j=0; j<3; j++)  L[i] += n[i][j]*n[i][j];
            L[i] = sqrt(L[i]);
        }

        //normalize arm vectors        
        for(i=0; i<3; i++)
            for(j=0;j<3;j++)
                n[i][j] /= L[i];
        


        // k dot r_i (source direction spacecraft locations)
        for(i=0; i<3; i++)
        {
            kdotr[i] = 0.0;
            for(j=0; j<3; j++)
                kdotr[i] += k[j+1] * x[i][j];
        }
        
        // k dot n_i (source direction w/ arm vectors)
        for(i=0; i<3; i++)
        {
            kdotn[i] = 0.0;
            for(j=0; j<3; j++)
                kdotn[i] += k[j+1] * n[i][j];
        }        
        
        //Antenna primitives
        for(i=0; i<3; i++)
        {
            dplus[i] = 0.0;
            dcross[i] = 0.0;
            for(j=0; j<3; j++)
            {
                for(int l=0; l<3; l++)
                {
                    dplus[i]  += (n[i][j] * n[i][l]) * eplus[j+1][l+1];
                    dcross[i] += (n[i][j] * n[i][l]) * ecross[j+1][l+1];
                }
            }
        }
        
        //Full Antenna patterns
        for(i=0; i<3; i++)
        {
            App[i] = 0.5 * dplus[i]  / (1.0 + kdotn[i]);
            Apm[i] = 0.5 * dplus[i]  / (1.0 - kdotn[i]);
            Acp[i] = 0.5 * dcross[i] / (1.0 + kdotn[i]);
            Acm[i] = 0.5 * dcross[i] / (1.0 - kdotn[i]);
        }
        
        // build X, Y, Z responses
        LISA_TDI_spline(X, Xf, 0, 1, 2, tarray, m, amp_spline, freq_spline, phase_spline, Aplus, Across, cos2psi, sin2psi, App, Apm, Acp, Acm, kdotr, L);
        LISA_TDI_spline(Y, Yf, 1, 2, 0, tarray, m, amp_spline, freq_spline, phase_spline, Aplus, Across, cos2psi, sin2psi, App, Apm, Acp, Acm, kdotr, L);
        LISA_TDI_spline(Z, Zf, 2, 0, 1, tarray, m, amp_spline, freq_spline, phase_spline, Aplus, Across, cos2psi, sin2psi, App, Apm, Acp, Acm, kdotr, L);

    }
    
    /*
    Separate TDI responses back into terms of phase and amplitude
    */

    //extract_amplitude_and_phase() is removing the carrier phase
    extract_amplitude_and_phase(N, tdi_amp->X, tdi_phase->X, R->X, Rf->X, phase_ref);
    extract_amplitude_and_phase(N, tdi_amp->Y, tdi_phase->Y, R->Y, Rf->Y, phase_ref);
    extract_amplitude_and_phase(N, tdi_amp->Z, tdi_phase->Z, R->Z, Rf->Z, phase_ref);
    
    // remove any phase wraps
    unwrap_phase(N, tdi_phase->X);
    unwrap_phase(N, tdi_phase->Y);
    unwrap_phase(N, tdi_phase->Z);
    
    free(App);
    free(Apm);
    free(Acp); 
    free(Acm);
    
    free_tdi(R);
    free_tdi(Rf);

}

/*
 adapted from LDC git
 https://gitlab.in2p3.fr/LISA/LDC/-/blob/master/ldc/lisa/noise/noise.py
 commit #491bf4b3
 (c) LISA 2019
*/
void get_noise_levels(char model[], double f, double *Spm, double *Sop)
{
    if (strcmp(model, "radler") == 0)
    {
      *Spm = 2.5e-48*(1.0 + pow(f/1e-4,-2))/f/f;
      *Sop = 1.8e-37*(0.5)*(0.5)*f*f;
    }
    else if (strcmp(model, "scirdv1") == 0)
    {
        
        double fonc = PI2*f/CLIGHT;
        
        double DSoms_d = 2.25e-22;//15.e-12*15.e-12;
        double DSa_a = 9.00e-30;//3.e-15*3.e-15;
      
      double Sa_a = DSa_a * (1.0 +(0.4e-3/f)*(0.4e-3/f)) * (1.0 + pow((f/8.e-3),4)); // in acceleration

      double Sa_d = Sa_a*pow(PI2*f,-4.); // in displacement
      
      *Spm = Sa_d*(fonc)*(fonc);
      
      double Soms_d  = DSoms_d*(1. + pow(2.e-3/f, 4.) );
      double Soms_nu = Soms_d*(fonc)*(fonc);
      *Sop = Soms_nu;
    }
    else if (strcmp(model, "sangria") == 0)
    {
        //https://gitlab.in2p3.fr/LISA/LDC/-/blob/2f2d54b8b5c2d23d6251850f443e786dd9f7780a/data_generation/sangria/config.yml
        *Spm = 5.76e-30 / (PI2*f*CLIGHT)/(PI2*f*CLIGHT) * (1.0 + pow(0.4e-3/f,2)) * (1.0 + pow(f/8.0e-3,4));
        *Sop = 1.28e-22 * (PI2*f/CLIGHT)*(PI2*f/CLIGHT) * (1.0 + pow(2.0e-3/f,4));
    }
    /* more else if clauses */
    else /* default: */
    {
        fprintf(stderr,"Unrecognized noise model %s\n",model);
        exit(1);
    }
}

double XYZnoise_FF(double L, double fstar, double f, double Spm, double Sop)
{
    double x = f/fstar;
    double cosx  = cos(x);

    return 16. * noise_transfer_function(x) * ( 2.*(1.0 + cosx*cosx)*Spm + Sop );
}

double XYZcross_FF(double L, double fstar, double f, double Spm, double Sop)
{
    double x = f/fstar;
    double cosx  = cos(x);

    return -8. * noise_transfer_function(x) * cosx * ( 4.*Spm + Sop );
}

double AEnoise_FF(double L, double fstar, double f, double Spm, double Sop)
{
    double x = f/fstar;
    
    double cosx  = cos(x);
    double cos2x = cos(2.*x);
    
    return  8. * noise_transfer_function(x) * ( 2.*Spm*(3. + 2.*cosx + cos2x) + Sop*(2. + cosx) );
}

double Tnoise_FF(double L, double fstar, double f, double Spm, double Sop)
{
    double x = f/fstar;
    
    double cosx  = cos(x);

    return 16.0 * Sop * (1.0 - cosx) * noise_transfer_function(x) + 128.0 * Spm * noise_transfer_function(x) * sin(0.5*x)*sin(0.5*x)*sin(0.5*x)*sin(0.5*x);
}

/*
 adapted from LDC git
 https://gitlab.in2p3.fr/LISA/LDC/-/blob/master/ldc/lisa/noise/noise.py
 commit #491bf4b3
 (c) LISA 2019
*/
double GBnoise_FF(double T, double fstar, double f)
{
    double x = f/fstar;
    double t = 4.*x*x*sin(x)*sin(x);
    
    double Ampl = 1.2826e-44;
    double alpha = 1.629667;
    double fr2 = 4.810781e-4;
    double af1 = -2.235e-1;
    double bf1 = -2.7040844;
    double afk = -3.60976122e-1;
    double bfk = -2.37822436;
    
    double fr1 = pow(10., af1*log10(T/YEAR) + bf1);
    double knee = pow(10., afk*log10(T/YEAR) + bfk);
    double SG_sense = Ampl * exp(-pow(f/fr1,alpha)) * pow(f,-7./3.) * 0.5 * (1.0 + tanh(-(f-knee)/fr2) );
    double SGXYZ = t*SG_sense;
    double SGAE  = 1.5*SGXYZ;
    return SGAE;
}

static double rednoise(double f)
{
    return 1. + 16.0*(pow((2.0e-5/f), 10.0)+ ipow(1.0e-4/f,2));
}

double AEnoise(double L, double fstar, double f)
{
    //Power spectral density of the detector noise and transfer frequency
    
    double fonfstar = f/fstar;
    double trans = ipow(sin(fonfstar),2.0);
    
    return (16.0/3.0)*trans*( (2.0+cos(fonfstar))*(SPS + SLOC) + 2.0*( 3.0 + 2.0*cos(fonfstar) + cos(2.0*fonfstar) ) * ( SLOC/2.0 + SACC/ipow(PI2*f,4)*rednoise(f) ) ) / ipow(2.0*L,2);
    
    
}

double XYZnoise(double L, double fstar, double f)
{
    //Power spectral density of the detector noise and transfer frequency
    
    double fonfstar = f/fstar;
    double trans = ipow(sin(fonfstar),2.0);
    
    return (4.0)*trans*( (4.0)*(SPS + SLOC) + 8.0*( 1.0 + ipow(cos(fonfstar),2) ) * ( SLOC/2.0 + SACC*(1./(ipow(PI2*f,4)))*rednoise(f) ) ) / ipow(2.0*L,2);
}

double GBnoise(double T, double f)
{
    /* Fits to confusion noise from Cornish and Robson https://arxiv.org/pdf/1703.09858.pdf */
    double A = 1.8e-44;
    double alpha;
    double beta;
    double kappa;
    double gamma;
    double fk;
    
    //map T to number of years
    double Tyear = T/YEAR;
    
    if(Tyear>3)
    {
        alpha = 0.138;
        beta  = -221.0;
        kappa = 512.0;
        gamma = 1680.0;
        fk    = 0.00113;
    }
    else if(Tyear>1.5)
    {
        alpha = 0.165;
        beta  = 299.;
        kappa = 611.;
        gamma = 1340.;
        fk    = 0.00173;
    }
    else if(Tyear>0.75)
    {
        alpha = 0.171;
        beta  = 292.0;
        kappa = 1020.;
        gamma = 1680.0;
        fk    = 0.00215;
    }
    else
    {
        alpha = 0.133;
        beta  = 243.0;
        kappa = 482.0;
        gamma = 917.0;
        fk    = 0.00258;
    }
    return A*pow(f,-7./3.)*exp(-pow(f,alpha) + beta*f*sin(kappa*f))*(1. + tanh(gamma*(fk-f)));
}

double noise_transfer_function(double x)
{
    double sinx = sin(x);
    return sinx*sinx;
}

void test_noise_model(struct Orbit *orbit)
{
    double Spm,Sop;
    double fstart = log10(1e-5);
    double fstop = log10(0.1);
    int Nf = 1000;
    double df=(fstop-fstart)/(double)Nf;
    FILE *psdfile = fopen("psd.dat","w");
    for(int nf=0; nf<Nf; nf++)
    {
        double ftemp = pow(10,fstart + nf*df);
        get_noise_levels("radler", ftemp, &Spm, &Sop);
        fprintf(psdfile,"%lg %lg %lg\n",ftemp,AEnoise_FF(orbit->L,orbit->fstar,ftemp,Spm,Sop),Tnoise_FF(orbit->L,orbit->fstar,ftemp,Spm,Sop));
    }
    fclose(psdfile);
}

void alloc_tdi(struct TDI *tdi, int N, int Nchannel)
{
    //Number of frequency bins (2*N samples)
    tdi->N = N;
    
    //Michelson
    tdi->X = calloc(tdi->N,sizeof(double));
    tdi->Y = calloc(tdi->N,sizeof(double));
    tdi->Z = calloc(tdi->N,sizeof(double));
    
    //Noise-orthogonal
    tdi->A = calloc(tdi->N,sizeof(double));
    tdi->E = calloc(tdi->N,sizeof(double));
    tdi->T = calloc(tdi->N,sizeof(double));
        
    //Number of TDI channels (X or A&E or maybe one day A,E,&T)
    tdi->Nchannel = Nchannel;
}

void copy_tdi(struct TDI *origin, struct TDI *copy)
{
    copy->N        = origin->N;
    copy->Nchannel = origin->Nchannel;
    
    memcpy(copy->A, origin->A, origin->N*sizeof(double));
    memcpy(copy->E, origin->E, origin->N*sizeof(double));
    memcpy(copy->T, origin->T, origin->N*sizeof(double));
    memcpy(copy->X, origin->X, origin->N*sizeof(double));
    memcpy(copy->Y, origin->Y, origin->N*sizeof(double));
    memcpy(copy->Z, origin->Z, origin->N*sizeof(double));
}

void copy_tdi_segment(struct TDI *origin, struct TDI *copy, int index, int N)
{
    copy->N        = origin->N;
    copy->Nchannel = origin->Nchannel;
    index*=2;
    memcpy(copy->A+index, origin->A+index, N*sizeof(double));
    memcpy(copy->E+index, origin->E+index, N*sizeof(double));
    memcpy(copy->T+index, origin->T+index, N*sizeof(double));
    memcpy(copy->X+index, origin->X+index, N*sizeof(double));
    memcpy(copy->Y+index, origin->Y+index, N*sizeof(double));
    memcpy(copy->Z+index, origin->Z+index, N*sizeof(double));
}

void free_tdi(struct TDI *tdi)
{
    free(tdi->X);
    free(tdi->Y);
    free(tdi->Z);
    free(tdi->A);
    free(tdi->E);
    free(tdi->T);
    
    free(tdi);
}

/* LDC HDF5 */

#define DATASET "/obs/tdi"

void LISA_Read_HDF5_LDC_TDI(struct TDI *tdi, char *fileName, const char *dataName)
{    
    /* LDC-formatted structure for compound HDF5 dataset */
    typedef struct tdi_dataset {
        double t;
        double X;
        double Y;
        double Z;
    } tdi_dataset;
    static tdi_dataset *s1;
    
    hid_t  file, dataset, dspace; /* identifiers */
    int ndims; /* dimension of dataset */
    
    /* Open an existing file. */
    file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /* Open an existing dataset. */
    dataset = H5Dopen(file, dataName, H5P_DEFAULT);

    /* Get size of dataset */
    dspace = H5Dget_space(dataset);
    ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int Nsamples = (int)dims[0];
    
    s1 = malloc(Nsamples*sizeof(struct tdi_dataset));
    alloc_tdi(tdi, Nsamples, 3);
    
    hid_t s1_tid; /* Memory datatype handle */
    
    s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(struct tdi_dataset));
    H5Tinsert(s1_tid, "t", HOFFSET(struct tdi_dataset, t), H5T_IEEE_F64LE);
    H5Tinsert(s1_tid, "X", HOFFSET(struct tdi_dataset, X), H5T_IEEE_F64LE);
    H5Tinsert(s1_tid, "Y", HOFFSET(struct tdi_dataset, Y), H5T_IEEE_F64LE);
    H5Tinsert(s1_tid, "Z", HOFFSET(struct tdi_dataset, Z), H5T_IEEE_F64LE);
    
    H5Dread(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, s1);
        
    
    /* Copy LDC-formatted structure into ldasoft format */
    for(int i=0; i<Nsamples; i++)
    {
        double X = s1[i].X;
        double Y = s1[i].Y;
        double Z = s1[i].Z;
        
        tdi->X[i] = X;
        tdi->Y[i] = Y;
        tdi->Z[i] = Z;
        
        /* LDC conventions for AET channels */
        XYZ2AET(X,Y,Z,&tdi->A[i],&tdi->E[i],&tdi->T[i]);
        
    }
    tdi->delta = s1[1].t - s1[0].t;
    
    /* Close the dataset. */
    H5Dclose(dataset);
    
    /* Close the file. */
    H5Fclose(file);
    
    /* Free up memory */
    free(s1);

}

void LISA_Read_HDF5_LDC_RADLER_TDI(struct TDI *tdi, char *fileName)
{
    hid_t  file, dataset, dataspace, memtype; /* identifiers */
    int ndims, Nrow, Ncol; /* dimension of dataset */
    double *data; /* array for dataset */

    /* Open an existing file. */
    file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

    /* Open an existing dataset. */
    dataset = H5Dopen(file, "/H5LISA/PreProcess/TDIdata", H5P_DEFAULT);
    memtype = H5Dget_type(dataset);

    /* Get dimension of dataset */
    dataspace = H5Dget_space(dataset);
    ndims = H5Sget_simple_extent_ndims(dataspace);
    hsize_t adims[ndims];
    H5Sget_simple_extent_dims(dataspace, adims, NULL);

    /* Allocate memory for reading dataset */
    Nrow = (int)adims[0];
    Ncol = (int)adims[1];
    data = malloc(Nrow*Ncol*sizeof(double));
    
    alloc_tdi(tdi, Nrow/2, 3);

    
    /* Read the data */
    H5Dread (dataset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);

    
    for(int i=0; i<Nrow; i++)
    {
        
        double X = data[Ncol*i+1];
        double Y = data[Ncol*i+2];
        double Z = data[Ncol*i+3];
        
        tdi->X[i] = X;
        tdi->Y[i] = Y;
        tdi->Z[i] = Z;
        
        /* LDC conventions for AET channels */
        XYZ2AET(X,Y,Z,&tdi->A[i],&tdi->E[i],&tdi->T[i]);
    }
    tdi->delta = data[Ncol] - data[0];

    free(data);

    /* Close the dataset. */
    H5Dclose(dataset);

    /* Close the file. */
    H5Fclose(file);

}

void LISA_detector_tensor(double L, double eplus[4][4], double ecross[4][4], double x[4], double y[4], double z[4], double k[4], double dplus[4][4], double dcross[4][4], double kdotr[4][4])
{

    /*   Spacecraft position and separation vector   */
    double r12[4]={0},r13[4]={0},r23[4]={0};

    //Unit separation vector from spacecrafts i to j
    r12[1] = (x[2] - x[1])/L;   r13[1] = (x[3] - x[1])/L;   r23[1] = (x[3] - x[2])/L;
    r12[2] = (y[2] - y[1])/L;   r13[2] = (y[3] - y[1])/L;   r23[2] = (y[3] - y[2])/L;
    r12[3] = (z[2] - z[1])/L;   r13[3] = (z[3] - z[1])/L;   r23[3] = (z[3] - z[2])/L;
    
    
    //Zero arrays to be summed
    dplus[1][2]  = dplus[1][3]  = dplus[2][1]  = dplus[2][3]  = dplus[3][1]  = dplus[3][2]  = 0.;
    dcross[1][2] = dcross[1][3] = dcross[2][1] = dcross[2][3] = dcross[3][1] = dcross[3][2] = 0.;
    
    //Convenient quantities d+ & dx
    for(int i=1; i<=3; i++)
    {
        for(int j=1; j<=3; j++)
        {
            dplus[1][2]  += r12[i]*r12[j]*eplus[i][j];   dcross[1][2] += r12[i]*r12[j]*ecross[i][j];
            dplus[2][3]  += r23[i]*r23[j]*eplus[i][j];   dcross[2][3] += r23[i]*r23[j]*ecross[i][j];
            dplus[1][3]  += r13[i]*r13[j]*eplus[i][j];   dcross[1][3] += r13[i]*r13[j]*ecross[i][j];
        }
    }
    
    //Make use of symmetry
    dplus[2][1] = dplus[1][2];  dcross[2][1] = dcross[1][2];
    dplus[3][2] = dplus[2][3];  dcross[3][2] = dcross[2][3];
    dplus[3][1] = dplus[1][3];  dcross[3][1] = dcross[1][3];
    
    //Zero arrays to be summed
    kdotr[1][2] = kdotr[1][3] = kdotr[2][1] = kdotr[2][3] = kdotr[3][1] = kdotr[3][2] = 0.;
    for(int i=1; i<=3; i++)
    {
        kdotr[1][2] += k[i]*r12[i];   kdotr[1][3] += k[i]*r13[i];   kdotr[2][3] += k[i]*r23[i];
    }
    
    //Make use of antisymmetry
    kdotr[2][1] = -kdotr[1][2];
    kdotr[3][1] = -kdotr[1][3];
    kdotr[3][2] = -kdotr[2][3];
}

void LISA_polarization_tensor(double costh, double phi, double eplus[4][4], double ecross[4][4], double k[4])
{
    /*   Gravitational Wave basis vectors   */
    double u[4],v[4];

    /*   Sky location angles   */
    double sinth = sqrt(1.0 - costh*costh);
    double cosph = cos(phi);
    double sinph = sin(phi);

    /*   Tensor construction for building slowly evolving LISA response   */
    //Gravitational Wave source basis vectors
    u[1] =  costh*cosph;  u[2] =  costh*sinph;  u[3] = -sinth;
    v[1] =  sinph;        v[2] = -cosph;        v[3] =  0.;
    k[1] = -sinth*cosph;  k[2] = -sinth*sinph;  k[3] = -costh;
    
    //GW polarization basis tensors
    /*
     * LDC convention:
     * https://gitlab.in2p3.fr/LISA/LDC/-/blob/develop/ldc/waveform/fastGB/GB.cc
     */
    for(int i=1;i<=3;i++)
    {
        for(int j=1;j<=3;j++)
        {
            eplus[i][j]  = v[i]*v[j] - u[i]*u[j];
            ecross[i][j] = u[i]*v[j] + v[i]*u[j];
        }
    }

}


