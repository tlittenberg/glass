/*
 * Copyright 2020 Tyson B. Littenberg
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

/*
static void printProgress (double percentage)
{
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf ("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush (stdout);
}*/

void printUsage(const char *program)
{
    fprintf( stdout,"\n");
    fprintf( stdout, "Gaussian Mixture Model (GMM):\n");
    fprintf( stdout, "  Expectation-Maximization algorithm to fit GMM to \n");
    fprintf( stdout, "  input list of data samples e.g. from MCMC\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "Usage: %s required [optional] \n", program );
    fprintf( stdout, "\n");
    fprintf( stdout, "  Required:\n");
    fprintf( stdout, "     -f, --file=FILE        filename for input chain\n" );
    fprintf( stdout, "     -n, --nparams=INT      number of model parameters\n" );
    fprintf( stdout, "                            extra columns are ignored\n");
    fprintf( stdout, "\n");
    fprintf( stdout, "  Optional:\n");
    fprintf( stdout, "    [-h, --help]            print this message and exit\n" );
    fprintf( stdout, "    [-l, --log=INT]         use log(params) in specified column number [0,m-1]\n" );
    fprintf( stdout, "                            multiple arguments can be given.\n");
    fprintf( stdout, "    [-m, --modes=INT]       number of modes in fit (2)\n");
    fprintf( stdout, "    [-s, --seed=LONG]       RNG seed\n");
    fprintf( stdout, "    [-t, --thin=INT]        downsample rate for chains (1)\n");
    fprintf( stdout,"\n");
}

void read_gmm_binary(struct GMM *gmm, char filename[])
{
    FILE *fptr = NULL;
    /* Read GMM results to binary for pick up by other processes */
    if( (fptr = fopen(filename,"rb"))!=NULL)
    {
        fread(&gmm->NMODE, sizeof gmm->NMODE, 1, fptr);
        
        gmm->modes = malloc(gmm->NMODE*sizeof(struct MVG*));
        for(size_t n=0; n<gmm->NMODE; n++)
        {
            gmm->modes[n] = malloc(sizeof(struct MVG));
            alloc_MVG(gmm->modes[n],gmm->NParams);
        }
        
        for(size_t n=0; n<gmm->NMODE; n++) read_MVG(gmm->modes[n],fptr);
        fclose(fptr);
    }
    else
    {
        fprintf(stderr,"Error reading %s\n",filename);
        fprintf(stderr,"Exiting to system\n");
        exit(1);
    }
}

void alloc_MVG(struct MVG *mode, size_t N)
{
    mode->size = N;

    mode->mu       = double_vector(mode->size);
    mode->C        = double_matrix(mode->size,mode->size);
    mode->L        = double_matrix(mode->size,mode->size);
    mode->Cinv     = double_matrix(mode->size,mode->size);
    mode->evalues  = double_vector(mode->size);
    mode->evectors = double_matrix(mode->size,mode->size);
    mode->minmax   = double_matrix(mode->size,2);
}

void free_MVG(struct MVG *mode)
{
    free_double_vector(mode->mu);
    free_double_matrix(mode->C,mode->size);
    free_double_matrix(mode->L,mode->size);
    free_double_matrix(mode->Cinv,mode->size);
    free_double_matrix(mode->evectors,mode->size);
    free_double_vector(mode->evalues);
    free_double_matrix(mode->minmax,2);
    free(mode);
}

static void glass_vector_fwrite(FILE *fptr, int n, double *v)
{
    fwrite(v, sizeof(double), (size_t)n, fptr);
}

static void glass_matrix_fwrite(FILE *fptr, int n, double **m)
{
    for(int i=0; i<n; i++) glass_vector_fwrite(fptr, n, m[i]);
}

static void glass_vector_fread(FILE *fptr, int n, double *v)
{
    fread(v, sizeof(double), (size_t)n, fptr);
}

static void glass_matrix_fread(FILE *fptr, int n, double **m)
{
    for(int i=0; i<n; i++) glass_vector_fread(fptr, n, m[i]);
}


void write_MVG(struct MVG *mode, FILE *fptr)
{
    //pack detC,p,and Neff into a vector to make life easier
    double *temp = double_vector(3);
    temp[0] = mode->detC;
    temp[1] = mode->p;
    temp[2] = mode->Neff;

    //write 'em!
    glass_vector_fwrite(fptr,mode->size,mode->mu);
    glass_matrix_fwrite(fptr,mode->size,mode->C);
    glass_matrix_fwrite(fptr,mode->size,mode->L);
    glass_matrix_fwrite(fptr,mode->size,mode->Cinv);
    glass_matrix_fwrite(fptr,mode->size,mode->evectors);
    glass_vector_fwrite(fptr,mode->size,mode->evalues);
    glass_vector_fwrite(fptr,3,temp);
    glass_matrix_fwrite(fptr,2,mode->minmax);
    free_double_vector(temp);
}

void read_MVG(struct MVG *mode, FILE *fptr)
{
    //vector for holding packed detC,p,and Neff
    double *temp = double_vector(3);

    //read 'em!
    glass_vector_fread(fptr,mode->size,mode->mu);
    glass_matrix_fread(fptr,mode->size,mode->C);
    glass_matrix_fread(fptr,mode->size,mode->L);
    glass_matrix_fread(fptr,mode->size,mode->Cinv);
    glass_matrix_fread(fptr,mode->size,mode->evectors);
    glass_vector_fread(fptr,mode->size,mode->evalues);
    glass_vector_fread(fptr,3,temp);
    glass_matrix_fread(fptr,2,mode->minmax);

    //unpack 'em!
    mode->detC = temp[0];
    mode->p    = temp[1];
    mode->Neff = temp[2];

    free_double_vector(temp);
}

static void glass_vector_memcpy(double *copy, double *origin, int n)
{
    memcpy(copy,origin,n*sizeof(double));
}
static void glass_matrix_memcpy(double **copy, double **origin, int n)
{
    for(int i=0; i<n; i++) glass_vector_memcpy(copy[i], origin[i], n);
}

void copy_MVG(struct MVG *origin, struct MVG *copy)
{
    copy->size = origin->size;
    glass_vector_memcpy(copy->mu, origin->mu, origin->size);
    glass_matrix_memcpy(copy->C, origin->C, origin->size);
    glass_matrix_memcpy(copy->L, origin->L, origin->size);
    glass_matrix_memcpy(copy->Cinv, origin->Cinv, origin->size);
    glass_matrix_memcpy(copy->evectors, origin->evectors, origin->size);
    glass_matrix_memcpy(copy->minmax, origin->minmax, origin->size);
    glass_vector_memcpy(copy->evalues, origin->evalues, origin->size);
    copy->detC = origin->detC;
    copy->p = origin->p;
    copy->Neff = origin->Neff;
}

double multivariate_gaussian(double *x, struct MVG *mvg, int N)
{
    
    double *dx = double_vector(N);
    double *mu = mvg->mu;
    double **Cinv = mvg->Cinv;
    double detC = mvg->detC;
    
    // x-mu
    for(size_t n=0; n<N; n++)
    {
        double xi = x[n];
        double mi = mu[n];
        double d  = xi-mi;
        
        dx[n]=d;
    }
    
    
    /*
     (x-mu)^T C^-1 (x-mu)
     */
    double CdotX;
    double chi2 = 0.0;
    for(size_t m=0; m<N; m++)
    {
        CdotX = 0.0;
        for(size_t n=0; n<N; n++)
        {
            CdotX += Cinv[m][n]*dx[n];
        }
        chi2 += CdotX*dx[m];
    }
    free_double_vector(dx);
    
    return exp(-0.5*chi2)/sqrt(pow(PI2,N)*detC);
}

double log_likelihood(struct MVG **modes, struct Sample **samples, int NMCMC, int NMODE)
{
    
    double logL = 0.0;
    for(size_t i=0; i<NMCMC; i++)
    {
        double P = PMIN;
        for(size_t k=0; k<NMODE; k++)
        {
            P += modes[k]->p*samples[i]->p[k];
        }
        logL += log(P);
    }
    return logL;
}

void print_1D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMODE, size_t NMCMC, char root[], size_t ix)
{
    char filename[128];
    sprintf(filename,"%s_%i.dat",root,(int)ix);
    FILE *fptr = fopen(filename,"w");
        
    //get original parameter boundaries
    double pmin = modes[0]->minmax[ix][0];
    double pmax = modes[0]->minmax[ix][1];
    
    double *xvec = double_vector(NMCMC);
    double xmin,xmax;
    double x0 =  1e60;
    double xf = -1e60;
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t i=0; i<NMCMC; i++) xvec[i] = samples[i]->x[ix];
        get_min_max(xvec, NMCMC, &xmin, &xmax);

        if(xmin<x0) x0 = xmin;
        if(xmax>xf) xf = xmax;
    }
    
    double p;
    double x;
    double dx = (xf-x0)/100.;
    for(int n=0; n<=100; n++)
    {
        p = 0.0;
        x = x0 + (double)n*dx;
        for(size_t k=0; k<NMODE; k++)
        {
            double mean = modes[k]->mu[ix];
            double var  = modes[k]->C[ix][ix];
            
            /*
             Probability density p is weight * normal / Jacobian
             */
            p += modes[k]->p * exp( -0.5*(x-mean)*(x-mean)/var )/sqrt(2*M_PI*var) / ((pmax-pmin)*exp(-x)/pow(1. + exp(-x),2));
        }
        
        
        fprintf(fptr,"%.16g %.16g\n",sigmoid(x,pmin,pmax),p);
    }
    
    fclose(fptr);
    
    
    free(xvec);
}

void print_2D_pdfs(struct MVG **modes, struct Sample **samples, size_t NMODE, size_t NMCMC, char root[], size_t ix, size_t iy)
{
    char filename[128];
    sprintf(filename,"%s_%i_%i.dat",root,(int)ix,(int)iy);
    FILE *fptr = fopen(filename,"w");
        
    //get original parameter boundaries
    double pmin_x = modes[0]->minmax[ix][0];
    double pmax_x = modes[0]->minmax[ix][1];
    double pmin_y = modes[0]->minmax[iy][0];
    double pmax_y = modes[0]->minmax[iy][1];

    double *xvec = double_vector(NMCMC);
    double *yvec = double_vector(NMCMC);
    double xmin,xmax;
    double ymin,ymax;
    double x0 =  1e60;
    double xf = -1e60;
    double y0 =  1e60;
    double yf = -1e60;
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t i=0; i<NMCMC; i++) xvec[i] = samples[i]->x[ix];
        get_min_max(xvec,NMCMC,&xmin,&xmax);
        if(xmin<x0) x0 = xmin;
        if(xmax>xf) xf = xmax;
        
        for(size_t i=0; i<NMCMC; i++) yvec[i] = samples[i]->x[iy];
        get_min_max(yvec,NMCMC,&ymin,&ymax);
        if(ymin<y0) y0 = ymin;
        if(ymax>yf) yf = ymax;
    }
    
    double p;
    double x,y;
    double dx = (pmax_x-pmin_x)/100.;
    double dy = (pmax_y-pmin_y)/100.;
    for(int n=1; n<100; n++)
    {
        x = pmin_x + (double)n*dx;
        x = logit(x,pmin_x,pmax_x);
        
        for(int m=1; m<100; m++)
        {
            p = 0.0;

            y = pmin_y + (double)m*dy;
            y = logit(y,pmin_y,pmax_y);
            
            for(size_t k=0; k<NMODE; k++)
            {
                double mean_x = modes[k]->mu[ix];
                double mean_y = modes[k]->mu[iy];
                if(n==0 && m==0 && ix==2 && iy==5)
                {
                    printf("means: %g %g -> %g %g\n",mean_x,mean_y, sigmoid(mean_x,pmin_x,pmax_x),sigmoid(mean_y,pmin_y,pmax_y));
                }
                double C_xx  = modes[k]->C[ix][ix];
                double C_xy  = modes[k]->C[ix][iy];
                double C_yy  = modes[k]->C[iy][iy];
                double detC = C_xx*C_yy - C_xy*C_xy;
                double iC_xx = C_yy/detC;
                double iC_yy = C_xx/detC;
                double iC_xy = -C_xy/detC;

                /*
                 Probability density p is weight * normal / Jacobian
                 */
                p += modes[k]->p * exp( -0.5*( (x-mean_x)*(x-mean_x)*iC_xx + (y-mean_y)*(y-mean_y)*iC_yy + 2*(x-mean_x)*(y-mean_y)*iC_xy) )  * 1./(2.*M_PI*sqrt(detC)) / ((pmax_x-pmin_x)*exp(-x)/pow(1. + exp(-x),2)) / ((pmax_y-pmin_y)*exp(-y)/pow(1. + exp(-y),2)) ;
            }
            
            
            fprintf(fptr,"%.16g %.16g %.16g\n",sigmoid(x,pmin_x,pmax_x),sigmoid(y,pmin_y,pmax_y),p);
        }
    }
    
    fclose(fptr);

    
    free_double_vector(xvec);
    free_double_vector(yvec);
}

void print_2D_contours(struct MVG **modes, size_t NMODE, char root[], size_t x1, size_t x2)
{
    char filename[128];
    char filename2[128];
    FILE *fptr = NULL;
    FILE *fptr2 = NULL;

    struct MVG **submodes = malloc(NMODE*sizeof(struct MVG*));
    for(size_t k=0; k<NMODE; k++)
    {
        submodes[k] = malloc(sizeof(struct MVG));
        alloc_MVG(submodes[k], 2);
    }
    
    //get original parameter boundaries
    double pmin[2] = {modes[0]->minmax[x1][0],modes[0]->minmax[x2][0]};
    double pmax[2] = {modes[0]->minmax[x1][1],modes[0]->minmax[x2][1]};

    
    //Pick parameters
    size_t X[2] = {x1,x2};
    
    for(size_t k=0; k<NMODE; k++)
    {
        for(size_t n=0; n<2; n++)
        {
            submodes[k]->mu[n] = modes[k]->mu[X[n]];
            for(size_t m=0; m<2; m++) submodes[k]->C[m][n] = modes[k]->C[X[m]][X[n]];
        }
        
        matrix_eigenstuff(submodes[k]->C, submodes[k]->evectors, submodes[k]->evalues, 2);
    }
    
    
    double x,y;
    for(size_t k=0; k<NMODE; k++)
    {
        sprintf(filename,"%s_%i_%i_%i.dat",root,(int)x1,(int)x2,(int)k);
        fptr=fopen(filename,"w");
        sprintf(filename2,"remapped_%s_%i_%i_%i.dat",root,(int)x1,(int)x2,(int)k);
        fptr2=fopen(filename2,"w");
        for(int n=0; n<=100; n++)
        {
            double theta = atan2( submodes[k]->evectors[0][1], submodes[k]->evectors[0][0]);
            double Rx = sqrt(submodes[k]->evalues[0]);
            double Ry = sqrt(submodes[k]->evalues[1]);
            double Cx = submodes[k]->mu[0];
            double Cy = submodes[k]->mu[1];
            
            double angle = n*(2.*M_PI/100.);
            x = 1.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 1.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 2.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 2.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            x = 3.*( Rx*cos(angle)*cos(theta) + Ry*sin(angle)*sin(theta) ) + Cx;
            y = 3.*(-Rx*cos(angle)*sin(theta) + Ry*sin(angle)*cos(theta) ) + Cy;
            fprintf(fptr2,"%.16lg %.16lg ",x,y);
            x = sigmoid(x,pmin[0],pmax[0]);
            y = sigmoid(y,pmin[1],pmax[1]);
            fprintf(fptr,"%.16lg %.16lg ",x,y);
            
            fprintf(fptr,"\n");
            fprintf(fptr2,"\n");
        }
        fclose(fptr);
        fclose(fptr2);
    }
    
    for(size_t k=0; k<NMODE; k++) free_MVG(submodes[k]);
    free(submodes);
}

void print_model(struct MVG **modes, struct Sample **samples, size_t NP, size_t NMODE, size_t NMCMC, double logL, double BIC, size_t step)
{
    char filename[128];
    
    sprintf(filename,"gmm.dat");
    FILE *fptr = fopen(filename,"a");
    fprintf(fptr,"%i %g %g\n",(int)step, logL, BIC);
    fclose(fptr);
    
    for(size_t m=0; m<NP; m++)
    {
        for(size_t n=m; n<NP; n++)
        {
            /* Get 1D PDFs for plotting */
            if(m==n)
            {
                sprintf(filename,"pdf_%i",(int)step);
                print_1D_pdfs(modes, samples, NMODE, NMCMC, filename, m);
            }
            
            /* Get Eigenvectors & Eigenvalues of Covariance matrix for plotting */
            else
            {
                sprintf(filename,"contours_%i",(int)step);
                print_2D_contours(modes, NMODE, filename, m, n);
                
                sprintf(filename,"2Dpdf_%i",(int)step);
                print_2D_pdfs(modes, samples, NMODE, NMCMC, filename, m, n);

            }
        }
    }
}


int expectation_maximization(struct Sample **samples, struct MVG **modes, size_t NP, size_t NMODE, size_t NMCMC, double *logL, double *BIC)
{
    // aliases to structures
    struct MVG *M = NULL;
    struct Sample *s = NULL;
    
    // helper quantities for building sums etc.
    double norm;
    double mu=0.0;
    double C=0.0;
    double p;
    
    /*
     E-step:
     compute probability for each sample to belong to each mode
     */
    
    /* compute p(x|mode) for each sample for each mode */
    for(size_t i=0; i<NMCMC; i++)
    {
        
        norm=0.0;
        s = samples[i];
        
        //loop over modes
        for(size_t k=0; k<NMODE; k++)
        {
            M = modes[k];
            
            //compute p(x|mode)
            p = multivariate_gaussian(s->x, M, NP);
            p = p > PMIN ? p : PMIN;

            s->p[k] = p;
            s->w[k] = M->p*p;
            norm += M->p*p;
        }
        for(size_t k=0; k<NMODE; k++) s->w[k]/=norm;
    }
    
    /* weigh the number of samples in each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        modes[k]->Neff = 0;
        for(size_t i=0; i<NMCMC; i++) modes[k]->Neff += samples[i]->w[k];
        if(modes[k]->Neff < 1.0 || modes[k]->Neff != modes[k]->Neff) return 1;
        modes[k]->p = modes[k]->Neff/(double)NMCMC;
    }
    
    /* check convergence with log likelihood & BIC */
    *logL = log_likelihood(modes, samples, (int)NMCMC, (int)NMODE);
    *BIC = -2.*(*logL) + (double)NMODE*((double)NP*((double)NP+3.)/2. + 1)*log((double)NMCMC);
    //printf(" logL = %g,  BIC = %g     ",*logL, *BIC);
    
    
    /*
     M-Step:
     recompute mean, covariance, and relative weight of newly weighted samples
     */
    
    /* compute weighted mean and covariance for each mode */
    for(size_t k=0; k<NMODE; k++)
    {
        M = modes[k];
        
        //get new means
        for(size_t n=0; n<NP; n++)
        {
            mu = 0.0;
            
            //mu is a weighted average for each mode
            for(size_t i=0; i<NMCMC; i++)
                mu += (samples[i]->w[k]/modes[k]->Neff) * samples[i]->x[n];
            
            M->mu[n] = mu;
        }
        
        //get new covariance
        for(size_t m=0; m<NP; m++)
        {
            //loop over parameters again (wasteful, not using symmetry
            for(size_t n=0; n<NP; n++)
            {
                C = 0.0;
                
                //loop over samples
                for(size_t i=0; i<NMCMC; i++)
                {
                    double dx_m = samples[i]->x[m] - M->mu[m];
                    double dx_n = samples[i]->x[n] - M->mu[n];
                    
                    C += (samples[i]->w[k]/modes[k]->Neff)*(dx_m)*(dx_n);
                }
                M->C[m][n] = C;
            }
        }
        
        //invert new matrix to evaluate the probabilities
        decompose_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, NP);
    }
    
    return 0;
}

int GMM_with_EM(struct MVG **modes, struct Sample **samples, size_t NP, size_t NMODE, size_t NMCMC, size_t NSTEP, unsigned int *r, double *logL, double *BIC)
{

    /* construct diagonal covariance matrix of full sample variances */
    double x_temp[NMCMC];
    double var_temp;
    for(size_t i=0; i<NP; i++)
    {
        for(size_t n=0; n<NMCMC; n++) x_temp[n] = samples[n]->x[i];
        var_temp  = get_variance(x_temp,NMCMC);
        
        //set diagonals of C
        for(size_t n=0; n<NMODE; n++) modes[n]->C[i][i] = var_temp;
        
    }
    
    /* place covariance matrices at random draws from the chain file */
    for(size_t k=0; k<NMODE; k++)
    {
        //pick a sample from the chain to anchor each covariance matrix
        int fair_draw = (int)(rand_r_U_0_1(r)*NMCMC);
        for(size_t n=0; n<NP; n++) modes[k]->mu[n] = samples[fair_draw]->x[n];
        
        //set priors for each model
        modes[k]->p = (double)1./(double)NMODE;
        
        //get inverse, determinant, etc.
        decompose_matrix(modes[k]->C, modes[k]->Cinv, modes[k]->L, &modes[k]->detC, NP);
    }
    
    /* EM Algorithm for Gaussian Mixture Models */
    size_t step=0;
    double BICmin = 1e60;
    while(step<NSTEP)
    {
        //printProgress((double)(step+1)/NSTEP);
        if(expectation_maximization(samples, modes, NP, NMODE, NMCMC, logL, BIC)) return 1;
        else
        {
            if(floor(*BIC) < floor(BICmin))
            {
                BICmin = *BIC;
                step=0;
            }
            step++;
        }
    }
    //printf("\n");
    return 0;
}

double logit(double x,double xmin,double xmax)
{
    return x;//log( (x-xmin)/(xmax-x) );
}

double sigmoid(double x,double xmin,double xmax)
{
    return x;//xmin + (1./(1. + exp(-x)))*(xmax - xmin);
}


double dsigmoid(double x, double xmin, double xmax)
{
    //double expyn = exp(-x);
    return 1.0;//(xmax-xmin) * expyn/((1.+expyn)*(1.+expyn));
}


void logit_mapping(double *x_vec, double *y_vec, double xmin, double xmax, int N)
{
    for(size_t n=0; n<N; n++)
    {
        double x = x_vec[n];
        if(x<xmin || x>xmax)
        {
            /*
             This is a safety check incase there are small differences
             between the prior setup by gb_catalog and the one used by
             the analysis.  This was happening for the fdot parameter
             which is a strong function of fmin, and inconsistencies
             in fmin between global_fit and gb_catalog runs.
             */
            if(x>xmax) x -= 2.*(x-xmax);
            if(x<xmin) x += 2.*(xmin-x);
        }
        double y = logit(x,xmin,xmax);//log( (x - xmin)/(xmax - x) );
        y_vec[n]=y;
    }
}

/* sigmoid function */
void sigmoid_mapping(double *x_vec, double *y_vec, double xmin, double xmax, int N)
{
    for(size_t n=0; n<N; n++)
    {
        double y = y_vec[n];
        double x = sigmoid(y,xmin,xmax);//xmin + (1./(1. + exp(-y)))*(xmax - xmin);
        x_vec[n] = x;
    }
}
