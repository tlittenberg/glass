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


#include <glass_utils.h>

#define RAD2ARCMIN 3437.75
#define RAD2DEGREE 57.295833313961

/* ============================  MAIN PROGRAM  ============================ */

double M_fdot(double f, double fdot)
{
  return pow( fdot*pow(f,-11./3.)*(5./96.)*pow(M_PI,-8./3.)  ,  3./5.)/TSUN;
}

double M_fddot(double f, double fddot)
{
  return pow( fddot*pow(f,-19./3.)*(25./33792.)*pow(M_PI,-16./3.)  ,  3./10.)/TSUN;
}

double M_fdot_fddot(double f, double fdot, double fddot)
{
  //Eq 174 in SVN:LISA-WDWD-Test/Mathematica/with-Tyson.nb.pdf
  return pow(5.,3./5.)*pow((-(fdot*sqrt(fddot))/(pow(f,19./6.)*(30.*sqrt(f*fddot) - 11.*sqrt(33.)*fdot))),3./5.)/8./pow(M_PI,8./5.)/TSUN;
}

double beta(double f, double fdot, double fddot)
{
  //Eq 175 in SVN:LISA-WDWD-Test/Mathematica/with-Tyson.nb.pdf
  return (11.*(-3*sqrt(f*fddot)+sqrt(33.)*fdot)*pow((sqrt(fddot)*fdot)/(-30.*sqrt(f*fddot) + 11.*sqrt(33.)*fdot),2./5.)*pow(5./M_PI,2./5.))/(12.*pow(f,11./10.)*sqrt(fddot));
}

double distance(double f0, double dfdt, double A)
{
  double f    = f0;//T;
  double fd = dfdt;//(T*T);
  double amp   = A;
  return ((5./48.)*(fd/(M_PI*M_PI*f*f*f*amp))*CLIGHT/PC); //seconds  !check notes on 02/28!
}

int main(int argc, char* argv[])
{

  if(argc<3)
  {
    fprintf(stdout,"Usage: ucb_chirpmass parameter_chain.dat.0 9 [alpha] [dalpha]\n");
    return 1;
  }
  
  FILE *ifile = fopen(argv[1],"r");
  
  int NP = atoi(argv[2]);
  
  int tideFlag   = 0;
  double alpha0  = 0.0;
  double dalpha0 = 0.0;
  double alpha   = 0.0;
  if(argc>3)
  {
    alpha0 = atof(argv[3]);
    dalpha0=alpha0*.1; //10% error
    tideFlag = 1;
  }
  if(argc>4)
  {
    dalpha0 = atof(argv[4]);
  }
  
  //set up RNG in case tidal params are used
  unsigned int r = 150914;
  
  
  char filename[128];
  sprintf(filename,"mchirp_%s",argv[1]);
  FILE *ofile = fopen(filename,"w");
  
  sprintf(filename,"frequencies_%s",argv[1]);
  FILE *ofile2 = fopen(filename,"w");

  sprintf(filename,"localization_%s",argv[1]);
  FILE *ofile3 = fopen(filename,"w");

  //parse file
  
  /*
   0)  index
   1)  logL
   2)  t0
   3)  tgap
   4)  f
   5)  fdot
   6)  A
   7)  phi
   8)  theta
   9)  cosi
   10) psi
   11) phase
   12) fddot
   */
  
  double A=0,phi=0.0,theta=0,cosi=0,psi=0,phase=0;
    double f = 0.0;
    double fdot = 0.0;
    double fddot = 0.0;
  double Mc,B;
    int check;
  while(!feof(ifile))
  {
    if(NP==9)
    {
        check = fscanf(ifile,"%lg%lg%lg%lg%lg%lg%lg%lg%lg",&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase,&fddot);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",argv[1]);
            exit(1);
        }
    }
    if(NP==8)
    {
      check = fscanf(ifile,"%lg%lg%lg%lg%lg%lg%lg%lg",&f,&fdot,&A,&phi,&theta,&cosi,&psi,&phase);
        if(!check)
        {
            fprintf(stderr,"Error reading %s\n",argv[1]);
            exit(1);
        }

      fddot = 11.0/3.0*fdot*fdot/f;
    }
    
    //get tides
    if(tideFlag)
    {
      do
      {
        alpha = alpha0+dalpha0*rand_r_N_0_1(&r);
      }while(alpha < 0.0);
    }
    else alpha = 0.0;
    
    
    if(fdot>0.0&&fddot>0.0)
    {
      Mc = M_fdot_fddot(f,fdot,fddot);
      B = beta(f,fdot,fddot);
      fprintf(ofile3,"%.12g %.12g %.12g %.12g %.12g\n",phi*RAD2ARCMIN,theta*RAD2ARCMIN,distance(f,fdot,A)/(1.+alpha),acos(cosi)*RAD2DEGREE,M_fdot(f,fdot));
      fprintf(ofile2,"%.12g %.12g %.12g\n",f,fdot,fddot);
      if(Mc==Mc && B==B) fprintf(ofile,"%.12g %.12g %.12g %.12g\n",M_fdot(f,fdot),M_fddot(f,fddot),M_fdot_fddot(f,fdot,fddot),beta(f,fdot,fddot));
    }
  }
  
  return 0;
}

