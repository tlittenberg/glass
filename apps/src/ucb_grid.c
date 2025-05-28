/*
 * Copyright 2024 Tyson B. Littenberg
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

#define UCB_PER_SEGMENT 15
#define MAX_UCB_FREQUENCY 0.012
#define SEGMENT_BANDWIDTH 1E-5
#define MAX_SEGMENT_BANDWIDTH 4*2*SEGMENT_BANDWIDTH

int main(int argc, char* argv[])
{
    FILE *cachefile = fopen("ucb_catalog.cache","r");
    char name[256];
    char path[256];
    double snr;
    double confidence;
    
    int count=0;
    double x;
    while(!feof(cachefile))
    {
        fscanf(cachefile,"%s %lg %lg %lg %s",name,&x,&snr,&confidence,path);
        count++;
    }
    count--;
    rewind(cachefile);
    
    printf("found %i sources in cache file\n",count);
    
    double *f = malloc(count*sizeof(double));
    for(int n=0; n<count; n++)
    {
        fscanf(cachefile,"%s %lg %lg %lg %s",name,&f[n],&snr,&confidence,path);
    }
    fclose(cachefile);
    
    double_sort(f,count);
    
    
    FILE *outfile = fopen("ucb_grid.dat","w");
    
    int n=0;
    int seg=0;
    
    if(f[0]<MAX_UCB_FREQUENCY)
    {
        double fstart = f[0]-(f[1]-f[0]);
        double fstop;
        
        if(n+UCB_PER_SEGMENT+1 < count)
             fstop = (f[n+UCB_PER_SEGMENT]+f[n+UCB_PER_SEGMENT+1])/2;
        else
            fstop = f[count-1];


	//check that the step isn't too large
        double segment_sub_bandwidth = fstop - fstart;
        if(segment_sub_bandwidth > MAX_SEGMENT_BANDWIDTH)
        {
            double overfactor = round((fstop-fstart)/(MAX_SEGMENT_BANDWIDTH));
            fprintf(stdout,"WARNING: Wide segment at %lg Hz by a factor of %lg (%i bins)\n",f[n], overfactor, (int)(segment_sub_bandwidth*YEAR));
            segment_sub_bandwidth = segment_sub_bandwidth / overfactor;
            for(int j=0; j<overfactor; j++)
            {
                fstop=fstart+segment_sub_bandwidth;
                fprintf(outfile,"%i %i %lg %lg\n",seg, (int)(segment_sub_bandwidth*YEAR), fstart, fstop);
                fstart=fstop;
                seg++;
            }

        }
        else
        {
            fprintf(outfile,"%i %i %lg %lg\n",seg, (int)(segment_sub_bandwidth*YEAR), fstart, fstop);
            seg++;
        }
        //n+=UCB_PER_SEGMENT;

        //fprintf(outfile,"%i %i %lg %lg\n",0, (int)((fstop-fstart)*YEAR), fstart,fstop);
        //seg++;
        n=UCB_PER_SEGMENT;
    }
    
    while(f[n]<MAX_UCB_FREQUENCY)
    {
        double fstart = (f[n]+f[n+1])/2.;
        double fstop;
        
        if(n+UCB_PER_SEGMENT+1 < count)
             fstop = (f[n+UCB_PER_SEGMENT]+f[n+UCB_PER_SEGMENT+1])/2;
        else
            fstop = f[count-1];

        //check that the step isn't too large
        double segment_sub_bandwidth = fstop - fstart;
        if(segment_sub_bandwidth > MAX_SEGMENT_BANDWIDTH)
        {
            double overfactor = round((fstop-fstart)/(MAX_SEGMENT_BANDWIDTH));
            fprintf(stdout,"WARNING: Wide segment at %lg Hz by a factor of %lg (%i bins)\n",f[n], overfactor, (int)(segment_sub_bandwidth*YEAR));
            segment_sub_bandwidth = segment_sub_bandwidth / overfactor;
            for(int j=0; j<overfactor; j++)
            {
                fstop=fstart+segment_sub_bandwidth;
                fprintf(outfile,"%i %i %lg %lg\n",seg, (int)(segment_sub_bandwidth*YEAR), fstart, fstop);
                fstart=fstop;
                seg++;
            }
            
        }
        else
        {
            fprintf(outfile,"%i %i %lg %lg\n",seg, (int)(segment_sub_bandwidth*YEAR), fstart, fstop);
            seg++;
        }
        n+=UCB_PER_SEGMENT;
        if(n>count-UCB_PER_SEGMENT) 
        {
            double fstart = (f[n]+f[n+1])/2.;
            double fstop  = f[count-1]+SEGMENT_BANDWIDTH/2.;
            fprintf(outfile,"%i %i %lg %lg\n",seg, (int)(segment_sub_bandwidth*YEAR), fstart, fstop);
           return 0;
        }
    }
    
    n++;
    while(n<count)
    {
        double fstart = f[n]-SEGMENT_BANDWIDTH;
        double fstop  = f[n]+SEGMENT_BANDWIDTH;
        
        
        while(n<count-1 && fstop>f[n+1])
        {
            fstop  = f[n+1]+SEGMENT_BANDWIDTH;
            n++;
        }

        fprintf(outfile,"%i %i %lg %lg\n",seg, (int)((fstop-fstart)*YEAR), fstart, fstop);
        n++;
        seg++;
    }
    
    
    
    
    return 0;
}

