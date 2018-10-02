
/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include "param.h"
#include "definitions.h"
#include "functions.h"

extern myfloat prob[nT][nB], probg[nT][nB][nQ];

#if grey==0
void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP, myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax)
{

    myfloat dummy,dummy1;
    //read planck-absorption coefficient
    FILE *fp = fopen("../distr/tables/planck-mean.txt","r");
    for(int i=0; i<nT; i++)
    {
 	   fscanf(fp,readf" "reade" "reade" "reade"\n",&Tnb[i],&kP[i],&dummy,&dummy1);
    }
    fclose(fp);


    // absorption coefficient at maximum temperature in the system
    int t = (int) ((Tmax - Tnb[0])/(Tnb[1]-Tnb[0]));
    *kappamax = (kP[t+1] - kP[t])/(Tnb[t+1]-Tnb[t])*
   	   	   	   (Tmax - Tnb[t])+kP[t];

    //read tables of narrow-bands and probability
    char *file[nB];
    for (int i=0; i<nB; i++)
    {
       file[i] = (char*)malloc(68*sizeof(char));
 	   sprintf(file[i],"../distr/tables/NarrBand%d.txt",i);
 	   fp = fopen(file[i],"r");
 	   fscanf(fp,reade"\t"reade"\t"reade"\n",&dummy,&dummy,&narrBand[i].wvc);
 	   for (int j=0; j<nT; j++)
 	   {
 		   fscanf(fp," "readf"\t",&dummy);
 		   fscanf(fp," "readf"\t",&dummy1);
		   fscanf(fp," "readf"\t",&narrBand[i].kavg);
 		   for(int g=0; g<nQ; g++)
 		   {
 			   fscanf(fp," "reade"\t",&narrBand[i].kq[j][g]);
 		   }
 		   fscanf(fp,"\n");
 	   }
 	   fclose(fp);
 	   narrBand[i].idx = i;
    }

    // read probability and define the host dynamic


    int dummy3;
    fp = fopen("../distr/tables/prob_new2.txt","r");
    int n=0;
    int m=0;
    for (int i=0; i<nT; i++)
    {
 	   for (int j=0; j<nB; j++)
 	   {
 		   fscanf(fp,readf" \t %d \t "reade" \t",&dummy,&dummy3,&prob[i][j]);
 		   prob_h[n] = prob[i][j];
 		   n += 1;

 		   for (int g=0; g<nQ; g++)
 		   {
 			   fscanf(fp,reade"\t",&probg[i][j][g]);
 	 		   probg_h[m] = probg[i][j][g];
 	 		   m += 1;

 		   }
 		   fscanf(fp,"\n");
 	   }
    }
    fclose(fp);

    for (int i=0; i<nB; i++)
    {
	free(file[i]);
    }
}

#elif grey==1
#include <math.h>
myfloat I_blackC2( myfloat T, myfloat nu)
{
    // way of calculating C1 and C2 (nu is in cm^(-1) while the output is in W/m^2
    myfloat C1 = 3.741771790075259e-16;
    myfloat C2 = 0.014387741858429;
    return 1.0 / pi * C1 * pow3(nu*100) / (expf(C2*nu*100/T)-1);
}

void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP, myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax)
{

        //Planck absorption coefficient is the grey absorption coefficient
	for(int i=0; i<nT; i++)
	{
		kP[i] = abscoeff;
		Tnb[i] = Tmins*0.9 + (Tmax*1.1 - Tmins*0.9) * ((myfloat)i)/((myfloat)nT);
	}

	// absorption coefficient at maximum temperature in the system
	*kappamax = abscoeff;

	//read tables of narrow-bands and probability
	myfloat wvl[nB],wvr[nB];
	wvl[0] = 0;
	for (int i=1; i<nB; i++)
	{
		wvl[i] = wvl[i-1]+wvmax/nB;
	}
	wvr[0] = wvmax/nB;
	for (int i=1; i<nB; i++)
	{
		wvr[i] = wvr[i-1]+wvmax/nB;
	}
	for (int i=0; i<nB; i++)
	{
		narrBand[i].wvc = (wvl[i]+wvr[i])/2.;
		for (int j=0; j<nT; j++)
		{
			for(int g=0; g<nQ; g++)
			{
				narrBand[i].kq[j][g] = abscoeff;
			}
		}
		narrBand[i].idx = i;
	}

	// read probability and define the host dynamic

	myfloat Ib[nT][nB];

	for(int t = 0; t<nT; t++ )
	{
		for(int nb = 0; nb<nB; nb++ )
		{
#if srt ==1
			Ib[t][nb] = I_blackC(Tnb[t],narrBand[narrBand[nb].idx].wvc);
#else
			Ib[t][nb] = I_blackC2(Tnb[t],narrBand[nb].wvc);
#endif
		}
	}
	int n=0;
	int m=0;
	for (int i=0; i<nT; i++)
	{
		prob[i][0] =  Ib[i][0]*(wvr[0]-wvl[0])*100;
		for (int j=1; j<nB; j++)
		{
			prob[i][j] = prob[i][j-1] + Ib[i][j]*(wvr[j]-wvl[j])*100;
		}
	}
	for (int i=0; i<nT; i++)
	{
		for (int j=0; j<nB; j++)
		{
			prob[i][j] = prob[i][j]/prob[i][nB-1];
			prob_h[n] = prob[i][j];
			n += 1;

			for (int g=0; g<nQ; g++)
			{
				probg[nT][nB][nQ] = g;
				probg_h[m] = probg[nT][nB][nQ];
				m += 1;
			}
		}
	}
}

#elif grey==2 

void readT(NarrowBand *narrBand, myfloat *Tnb, myfloat *kP, myfloat *prob_h, myfloat *probg_h, myfloat Tmax, myfloat *kappamax)
{

    //read planck-absorption coefficient
    FILE *fp = fopen("tables/particles/planck-mean.txt","r");
    for(int i=0; i<nT; i++)
    {
 	   fscanf(fp,readf" "reade"\n",&Tnb[i],&kP[i]);
    }
    fclose(fp);

    // absorption coefficient at maximum temperature in the system
    int t = (int) ((Tmax - Tnb[0])/(Tnb[1]-Tnb[0]));
    *kappamax = (kP[t+1] - kP[t])/(Tnb[t+1]-Tnb[t])*
   	   	   	   (Tmax - Tnb[t])+kP[t];

    //read tables of narrow-bands and probability
    char *file[nB];
    for (int i=0; i<nB; i++)
    {
      	   file[i] = (char*)malloc(68*sizeof(char));
 	   sprintf(file[i],"../distr/tables/particles/NarrBand%d.txt",i);
 	   fp = fopen(file[i],"r");
 	   fscanf(fp,reade"\t"reade"\t"reade"\n",&dummy,&dummy,&narrBand[i].wvc,narrband[i].kq[0][0]);
	   for (int j=0; j<nT; j++)
                 for(int g=0; g<nQ, g++)           
                     narrBand[i].kq[j][g] = narrBand[i].kq[0][0];
																			                           }
 	   fclose(fp);
 	   narrBand[i].idx = i;
    }

    int dummy3;
    fp = fopen("../distr/tables/particles/prob.txt","r");
    int n=0;
    int m=0;
    for (int i=0; i<nT; i++)
    {
 	   for (int j=0; j<nB; j++)
 	   {
 		   fscanf(fp,readf" \t %d \t "reade" \t",&dummy,&dummy3,&prob[i][j]);
 		   prob_h[n] = prob[i][j];
 		   n += 1;

 		   for (int g=0; g<nQ; g++)
 		   {
 			   probg_h[m] = 1;
 	 		   m += 1;
 		   }
 		   fscanf(fp,"\n");
 	   }
    }
    fclose(fp);
}


#endif
