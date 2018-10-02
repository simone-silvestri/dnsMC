
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
