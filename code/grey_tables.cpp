
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
#include <math.h>
#include "param.h"
#include "definitions.h"
#include "functions.h"

extern myfloat prob[nT][nB], probg[nT][nB][nQ];

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

