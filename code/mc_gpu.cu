#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <cuda.h>
#include <curand_kernel.h>
#include "param.h"
#include "definitions.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include "functions.h"
#include "device_functions.h"
#include "memory.h"

using namespace std;
const int Ntot = imax*jmax*kmax;
const int N = Ntot/num_streams/num_gpu/p_row;
myfloat kP[nT] ,prob[nT][nB] ,probg[nT][nB][nQ];
myfloat kPp[nT],probp[nT][nB];

// global declaration of 2D float texture (visible for host and device code)
texture<myfloat, cudaTextureType2D, cudaReadModeElementType> tex_Ib;
__device__ __constant__ myfloat probd[2][nB];
__device__ myfloat solution[num_streams][N];
__device__ myfloat variance[num_streams][N];
#if srt==1
void sort_idx(NarrowBand *narrBand);
#endif

// pointer definitions
Gridn *gridGPU[num_gpu];
EmissSpec *Ibw_d[num_gpu];
cudaArray *cuArray[num_gpu];
cudaTextureObject_t *tex_tempf_d[num_gpu];
int *idx_d[num_gpu];
myfloat *wvc_d[num_gpu];
myfloat *Tnb_d[num_gpu];
cudaTextureObject_t *tex_d[num_gpu];
cudaTextureObject_t *tex_prob_d[num_gpu];


__global__ void kernel_fluid(Gridn *my_grid, myfloat *wvc_d, myfloat *Tnb_t, int n, int ns, int stream, myfloat kappamax, myfloat Tmax, int *idx_nb,
		cudaTextureObject_t *tex_Tf, cudaTextureObject_t *tex, cudaTextureObject_t *tex_prob,
		int gpu, int ystart, EmissSpec *Ibw)
{

	Count cnt;

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	volatile int tx = threadIdx.x;

	// for the moment ray stays in the register memory, if it spills,
	// we have to move it to the shared memory (it should fit in register memory)
	Beam ray;
	// we start defining the ray-Ibmax that is equal for every position

	__const__ myfloat Tnb[] = {Tnb_t[1] - Tnb_t[0], (Tnb_t[0])/(Tnb_t[1] - Tnb_t[0])};
	myfloat abscotemp1,abscotemp2;
	myfloat ratioI,invTi4;
	myfloat pos;
	myfloat Ibmax = 4.0f * kappamax * pow4(Tmax) * stefan / photmax;


	//define a narrow band
	myfloat wvc;
	__shared__ curandState_t state[blockdim];
	// this is basically the loop over the indices in the CPU.
	// the indices are found in relation with the stream-block-Thread design.
	// this is a grid-stride loop, hence every thread will do multiple cells
	for (int idx = tid; idx < n; idx += blockDim.x * gridDim.x)
	{
		solution[stream][idx] = 0;
		ray.i = (idx / (kmax*jmax/num_gpu/p_row) + 1) + (stream)*imax/ns;
		ray.j = (idx / (kmax/num_gpu) + 1 - (ray.i-1-(stream)*imax/ns)*jmax/p_row);
		ray.k = (idx - (kmax/num_gpu) * (ray.j-1 + (ray.i-1-(stream)*imax/ns)*jmax/p_row) + 1) + gpu*kmax/num_gpu;
		ray.j+= ystart-1;

		ray.Ti = tex3D<myfloat>(tex_Tf[0], ray.i+0.5f, ray.j+0.5f, ray.k+0.5f);

		myfloat De_var[nVar];

		for (int v = 0; v<nVar; v++)
		{
			De_var[v] = 0;
			int g=0;
			int nb=0;

			curand_init(1234, (idx+1)*(stream+1)*(v+1)*(1+gpu)*ystart, 0, &state[tx]);  // 	Initialize CURAND
#if srt == 1
			kernel_find(&cnt, Tnb_t, Tmax, &state[tx], tex_prob);
#endif
			int countnb=0;
			int countg=0;

			myfloat De_OERMc=0;
			/***********  ENTERING THE PHOTON LOOP, DEVIDED INTO SUB-BUNDLES *****************/

			for (int h=0; h<photmax; h++)
			{
				ray.xp = my_grid[0].x[ray.i];
				ray.yp = my_grid[0].y[ray.j];
				ray.zp = my_grid[0].z[ray.k];
				ray.ic = ray.i;
				ray.jc = ray.j;
				ray.kc = ray.k;
				int flag[3];
				// Now we have to deal with random numbers, selecting the angles
				// define the scattering function
				emiss_ang(&ray, &state[tx]);

				if(ray.sx == 0)
					ray.sx = 1e-10;
				if(ray.sy == 0)
					ray.sy = 1e-10;
				if(ray.sz == 0)
					ray.sz = 1e-10;
				flag[0] = (int)(ray.sx<0);
				flag[1] = (int)(ray.sy<0);
				flag[2] = (int)(ray.sz<0);

				//Find the narrowband index and the g index, based on random or based on the previous count
#if srt == 1
				wave_find(&nb,&g,&countnb,&countg,idx_nb,&cnt);
#else
				// non sorted monte carlo
				int tm =  (int) ((Tmax - Tnb_t[0])/(Tnb_t[1]-Tnb_t[0]));
				find_band(Tnb_t,Tmax,&g,&nb,&state[tx],tm, tex_prob);
#endif
				// the ray emanating from cell i,j,k is completely defined now
				// now we have to march it
				// initializing transmissivity

				ray.tra = 1.0f;
				wvc  =  wvc_d[nb];

				Emission emiss;

				// wall emissions
				emiss.west  =  Ibw->west[nb];
				emiss.east  =  Ibw->east[nb];
				//    			emiss.north =  Ibw->north[nb];
				//    			emiss.south =  Ibw->south[nb];
				//    			emiss.top   =  Ibw->top[nb];
				//    			emiss.bot   =  Ibw->bot[nb];

				// Calculate parameters for calculation of De_OERMc in march_ray
				pos = ((ray.Ti - Tnb_t[0])/(Tnb_t[1]-Tnb_t[0]));
				abscotemp2 = tex3D<myfloat>(*tex, nb+0.5f, pos+0.5f, g+0.5f);
				invTi4 = 1.0f / tex2D(tex_Ib,nb+0.5f ,pos+0.5f);

				pos = ((Tmax - Tnb_t[0])/(Tnb_t[1]-Tnb_t[0]));
				abscotemp1 = tex3D<myfloat>(*tex, nb+0.5f, pos+0.5f, g+0.5f);
				ratioI = (1.0f / invTi4) / tex2D(tex_Ib,nb+0.5f ,pos+0.5f) * abscotemp2 / abscotemp1 ;

				// Adjust the cells according to the scattering of the ray
				// if the ray originates from the negative cell boundary and it
				// is going in the negative direction, shift cells
				if ( flag[0] && ray.xp == my_grid[0].xu[ray.ic - 1] )
					ray.ic = ray.ic - 1;
				if ( flag[1] && ray.yp == my_grid[0].yu[ray.jc - 1] )
					ray.jc = ray.jc - 1;
				if ( flag[2] && ray.zp == my_grid[0].zu[ray.kc - 1] )
					ray.kc = ray.kc - 1;
				// Do the same thing if it is happening in the positive direction
				if ( !flag[0] && ray.xp == my_grid[0].xu[ray.ic] )
					ray.ic = ray.ic + 1;
				if ( !flag[1] && ray.yp == my_grid[0].yu[ray.jc] )
					ray.jc = ray.jc + 1;
				if ( !flag[2] && ray.zp == my_grid[0].zu[ray.kc] )
					ray.kc = ray.kc + 1;

				/**************************************************************************************/
				/************************** ENTERING THE MARCHING LOOP ********************************/
				/**************************************************************************************/

				//loop on all grids with maximum counter my_grid[grd].sm
				for(int grd = 0; grd < grid_num; grd++)
				{
					march_ray(&ray, nb, g, flag, &De_OERMc, Ibmax, my_grid[grd], ratioI, invTi4,
							tex_Tf[grd], emiss, wvc, Tnb, tex);

					if(grd<grid_num-1) {
						//adapting grid
						ray.ic += my_grid[grd].im/my_grid[grd+1].im-1;
						ray.jc += my_grid[grd].jm/my_grid[grd+1].jm-1;
						ray.kc += my_grid[grd].km/my_grid[grd+1].km-1;
						ray.ic /= my_grid[grd].im/my_grid[grd+1].im;
						ray.jc /= my_grid[grd].jm/my_grid[grd+1].jm;
						ray.kc /= my_grid[grd].km/my_grid[grd+1].km;
					}
				};

			}

			/***********  OUT OF THE PHOTON LOOP ****************************/
			solution[stream][idx] += De_OERMc;
			De_var[v] = De_OERMc;
		}
		/***********  OUT OF THE VARIANCE LOOP ****************************/

		solution[stream][idx] /= nVar;
		variance[stream][idx] = 0;
		for (int v = 0; v<nVar; v++)
		{
			variance[stream][idx] += 1./(nVar-0.5) * 1./(nVar-1) * (De_var[v] - solution[stream][idx]) * (De_var[v] - solution[stream][idx]);
		}
	}
	/***********  OUT OF THE GRID-STRIDE LOOP ****************************/

}

__device__ __forceinline__ void march_ray(Beam *ray, int nb, int g, int *flag, myfloat *De_OERMc,
		myfloat Ibmax, Gridn grid, myfloat ratioI, myfloat invTi4, cudaTextureObject_t tex_Tf,
		Emission emiss, myfloat wvc, const myfloat Tnb[], cudaTextureObject_t *tex)
{

	int counter = 0;
	while ( (ray->tra > toll) && (counter<grid.sm) )
	{

		// find the distance to cell boundary x,y,z -> the minimal will be the crossing
		myfloat dsx, dsy, dsz;
		dsx = ( grid.xu[ray->ic-flag[0]] - ray->xp ) / ray->sx;
		dsy = ( grid.yu[ray->jc-flag[1]] - ray->yp ) / ray->sy;
		dsz = ( grid.zu[ray->kc-flag[2]] - ray->zp ) / ray->sz;

		// interpolate to find the temperature of the cell (particle and fluid)
		myfloat tf = tex3D<myfloat>(tex_Tf,ray->ic+0.5f,ray->jc+0.5f,ray->kc+0.5f);

		// interpolate to find the absorption of the cell (particle and fluid)
		myfloat pos = (tf/Tnb[0] - Tnb[1]);
		myfloat absco =  tex3D<myfloat>(*tex, nb+0.5f, pos+0.5f, g+0.5f);

		// black-body radiation of the cell (particle and fluid)
		myfloat blackpow = tex2D(tex_Ib,nb+0.5f ,pos+0.5f);
		myfloat ds = dsx;
		(void) ( (ds > dsy) && (ds = dsy) );
		(void) ( (ds > dsz) && (ds = dsz) );

		// update ray position and scattering length
		(void) ( (ds==dsx) && (ray->xp = grid.xu[ray->ic-flag[0]]) );
		(void) ( (ds!=dsx) && (ray->xp =   ray->xp + ds * ray->sx) );
		(void) ( (ds==dsy) && (ray->yp = grid.yu[ray->jc-flag[1]]) );
		(void) ( (ds!=dsy) && (ray->yp =   ray->yp + ds * ray->sy) );
		(void) ( (ds==dsz) && (ray->zp = grid.zu[ray->kc-flag[2]]) );
		(void) ( (ds!=dsz) && (ray->zp =   ray->zp + ds * ray->sz) );

		myfloat alpha  = 1.0f - __expf(-ds * (absco));

		*De_OERMc -= Ibmax * ray->tra * alpha * ratioI *
				( invTi4 * blackpow - 1.0f );

		// update transmissivity of the ray and total distance travelled
		ray->tra = ray->tra * (1-alpha);
		if ( ray->tra < toll)
		{
			*De_OERMc -= Ibmax * ray->tra * ratioI *
					( invTi4 * blackpow - 1.0f );
			ray->tra = 0;
		}

		// Updating cell indices and Boundary conditions if end is reached
		// efficient ray tracing method?

		if ( dsx<dsy )
		{
			if ( dsx<dsz )
			{
				if ( flag[0] )
				{
					ray->ic = ray->ic-1;
					if ( ray->ic == 0 )
					{
						if(bdw==1)
						{
							ray->ic = grid.im;
							ray->xp = Lx;
						}
						else if(bdw==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.west - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
				else
				{
					ray->ic = ray->ic+1;
					if ( ray->ic == grid.im+1 )
					{
						if(bde==1)
						{
							ray->ic = 1;
							ray->xp = 0;
						}
						else if(bde==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.east - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
			}
			else
			{
				if ( flag[2] )
				{
					ray->kc = ray->kc-1;
					if ( ray->kc == 0 )
					{
						if(bdb==1)
						{
							ray->kc = grid.km;
							ray->zp = Lz;
						}
						else if(bdb==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.bot - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
				else
				{
					ray->kc = ray->kc+1;
					if ( ray->kc == grid.km+1 )
					{
						if(bdt==1)
						{
							ray->kc = 1;
							ray->zp = 0;
						}
						else if(bdt==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.top - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
			}
		}
		else
		{
			if ( dsy<dsz )
			{
				if ( flag[1] )
				{
					ray->jc = ray->jc-1;
					if ( ray->jc == 0 )
					{
						if(bds==1)
						{
							ray->jc = grid.jm;
							ray->yp = Ly;
						}
						else if(bds==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.south - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
				else
				{
					ray->jc = ray->jc+1;
					if ( ray->jc == grid.jm+1 )
					{
						if(bdn==1)
						{
							ray->jc = 1;
							ray->yp = 0;
						}
						else if(bdn==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.north - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
			}
			else
			{
				if ( flag[2] )
				{
					ray->kc = ray->kc-1;
					if ( ray->kc == 0 )
					{
						if(bdb==1)
						{
							ray->kc = grid.km;
							ray->zp = Lz;
						}
						else if(bdb==2)
						{
							*De_OERMc -= Ibmax * ray->tra*( invTi4 * emiss.bot - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
				else
				{
					ray->kc = ray->kc+1;
					if ( ray->kc == grid.km+1 )
					{
						if(bdt==1)
						{
							ray->kc = 1;
							ray->zp = 0;
						}
						else if(bdt==2)
						{
							*De_OERMc -= Ibmax * ray->tra *( invTi4 * emiss.top - 1.0f )* ratioI;
							ray->tra = 0;
						}
					}
				}
			}
		}
		counter+=1;

	}  // closing the while statement
}

extern "C" void mc_gpu_(myfloatF *Tfort, int *ystart)
{

	myfloat kappamax,Tmax;

	/**************************************************************/
	/********** CREATING GRID AND TEMPERATURE FIELD ***************/
	/**************************************************************/

	Var_CPU varCPU[grid_num];
	Gridn gridCPU[grid_num];
	for (int grd = 0; grd < grid_num; grd++)
	{
		varCPU[grd].mk_grid(maxi[grd],maxj[grd],maxk[grd],maxi[0]);
		gridCPU[grd].mk_grid(maxi[grd],maxj[grd],maxk[grd],maxi[0],maxs[grd]);
	}
	for (int k = 0; k < (kmax+2); k++)
	{
		for (int j = 0; j < (jmax+2); j++)
		{
			for (int i = 0; i < (imax+2); i++)
			{
				varCPU[0].T[idx_T(i,j,k,imax,jmax)] = (myfloat) Tfort[idx_T(i,j,k,imax,jmax)];
			}
		}
	}

	// interpolating temperature and finding new concentration on coarser grid
	interp3D(varCPU);

	Tmax = 0;
	for (int k=1; k<kmax+1; k++)
	{
		for (int j=1; j<jmax+1; j++)
		{
			for (int i=1; i<imax+1; i++)
			{
				Tmax = MAX(Tmax,varCPU[0].T[idx_T(i,j,k,maxi[0],maxj[0])]);
			}
		}
	}
	Tmax = MAX(Tmax,Tww);
	Tmax = MAX(Tmax,Twe);
	Tmax = MAX(Tmax,Twn);
	Tmax = MAX(Tmax,Tws);
	Tmax = MAX(Tmax,Twt);
	Tmax = MAX(Tmax,Twb);

	/**************************************************************/
	/***********  READ THE TABLES *********************************/
	/**************************************************************/

	NarrowBand *narrBand;
	myfloat *Tnb;
	myfloat *prob_h, *probg_h, *prob_h2;

	prob_h   = (myfloat*)malloc(nB*nT*   sizeof(myfloat));
	probg_h  = (myfloat*)malloc(nB*nT*nQ*sizeof(myfloat));
	prob_h2  = (myfloat*)malloc(nB*2*sizeof(myfloat));

	narrBand = (NarrowBand*)malloc(nB *  sizeof(NarrowBand));
	Tnb = (myfloat*)malloc(nT*sizeof(myfloat));

	readT(narrBand, Tnb, kP, prob_h, probg_h, Tmax, &kappamax);

	/**************************************************************/
	/***********  FINISHED READING ********************************/
	/**************************************************************/

	//Sorting NarrowBands based on kavg of the band
#if srt == 1
	sort_idx(narrBand);
#endif

	/**************************************************************/
	/************* MEMORY COPY TO THE GPU's ***********************/
	/**************************************************************/

	cudaStream_t streams0[num_streams];
	cudaStream_t streams1[num_streams];

	for(int gpu = 0; gpu < num_gpu; gpu++) {
		cudaSetDevice(gpu);

		for (int i=0; i<num_streams; i++)
		{
			if(gpu==0) {
				cudaStreamCreate(&streams0[i]);
			}
			else {
				cudaStreamCreate(&streams1[i]);
			}
		}

		// grid copy to GPU
		cudaMalloc((void**)&gridGPU[gpu], grid_num * sizeof(Gridn));
		grid_copy(gridCPU, gridGPU[gpu]);
		// Wall emission copy to GPU
		cudaMalloc((void**)&Ibw_d[gpu], nB * sizeof(EmissSpec));
		black_copy(Ibw_d[gpu],narrBand);

		myfloat Ib[nT][nB];

		for(int t = 0; t<nT; t++ )
			for(int nb = 0; nb<nB; nb++ )
			{
#if srt ==1
		Ib[t][nb] = I_blackC(Tnb[t],narrBand[narrBand[nb].idx].wvc);
#else
		Ib[t][nb] = I_blackC(Tnb[t],narrBand[nb].wvc);
#endif
			}

		// Create explicit channel description (could use an implicit as well)
		cudaChannelFormatDesc DescIb = cudaCreateChannelDesc<myfloat>();
		cudaMallocArray(&cuArray[gpu], &DescIb, nB, nT);
		cudaMemcpyToArray(cuArray[gpu], 0, 0, Ib, nB*nT*sizeof(myfloat), cudaMemcpyHostToDevice);
		tex_Ib.addressMode[0] = cudaAddressModeClamp;
		tex_Ib.addressMode[1] = cudaAddressModeClamp;
		tex_Ib.filterMode = cudaFilterModeLinear;
		tex_Ib.normalized = false;
		cudaBindTextureToArray(tex_Ib, cuArray[gpu], DescIb);

		//textured memory copy of interpolated temperature
		cudaMalloc((void**)&tex_tempf_d[gpu], grid_num*sizeof(cudaTextureObject_t) );
		temp_fluid_copy(tex_tempf_d[gpu], varCPU);

		//memory copy of -> sorted index 			idx_d
		//				 -> central wavenumber 		wvc_d
		//				 -> discrete temperature	Tnb_d
		//				 -> textured asb coeff		tex_d
		//				 -> textured emiss prob		tex_prob_d
		//				 -> phase function prob		prob_A_d
		//				 -> textured part prob		tex_part_d
		//				 -> particle scatt coeff	Csca_d
		//				 -> particle abs coeff		Cabs_d
		cudaMalloc((void**)&idx_d[gpu],nB * sizeof(int));
		cudaMalloc((void**)&wvc_d[gpu],nB * sizeof(myfloat));
		cudaMalloc((void**)&Tnb_d[gpu],nT*sizeof(myfloat));
		cudaMalloc((void**)&tex_d[gpu], sizeof(cudaTextureObject_t) );
		cudaMalloc((void**)&tex_prob_d[gpu], sizeof(cudaTextureObject_t) );
		narrowband_copy(narrBand, wvc_d[gpu], idx_d[gpu], tex_d[gpu], tex_prob_d[gpu], Tnb_d[gpu], Tnb);

		// CUDA memory allocation
		int tm = (int) ((Tmax - Tnb[0])/(Tnb[1]-Tnb[0]));
		for(int j = 0; j < 2; j++)
			for(int i = 0; i < nB; i++ )
			{
				prob_h2[idx_p(j,i)] = prob[tm+j][i];
			}
		cudaMemcpyToSymbol(probd, prob_h2 , nB*2*sizeof(myfloat) ,0,cudaMemcpyHostToDevice);
		cudaCheckErrors("Malloc fail");

		/**************************************************************/
		/***********  STARTING CUDA ROUTINES **************************/
		/**************************************************************/

		/**************************************************************/
		/***************** FLUID MONTE CARLO **************************/
		/**************************************************************/



		for (int i=0; i<num_streams; i++)
		{
			if(gpu==0) {
				// launch one worker kernel per stream
				kernel_fluid<<<nblocks, blockdim, 0, streams0[i]>>>(gridGPU[gpu], wvc_d[gpu], Tnb_d[gpu], N, num_streams, i, kappamax, Tmax, idx_d[gpu],
						tex_tempf_d[gpu], tex_d[gpu], tex_prob_d[gpu], gpu, *ystart, Ibw_d[gpu]);
			}
			else {
				kernel_fluid<<<nblocks, blockdim, 0, streams1[i]>>>(gridGPU[gpu], wvc_d[gpu], Tnb_d[gpu], N, num_streams, i, kappamax, Tmax, idx_d[gpu],
						tex_tempf_d[gpu], tex_d[gpu], tex_prob_d[gpu], gpu, *ystart, Ibw_d[gpu]);
			}
		}
		cudaCheckErrors("Failed kernel execution");
	}

	//freeing all the CPU used variables (GPU is freed automatically by cudaDeviceReset(); )
	for(int grd=0; grd<grid_num; grd++)
	{
		varCPU[grd].destroyVar();
		gridCPU[grd].destroyVar();
	}
	free(prob_h);
	free(probg_h);
	free(prob_h2);
	free(narrBand);
	free(Tnb);

}

extern "C" void get_results_(myfloatF resfort[(imax+2)*(jmax/p_row+2)*(kmax+2)], myfloatF varfort[(imax+2)*(jmax/p_row+2)*(kmax+2)])
{
	// clock_t start,end;
	// start = clock();
	myfloat *host[num_streams][num_gpu];
	myfloat *varh[num_streams][num_gpu];
	for (int i=0; i<num_streams; i++) {
		for(int gpu = 0; gpu < num_gpu; gpu++) {
			host[i][gpu] = (myfloat*)malloc(N*sizeof(myfloat));
			varh[i][gpu] = (myfloat*)malloc(N*sizeof(myfloat));
		}
	}
	myfloat *device0[num_streams];
	myfloat *device1[num_streams];
	myfloat *vard0[num_streams];
	myfloat *vard1[num_streams];
	cudaStream_t streams0[num_streams];
	cudaStream_t streams1[num_streams];

	for(int gpu = 0; gpu < num_gpu; gpu++) {
		cudaSetDevice(gpu);
		cudaDeviceSynchronize();
		for (int i=0; i<num_streams; i++)
		{
			if(gpu==0) {
				cudaMalloc((void**)&device0[i]  ,N * sizeof(myfloat));
				cudaCheckErrors("Malloc fail device");
				cudaMalloc((void**)&vard0[i]  ,N * sizeof(myfloat));
				cudaCheckErrors("Malloc fail device");
				cudaStreamCreate(&streams0[i]);
			}
			else {
				cudaMalloc((void**)&device1[i]  ,N * sizeof(myfloat));
				cudaCheckErrors("Malloc fail device");
				cudaMalloc((void**)&vard1[i]  ,N * sizeof(myfloat));
				cudaCheckErrors("Malloc fail device");
				cudaStreamCreate(&streams1[i]);
			}
		}


		for (int i=0; i<num_streams; i++)
		{
			if(gpu==0) {
				// launch one worker kernel per stream
				kernel_results<<<nblocks, blockdim, 0, streams0[i]>>>(vard0[i],device0[i], N, num_streams, i);
			}
			else {
				kernel_results<<<nblocks, blockdim, 0, streams1[i]>>>(vard1[i],device1[i], N, num_streams, i);
			}
		}
		cudaCheckErrors("Failed kernel execution");
	}
	for(int gpu = 0; gpu < num_gpu; gpu++) {
		cudaSetDevice(gpu);
		for (int i = 0; i < num_streams; i++)
		{
			if(gpu==0) {
				cudaMemcpyAsync(host[i][gpu],device0[i],N*sizeof(myfloat),cudaMemcpyDeviceToHost,streams0[i]);
				cudaCheckErrors("Cuda memory copy asynchronous, device to host");
				cudaMemcpyAsync(varh[i][gpu],vard0[i],N*sizeof(myfloat),cudaMemcpyDeviceToHost,streams0[i]);
				cudaCheckErrors("Cuda memory copy asynchronous, device to host");
			}
			else {
				cudaMemcpyAsync(host[i][gpu],device1[i],N*sizeof(myfloat),cudaMemcpyDeviceToHost,streams1[i]);
				cudaCheckErrors("Cuda memory copy asynchronous, device to host");
				cudaMemcpyAsync(varh[i][gpu],vard1[i],N*sizeof(myfloat),cudaMemcpyDeviceToHost,streams0[i]);
				cudaCheckErrors("Cuda memory copy asynchronous, device to host");
			}
		}
		cudaCheckErrors("Copying to host fail");
	}

	cudaCheckErrors("unbind and/or free fail");

	/**************************************************************/
	/***************** RESETTING DEVICE MEMORY ********************/
	/**************************************************************/


	for(int gpu = 0; gpu < num_gpu; gpu++)
	{
		cudaSetDevice(gpu);

#if (bdw==2) || (bde==2) || (bdn==2) || (bds==2) || (bdb==2) || (bdt==2)
		cudaFree(Ibw_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
#endif
		cudaFree(gridGPU[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(tex_tempf_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(idx_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(wvc_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(Tnb_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(tex_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFree(tex_prob_d[gpu]);
		cudaCheckErrors("unbind and/or free fail");
		cudaFreeArray(cuArray[gpu]);
		cudaCheckErrors("unbind and/or free fail");

		for(int i=0; i<num_streams; i++) {
			if(gpu==0) {
				cudaFree(device0[i]);
				cudaFree(vard0[i]);
			}
			else {
				cudaFree(device1[i]);
				cudaFree(vard1[i]);
			}
		}
		cudaDeviceReset();
	}

	/**************************************************************/
	/************ RETURNING RESULTS IN A 3D FASHION ***************/
	/**************************************************************/

	myfloat result[imax+2][jmax/p_row+2][kmax+2];

	int i, j, k;
	for(int stream = 0; stream < num_streams; stream++) {
		for(int n = 0; n < N; n++) {
			for(int gpu = 0; gpu < num_gpu; gpu++) {
				i = (n / (kmax*jmax/num_gpu/p_row) + 1) + (stream)*imax/num_streams;
				j = (n / (kmax/num_gpu) + 1 - (i-1-(stream)*imax/num_streams)*jmax/p_row);
				k = (n - (kmax/num_gpu) * (j-1 + (i-1-(stream)*imax/num_streams)*jmax/p_row) + 1) + gpu*kmax/num_gpu;
				result[i][j][k] = host[stream][gpu][n];
				resfort[idx_F(i,j,k)] = (myfloatF) result[i][j][k];
				varfort[idx_F(i,j,k)] = (myfloatF) powf(varh[stream][gpu][n],0.5);
			}
		}
	}
	for (int i=0; i<num_streams; i++) {
		for(int gpu = 0; gpu < num_gpu; gpu++) {
			free(host[i][gpu]);
			free(varh[i][gpu]);
		}
	}

}

__device__ __forceinline__ void kernel_find(Count *count, myfloat *Tnb, myfloat Tmax, curandState_t *state, cudaTextureObject_t *tex_prob)
{
	int tm = (int) ((Tmax - Tnb[0])/(Tnb[1]-Tnb[0]));

	int nb = 0;
	int g = 0;
	for (int h = 0; h<nB; h++)
	{
		count->nb_cnt[h] = 0;
		for (int f = 0; f<nQ; f++)
		{
			count->g_cnt[f][h] = 0;
		}
	}
	/***********  ENTERING THE PHOTON LOOP ****************************/
	for (int h=0; h<photmax; h++)
	{
		// now we have to define the absorption narrow band wavenumber and the
		// quadrature point from prob and probg
		find_band(Tnb,Tmax,&g,&nb,state,tm,tex_prob);
		count->nb_cnt[nb] += 1;
		count->g_cnt[g][nb] +=1;

		/***********  OUT OF THE VARIANCE LOOP ****************************/
	}

}
__device__ __forceinline__ void emiss_ang(Beam *ray, curandState_t *state)
{
	myfloat phi   = curand_uniform(state)*2*pi;
	myfloat theta = acosf( 1 - 2*curand_uniform(state) );
	ray->sx = __cosf(theta);
	ray->sy = __sinf(theta)*__cosf(phi);
	ray->sz = __sinf(theta)*__sinf(phi);
}
__device__ __forceinline__ void find_band(myfloat *Tnb, myfloat Tmax, int *g, int *nb, curandState_t *state, int tm, cudaTextureObject_t *tex_prob)
{

	myfloat Rwave = curand_uniform(state);
	// find index of temperature
	//temperature index is t and t+1, now search for R on t and t+1
	int t;
	if( (Tmax - Tnb[tm]) < (Tnb[tm+1] - Tmax) )
	{
		t = tm;
		int nb1 = 0;
		int w = nB;
		while (w > nb1+1)
		{
			int d = (nb1+w)/2;
			(void)((Rwave <  probd[t-tm][d] ) && (w=d));
			(void)((Rwave >= probd[t-tm][d] ) && (nb1=d));
		};
		nb1 = nb1+1;
		if( Rwave <  probd[t-tm][0] )
			nb1 = 0;
		if(Rwave == 1 )
			nb1 = nB -1;

		*nb = nb1;
	}
	else
	{
		t = tm+1;
		int nb2 = 0;
		int w = nB;
		while (w > nb2+1)
		{
			int d = (nb2+w)/2;
			(void)((Rwave <  probd[t-tm][d] ) && (w=d));
			(void)((Rwave >= probd[t-tm][d] ) && (nb2=d));
		};
		nb2 = nb2+1;
		if(Rwave <  probd[t-tm][0] )
			nb2 = 0;
		if(Rwave == 1 )
			nb2 = nB -1;
		*nb = nb2;
	}

#if grey == 1
	*g=0;
#else
	myfloat Rwave2 = curand_uniform(state);
	myfloat prob;

	*g =0;
	int w = nQ;
	while (w > *g+1)
	{
		int d = (*g+w)/2;
		prob = tex3D<myfloat>(*tex_prob, d+0.5f, *nb+0.5f, t+0.5f);
		(void)((Rwave2 < prob ) && ( w=d));
		(void)((Rwave2 >= prob ) && (*g=d));
	};
	*g = *g+1;
	prob = tex3D<myfloat>(*tex_prob, 0.5f, *nb+0.5f, t+0.5f);
	(void)((Rwave2 < prob ) && (*g = 0));
	(void)((Rwave2 == 1)    && (*g = nQ-1));
#endif

}
__device__ __forceinline__ void wave_find(int *nb, int *g, int *countnb, int *countg, int *idx_nb, Count *cnt)
{
	int tmp=1;
	while(tmp)
	{
		if(*countnb < cnt->nb_cnt[idx_nb[*nb]])
		{
			*countnb+=1;
			tmp=0;
		}
		else
		{
			*nb +=1;
			*countnb=0;
			*countg=0;
			*g=0;
			tmp=1;
		}
	};
	while(!tmp)
	{
		if(*countg < cnt->g_cnt[*g][idx_nb[*nb]])
		{
			*countg+=1;
			tmp=1;
		}
		else
		{
			*g +=1;
			*countg=0;
			tmp=0;
		}
	};
}


struct PAIR
{
	int pos;
	myfloat val;
};

bool compare(PAIR p1, PAIR p2) {return p1.val < p2.val;}

void sort_idx(NarrowBand *narrBand)
{
	vector<PAIR> p(nB);

	for (int nb = 0; nb<nB; nb++)
	{
		p[nb].pos = nb;
		p[nb].val = narrBand[nb].kavg;
	}

	sort( p.begin(), p.end(), compare );

	for (int nb = 0; nb<nB; nb++)
	{
		narrBand[nb].idx = p[nb].pos;
	}
}

__global__ void kernel_results(myfloat *var, myfloat *device, int n, int ns, int stream)
{

	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	volatile int tx = threadIdx.x;

	for (int idx = tid; idx < n; idx += blockDim.x * gridDim.x)
	{
		device[idx] = solution[stream][idx];
		var[idx]    = variance[stream][idx];
	}
}
