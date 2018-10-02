// NARROWBAND MONTE CARLO OERM GPU

//simulation parameters
#define  num_gpu 1    //2
#define  p_row   1    //32 
#define	 imax	 64   //190  // REMEMBER!!!! imax must be a multiple of Num_streams!!!!! otherwise seg_fault
#define	 jmax	 64   //192
#define	 kmax	 64   //192
#define	 Lx	 1.f  //2.f
#define  Ly      1.f  //6.2831853f
#define  Lz      1.f  //12.566371f
//variance block calculator and photon per block
#define  nVar 	 3 
#define	 photmax 20000
#define  photmaS 1

//number of grid used (maximum 4, triggers preprocessor #if statements)
#define grid_num 5 //6 

//const int maxi[] = {imax, imax  , imax  , imax/5, imax/5 , imax/10};
const int maxi[] = {imax, imax/2, imax/4, imax/8, imax/16, imax/32};
const int maxj[] = {jmax, jmax/2, jmax/4, jmax/8, jmax/16, jmax/32};
const int maxk[] = {kmax, kmax/2, kmax/4, kmax/8, kmax/16, kmax/32};
//maximum number of steps in the grids
#if grid_num==1
const int maxs[] = {100000};
#elif grid_num==2
const int maxs[] = {5, 100000};
#elif grid_num==3
const int maxs[] = {5, 5, 100000};
#elif grid_num==4
const int maxs[] = {5, 5, 5, 100000};
#elif grid_num==5
const int maxs[] = {5, 5, 5, 5, 100000};
#elif grid_num==6
const int maxs[] = {5, 5, 5, 5, 5, 100000};
#endif

//GPU details
#define num_streams     16     //10
#define blockdim 	128
#define nblocks 	32
#define blockdimS       256
#define nblocksS	128

//sorting or not sorting the narrow bands (1 -> yes, other int -> no)
#define srt 1
//Spectral radiation models (0 -> nbcK (fluid), 1 -> grey, 2 -> nb (particles))
#define grey 2 
//Angular distribution (0 -> Uniform, other int -> Sobol)
#define sobol 1
// 1 -> double precision, 0 -> single precision
#define floatingpoint 0

//define the boundary values
#define	 epsw	 1.f
#define	 epse	 1.f
#define	 epss	 1.f
#define	 epsn	 1.f
#define	 epsb	 1.f
#define	 epst	 1.f

#define	 Tww	 500 //955 
#define	 Twe	 500 //573
#define	 Tws	 0
#define	 Twn	 0
#define	 Twb	 0
#define	 Twt	 0

#define	 stefan	 5.670373e-08
#define	 scatco	 0.0f
#define	 toll	 1.0e-3
#define	 pi	 	 3.14159265358979323846

// 1 -> periodic walls
// 2 -> black walls
// 3 -> diffuse grey walls (not implemented in GPU as of now)
#define bdw 2
#define bde 2
#define bds 1
#define bdn 1
#define bdb 1
#define bdt 1

#define cudaCheckErrors(msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s (%s at %s:%d)\n", \
                msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0)


#define idx_p(t,b) \
		({ t*nB + b; }) //for values going from 0 to max-1, 2-D
#define idx_I(t,b) \
		({ t + b*nT; }) //for values going from 0 to max-1, 2-D
#define idx_pg(t,b,g) \
		({ g + nQ*( b + t*nB ) ; }) //for values going from 0 to max-1, 3-D
#define idx_T(i,j,k,im,jm) \
		({ i + (im+2)*( j + k*(jm+2) ) ; }) //for values going from 0 to max-1, 3-D
#define idx_S(i,k) \
		({ (k-1) + (i-1) * kmax ;} ) //for values going from 1 to max 2-D
#define idx_D(i,j,k) \
        ({ (k-1)+kmax*((j-1)+(i-1)*jmax); }) //for values going from 1 to max 3-D
#define idx_F(i,j,k) \
		({ i + (imax+2)*( j + k*(jmax/p_row+2) ) ; }) //for values going from 0 to max-1, 3-D

#define pow4(a) \
        ({ a*a*a*a ; })

#define pow3(a) \
        ({ a*a*a ; })

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a >= _b ? _a : _b; })
#define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a <= _b ? _a : _b; })

// modifying accordingly the read statements

#if floatingpoint == 1
#define myfloat double
#define reade "%le"
#define readf "%lf"
#elif floatingpoint == 0
#define myfloat float
#define reade "%e"
#define readf "%f"
#endif

#define myfloatF float 
