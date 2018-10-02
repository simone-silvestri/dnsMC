/* Program that calculates the mie coefficients for different particle clouds,
 * based on the particle radius, the particle coefficients and the volume
 * fraction of the cloud,
 *
 * SIMONE SILVESTRI - 28 / 03 / 2017
 */


#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <complex.h>
#include <deque>

using namespace std;

#define PI  3.14159265358979323846
#define wmax 530
#define thetamax 180
const double rad = 1.0e-1;
const double fv  = 5.e-7;

int v1;

class RefIdx
{
public:
   deque<double> wv,kv,nv;
};

void Refraction_index(RefIdx &refidx, char *namefile);
void rjbes(double x, double *sol, double *dsol, int NMAX, int NNMAX);
void bess(double *j, double *y, double *jr, double *ji, double complex *h, double complex *dh, 
          double *dj, double *dy, double *djr, double *dji, double x, double xr, double xi, int NMAX, int NMAX2);
void rybes(double x, double *sol, double *dsol, int NMAX);
void cjbes(double xr, double xi, double *solr, double *soli, double *dsolr,
		   double *dsoli, int NMAX, int NNMAX);
void scattering(double x, double Csca, double complex *an, double complex *bn, double *phi, int NMAX);


int main()
{

	/**************************************************************/
	/****************** READ REFRACTION INDEX *********************/
	/**************************************************************/

	char namefile[] ="DATA/AG";

	RefIdx r;
	Refraction_index(r, namefile);

	int NMAX,NMAX2;

	/**************************************************************/
	/****************** LOOP OVER WAVENUMBERS *********************/
	/**************************************************************/

	double Csca[r.wv.size()],Cscar[r.wv.size()];
	double Cext[r.wv.size()],Cextr[r.wv.size()];
	double phi[r.wv.size()][thetamax];

	FILE *fp  = fopen("coeff_mie","w+");
	FILE *fp4 = fopen("eff_mie","w+");
	FILE *fp2 = fopen("coeff_ray","w+");
	FILE *fp3 = fopen("coeff_nano","w+");
	FILE *fp1= fopen("scatt","w+");
	FILE *fp5 = fopen("scatt_to_cpy","w+");

	for(unsigned int v = 0; v<r.wv.size(); v++)
	{
		v1 = v;
		Csca[v] = 0;
		Cext[v] = 0;


		//calculate the size parameter where the radius and the wavenumber are in cm and cm-1
		double x = 2.0 * PI * rad * r.wv[v];
	
                if(x>30) {	
	                NMAX=1000;
        	        NMAX2=5000;
		}
		else if(x<=30 && x>10) {
                        NMAX=1000;
 			NMAX2=5000;
		}
                else {
                    	NMAX=500;
			NMAX2=1000;
                }



		double xr= r.nv[v];
		double xi= r.kv[v];
		double complex m = xr + xi * I;
		double complex c = m * x ;
		double Conc = fv / ( rad * rad * rad * 1.e-6 * M_PI * 4./3. );

		double complex an[NMAX+1];
		double complex bn[NMAX+1];
		double j[NMAX+1],y[NMAX+1],jr[NMAX+1],ji[NMAX+1];
		double complex h[NMAX+1],dh[NMAX+1];
		double complex jc[NMAX+1], djc[NMAX+1];
		double dj[NMAX+1],dy[NMAX+1],djr[NMAX+1],dji[NMAX+1];

		bess(j,y,jr,ji,h,dh,dj,dy,djr,dji,x,creal(c),cimag(c),NMAX,NMAX2);
		
		for(int n=1; n<NMAX+1; n++ )
		{
			h[n]   = x * (j[n] + I * y[n]);
			dj[n]  = x * dj[n];
			dy[n]  = x * dy[n];
			dh[n]  = (dj[n] + I * dy[n]) ;

			j[n]   = x * j[n];
			y[n]   = x * y[n];

			jc[n]  = c * (jr[n]  + ji[n] * I);
			djc[n] = c * (djr[n] + dji[n]* I);

			an[n] = (djc[n] * j[n] - m * jc[n] * dj[n]);
			an[n]/= (djc[n] * h[n] - m * jc[n] * dh[n]);
			bn[n] = (m * djc[n] * j[n] - jc[n] * dj[n]);
			bn[n]/= (m * djc[n] * h[n] - jc[n] * dh[n]);

 			if((isnan(creal(an[n])))||isnan(cimag(an[n])))
				an[n] = 0 + 0*I;
 			if((isnan(creal(bn[n])))||isnan(cimag(bn[n])))
				bn[n] = 0 + 0*I;
			Csca[v] += 2/(x*x) * (2*n+1) * (cabsl(an[n])*cabsl(an[n]) + cabsl(bn[n])*cabsl(bn[n])) * M_PI * (rad * rad * 1.e-04) ;
			Cext[v] += 2/(x*x) * (2*n+1) * creall(an[n] + bn[n]) * M_PI * (rad * rad * 1.e-04) ;

		}
		Cscar[v] = 8.0L/3.0L * cabsl((m*m - 1)/(m*m + 2)) * cabsl((m*m - 1)/(m*m + 2)) * powl(x,4) * M_PI * (rad * rad * 1.e-04);
		Cextr[v] = 4.0L * cimagl((m*m - 1)/(m*m + 2)) * x * M_PI * (rad * rad * 1.e-04);
           
		double phit[thetamax];
		double integral;
		for (int t=0; t<thetamax; t++)
		{
			phit[t] = 0;
			phi[v][t] = 0;
		}

		NMAX=100;
		scattering(x,Csca[v],an,bn,phit,NMAX);
		for (int t=0; t<thetamax; t++)
		{
			double theta = t*PI/(thetamax*1.0);
			phi[v][t] = phit[t] * 2 /(x*x*Csca[v] / (M_PI * (rad * rad)) ) ;
			fprintf(fp5,"%le %le %le\n",r.wv[v],theta,phi[v][t]);
		}
		for (int t=1; t<thetamax; t++)
		{
			double theta = t*M_PI/(thetamax*1.0);
			double thetap = (t-1)*M_PI/(thetamax*1.0);
			integral += (theta - thetap)*(phi[v][t]+phi[v][t-1])/2.0L * (sinl(theta)+sinl(thetap))/2.0L;
		}
		if(v>5)
			if( ((Cext[v]-Cext[v-1])/Cext[v-1]>0.1) || ((Cext[v]-Cext[v-1])/Cext[v-1]<-0.1) ){
				Cext[v] = (Cext[v-1]*2.0L-Cext[v-2]) ;
				Csca[v] = (Csca[v-1]*2.0L-Csca[v-2]) ;
			}
		fprintf(fp, "%d %le %le %le %le %le\n",v,r.wv[v] , Csca[v] , Cext[v] , Cscar[v], Cextr[v]);
		fprintf(fp4,"%le %le %le %le\n",x,r.wv[v], Csca[v]/(M_PI * (rad * rad * 1.e-04)), Cext[v]/(M_PI * (rad * rad * 1.e-04)));
		fprintf(fp2,"%d %le %le %le %le %le\n",v,r.wv[v], Cscar[v], Cextr[v], Csca[v] , Cext[v]);
		fprintf(fp3,"%d %le %le %le %le %le\n",v,r.wv[v], Cscar[v]*Conc, Cextr[v]*Conc,Csca[v], Cext[v]);
	}
	fclose(fp);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);
	fclose(fp5);
	for (int t=0; t<thetamax; t++)
	{
		double theta = t*180.0L/(thetamax*1.0);

		fprintf(fp1,"%le\t",theta);

		for(unsigned int v = 0; v<r.wv.size(); v++)
		{
			fprintf(fp1,"%le\t",phi[v][t]);
		}
		fprintf(fp1,"\n");
	}
	fclose(fp1);

	printf("Ciao\n");

	return 0;
}


void Refraction_index(RefIdx &r, char *namefile)
{
	double wv,kv,nv;

	FILE *fp = fopen(namefile,"r");

	int i = 0, n = 0;
	while (fscanf(fp,"%le %le %le\n",&wv,&nv,&kv) != EOF)
	{
		if(wv>0.1)
		{
			r.wv.push_front(1.e7/wv/1000);
			r.nv.push_front(nv);
			r.kv.push_front(kv);
			printf("%lf %lf %lf\n",1.e7/wv/1000,nv,kv);
		}
	};
	fclose(fp);
	printf("size of refraction index: %d\n",r.wv.size());
}
 
void bess(double *j, double *y, double *jr, double *ji, double complex *h, double complex *dh, 
          double *dj, double *dy, double *djr, double *dji, double x, double xr, double xi, int NMAX, int NMAX2)
{
	double aj[NMAX+1] ,ay[NMAX+1] ,ajr[NMAX+1] ,aji[NMAX+1];
	double adj[NMAX+1],ady[NMAX+1],adjr[NMAX+1],adji[NMAX+1];
	double xx;
	xx = x;
	rjbes(xx,aj,adj,NMAX,NMAX2);
	rybes(xx,ay,ady,NMAX);
	double xxr = xr;
	double xxi = xi;
	cjbes(xxr,xxi,ajr,aji,adjr,adji,NMAX,NMAX2);
	for(int n=0; n<NMAX+1; n++)
	{
		j[n]   = aj[n];
		y[n]   = ay[n];
		jr[n]  = ajr[n];
		ji[n]  = aji[n];
		dj[n]  = adj[n];
		dy[n]  = ady[n];
		djr[n] = adjr[n];
		dji[n] = adji[n];
	}
}
void rjbes(double x, double *sol, double *dsol, int NMAX, int NNMAX)
{
	double z[NMAX+NNMAX+10];
	int l = NMAX+NNMAX+1;
	double xx = 1.0/x;
	z[l] = 1.0/((double)(2*l+1)*xx);
	int l1 = l-1;
	for(int i=1; i<l1+1; i++)
	{
		int i1 = l-i;
		z[i1] = 1.0/((double)(2*i1+1)*xx - z[i1+1]);
	}
	double z0 = 1.0/(xx-z[1]);
	double y0 = z0*cos(x)*xx;
	double y1 = y0 * z[1];
	dsol[1] = y0 - y1*xx;
	sol[1] = y1;
	double yi1,yi;
	for(int i=2; i<NMAX+1; i++)
	{
		yi1 = sol[i-1];
		yi = yi1*z[i];
		sol[i]  = yi;
		dsol[i] = yi1 - (double) (i) * yi * xx;
	}
}
void rybes(double x, double *sol, double *dsol, int NMAX)
{
	double c  = cos(x);
	double s  = sin(x);
	double x1 = 1.0/x;
	double x2 = x1*x1;
	double x3 = x2*x1;
	double y1 = - c * x2 - s * x1;
	sol[1] = y1;
	sol[2] = (-3.0*x3 + x1)*c - 3.0*x2*s;
	for(int i=2; i<NMAX+1; i++)
		sol[i+1] = (double)(2*i+1)*x1*sol[i] - sol[i-1];
	dsol[1] = -x1*(c+y1);
	for(int i=2; i<NMAX+1; i++)
		dsol[i] = sol[i-1] - (double)(i)*x1*sol[i];
}
void cjbes(double xr, double xi, double *solr, double *soli, double *dsolr,
		   double *dsoli, int NMAX, int NNMAX)
{
	double cyr[NMAX+1],cyi[NMAX+1],czr[NMAX+NNMAX+10],czi[NMAX+NNMAX+10],cur[NMAX+1],cui[NMAX+1];
	double ar,ai,ari;
	int l = NMAX+NNMAX+1;
	double xrxi = 1.0/(xr*xr+xi*xi);
	double cxxr = xr*xrxi;
	double cxxi = -xi*xrxi;
	double qf = 1.0 / (double)(2*l+1);
	czr[l] = xr*qf;
	czi[l] = xi*qf;
	int l1=l-1;
	for(int i=1; i<l1+1; i++)
	{
		int i1 = l-i;
		qf = (double)(2*i1+1);
		ar = qf*cxxr - czr[i1+1];
		ai = qf*cxxi - czi[i1+1];
		ari= 1.0 / (ar*ar + ai*ai);
		czr[i1] =  ar*ari;
		czi[i1] = -ai*ari;
	}
	ar = cxxr - czr[1];
	ai = cxxi - czi[1];
	ari= 1.0 / (ar*ar + ai*ai);

	double cz0r =  ar*ari;
	double cz0i = -ai*ari;
	double cr =  cos(xr)*cosh(xi);
	double ci = -sin(xr)*sinh(xi);
	ar = cz0r*cr-cz0i*ci;
	ai = cz0i*cr+cz0r*ci;

	double cy0r = ar*cxxr-ai*cxxi;
	double cy0i = ai*cxxr+ar*cxxi;
	double cy1r = cy0r*czr[1]-cy0i*czi[1];
	double cy1i = cy0i*czr[1]+cy0r*czi[1];
	ar = cy1r*cxxr-cy1i*cxxi;
	ai = cy1i*cxxr+cy1r*cxxi;

	double cu1r = cy0r - ar;
	double cu1i = cy0i - ai;
	cyr[1]  =cy1r;
	cyi[1]  =cy1i;
	cur[1]  =cu1r;
	cui[1]  =cu1i;
	solr[1] =cy1r;
	soli[1] =cy1i;
	dsolr[1]=cu1r;
	dsoli[1]=cu1i;


	double cyi1r,cyi1i,cyir,cyii,cuir,cuii;
	for(int i=2; i<NMAX+1; i++)
	{
		double qi = (double) i;
		cyi1r = cyr[i-1];
		cyi1i = cyi[i-1];
		cyir = cyi1r*czr[i] - cyi1i*czi[i];
		cyii = cyi1i*czr[i] + cyi1r*czi[i];
		ar = cyir*cxxr - cyii*cxxi;
		ai = cyii*cxxr + cyir*cxxi;
		cuir = cyi1r - qi*ar;
		cuii = cyi1i - qi*ai;

		cyr[i]  =cyir;
		cyi[i]  =cyii;
		cur[i]  =cuir;
		cui[i]  =cuii;
		solr[i] =cyir;
		soli[i] =cyii;
		dsolr[i]=cuir;
		dsoli[i]=cuii;
	}
}
void scattering(double x, double Csca, double complex *an, double complex *bn, double *phi, int NMAX)
{
	double complex S1[thetamax],S2[thetamax];
	double pn[NMAX+1], tn[NMAX+1];
	double p1,p2,t1,t2;
	for (int t=0; t<thetamax; t++)
	{
		double theta = t*PI/(double)(thetamax);
		double mu = cosl(theta);
		pn[0] = 0;
		pn[1] = 1;
		tn[1] = mu*pn[1] - 2*pn[0];
		for (int n=2; n<NMAX+1; n++)
		{
			p1 = pn[n-1];
			p2 = pn[n-2];
			pn[n] = ((2.0L*n-1.0L)*mu*p1 - n*p2)/(n-1.0L);
			t1 = n * mu * pn[n];
			t2 = (n+1) * pn[n-1];
			tn[n] = t1-t2;
		}
		S1[t] = 0.0L + 0.0L * I;
		S2[t] = 0.0L + 0.0L * I;
		double temp;
		for (int n=1; n<NMAX+1; n++)
		{
			temp = (2.0L*n+1.0L)/(n*(n+1.0L));
			S1[t] += temp*(an[n]*pn[n] + bn[n]*tn[n]);
			S2[t] += temp*(bn[n]*pn[n] + an[n]*tn[n]);
		}
		phi[t]=creall(conjl(S1[t])*S1[t])+creall(conjl(S2[t])*S2[t]);
	}
}
