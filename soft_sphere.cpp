#include <iostream>
#include <fstream>
#include <math.h>
#include <time.h>
#include <iomanip>

using namespace std;

// double BOX=2*7.2005488787721159;
// double BOX=2*26.2107109945647245;
double BOX;
int N=1000;
double epsilon=1.0;
double magnitude=0.;
double p_applied=1e-2;
double STRAIN=0.0;
double tilt=STRAIN*BOX;

double *ATOMx = new  double[N];
double *ATOMy = new  double[N];
// double *ATOMz = new  double[N];
double *RAD   = new  double[N];
//static  double *DIST= new  double[N];
static int *VLIST= new int [N*(N+1)];

static  double *tATOMx = new  double[N];
static  double *tATOMy = new  double[N];
//static  double *tATOMz = new  double[N];
// gradient is calculated and stored in this 

//double *CONJ_D= new  double[2*N+1];
//double *G= new  double[2*N+1];
double *X= new  double[2*N+1];
static double *tX= new  double[2*N+1];
double deltadrv=0.5;
double deltadrvsq=deltadrv*deltadrv;

void make_list()
{
	for(int i=0; i<2*N; i++)
	{
		tX[i]=-999;
		//tATOMy[i]=-999;
//		tATOMz[i]=ATOMz[i];
	}
	//cout<<tATOMx[0]<<"\t"<<ATOMx[0]<<" htis \n";
}

double energy_force( double X[], double GRAD [], int ndim)
{
	
	//BOX=X[ndim-1];
	clock_t begin=clock();
	double delta_rx,delta_ry,delta_rz;
	double delta_r;
	double pot_e=0;
    double skin_depth=1.4+deltadrv;
    double skin_depth_sq=skin_depth*skin_depth;
    int nnei=0;
    int nlistbeg=0;
	double iBOX=1./BOX;
	
	double maxdisp=0.;
	double maxdisp2=0.;
	//cout<<tATOMx[0]<<"\t"<<X[0]<<"\n";
	for(int i=0; i<N; i++)
	{
		ATOMx[i]=X[i*2];
		ATOMy[i]=X[i*2+1];
		GRAD[i*2] = 0.;
		GRAD[i*2+1] = 0.;
//		ATOMz[i]=X[i*3+2];

		delta_rx=(tX[2*i]-X[i*2]);
		delta_ry=(tX[2*i+1]-X[i*2+1]);
		//cout<<i<<"\t"<<delta_rx<<"\n";
//		delta_rz=tATOMz[i]-ATOMz[i];
		    
		if(tilt)
	    	delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
	    delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
  	    delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
  //	    delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));
		
		delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
		
		delta_r=sqrt(delta_r);

		if(delta_r > maxdisp)
		{
			maxdisp2=maxdisp;
			maxdisp=delta_r;
			
		} 
		else if(delta_r > maxdisp2)
		{
			maxdisp2=delta_r;
		}
		//cout<<ATOMx[i]<<"\n";
	}	
	//cout<<"( energyforce ) \t"<<maxdisp<<"\t"<<maxdisp2<<"\t"<<maxdisp+maxdisp2<<"\t"<<deltadrv<<"\n" ;
	if(maxdisp+maxdisp2 > deltadrv )
	{
		//cout<<"inside ( energyforce ) \t"<<maxdisp<<"\t"<<maxdisp2<<"\t"<<maxdisp+maxdisp2<<"\t"<<deltadrv<<"\n" ;
		//cout<<tATOMx[0]<<"\t"<<X[0]<<"\n";
		//make_list();
		for(int i=0; i<N; i++)
		{
			nnei=0;
			
			for(int j=0 ; j<N; j++)
			{
				if ( i != j )
				{
					delta_rx=(X[2*i]-X[2*j]);
					delta_ry=(X[2*i+1]-X[2*j+1]);
					//delta_rz=ATOMz[i]-ATOMz[j];

					if(tilt)
						delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
					delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
					delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
					//delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));

					delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
					
					//delta_r=sqrt(delta_r);

					if(delta_r< skin_depth_sq)
					{
						nnei=nnei+1;
						VLIST[nlistbeg+nnei]=j;
					}
				}
			}
			VLIST[nlistbeg]=nnei;
			nlistbeg=nlistbeg+nnei+1;
		}
		for(int i=0; i<2*N; i++)
		{
			tX[i]=X[i];
		}
	}
	nlistbeg=0;
	for(int i=0; i<N; i++)
	{
		int jj;
		//cout<<i<<"\n";
		for(jj=nlistbeg+1; jj<nlistbeg+VLIST[nlistbeg]+1; jj++)
		{
			int	j=VLIST[jj];
			if(j>i)
			{
				delta_rx=ATOMx[i]-ATOMx[j];
				delta_ry=ATOMy[i]-ATOMy[j];
//				delta_rz=ATOMz[i]-ATOMz[j];

				if(tilt)
					delta_rx=(delta_rx-(tilt*lroundl(delta_ry*iBOX)));
				delta_rx=(delta_rx-(BOX*lroundl(delta_rx*iBOX)));
				delta_ry=(delta_ry-(BOX*lroundl(delta_ry*iBOX)));
//					delta_rz=(delta_rz-(BOX*lroundl(delta_rz/BOX)));

				delta_r=delta_rx*delta_rx+delta_ry*delta_ry;//+delta_rz*delta_rz;
				
				//delta_r=pow(delta_r,0.5);
				double sigma=RAD[i]+RAD[j];
				double sigmasq=sigma*sigma;
				//cout<<i<<"\t"<<j<<"\t"<<delta_r<<"\n";
				if(delta_r<sigmasq)
				{
					//pot_e=pot_e+U_ij(delta_r,RAD[i]+RAD[j]);
					double isigma=1./sigma;
					delta_r=sqrt(delta_r);
					double vh=(1-delta_r*isigma);
					pot_e=pot_e+epsilon*vh*vh;
					//dwreturn epsilon*2.*(1-r/sigma)*(-1./sigma);
					double U=epsilon*2.*vh*(-1.*isigma)/delta_r;
					GRAD[i*2]   = GRAD[i*2]   + U*delta_rx;
					GRAD[i*2+1] = GRAD[i*2+1] + U*delta_ry;

					GRAD[j*2]   = GRAD[j*2]   - U*delta_rx;
					GRAD[j*2+1] = GRAD[j*2+1] - U*delta_ry;
				}
			}
		}
		nlistbeg=jj;//nlistbeg+jj+1;
	}
	clock_t end=clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout<<time_spent<<" in energy \n";
////for(int i=0; i<N; i++)
////{
////	GRAD[2*i]=GRAD[2*i]*BOX;
////	GRAD[2*i+1]=GRAD[2*i+1]*BOX;
////	GRAD[ndim-1]=GRAD[ndim-1]+(GRAD[2*i]/BOX*ATOMx[i]+GRAD[2*i+1]/BOX*ATOMy[i]);
////}
////GRAD[ndim-1]=GRAD[ndim-1]/BOX+2*BOX*p_applied;
	return pot_e;
}
void write_config(double k,  double X[], char name[])
{
////for(int i=0;i <N; i++)
////{
////	ATOMx[i]=X[i*3];
////	ATOMy[i]=X[i*3+1];
////	ATOMz[i]=X[i*3+2];
////}	
	fstream input(name+std::to_string(k)+".xyz", fstream::out);
	input<<N<<"\n";
	input<<std::setprecision(16);
	input<<BOX<<"\n";
	for(int i=0; i<N; i++)
	{
		if(tilt)
			X[i]=(X[i]-(tilt*lroundf(X[i+N]/BOX)));
		X[i]=(X[i]-(BOX*lroundf(X[i]/BOX)));
		X[i+N]=(X[i+N]-(BOX*lroundf(X[i+N]/BOX)));

		input<<i+1<<"\t"<<X[2*i]<<"\t"<<X[2*i+1]<<"\t"<<RAD[i]<<"\n";
	}
}
