#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "conj_grad.h"
#include "soft_sphere.h"

using namespace std;

int main(int argc, char *argv[])
{
    int rand();
    char buffer[64];
    fstream input;
    input.open(argv[1]);
    BOX=2.*stod(argv[2]);
    double a,b,c,d,e;
    while(input>>a>>b>>c>>d)
    {
        ATOMx[int(a)-1]=b;
        ATOMy[int(a)-1]=c;
        //ATOMz[int(a)-1]=e;
        //if(b==1)
            RAD[int(a)-1]=d;
        //if(b==2)
          //  RAD[int(a)-1]=0.7;
      	 //cout<<ATOMx[int(a)-1]<<"\t"<<ATOMy[int(a)-1]<<"\t"<<RAD[int(a)-1]<<"\n";
    }
	
	int ndim=2*N;
	for(int i=0;i <N; i++)
	{
		X[i*2] =ATOMx[i];
		X[i*2+1]=ATOMy[i];
		//X[i*3+2]=ATOMz[i];
	}	
	//X[ndim-1]=BOX;
	make_list();
////cout<<"this "<<energy_force(X,G,ndim)/2000<<"\n";
//////calculate_gradient(X,G,ndim);
////double XX=0.;
////for(int i=0; i<N; i++)
////{
////	//cout<<i<<"\t"<<G[2*i+1]<<"\n";
////	XX=XX+G[i]*G[i];
////	//G[i]=GRAD[i];
////}
////cout<<XX<<"\n";
	minimize(ndim,X,100000,energy_force,write_config);
	int compress=1;
	double dphi=0.0001;
	double dphi_accum=0.;
	double fact=N/2.*3.14159265359*(0.5*0.5+0.7*0.7);
	double phi;
	double phi_new;

	double BOX_new;
	double rat;
	//return 0;
	fstream energy("energy_compress.dat", fstream::out);
	while(1)
	{
		phi=fact/(BOX*BOX);
		cout<<phi<<"\t"<<BOX<<"\t"<<compress<<"\n";
		phi_new=phi+compress*dphi;
		dphi_accum=dphi_accum+compress*dphi;

		BOX_new=sqrt(fact/phi_new);
		rat=(BOX_new/BOX);
		BOX=BOX_new;
		for(int i=0;i <N; i++)
		{
			X[i*2]  = rat*X[i*2];
			X[i*2+1]= rat*X[i*2+1];
			//X[i*3+2]=ATOMz[i];
		}	
		//print('%.16f'%(BOX_new/2.))
		//for i in range (0, 1000):
		double e=minimize(ndim,X,100000,energy_force,write_config);
		energy<<std::setprecision(16);
		energy<<phi_new<<"\t"<<e<<"\n";
		energy<<std::flush;
		if(e>1.e-7)
		{
			compress=-1;
			dphi=0.00001;
			//dphi=dphi/10.;
			//break;
		}
	        if(compress == -1 && e<1.e-16)
	        {
			write_config(phi_new,X,"compressed");
	        	cout<<"phi_J ="<<phi_new<<"\t energy= "<<e<<"\n";
	        	break;
	        }

	}
}
