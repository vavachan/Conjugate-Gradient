/* CONSTANT PRESSURE ENERGY MINIMIZATION USING CONJUGATE GRADIENT 
 *
 * This code should minimize the enthalpy U+pV  where p is the applied pressure 
 * and V is the volume of the system
 * 
 * let x_0 be the initial vector, consiting of the inital atomic configuration
 * 
 *
 * varghese babu 20/03/2019 
 */
using namespace std;
#include<iostream>
#include<iomanip>
//#include<fstream>
#include<math.h>
//#include <cstdlib>
#include <limits>
#include <time.h>

//#include "soft_sphere_hip.h"

double *TMP;
double *TMP_VEC;
double function1dim(double alpha,double X[],double CONJ_D[], double (*f)(double *, int),int ndim )
{
	//cout<<"how many\n";
	//cout<<alpha<<" alpha\n";
	for(int i=0; i<ndim; i++)
	{
		TMP[i]=X[i]+alpha*CONJ_D[i];
	}
	//delete [] TMP;
	return f(TMP,ndim);
}
double dfunction1dim(double alpha,double X[],double CONJ_D[], double (*fdf)(double *, double*, int) ,int ndim )
{
	//cout<<"how many\n";
	for(int i=0; i<ndim; i++)
	{
		TMP[i]=X[i]+alpha*CONJ_D[i];
	}
	//delete [] TMP;
	double e=fdf(TMP,TMP_VEC,ndim);	

	double der=0.;
	for(int i=0; i<ndim; i++)
	{
		der=der+TMP_VEC[i]*CONJ_D[i];
		//TMP[i]=X[i]+alpha*CONJ_D[i];
	}
	return der;	
	//return f(TMP,ndim);
}
double fdfunction1dim(double alpha,double X[],double CONJ_D[], double (*fdf)(double *, double*, int) ,int ndim, double& energy)
{
	//cout<<"how many\n";
	//cout<<alpha<<" alpha fdf\n";
	for(int i=0; i<ndim; i++)
	{
		TMP[i]=X[i]+alpha*CONJ_D[i];
	}
	//delete [] TMP;
	energy=fdf(TMP,TMP_VEC,ndim);	

	double der=0.;
	for(int i=0; i<ndim; i++)
	{
		der=der+TMP_VEC[i]*CONJ_D[i];
		//TMP[i]=X[i]+alpha*CONJ_D[i];
	}
	return der;	
	//return f(TMP,ndim);
}
void bracket(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc, double X[],double CONJ_D[], double (*f) (double *,int),  int ndim)
{
/* In this function bracketing the minima is carried out 
 * This is copied from numerical recipes, because i am a useless fuck.
 */
	double GOLD=1.618034;
	double GLIMIT=100;
	double ulim;
	fa=function1dim(ax,X,CONJ_D,f,ndim);		
	fb=function1dim(bx,X,CONJ_D,f,ndim);		
	double fu;
	double u;
	//Here we swap the two data points such that fb<fa always
	if(fb>fa)
	{
		double tmp;	
		tmp=ax;
		ax=bx;
		bx=tmp;

		tmp=fa;
		fa=fb;
		fb=tmp;
		
	}
	
	cx=bx+GOLD*(bx-ax); // first guess of c using the golden section rule. 
	fc=function1dim(cx,X,CONJ_D,f,ndim);		
	
	while(1)
	{
		// We will have to come here till the minima is bracketed 
		//cout << ax << "\t" <<bx <<"\t" << cx <<" this  guy\n";
		if(fb > fc) // Then the minima is not bracketed
		{
			double r=(bx-ax)*(fb-fc);			
			double q=(bx-cx)*(fb-fa);			
			u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*(q-r)); // This is the parabolic approximation, u is the minima of the 
												  // parabola with goes through a,b and c.
			double ulim=bx+GLIMIT*(cx-bx);
			if( (bx-u)*(u-cx) > 0. ) // u is between b and c
			{
				fu=function1dim(u,X,CONJ_D,f,ndim);		
				if( fu < fc ) // b u c brackets the minima 
				{
					ax=bx;
					fa=fb;
					bx=u;
					fb=fu;
					return ;
				}
				else if ( fu > fb ) // a b u brackets the minima 
				{
					cx=u;
					fc=fu;
					return;
				}
				
				u=cx+GOLD*(cx-bx);			
				fu=function1dim(u,X,CONJ_D,f,ndim);		
			}
			else if ( (cx-u)*(u-ulim) > 0. ) // parabolif fit is yadfaya
			{
				fu=function1dim(u,X,CONJ_D,f,ndim);		
				if(fu < fc ) // I do not know what this means ! Now we consider a, c and u as the triplet of points
				{
					bx=cx;	
					cx=u;
					u=cx+GOLD*(cx-bx);
					fb=fc;
					fc=fu;
					fu=function1dim(u,X,CONJ_D,f,ndim);
				}
			}
			else if ( (u-ulim)*(ulim-cx) > 0. ) // Then we set u to ulim and call it a day ?
			{
				u=ulim;
				fu=function1dim(u,X,CONJ_D,f,ndim);
			}
			else 
			{
				u=cx+GOLD*(cx-bx);
				fu=function1dim(u,X,CONJ_D,f,ndim);
			}
			ax=bx;
			bx=cx;
			cx=u;
			fa=fb;
			fb=fc;
			fc=fu;
		}	
		else
		{
			//cout<<"has this been happening >>>>???\n";
			return ;
		}
	}	

}

double brent(double ax, double bx, double cx, double fa, double fb, double fc,double &fmin, double X[], double CONJ_D[], double (*f)(double *, int), int ndim)
{
/* Given three variables with ax,bx and cx with f(bx) < f(ax) and f(bx) < f(cx) 
 * we isolate the minimum to some precision. 
 */
	double a,b,c;
	double v,w,x;	
	
	double xm;

	double fv,fw,fx;	

	const double ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
	const double CGOLD=0.3819660;
	const double tol=3.0e-14;
	double tol2=0.;

	double e=0.,d=0.,etemp;
	
	double p,q,r;

	double xmin;
	//double fmin;

	double u,fu;

	if(ax<cx)
	{
		a=ax;
		b=cx;
	}
	if(ax>cx)
	{
		a=cx;
		b=ax;
	}
	
	x=w=v=bx;
	fx=fw=fv=function1dim(x,X,CONJ_D,f,ndim);

	for(int i=0; i<100; i++)
	{
		//cout<<i<<"\n";
		xm=0.5*(a+b);	
		double tol1=tol*fabsl(x)+ZEPS;
		tol2=2.0*(tol1);
		//cout<<i<<"\t"<<b-a<<"\t"<<e<<"\n";;
		if (fabsl(x-xm) <= (tol2-0.5*(b-a)))
		{
			fmin=fx;
			return xmin=x;
		}
		if(fabsl(e) > tol1) 
		{
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if(q>0.) p=-p;
			etemp=e;
			e=d;
			
			if(fabsl(p) >= fabsl(0.5*q*etemp) || p <= q*(a-x) || p>= q*(b-x))
			{
				if( x > xm )
				{
					e=a-x;
					d=CGOLD*e;
				}
				else 
				{
					e=b-x;
					d=CGOLD*e;
				}
			}
			else
			{
				d=p/q;
				u=x+d;
				if ( u-a < tol2 || b-u < tol2 )
				{
					if( xm-x < 0. )
						d=-1.*tol1;
					else 
						d=tol1;
				}
			}
		}
		else
		{
			if( x >= xm )
			{
				e=a-x;
				d=CGOLD*e;
			}
			else 
			{
				e=b-x;
				d=CGOLD*e;
			}
		}
		if(fabsl(d) >= tol1)
		{
			u=x+d;
		}
		else 
		{
			if( d < 0. )
				d=-1.*tol1;
			else 
				d=tol1;
			u=x+d;
				
		}
		fu=function1dim(u,X,CONJ_D,f,ndim);
		
		if( fu <= fx ) 
		{
			if( u >= x ) a=x; else b=x;	
			v=w;
			w=x;
			x=u;

			fv=fw;
			fw=fx;
			fx=fu;

		}
		else 
		{
			if( u<=x ) a=u; else b=u;
			if(fu  <= fw || w==x )
			{
				w=u;
				v=w;
				fw=fu;
				fv=fw;
			}
			else if ( fu <= fv || v==x || v==w )
			{
				v=u;
				fv=fu;
			}
		}
	}
}

double dbrent(double ax, double bx, double cx, double fa, double fb, double fc,double &fmin, double X[], double CONJ_D[], double (*f)(double *, int), double (*fdf)(double*, double*,int), int ndim)
{
/* Given three variables with ax,bx and cx with f(bx) < f(ax) and f(bx) < f(cx) 
 * we isolate the minimum to some precision. 
 */
	double a,b,c;
	double v,w,x;	
	
	double d1,d2;
	double u1,u2;
	double ok1,ok2;

	double xm;

	double fv,fw,fx;	
	double dv,dw,dx;	

	const double ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
	const double CGOLD=0.3819660;
	const double tol=3.0e-14;
	double tol2=0.;

	double e=0.,d=0.,etemp,olde;
	
	double p,q,r;

	double xmin;
	//double fmin;

	double u,fu,du;

	if(ax<cx)
	{
		a=ax;
		b=cx;
	}
	if(ax>cx)
	{
		a=cx;
		b=ax;
	}
	
	x=w=v=bx;
	//fx=fw=fv=function1dim(x,X,CONJ_D,f,ndim);
	//dx=dw=dv=dfunction1dim(x,X,CONJ_D,fdf,ndim);

	dx=fdfunction1dim(x,X,CONJ_D,fdf,ndim,fx);
	dw=dv=dx;
	fw=fv=fx;

	for(int i=0; i<100; i++)
	{
		//cout<<i<<"\n";
		xm=0.5*(a+b);	
		double tol1=tol*fabsl(x)+ZEPS;
		tol2=2.0*(tol1);
		//cout<<i<<"\t"<<b-a<<"\t"<<e<<"\n";;
		if (fabsl(x-xm) <= (tol2-0.5*(b-a)))
		{
			fmin=fx;
			return xmin=x;
		}
		if(fabsl(e) > tol1) 
		{
			d1=2.*(b-a);
			d2=d1;
			if(dw != dx) d1=(w-x)*dx/(dx-dw); // secant method
			if(dv != dx) d2=(v-x)*dx/(dx-dv);
			
			u1=x+d1;
			u2=x+d2;
			
			ok1= (a-u1)*(u1-b) > 0. && dx*d1 <= 0. ;
			ok2= (a-u2)*(u2-b) > 0. && dx*d2 <= 0. ;

			olde=e;
			e=d;

			if( ok1 || ok2 )
			{
				if( ok1 && ok2 )
				{
					if( fabsl(d1) < fabsl(d2) )
					{
						d=d1;
					}
					else
					{
						d=d2;
					}
							
				}
				else if (ok1)
				{
					d=d1;
				}
				else if (ok2)
				{
					d=d2;
				}
				if( fabsl(d) <= fabsl( 0.5*olde) ) 
				{
					u=x+d;
					if ( u-a < tol2 || b-u < tol2 )
					{
						if( xm-x < 0. )
							d=-1.*tol1;
						else 
							d=tol1;
					}
				}
				else
				{
					if ( dx >= 0. )
					{
						e=a-x;
					}
					else 
					{
						e=b-x;
					}
					d=0.5*e;
				}
					
			}
			else
			{
				if ( dx >= 0. )
				{
					e=a-x;
				}
				else 
				{
					e=b-x;
				}
				d=0.5*e;
			}
		}
		else 
		{
			if ( dx >= 0. )
			{
				e=a-x;
			}
			else 
			{
				e=b-x;
			}
			d=0.5*e;
		}
		if(fabsl(d) >= tol1)
		{
			u=x+d;
			//fu=function1dim(u,X,CONJ_D,f,ndim);
			du=fdfunction1dim(u,X,CONJ_D,fdf,ndim,fu);
			
		}
		else 
		{
			if( d < 0. )
				d=-1.*tol1;
			else 
				d=tol1;
			u=x+d;
			//fu=function1dim(u,X,CONJ_D,f,ndim);
			du=fdfunction1dim(u,X,CONJ_D,fdf,ndim,fu);
			if( fu > fx )
			{
				xmin=x;
				fmin=fx;
				return xmin;
			}
		}
		//du=dfunction1dim(u,X,CONJ_D,fdf,ndim);
		if( fu <= fx ) 
		{
			if( u >= x ) a=x; else b=x;	
			v=w;
			w=x;
			x=u;

			fv=fw;
			fw=fx;
			fx=fu;

			dv=dw;
			dw=dx;
			dx=du;

		}
		else 
		{
			if( u<=x ) a=u; else b=u;
			if(fu  <= fw || w==x )
			{
				w=u;
				v=w;
				fw=fu;
				fv=fw;
				dw=du;
				dv=dw;
			}
			else if ( fu <= fv || v==x || v==w )
			{
				v=u;
				fv=fu;
				dv=du;
			}
		}
	}
}
double cubic_approx(double &fmin, double x0, double x1, double f0, double f1, double g0, double g1)
{
	double h = x1 - x0 ;	
	double F = f1 - f0 ;
	double G = (g1 - g0)*h ;
	double c = G - 2.*(F - g0*h);

	double gamma1 = -2.*g0*h/(( G- 3.*c)  + sqrt((G-3.*c)*(G-3.*c)-12.*c*g0*h));
	double gamma2 = -2.*g0*h/(( G- 3.*c)  - sqrt((G-3.*c)*(G-3.*c)-12.*c*g0*h));

	double b1 = (G - 3.*c*(1-2.*gamma1))/2.;
	double b2 = (G - 3.*c*(1-2.*gamma2))/2.;
	
	//cout<<gamma1*h<<"\t"<<gamma2*h<<"\t"<<2.*b1/(h*h)<<"\t"<<2.*b2/(h*h)<<"\n";	
	if ( 2.*b1/h*h > 0. )
	{
		
		double alpha = h*gamma1;
		double xmin = alpha+x0;
	//double b = (G - 3.*c*(1-2.*gamma))/2.;

		return xmin;
	}
	else 
	{
		double alpha = h*gamma2;
		double xmin = alpha+x0;
	//double b = (G - 3.*c*(1-2.*gamma))/2.;

		return xmin;

	}

}
void line_search(double &fmin, double X[],double CONJ_D[], double (*fdf)(double *,double *, int), int ndim, double f0, double* GRAD)
{
	//cout<<"in_linesearch\n";
	clock_t begin=clock();
	double ax,cx,fa,fb,fc;
	double ga,gb,gc;
	double xmin=0.;
	//double fmin=0.;
    double u=0.0;
  //cout<<X[0]<<"\t"<<CONJ_D[0]<<"\n";;
	cout<<std::setprecision(5);
  //for(int i=0;i < 100; i++)
  //{
  //	double tmp=fdfunction1dim(u,X,CONJ_D,fdf,ndim,fc);
  //	cout<<u<<"\t"<<fc/2000<<"\n";
  //	u=u+1./100000;
  //}
	//cout<<X[2*8]<<"\t"<<CONJ_D[2*8]<<"\n";;
	static double bx=0.2;
	ax=0.;
	cx=0.;
	//fa = function1dim(ax,X,CONJ_D,f,ndim);
    fa = f0;
////
    double der=0.;
    for(int i=0; i<ndim; i++)
    {
    	der=der+GRAD[i]*CONJ_D[i];
    	//TMP[i]=X[i]+alpha*CONJ_D[i];
    }

    ga = der;
	gb = fdfunction1dim(bx,X,CONJ_D,fdf,ndim,fb);
	//gb = dfunction1dim(bx,X,CONJ_D, df, ndim);
	//cout<<ax<<"\t"<<bx<<"\t"<<fa/2000<<"\t"<<fb/2000<<"\t"<<ga/2000<<"\t"<<gb/2000<<"\n";
	cout<<std::setprecision(16);
	double gsol,xsol,fsol,fprev;
	int index=0;
	while(1)
	{
		//cout<<index+1<<"\t"<<ax<<"\t"<<bx<<"\t"<<fa/2000<<"\t"<<fb/2000<<"\t"<<ga/2000<<"\t"<<gb/2000<<"\t"<<fprev/2000<<"\n";
		if(fabsl(gb/2000) < 1.0e-16)
		{
			fmin=fb;
			xmin=bx;
			break;
		}
		if ( ga*gb < 0. )
		{
			//cout<<"here\n";
			xsol=cubic_approx(fmin,ax,bx,fa,fb,ga,gb);
			gsol = fdfunction1dim(xsol,X,CONJ_D,fdf,ndim,fsol);
			//cout<<xsol<<"\t"<<fsol<<" which fuc\n";
			//return 0.;
			//cout<<fabsl((fsol - fprev )/fsol)<<"\n";
			if( fabsl((fsol - fprev )/fsol) < 1.0e-12)
			{
				//cout<<"anytime here\t"<<fprev/2000<<"\t"<<fsol/2000<<"\n";;
				fmin=fsol;
				xmin=xsol;
				break;
				//return xmin;
			}
			if(gsol*ga < 0.)
			{
				gb=gsol;
				bx=xsol;
				fb=fsol;
			}
			else 
			{
				ga=gsol;
				ax=xsol;
				fa=fsol;
			}
			if( fb < fa )
			{
				fprev=fb;	
			}
			else
			{
				fprev=fa;
			}
		}
		else 
		{
			//cout<<gb/2000<<"\t"<<ga/2000<<"\n";
			double correc= (gb-ga)/(bx-ax);	
			double adx = fabsl(bx-ax);
			//cout<<correc<<"\t"<<adx<<"\n";
			if(correc>0.)
			{
				xsol = ax - ga/correc;
				//cout<<xsol<<"\t"<<ga<<"\n";
				if( fabsl( xsol - ax ) > 3.*adx)
				{
					xsol = ax + 3.*adx;
				}
				else 
				{
					xsol = ax + 1.5*adx;
				}
				fsol=0.;
				gsol=fdfunction1dim(xsol,X,CONJ_D,fdf,ndim,fsol);
			}
			if( fabsl((fsol - fprev )/fsol) < 1.0e-12)
			{
				//cout<<"anytime here\t"<<fprev/2000<<"\t"<<fsol/2000<<"\n";;
				fmin=fsol;
				xmin=xsol;
				break;
				//return xmin;
			}
			if(fb<fa)
			{
				fa=fb;
				ax=bx;
				ga=gb;
			}
			bx=xsol;
			fb=fsol;
			gb=gsol;
			if( fb < fa )
			{
				fprev=fb;	
			}
			else
			{
				fprev=fa;
			}
		}
		index++;
	}
	//cout<<index<<"\n";
	bx=xmin;
	//cout<<xmin<<"\t"<<fmin/2000<<" cubic\n";
	//cout<<function1dim(bx,X,CONJ_D,f,ndim)/2000<<"\n";
	//cout<<function1dim(bx,X,CONJ_D,f,ndim)/2000<<"\n";
	//cout<<function1dim(bx,X,CONJ_D,f,ndim)/2000<<"\n";
	//cout<<function1dim(bx,X,CONJ_D,f,ndim)/2000<<"\n";
	//bracket(ax,bx,cx,fa,fb,fc,X,CONJ_D,f,ndim);	
	//cout<<fmin/2000<<"\t";
	//xmin=brent(ax,bx,cx,fa,fb,fc,fmin,X,CONJ_D,f,ndim);	
	//cout<<xmin<<"\t"<<fmin/2000<<" brent \n";
    //xmin=dbrent(ax,bx,cx,fa,fb,fc,fmin,X,CONJ_D,f,fdf,ndim);	
    //cout<<xmin<<"\t"<<fmin/2000<<" dbrent \n";
	clock_t end=clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	//cout<<time_spent<<"\n";
	cout<<std::setprecision(16);

	//cout<<xmin<<"\t"<<fmin/2000<<"\t"<<gsol/2000<<"\n";
	//return xmin;
	for(int i=0; i<ndim; i++)
	{
		X[i]=X[i]+xmin*CONJ_D[i];
		GRAD[i]=TMP_VEC[i];
		//cout<<i<<"\t"<<xmin*CONJ_D[i]<<"\n";
	}
	
	//cout<<"end_linesearch\n";
}
double minimize ( int ndim, double X[], int ITMAX, double (*fdf)(double *, double *, int), void (*ws)(double, double*, char*) )
{
	double *GRAD   = new (nothrow) double[ndim];
	double *CONJ_D = new (nothrow) double[ndim];
	double *G	    = new (nothrow) double[ndim];
	TMP	= new (nothrow) double[ndim];
	TMP_VEC = new (nothrow) double[ndim];

	
	double e=fdf(X,GRAD,ndim);
	cout<<e/2000<<"= initial energy \n";
	//return;

	//df(X,GRAD,ndim);

	for(int i=0; i<ndim; i++)
	{
		G[i]=-GRAD[i];	
		CONJ_D[i]=G[i];
	}
	
  	double fmin=e;
	for(int i=0; i<ITMAX; i++)
	{
		double fprev=fmin;
  		line_search(fmin,X,CONJ_D,fdf,ndim,fprev,GRAD);
		//cout<<fmin<<"\n";
		//double fnew=fmin;
		//fnew=fdf(X,GRAD,ndim);
  		double GG=0.;
  		double dGG=0.;
  		double XX=0.;
  		for(int i=0; i<ndim; i++)
  		{
  			dGG=dGG+(GRAD[i]+G[i])*GRAD[i];
  			GG=GG+G[i]*G[i];
  			XX=XX+GRAD[i]*GRAD[i];
  			//G[i]=GRAD[i];
  		}
		if(fabsl(XX/2000) < 1e-16)
		{
			cout<<i<<"\t"<<fmin/2000<<"\t"<<(fprev-fmin)/fprev<<"\t"<<XX/2000<<"\n";
			//ws(i,X);
			return fmin/2000;
			//break;
		}
	//  if((i+1)%100==0)
	//  {
	//  	cout<<i+1<<"\t"<<fmin/2000<<"\t"<<(fprev-fmin)/fprev<<"\t"<<XX/2000<<"\n";
	//  	//ws(i,X);
	//  }
  		
		//cout<<sqrt(XX)<<"\n";
		//cout<<dGG<<"\t"<<GG<<"\n";
  		double gamma=dGG/GG;
		//cout<<gamma<<"\n";
  
  		for(int i=0; i<ndim; i++)
  		{
  			G[i]=-GRAD[i];
  			CONJ_D[i]=G[i]+gamma*CONJ_D[i];
  		}	
	} 
}
