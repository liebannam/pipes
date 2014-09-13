/***********
 *Class delcarations for Channel, Junction1, Junction2, and Junction3 (definitions in channel.cpp)
 *Class definitions for fRI and dfRI (for Riemann invariant calculation at boundary)
 *Random function prototypes
 *typedef for numerical flux function
 ***********/


#ifndef CHANNEL_H
#define CHANNEL_H


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include "globals.h"

#include <valarray>

#include "real_def.h"
#include "ridders.h"

#include "newton.h"
#include <vector>
using namespace std;

using std::max;
using std::min;
///////
//horrifying macros for easily changing numerical flux routines
//numFlux choices are numFluxHLL and numFluxExact (THEY'RE BOTH IDIOTIC SO DON'T WORRY TOO MUCH HERE)
//if you choose HLL you also need to define speeds-- choices are speedsHLL and speedsROE (also both dumb)
/////

#define numFlux(q1m, q1p, q2m, q2p, flux, Pm, Pp)  numFluxHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 
//#define numFlux(q1m, q1p, q2m, q2p, flux, Pm, Pp)  numFluxExact(q1m, q1p, q2m, q2p, flux, Pm, Pp) 
#define speeds(q1m, q1p, q2m, q2p, flux, Pm, Pp)   speedsHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 
//#define speeds speedsRoe

//////
//various prototype detritus
//////

//so that you can pass objective function derivative as a function easily...pathetically restrictive at present...
//typedef double(*drdA_t )(double Am, double A, double Ap, int j, double w, double dx, int N);


inline double gettheta(double A, double D);
inline double powphi(double t, double D);
 double ePhi  (double A, double D, double At, double Ts, bool P);
 double eCgrav(double A, double D, double At, double Ts, bool P);
 double eEta  (double A, double D, double At, double Ts, bool P);

/////
///Class for channel data 
/////////

//so that you can pass numflux as a function easily--screw it, too hard :(
//typedef void(Channel::*numFlux_t)(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp);
//typedef void (Channel::*speeds_t)(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);



class Channel
{

	public:
		int channeltype;                           //0 for uniform cross section, 1 for Preissman slot
		const double kn;                           // manning equation unit coeff (1.0 in MKS)
/****basic geometry*/
		const double w;                            // channel width
		const double L;                            // pipe length
		const int N;                               // number of grid points
		const int M;				   // number of time steps (each channel instance must specify!)
		int n;					   // what time step we're at...		
		double dx;                                 // grid spacing
		double At, Af, a, Ts;				   //Preissman parameters (won't get initialized or used for uniform cross section
/**** slope and friction terms*/
		double S0;                                 // bed slope (-dz/dx)
		double Mr;                                 // Manning Roughness coefficient
		double cmax;				   //maximum wave speed
/****dynamic variables and boundary info */
		double *q0, *q, *qhat;                     // initial and final states q0 and q = [q(0,0), q(0,1)...,q(0,N-1), q(1,0),....q(1,N-1)] where q(0,:) is area A and q(1,:) is discharge Q
		double *q_hist;             		   // history of states q0; laid out as [q(t=0), q(t=dt), ....q(t=dt*(M-1))] and each q is a row vector arranged as above.
		vector<bool> P;   			   //if P =[P[0], P[1], ...P[N+1]], P[0] and P[N+1] are ghost cell values; if P(i) = 1, cell i is pressurized--details in Bourdarias 2007, pg 122
		bool Pnow; 				   // this is a hilariously bad idea!
		int idx(int i_in, int j_in){return (N*i_in+j_in);}//indexing function to access q(i,j)  where i =0,1 and j= 0,1...N-1
		int idx_t(int i_in, int j_in, int n_in){return (2*(N+2)*n_in+(N+2)*i_in+j_in);}//indexing function to accessq^n(i,j)with i,j as above and n = 0,...M-1 (n*dt = t at which this slice is taken)
		int pj(int i){return i+1;} 		    //indexing for pressurization states vector
		double bfluxleft[2], bfluxright[2];         // left and right boundary fluxes
/*****Methods*/
		Channel (int Nin, double win, double Lin, int Min); // constructor
		~Channel();
		void showVals(int Iwantq);                 //show values of q		
		virtual void showGeom() =0;		   //show geometry paramters	
		
		void setq0(double A0, double Q0);          // set q0 to (A0,Q0)  = const
		void setq(double A0, double Q0);           // set q0 to (A0,Q0)  = const
		void setq0(double *A0, double *Q0);	   // set q to nonconst value	
		void setq(vector<double>A0, vector<double>Q0);	   // set q and q0 to nonconst value	
		
		void stepEuler(double dt);          // Take an Euler step
	//	numFlux_t numFlux;
	//	speeds_t speeds;
		void numFluxHLL(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp);      //HLL Flux (need to define speed as Roe or other)
		void numFluxExact(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp);    //Exact RS flux for Preissman slot-not working at present
		

		void physFlux(double q1, double q2, double *flux, bool P);
		virtual void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
		virtual void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
		virtual void updateExactRS(double q1m, double q1p, double q2m, double q2p, double *qnew, bool Pl, bool Pr, bool Px)=0;
//Specify geometry of specific cross section in a derived class:
		virtual void showp()=0;
		virtual double pbar(double A, bool P) = 0;		  		 //pbar = average hydrostatic pressure term 
		virtual double hofA(double A)=0;		 	 	 //height as a function of cross sectional area
		virtual double fakehofA(double A, bool P)=0;		 	 //for ``negative Preissman'' model, height as a function of cross sectional area
		virtual double fakeAofh(double h, bool P) =0;
		virtual double Aofh(double h)=0;			  	//crosssectional area as a function of height
		virtual double cgrav(double h)=0;                   		//gravity wavespeed actually redudant, I think!?!?!?
		virtual double getHydRad(double A)=0;               		// Hydraulic radius =(Area/Wetted Perimeter)	
		virtual double intPhi(double A) = 0;				//integral (c(x)/x)dx
		virtual double Aofphi(double phi) =0;				//solve A = phi(x) for x
//random detritus
		double min3(double a, double b, double c);
		double max3(double a, double b, double c);
//source terms
		void stepSourceTerms(double dt);
		double getSourceTerms(double A, double Q); // evaluate source terms for given A and Q
		
		double getTheGoddamnVolume();
		double getAveGradH(int i);  //returns int_0^L(dh/dx) dx at time ti (probably not accurate enough...) 

//file write
		int writeqToFile(int Mi, double dt);
		int writeRItofile(double dt, int sign);	
};



//uniform cross section with width w
class Cuniform: public Channel
{
	public:
		Cuniform(int Nin, double win, double Lin, int Min):Channel(Nin, win, Lin, Min)
		{
			channeltype = 0;
		}	
		void showp();
		double pbar(double A, bool P){return G/(2.*w)*A*A;}
		double hofA(double A){return A/w;}
		double fakehofA(double A, bool P){return A/w;}
		double fakeAofh(double h, bool P){return h*w;}
		double Aofh(double h){return h*w;}
		double cgrav(double h){return sqrt(G*h);}			
		double getHydRad(double A){return A/w*2.+w;} 

		void showGeom();
		void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
		void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp){return;};
		void updateExactRS(double q1m, double q1p, double q2m, double q2p, double *qnew, bool Pl, bool Pr, bool Px){return;};
		double intPhi(double A){return 2.*sqrt(G*A/w);}
		double Aofphi(double phi){return (phi/2.)*(phi/2.)*w/G;} //invert the phi thing that's such a pain in the ass...
};


extern double Ptheta(double A);//mumblemumle plzwork??
//preissmann slot geometry
class Cpreiss: public Channel{
	public:
		
		double D, yt, tt;  //Preissman parameters
		void setGeom();
	    		
		Cpreiss(int Nin, double win, double Lin, int Min):Channel(Nin, win, Lin, Min)
		{
			D = win;
			channeltype = 1;
			setGeom();
		}
		
		void showp();	
		double pbar(double A, bool P);
		double hofA(double A);
		double fakehofA(double A, bool P);
		double fakeAofh(double h, bool P);
		double Aofh(double h);
		double thetaofA(double A);
		double cgrav(double h);		
		double getHydRad(double A);
		void showGeom();
		double findOmega(double Astar, double Ak, bool Ps, bool Pk);
		double intPhi(double A);
		double Aofphi(double phi); //invert the phi thing that's such a pain in the ass...
		void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
		void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
		void updateExactRS(double q1m, double q1p, double q2m, double q2p, double *qnew, bool Pl, bool Pr, bool Px);
};
///////
//Junction class definitions 
///////////

//////
///Junction1: impose either Q or A according to external specification
//////

class Junction1
{
	public:
		Channel &ch0;
		int N;					
		double w; 			//width of associated channel
		int bvaltype;			//type of externally specified quantity.(bvaltype =0 specifies A; bvaltype  =1  specifies Q) 
		double *bval;			//value of externally specified quantity 
		int whichend;			//which end of pipe connects to this junction; 0 corresponds to left, 1 corresponds to right
	//	Junction1(Channel &ch0_, int which);
		Junction1(Channel &a_ch0, int a_which, double a_bval, int a_bvaltype);
		void setbVal(double bvalnew);		
		void setbVal(valarray<Real> x);
		void setbVal(double *x);
		void setbVal(vector<Real> x);
		~Junction1();
		void boundaryFluxes();//	
	//	int idx(int i_in, int j_in){return (N*i_in+j_in);}
		int reflect;//=1 means do (A,-Q) at boundary rather than try RI thing
};

///////
///Junction2: class for two channels connected in serial
//////
class Junction2   
{
	public:
		Channel &ch0, &ch1;
		int N0, N1, Ns0, Ns1; //to find and store boundary values correctly depending on which end we're at
		double valveopen; // set to 1 for open, 0 for closed (default is open)
		double offset; //change in elevation as you go from pipe 1 to pipe 0. positive value means pipe 0 is above pipe 1. initialized to zero.
		int whichend0, whichend1;
		Junction2(Channel &a_ch0, Channel &a_ch1, int a_which1, int a_which2, double a_valveopen);
		void boundaryFluxes(); // assign boundary fluxes to last cell of pipeleft and first cell of pipe left	
};

//////
///Junction 3: class for junction where 3 channels come together and !&^# gets real
////
class Junction3///time for excitement...
{
	public:
		Channel &ch0, &ch1, &ch2;
		Junction2 j2_01, j2_12, j2_21, j2_20;
		int Ns[3], whichend[3];                                                    // number of cells and which end connects to junction, for each pipe respectively
		Junction3(Channel &ch0_, Channel &ch1_, Channel &ch2_, int which1, int which2, int which3);
				  
		//offsets = [off01, off12, off20] is offsets between each pair. there should be some compatibility here or it will be whack.

		void boundaryFluxes(); // assign all the bloody boundary fluxes according to some voodoo magic	

};



//for writing data to binary file 
//int binaryWrite(std::vector<double> x, int n, char* filename);

//Function prototypes for available numerical fluxes (only HLL at the moment)
void fluxHLLuniformxc(double q1m, double q1p, double q2m, double q2p, double *flux, double w);
void fluxRoe(double q1m, double q1p, double q2m, double q2p, double *flux, double w);
class fRI   //function we need to find zero of to deal with Riemman Invariants at boundaries
{
	public:
		double c1,c2,c3;
		fRI(double a_c1, double a_c2, double a_c3)
		{
			c1 = a_c1;
			c2 = a_c2;
			c3 = a_c3;
		//	printf("solving %f- %f/x-%f *sqrt(x) =0\n", c1,c3,c2);
		}
		
		double operator()(double x) 
		{
			double y;
			if (x>0.){y = c1-c2*sqrt(x)-c3/x;}
			else { y = c1; }
			return y;
		}		
};


class dfRI  //derivative for Newton sovler
{
	public:
		double c2,c3;
		dfRI( double a_c2, double a_c3)
		{
			c2 = a_c2;
			c3 = a_c3;
		}
		double 	operator()(double x){ return -c2/(2.*sqrt(x))+c3/(x*x);}
};




//function for setting up geometry info for Preissman slot
class ftheta {
public:
	double A, D;
	ftheta(double a_,double D_) {
		A =a_; 
		D=D_; 
	}
	~ftheta(){
	}
	double operator()(double x) { return D*D/8.*(x-sin(x))-A; }
};

//for refining theta estimates with a few steps of Newton
class dftheta {
public:
	double D;
	dftheta(double D_) {
		D=D_; 
	}
	~dftheta(){
	}
	double operator()(double x) { return D*D/8.*(1-cos(x)); }
};


class fphimlhs {       //evaluate phi(x)-lhs 
public:
	double D, At, Ts, lhs;
	fphimlhs(double D_,double At_,double Ts_, double lhs_):D(D_), At(At_),Ts(Ts_), lhs(lhs_)
	{
	}
	~fphimlhs(){
	}

	double operator()(double x) 
	{
		double Phi;	
		double t;
		if(x<0){return lhs;}
		else if (x<=At)
		{
			t = gettheta(x,D);
			Phi = powphi(t,D);
		}
		else
		{
			t  = gettheta(At,D);
			Phi = powphi(t,D);		
			Phi += 2*sqrt(x/Ts)-2*sqrt(At/Ts);
		}

		//cout<<"x = "<<x<<"phi = "<<Phi<<endl;
		return -Phi+lhs;
	}


};




class flr{
public:
	double D, Ts, At, Ak;
	bool Pk, Px;
	flr(double D_, double Ts_, double At_, double Ak_, bool Pk_, bool Px_):D(D_), Ts(Ts_), At(At_), Ak(Ak_), Pk(Pk_), Px(Px_)
	{
	}
	~flr(){
	}
	double eps = 1e-8;
	double operator()(double x)
	{
		if(x<=Ak)
		{
		//	cout<<"Ak="<<Ak<<"c(Ak) = "<<eCgrav(Ak,D, At, Ts, Px)<<endl;
			return ePhi(x,D,At,Ts,Px)-ePhi(Ak,D,At,Ts,Pk);
		}
		else
		{
		//	printf("Ak = %.15f\nx =  %.15f, D = %f, At = %f, Ts = %f\n", Ak,x,D, At, Ts);
		//	printf("Etak= %.15f\n",eEta(Ak,D,At,Ts,Pk));
		//	printf("Etax= %.15f\n",eEta(x,D,At,Ts,Px));
		//	printf("Ak-x= %.15f\n",(Ak-x));

			if(Ak==0){return 0.;}
			else if(x<=Ak+eps)//if Ak-x is super small, the usual routine will fuck up
			{//so use Ak-eps small and taylor expand to evaluate Eta(Ak)-Eta(x) ~ dEta/dA(x)(x-Ak)
			//this is easy since c^2 = GdEta/dA, lol. 
				double c = eCgrav((x+Ak)/2.,D, At, Ts, Px);
				return c*(x-Ak)*sqrt(1/(x*Ak)); 
			//	printf("c = %f\n", c);	
			}
			else{
				return sqrt((eEta(Ak,D,At,Ts,Pk)-eEta(x,D,At,Ts,Px))*(Ak-x)/(Ak*x));
			}
		}
	}
};

class f_exactRS{

public:
	double D, At, Ts, ul, ur;
	flr fl, fr;
	f_exactRS(flr fl_, flr fr_, double ul_, double ur_):fl(fl_),fr(fr_),ul(ul_), ur(ur_)
	{
	}

	~f_exactRS(){
	}
	double operator() (double x)
	{
	//	cout<<"ul-ur="<<ul-ur<<endl;
		double f1,f2;
		f1 = fl(x);
		f2 = fr(x);
	//	cout<<"x= "<<x<<" fl+fr="<<f1+f2<<endl;
		return f1+f2+ur-ul;
	}
	

};
/*
class fshock{
	//find A2 to satisfy (Q1-Q2)^2/(A1-A2) = Q1^2/A1 + Q2^2/A2 +gI(A1)-I(A2))
	double D, At, Ts, lhs, Q1, Q2, A1;
	double operator()(double x)
	{
		double I1 = ;
		double I2 = ;
		return Q1*Q1/A1 +Q2*Q2/A2 
	}
}
*/

class fallpurpose{       //evaluate lhs - (cq*Qext/x +sign*phi(x)+cc*c(x))  c is wavespeed 
public:
	double D, At, Ts, lhs, Q, cq, cc;
	bool P;
	int sign;
	fallpurpose(double D_,double At_,double Ts_, double lhs_, double Q_, int sign_, double cq_, double cc_, bool P_):
		D(D_), At(At_),Ts(Ts_), lhs(lhs_), sign(sign_),Q(Q_), cq(cq_), cc(cc_), P(P_)
	{
	}
	~fallpurpose(){
	}

	double operator()(double x) 
	{
		/*
		double Phi,t,cgrav;
	//	cout<<"x = "<<x<<endl;
	//	if(x<=1e-8){return lhs;}
		if (x<=At &&(!P))
		{
			t = gettheta(x,D);
			Phi = powphi(t,D);
			cgrav = sqrt(G*x/(D*sin(t/2.)));	
		}
		else
		{
			t  = gettheta(At,D);
			//cout<<"At = "<<At<<"Ts = "<<Ts<<endl;
			Phi = powphi(t,D);		
			Phi += 2*sqrt(x/Ts)-2*sqrt(At/Ts);
			cgrav = sqrt(G*(x-At)/Ts);

		}
		*/
		//cout<<"x = "<<x<<"phi = "<<Phi<<"ans = "<<lhs-(Q/x+sign*Phi)<<endl;
		//cout<<"lhs = "<<lhs<<" cgrav = "<<cgrav<<"  Phi = "<<Phi<<" x ="<<x<<" ans ="<<cc*cgrav+sign*Phi<<endl;
	//	cout<<"x "<<x<<" return value "<<lhs-((x>0?cq*Q/x:0)+sign*Phi+cc*cgrav)<<endl;
	// the following line is wrong but what I had before (noted on 6/23)
	//return lhs-((x>0?cq*Q/x:0)+sign*Phi+cc*cgrav);
		return x*lhs-(cq*Q+x*(sign*ePhi(x,D,At,Ts,P)+cc*eCgrav(x,D,At,Ts,P)));
	}


};
/*
class fcplusphi {       //evaluate lhs - (-c(x) +/- phi(x)) 
public:
	double D, At, Ts, lhs;
	int sign;
	fcplusphi(double D_,double At_,double Ts_, double lhs_, int sign_):D(D_), At(At_),Ts(Ts_), lhs(lhs_), sign(sign_)
	{
	}
	~fcplusphi(){
	}
	double gettheta(double A)
	{
		int count;
		ftheta th(A,D);
		//double theta =::ridders(th, -.1, 2*PI+1, &count,1e-10, 1e-12);
		double theta = Ptheta(A/(D*D));
		//cout<<theta<<" "<<A<<endl;
		theta = (A>0.)? (theta -(theta-sin(theta)- 8.*A/(D*D))/(1.-cos(theta))):0.;
		theta = (A>0.)? (theta -(theta-sin(theta)- 8.*A/(D*D))/(1.-cos(theta))):0.;	
		//cout<<theta<<endl;
		return theta;
	}

	double operator()(double x) 
	{
		double Phi,c;	
		double t;
	//	cout<<"x = "<<x<<endl;
		if(x<0){return lhs;}
		else if (x<=At)
		{
			t = gettheta(x);
			Phi = sqrt(G*3.*D/8.)*(t-1/80.*pow(t,3)+19./448000.*pow(t,5)+1./10035200.*pow(t,7)+491./(27.*70647808000.)*pow(t,9));
			c = (x>0)? sqrt(G*x/(D*sin(t/2.))):0.;		}
		else
		{
			t  = gettheta(At);
			//cout<<"At = "<<At<<"Ts = "<<Ts<<endl;
			Phi = sqrt(G*3.*D/8.)*(t-1/80.*pow(t,3)+19./448000.*pow(t,5)+1./10035200.*pow(t,7)+491./(27.*70647808000.)*pow(t,9));
			Phi += 2*sqrt(x/Ts)-2*sqrt(At/Ts);
		//	cout<<"theta = "<<t<<" Phi = "<<Phi<<endl;
			c = sqrt(G*(x-At)/Ts);
		}
		
	//	cout<<"c = "<<c<<"phi = "<<Phi<<"ans = "<<lhs-sign*(c+Phi)<<endl;
		return (lhs-sign*(c+Phi));
	}


};
*/

#endif
