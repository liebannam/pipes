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
#include "chebyshevlite.h"

using namespace std;

using std::max;
using std::min;

//
//WTF defined as true prints lots of miscellanious debug info...
#define WTF false
//#define WTF true

///////
//horrifying macros for numerical flux routines--previously allowed for an exact riemann solver that was slow and crappy)
//numFluxHLL requires you to define speeds-- 
//choices are speedsHLL and speedsROE (both dumb)
/////

#define numFlux(q1m, q1p, q2m, q2p, flux, Pm, Pp)  numFluxHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 

#define speeds(q1m, q1p, q2m, q2p, flux, Pm, Pp)   speedsHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 
//#define speeds speedsRoe


//Templating magic for streaming into vectors (THANKS ROB!!!!!)
template<typename T>
struct AppendToVector
{
    std::vector<T>& vec;
    AppendToVector(std::vector<T>& vec) : vec(vec) {}
};

template<typename T>
std::istream& operator>>(std::istream& s, const AppendToVector<T>& app)
{
    T val;
    if (s >> val)  //checks to see if this string can be interpretted as type T
        app.vec.push_back(val);
    return s;
}

template<typename T>
AppendToVector<T> appendTo(std::vector<T>& vec)
{
    return AppendToVector<T>(vec);
}


//////
//various prototype detritus
//////

//so that you can pass objective function derivative as a function easily...pathetically restrictive at present...
//typedef double(*drdA_t )(double Am, double A, double Ap, int j, double w, double dx, int N);


//inline double getTheta(double A, double D);
//inline double powPhi(double t, double D);
//
double HofA   (double A, double D, double At, double Ts, bool P);
double PhiofA (double A, double D, double At, double Ts, bool P);
double AofPhi (double phi, double D, double At, double Ts, bool P);
double Cgrav  (double A, double D, double At, double Ts, bool P);
double Eta    (double A, double D, double At, double Ts, bool P);

///
/////
///Class for channel data 
/////////
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
	//	const int Mi;				   // number of time steps between recording data in q_hist
		int n;					   // what time step we're at...		
		double dx;                                 // grid spacing
		double At, Af, a, Ts;		           // Preissman parameters (won't get initialized or used for uniform cross section)
/**** slope and friction terms*/
		double S0;                                 // bed slope (-dz/dx)
		double Mr;                                 // Manning Roughness coefficient
		double cmax;				   // maximum wave speed

/****dynamic variables and boundary info
* initial and final states q0 and q = [q(0,0), q(0,1)...,q(0,N-1), q(1,0),....q(1,N-1)] where q(0,:) is area A andq(1,:) is discharge Q
* history of states is q_hist laid out as [q(t=0), q(t=dt), ....q(t=dt*(M/Mi-1))] and each q is a row vector arranged as above.
* associated indexing functions so you don't have to recall how this is laid out
*** idx: //indexing function to access q(i,j)  where i =0,1 and j= 0,1...N-1
*** idx_t: indexing function to accessq^n(i,j)with i,j as above and n = 0,...M-1 (n*dt = t at which this slice is taken)
* pressurization states of each cell are in vector P 
*if P =[P[0], P[1], ...P[N+1]], P[0] and P[N+1] are ghost cell values; if P(i) = 1, cell i is pressurized--details in Bourdarias 2007, pg 122
*associate indexing function pj keeps track of this shift*/

		double *q0, *q, *qhat;          // previous, current, and temporary dynamical variables	
		double *q_hist;           		// history of dynamical variables	  		  			
		bool *p_hist;					//history of pressurization states
		vector<bool> P;   			    // pressurization states = [p[leftend], p[0], p[1]....p[N-1], p[rightend]], size is N+2x1
		bool Pnow; 				        // ''current pressurization'' at cell under consideration-- this is a hilariously bad idea!

		
		//indexing functions		
		int idx(int i_in, int j_in){return (N*i_in+j_in);}    //access q(i,j)  where i =0,1 and j= 0,1...N-1
		int idx_t(int i_in, int j_in, int n_in){return (2*(N+2)*n_in+(N+2)*i_in+j_in);} //accessq^n(i,j)with i,j as above and n = 0,...M-1 (n*dt = t at which this slice is taken)

		int pj(int i){return i+1;} 		    // indexing for pressurization states vector 
		int pj_t(int i,int n){return (N+2)*n+i;}; //indexing for history of pressurization states, i = 0, [1,...N] N+1 
		double *bfluxleft, *bfluxright;         // left and right boundary fluxes
		
/*****Methods*/
		Channel (int Nin, double win, double Lin, int Min, double a); // constructor
		~Channel();
		void showVals(int Iwantq);                 // show values of q		
		virtual void showGeom() =0;		   // show geometry paramters	
		
		void setq0(double A0, double Q0);          // set q0 to (A0,Q0)  = const
		void setq(double A0, double Q0);           // set q0 to (A0,Q0)  = const
		void setq0(double *A0, double *Q0);	   // set q to nonconst value	
		void setq(vector<double>A0, vector<double>Q0);	   // set q and q0 to nonconst value	
		
		void stepEuler(double dt);                          // Take an Euler step
		void numFluxHLL(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp);      //HLL Flux (need to define speeds)
		

		void physFlux(double q1, double q2, double *flux, bool P);
		virtual void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
		virtual void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
/**Specify geometry of specific cross section in a derived class:*/
		virtual void showp()=0;
		virtual double pbar(double A, bool p) = 0;		         // pbar = average hydrostatic pressure term 
		virtual double HofA(double A, bool p)=0;		 	 //height as a function of cross sectional area
		virtual double fakehofA(double A, bool p)=0;		 	 //for ``negative Preissman'' model, height as a function of cross sectional area
		virtual double fakeAofh(double h, bool p) =0;
		virtual double AofH(double h, bool p)=0;			  	//crosssectional area as a function of height
		virtual double Cgrav(double A, bool p)=0;                   		//gravity wavespeed actually redudant, I think!?!?!?
		virtual double getHydRad(double A)=0;               		// Hydraulic radius =(Area/Wetted Perimeter)	
		virtual double PhiofA(double A, bool p) = 0;				//integral (c(x)/x)dx
		virtual double AofPhi(double phi, bool p) =0;				//solve A = phi(x) for x
//random detritus
		double min3(double a, double b, double c);
		double max3(double a, double b, double c);
//source terms
		void stepSourceTerms(double dt);
		double getSourceTerms(double A, double Q); // evaluate source terms for given A and Q
/** Other quantities of interest*/		
		double getTheGoddamnVolume();
		double getAveGradH(int i);  //returns int_0^L(dh/dx) dx at time ti (probably not accurate enough...) 
		double getKE(int i);
		double getPE(int i);
/**Write stuff */
		int writeqToFile(int Mi, double dt);    //write all information to file at intervals of dt*Mi/M 
		int writeRItofile(double dt, int sign); //Write Riemann Invariants to file
		void quickWrite(double *where, int *which, int K, double T, int skip);        //quickly output end results at location x = places[i] for i = 1...length(places)	
};



//uniform cross section with width w
class Cuniform: public Channel
{
	public:
		Cuniform(int Nin, double win, double Lin, int Min, double a):Channel(Nin, win, Lin, Min,a)
		{
			channeltype = 0;
		}	
		void showp();
		double pbar(double A, bool p){return G/(2.*w)*A*A;}
		double HofA(double A, bool p){return A/w;}
		double fakehofA(double A, bool p){return A/w;}
		double fakeAofh(double h, bool p){return h*w;}
		double AofH(double h, bool p){return h*w;}
		double Cgrav(double A, bool p){return sqrt(G*A/w);}			
		double getHydRad(double A){return A/w*2.+w;} 

		void showGeom();
		void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
		void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp){return;};
		void updateExactRS(double q1m, double q1p, double q2m, double q2p, double *qnew, bool Pl, bool Pr, bool Px){return;};
		double PhiofA(double A,bool p){return 2.*sqrt(G*A/w);}
		double AofPhi(double phi, bool p){return (phi/2.)*(phi/2.)*w/G;} //invert the phi thing that's such a pain in the ass...
};


//preissmann slot geometry
class Cpreiss: public Channel{
	public:
		
		double D, yt, tt;  //Preissman parameters
		void setGeom(double a_);	
		Cpreiss(int Nin, double win, double Lin, int Min, double a =1200.):Channel(Nin, win, Lin, Min,a)
		{
			D = win;
			channeltype = 1;
			setGeom(a);
		}
		
		void showp();	
		double pbar(double A, bool p);
		double HofA(double A,bool p){return ::HofA(A, D, At, Ts, p);}
	        double hofAold(double A);	
		double fakehofA(double A, bool p);
		double fakeAofh(double h, bool p);
		double AofH(double h, bool p);
	//	double thetaofA(double A);
		double Cgrav(double A, bool p){return ::Cgrav(A, w, At, Ts,p);}	
		double getHydRad(double A);
		void showGeom();
		double findOmega(double Astar, double Ak, bool Ps, bool Pk);
		double PhiofA(double A,bool p){return ::PhiofA(A, D, At, Ts, p);}
		double AofPhi(double phi,bool p){return ::AofPhi(phi, D, At, Ts, p);}
		void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
		void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp);
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
		double w; 			    //width of associated channel
		int bvaltype;			//type of externally specified quantity.(bvaltype =0 specifies A; bvaltype  =1  specifies Q) 
		double *bval;			//value of externally specified quantity 
		int whichend;			//which end of pipe connects to this junction; 0 corresponds to left, 1 corresponds to right
		Junction1(Channel &a_ch0, int a_which, double a_bval, int a_bvaltype);
		void setbVal(double bvalnew);		
		void setbVal(valarray<Real> x);
		void setbVal(double *x);
		void setbVal(vector<Real> x);
		~Junction1();
		void boundaryFluxes();//	
		double getFlowThrough(double dt);//compute total flow from left to right 	
		int reflect;        //reflect=-/+1 means do (A,+/-Q) at boundary rather than try RI thing
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
		vector<double>valvetimes;//time series of valve timings
		double offset; //change in elevation as you go from pipe 1 to pipe 0. positive value means pipe 0 is above pipe 1. initialized to zero.
		int whichend0, whichend1;
		Junction2(Channel &a_ch0, Channel &a_ch1, int a_which1, int a_which2, double a_valveopen);
		void setValveTimes(valarray<Real>x);
		void setValveTimes(vector<Real>x);
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
		int Ns[3], whichend[3];               // number of cells and which end connects to junction, for each pipe respectively
		Junction3(Channel &ch0_, Channel &ch1_, Channel &ch2_, int which1, int which2, int which3);
				  
		//offsets = [off01, off12, off20] is offsets between each pair. there should be some compatibility here or it will be whack.

		void boundaryFluxes(); // assign all the bloody boundary fluxes according to some voodoo magic	

};



//for writing data to binary file 
//int binaryWrite(std::vector<double> x, int n, char* filename);
//
/**function we need to find zero of to deal with uniform cross section Riemman invariants at boundaries*/
class fRI   
{
	public:
		double c1,c2,c3;
		fRI(double a_c1, double a_c2, double a_c3)
		{
			c1 = a_c1;
			c2 = a_c2;
			c3 = a_c3;
		}
		
		double operator()(double x) 
		{
			double y;
			if (x>0.){y = c1-c2*sqrt(x)-c3/x;}
			else { y = c1; }
			return y;
		}		
};

/**derivative of fRI (for Newton sovler)*/

class dfRI  
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




/**function for setting up geometry info for Preissman slot*/
class ftheta 
{
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

/**for refining theta estimates with a few steps of Newton*/
class dftheta 
{
public:
	double D;
	dftheta(double D_) {
		D=D_; 
	}
	~dftheta(){
	}
	double operator()(double x) { return D*D/8.*(1-cos(x)); }
};

template <typename T> int sgn(T val)
{
		return (T(0)<val)-(val<T(0));
}
/**all purpose function to evaluate things that have zeros of Riemann invariants
 *
 * evaluate lhs - (cq*Qext/x +sign*phi(x)+cc*c(x))  c is wavespeed */
class fallpurpose{       
public:
	double D, At, Ts, lhs, Q, cq, cc;
	bool P;
	int sign;
	int s;
	fallpurpose(double D_,double At_,double Ts_, double lhs_, double Q_, int sign_, double cq_, double cc_, bool P_):
		D(D_), At(At_),Ts(Ts_), lhs(lhs_), Q(Q_), cq(cq_), cc(cc_), P(P_),sign(sign_)
	{
	}
	~fallpurpose(){
	}
	double operator()(double x) 
	{
		return x*lhs-(cq*Q+x*(sign*PhiofA(x,D,At,Ts,P)+cc*Cgrav(x,D,At,Ts,P)));
	}


};


#endif
