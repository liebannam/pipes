/**
 * \file channel.h
 * \brief channel.h file documentation
 * 
 * Contains:
 * class delcarations for Channel and derived Channel classes (definitions in channel.cpp)
 * class declarations for Junction1, Junction2, and Junction3 (definitions in channel.cpp)
 * class definitions for fRI and dfRI (for Riemann invariant calculation at boundary)
 * random function prototypes, horrifying macros, and other fun detritus
 **/


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


using std::max;
using std::min;

/*\def flag for printing miscellaneous debug info. Set to false to supress printout.*/
#define WTF false
/*\def flag for tracking chlorine (for water quality sim). Set to false to not bother*/
#define Clplease true

/*
 * \def horrifying macros for numerical flux routines
 * numFluxHLL requires you to define speeds. 
 * Current choices are:
 *	--speedsHLL
 *	--speedsROE
*/
#define numFlux(q1m, q1p, q2m, q2p, flux, Pm, Pp)  numFluxHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 
#define speeds(q1m, q1p, q2m, q2p, flux, Pm, Pp)   speedsHLL(q1m, q1p, q2m, q2p, flux, Pm, Pp) 


/**\brief Templating magic for streaming into vectors 
 * (Thanks to R. Saye for help with this)*/
template<typename T>
struct AppendToVector
{
    std::vector<T>& vec;
    AppendToVector(std::vector<T>& vec) : vec(vec) {}
};
/** streaming operator for use with AppendToVector*/
template<typename T>
std::istream& operator>>(std::istream& s, const AppendToVector<T>& app)
{
    T val;
	//check to see if this string can be interpretted as type T
    if (s >> val) 
        app.vec.push_back(val);
    return s;
}
/**call for AppendtoVector*/
template<typename T>
AppendToVector<T> appendTo(std::vector<T>& vec)
{
    return AppendToVector<T>(vec);
}


//////
//various prototype detritus
//////

/**height h as a function area A (Preissman slot geometry) */
double HofA   (double A, double D, double At, double Ts, bool P);
/**Reimann invariant thing \phi as a function of A (Preissman slot geometry)*/ 
double PhiofA (double A, double D, double At, double Ts, bool P);
/**Inverse of Reimann invariant thing \phi as a function of A (Preissman slot geometry)*/ 
double AofPhi (double phi, double D, double At, double Ts, bool P);
/**Gravity wave speed as a function of A (Preissman slot geometry)*/ 
double Cgrav  (double A, double D, double At, double Ts, bool P);
/**\eta=A pbar/rho, where pbar is the depth-averaged hydrostatic pressure and rho is the density(Preissman slot geometry)*/ 
double Eta    (double A, double D, double At, double Ts, bool P);

//I'm keeping these for the time being so I can compare other people's power series with ours
//inline double getTheta(double A, double D);
//inline double powPhi(double t, double D);
//



/**\brief Class for channel data **/
class Channel
{

	public:
	//Background info	
		/** Channel cross section type. Current choices are 0 (uniform) and 1 (Preissman slot)*/
		int channeltype;                           
	//Basic geometry
		/** Channel width, (m) (for Preissman slot, this is pipe diameter*/
		const double w;                            
		/** Channel length (m) */
		const double L;                            
		/** Number of grid points */
		const int N;                               
		/** Number of time steps (each channel instance must specify!)*/
		const int M;				   
		/** What time step we are currently at*/		
		int n;					      
		/** Grid spacing*/
		double dx;                                 
		/** Preissman parameters (won't get initialized or used for uniform cross section)*/
		double At, Af, a, Ts;		           
	//Slope and friction terms
		/** bed slope (-dz/L), where dz is the elevation difference between x=0 and x=L*/
		double S0;                                 
		/** Manning roughness coefficient*/
		double Mr;                                 
		/** Manning equation unit conversion coefficient (1.0 in MKS units)*/
		const double kn;                           
		/** Maximum wave speed (m/s)*/
		double cmax;				   
                /**Chlorine decay constant*/
                double KCl;
	//dynamic variables and boundary info
		/** Current state q*/
		double *q0;
		/** Initial state q0*/
		double *q;
		/** Temporary state (for RK time steppingi)*/	
	        double *qhat;          
		/** History of states q_hist laid out as [q(t=0), q(t=dt), ....q(t=dt*(M/Mi-1))] 
		 * where each q(t) is a row vector arranged as above.*/
		double *q_hist;           	 		  			
        /**concentration of chlorine we're tracking (not to be confused with c, a wavespeed)*/
		double *Cl, *Cl0, *Clhat;  
        /** history of chlorine concentration*/
        double *Cl_hist;        
        /** Indexing function to access q. obtain q(i,j)  where i =0,1 and j= 0,1...N-1*/
		int idx(int i_in, int j_in){return (N*i_in+j_in);}    //access q(i,j)  where i =0,1 and j= 0,1...N-1
		/** Indexing function to access qhist. obtain q^n(i,j) with i,j as above and n = 0,...M-1 (n*dt = t at which this slice is taken)*/
		int idx_t(int i_in, int j_in, int n_in){return (2*(N+2)*n_in+(N+2)*i_in+j_in);} 
		/** Pressurization states of each cell [p[leftend], p[0], p[1]....p[N-1], p[rightend]], size is N+2x1
		 *P =[P[0], P[1], ...P[N+1]], P[0] and P[N+1] are ghost cell values; if P(i) = 1, cell i is pressurized--details in Bourdarias 2007, pg 122*/
		vector<bool> P;   			    
		/** Indexing function pj to acceess elements of P correctly*/
		int pj(int i){return i+1;} 		     
		/** History of pressurization states*/
		bool *p_hist;				
		/**indexing for history of pressurization states, i = 0, [1,...N] N+1 */
		int pj_t(int i,int n){return (N+2)*n+i;}; 
		/** Current pressurization at cell under consideration-- this is a hilariously bad idea*/
		bool Pnow; 				        
		/** Left and right boundary fluxes*/
		double *bfluxleft, *bfluxright;         
	        /**Left and right boundary values for Cl*/
                double bCll, bClr;        
	//Various methods
		/** Constructor */
		Channel (int Nin, double win, double Lin, int Min, double a, double KClin=0.); 
		/** Destructor*/
		~Channel();
	   	/** Show values of q (1) or q0 (0)*/		
		void showVals(int Iwantq);                
	   	/** Show geometry paramters	*/
		virtual void showGeom() =0;		  
		/**set q0 to (A0,Q0)  = const*/
		void setq0(double A0, double Q0);          
		/**set q to (A0,Q0)  = const*/
		void setq(double A0, double Q0);           
		/** set q0 to nonconst value*/	
		void setq0(double *A0, double *Q0);	      
		/** Set q and q0 to nonconst value*/	
		void setq(vector<double>A0, vector<double>Q0);	   
		/** Take an Euler step */
		void stepEuler(double dt);
		/** Compute numerical HLL flux (requires the function speeds to be defined)*/
		void numFluxHLL(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp);      
		/** Physical flux (i.e. the actual conservation law flux)*/
		void physFlux(double q1, double q2, double *flux, bool P);
		/**HLL speeds--defined for given cross sectional geometry in derived class*/
		virtual void speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
		/**Roe speeds--defined for given cross sectional geometry in derived class*/
		virtual void speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)=0;
	//Functions specific to cross-sectional geometry, definitions in the derived class*/
		virtual void showp()=0;
		/** Pressure term in conservation law. eta= A \bar{p}/\rho where
		 * \bar{p} is average hydrostatic pressure and \rho is density*/ 
		virtual double Eta(double A, bool p) = 0;		         
		/** Height as a function of cross sectional area*/
		virtual double HofA(double A, bool p)=0;		 	 
		/** Average hydrostatic pressure reported as pressure head (units are m)*/
		virtual double pbar(double A,bool p) = 0;
		/** For ``negative Preissman'' model, height as a function of cross sectional area (not supported at present) */
		virtual double fakehofA(double A, bool p)=0;		 	 
		/** For ``negative Preissman'' model, cross sectional area as a function of height (not supported at present) */
		virtual double fakeAofh(double h, bool p) =0;
		/** Cross sectional area as a function of height*/
		virtual double AofH(double h, bool p)=0;			  	
		/** Gravity wavespeed c=sqrt(gA/l(A)) where l(A) is the width of the free surface corresponding to height A */
		virtual double Cgrav(double A, bool p)=0;                   		
		/** Hydraulic radius =(Area/Wetted Perimeter)*/	
		virtual double getHydRad(double A)=0;               		
		/** From Reimann invariant expression, \phi = \int_0^A (c(x)/x)dx*/
		virtual double PhiofA(double A, bool p) = 0;				
		/** Inverse of PhiofA, i.e. solve A = phi(x) for x*/
		virtual double AofPhi(double phi, bool p) =0;			
		/** Take an RK2 step to update with source terms*/
		void stepSourceTerms(double dt);
		/** Evaluate source terms for given A and Q */
		double getSourceTerms(double A, double Q); 
		/**update chlorine terms*/
                void stepTransportTerms(double dt);
                /** Get the volume of water in the pipe*/		
		double getTheGoddamnVolume();
		/** Get average gradient of pressure head h. Returns int_0^L(dh/dx) dx at time t_i**/
		double getAveGradH(int i);  
		/** Get total kinetic energy at time t_i*/
		double getKE(int i);
		/** Get total potential energy at time t_i*/
		double getPE(int i);
	//Functions to write output
	   	/** Write all information to file at intervals of dt*Mi/M */
		int writeqToFile(int Mi, double dt);   
		/** Write Riemann invariants to file*/
		int writeRItofile(double dt, int sign); 
	   	/** Output results at either locations x_i or time t_i to terminal. 
		 * x_i (or t_i) = where[i] for i = 1...length(where). which[i] =0/1 means it's a time/place*/
		void quickWrite(double *where, int *which, int K, double T, int skip);       
        /** Set chlorine decay constant*/
        void setKCl(double KCl_){KCl=KCl_;}
		/** min of three numbers*/
		double min3(double a, double b, double c);
		/** max of three numbers*/
		double max3(double a, double b, double c);
};



/** \brief Derived channel class for uniform cross section with width w */
class Cuniform: public Channel
{
	public:
		Cuniform(int Nin, double win, double Lin, int Min, double a):Channel(Nin, win, Lin, Min,a)
		{
			channeltype = 0;
		}	
		void showp();
		double Eta(double A, bool p){return G/(2.*w)*A*A;}
		double pbar(double A, bool p){return A/(2.*w);}
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


/** \brief Derived channel class for Preissmann slot geometry */
class Cpreiss: public Channel{
	public:
		/** Preissman parameters*/
		double D, yt, tt;  
		/** Initialize parameters and set up geometry*/
		void setGeom(double a_);	
		/** Constructor */
		Cpreiss(int Nin, double win, double Lin, int Min, double a =1200.):Channel(Nin, win, Lin, Min,a)
		{
			D = win;
			channeltype = 1;
			setGeom(a);
		}
		
		void showp();	
		double Eta(double A, bool p);
		double HofA(double A,bool p){return ::HofA(A, D, At, Ts, p);}
		double pbar(double A, bool p){return (A>0? Eta(A,p)/(G*A)  : 0.);}
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

/**\brief Class for applying boundary conditions at
 * external boundary of single pipe. */
class Junction1
{
	public:
		/** The channel this is the boundary of */
		Channel &ch0;
		/** Number of gridpoints in this channel */
		int N;					
		/** Channel width*/
		double w; 			    
		/** Type of externally specified quantity.(bvaltype =0 specifies A; bvaltype  =1  specifies Q)*/ 
		int bvaltype;			
		/** Instruction to use reflect or extrapolate: -/+1 means do (A,+/-Q) at boundary rather. 0 means use specified quantity.*/
		int reflect;        
		/** Values of externally specified quantity at each of the M+1 time steps. */
		double *bval;			
		/** Which end of pipe connects to this junction; 0 corresponds to left, 1 corresponds to right */
		int whichend;
        /** Boundary Chlorine value (initialized to zeros; only used if Clplease flag set to true*/
        double *Clbval;
		/** Constructor*/		
		Junction1(Channel &a_ch0, int a_which, double a_bval, int a_bvaltype);
		/**Set boundary value to constant*/
		void setbVal(double bvalnew);		
		/**Set boundary value to nonconstant in valarray<Real>*/
		void setbVal(valarray<Real> x);
		/**Set boundary value to nonconstant in array of doubles*/
		void setbVal(double *x);
		/**Set boundary value to nonconstant in vector of Reals*/
		void setbVal(vector<Real> x);
		/**Set boundary chlorine value time series to constant */
        void setClbval(double bClvalnew);
        /**Set boundary chlorine value time series to nonconstant*/
        void setClbval(double *Clbvalnew); 
        /** Destructor */
		~Junction1();
		/**\brief Apply numerical flux function to get boundary fluxes for ch0*/
		void boundaryFluxes();//	
		/** Compute total flow from left to right */	
		double getFlowThrough(double dt);
	};

/** \brief Class for applying boundary conditions at a junction 
 *  with two channels connected in serial*/
class Junction2   
{
	public:
		/** The channels connected here*/
		Channel &ch0, &ch1;
		/** Indices to find and store boundary values correctly depending on which end we're at*/
		int N0, N1, Ns0, Ns1; 
		/** Set to 1 for open, 0 for closed (default is open)*/
		double valveopen; 
		/** Time series of valve timings*/
		vector<double>valvetimes;
		/** Change in elevation as you go from pipe 1 to pipe 0. positive value means pipe 0 is above pipe 1. initialized to zero.*/
		double offset; 
		/** Which end of each channel we are at.*/
		int whichend0, whichend1;
		/** Constructor */
		Junction2(Channel &a_ch0, Channel &a_ch1, int a_which1, int a_which2, double a_valveopen);
		/** Set valve times to valarray of Reals*/
		void setValveTimes(valarray<Real>x);
		/** Set valve times to vector of Reals*/
		void setValveTimes(vector<Real>x);
		/**assign boundary fluxes to last cell of ch0 and first cell of ch1*/
		void boundaryFluxes(); 
};

/**\brief Class for applying boundary conditions at a junction where 3 channels come together and !&^# gets real*/
class Junction3
{
	public:
		/** The channels involved*/
		Channel &ch0, &ch1, &ch2;
		/** The junction2s between each pair*/
		Junction2 j2_01, j2_12, j2_21, j2_20;
		/** Number of cells and which end connects to the junction, for each pipe respectively*/
		int Ns[3], whichend[3];               
		/** Constructor*/
		Junction3(Channel &ch0_, Channel &ch1_, Channel &ch2_, int which1, int which2, int which3);	  
		//offsets = [off01, off12, off20] is offsets between each pair. there should be some compatibility here or it will be whack.
		/**assign all the bloody boundary fluxes according to some voodoo magic	*/
		void boundaryFluxes(); 

};



//for writing data to binary file 
//int binaryWrite(std::vector<double> x, int n, char* filename);
//

/**\brief Class for function we need to use to rootfind when dealing with uniform cross section Riemman invariants at boundaries*/
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

/**\brief Class for derivative of function fRI (for Newton sovler)*/
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

/**\brief Class for function used in setting up geometry info for Preissman slot*/
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

/**\brief Class for refining theta estimates with a few steps of Newton*/
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
/**
 *\brief Class for "all-purpose" function to evaluate things that
 * have to do with zeros of Riemann invariants
 * evaluates lhs - (cq*Qext/x +sign*phi(x)+cc*c(x)), where  c is wavespeed */
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
