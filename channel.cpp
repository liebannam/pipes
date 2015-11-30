


/**
 * \file Channel.cpp
 * \brief Channel.cpp documentation
 * 
 * Definitions of functions for class Channel*
 * Definition of functions for classes Junction1, Junction2, and Junction3
 **/

/*This file is part of Pipes.

    Pipes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pipes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pipes.  If not, see <http://www.gnu.org/licenses/>.
*/



#include "channel.h"

/** Number of Chebyshev points*/
const int Ncheb = 20;

/** Chebyshev points on [-1,1]*/
const double xx []= {-1.0000000000000000,  -0.9876883405951378,  -0.9510565162951535,  -0.8910065241883679,  
		     -0.8090169943749475,  -0.7071067811865476,  -0.5877852522924732,  -0.4539904997395468,  
		     -0.3090169943749475,  -0.1564344650402309,  -0.0000000000000000,  0.1564344650402306,  
		      0.3090169943749473,  0.4539904997395467,  0.5877852522924730,  0.7071067811865475,  
		      0.8090169943749473,  0.8910065241883678,  0.9510565162951535,  0.9876883405951377, 1.0000000000000000, };

/** Chebyshev coefficients for h(A)*/
const double coeffs_h[]={2.412767243157768248135319100059554e-01, 2.490551009703833796234839617091749e-01, 8.588320472482467973567819078662384e-03,
 			 9.228793157998457405435172914627134e-04, 1.310660974448734367701266018026398e-04, 2.129558025374983921564763711749651e-05,
			 3.749094616467203989458105646535410e-06, 6.962795341462423230912527993610194e-07, 1.343547052451044239960227996180368e-07,
		  	 2.668125626637029475349687868096361e-08, 5.418626817526250261028102778989540e-09, 1.120393624965683396198531943326663e-09,
		  	 2.350950901068854991515379247450580e-10, 4.994029090156534158917749943719349e-11, 1.071960458760750008323825025497721e-11,   
		     	 2.321591726189518874185799075341815e-12, 5.067116805835826216315004584447364e-13, 1.113605910211862088384776819368985e-13,
		  	 2.467624478577450438100243851682381e-14, 5.743999679585741830689301264332563e-15, 1.221321604072796243118941605439226e-15};

/** Chebyshev coefficients for phi(A)  when 0<A<=pi/8  */
const double coeffs_p1[] = {2.6099703290809577, 2.6418910954954722, 0.0407198075367772, 0.0101962412156085, 0.0017585298273079,
       			0.0004510739954107, 0.0001148009501627, 0.0000322213950326, 0.0000093014159914, 0.0000027886790457, 
			0.0000008553527957, 0.0000002678904977, 0.0000000852610460, 0.0000000275113722, 0.0000000089810480, 
			0.0000000029602467, 0.0000000009847978, 0.0000000003306886, 0.0000000001120956, 0.0000000000423243, 0.0000000000127085 };
/** Chebyshev coefficiencts for phi(A) when pi/8<=A<=pi*/

const double coeffs_p2[]= {6.4670122383491355, -0.8635238603246120, -0.2591109614939733, -0.0292512185543841, -0.0076221992053693,
       			-0.0017389019022193, -0.0004503493365253, -0.0001208472298793, -0.0000332374280799, -0.0000093988438938, 
			-0.0000027068870730, -0.0000007850335833, -0.0000002345070328, -0.0000000682526355, -0.0000000214757085, 
			-0.0000000059227336, -0.0000000021646739, -0.0000000004279834, -0.0000000002981759, 0.0000000000305111, -0.0000000000597102};
/** Chebyshev coefficiencts for A(phi) when phi<phi(pi/8)*/
const double coeffs_a1[] = {0.1262575316112871, 0.1873031261944491, 0.0712845714324142, 0.0092391530562048, -0.0011662365223237,
       			-0.0001871946267119, -0.0000251074200750, -0.0000052509280877, -0.0000011451817876, -0.0000002739189119,
		       	-0.0000000680431439, -0.0000000175632951, -0.0000000046503898, -0.0000000012589959, -0.0000000003467156, 
			-0.0000000000969788, -0.0000000000276160, -0.0000000000079714, -0.0000000000024471, -0.0000000000003386, 0.0000000000001585};
//for A(phi) when phi(pi/8)<phi
const double coeffs_a2[] = {0.6277322274489641, -0.2023497520119139, -0.0391435613898155, 0.0060556625096610, 0.0004614027183724,
       			-0.0000553007320124, -0.0000015358044760, -0.0000001139988435, 0.0000000727565380, -0.0000000282843815,
		       	0.0000000124134323, -0.0000000058741562, 0.0000000029710635, -0.0000000015906890, 0.0000000008960590, 
			-0.0000000005293791, 0.0000000003289447, -0.0000000002183606, 0.0000000001557247, -0.0000000001192869, 0.0000000000532795};

/*powers for scaling x points before evaluation (determined via python-automated trial and error)*/
const double ph = 2./3.;
const double pA1 = 1./3.;
const double pA2 = 5./12.;
const double pPhi1 = 1.;
const double pPhi2 = 3./5.;

/**
 * Height as a function of cross sectional area
 * \param [in] A cross sectional area (m^2)
 * \param [in] D pipe diameter (m)
 * \param [in] At transition cross-sectional area (m^2) (above this, we're in Preissman slot)
 * \param [in] Ts preissman slot width (m)
 * \param [in] P pressurization (bool; true or false)
 **/
double HofA(double A, double D, double At, double Ts, bool P)
{
	double y; 
	if(A<1e-15){return 0.;}
//	if(A<=At &&(!P))  //below slot
	if(A<=At)  //below slot
	{
		A = A/(D*D);//normalize by full area;
		if (A<=PI/8.){
			double Ahat = 2*pow((A*8./PI),ph)-1.;
			y = D*ChebEval(&coeffs_h[0],Ncheb, Ahat);
		}
		else{
			double Ahat = 2*pow(8./PI*(PI/4.-A),ph)-1.;
			y = D*(1.- ChebEval(&coeffs_h[0],Ncheb, Ahat));

		}
	}
	else      //in Preissman Slot
	{
		y = D+(A-At)/Ts;
	}
	return y;
}


/** Phi(A) = \int_0^A (c(z)/z)dz, where c(z) is the gravity wavespeed as a function of A
 * \param [in] A  cross sectional area (m^2)
 * \param [in] D  pipe diameter (m)
 * \param [in] At transition cross-sectional area (m^2) (above this, we're in Preissman slot)
 * \param [in] Ts Preissman slot width (m)
 * \param [in] P  pressurization (bool; true or false)
 **/
double PhiofA(double A, double D, double At, double Ts, bool P)
{
	if(A<1e-15){return 0.;}
	double phi = 0.;
	if(A<=At)
	{	
		A = A/(D*D);
		if(A<PI/8.)
		{	
			double Ahat = 2*pow((A*8./PI),pA1)-1.;
			phi = D*ChebEval(&coeffs_p1[0],Ncheb, Ahat);
		}
		else
		{
			double Ahat = 2.*pow((PI/4.-A)*8/PI,pA2)-1;
			phi = D*ChebEval(&coeffs_p2[0], Ncheb, Ahat);
		}
	}
	else
	{
		phi = D*ChebEval(&coeffs_p2[0], Ncheb, -1.)+2.*sqrt(G)*(sqrt(A/Ts)-sqrt(At/Ts));

	}
	return phi;
}

/** Inverse of Phi(A)
 * \param [in] A  cross sectional area (m^2)
 * \param [in] D  pipe diameter (m)
 * \param [in] At transition cross-sectional area (m^2) (above this, we're in Preissman slot)
 * \param [in] Ts Preissman slot width (m)
 * \param [in] P  pressurization (bool; true or false)
 * */
double AofPhi(double phi, double D, double At, double Ts, bool P)
{
	if(phi<1e-15){return 0.;}
	double A = 0.;
	double phi1m = D*ChebEval(&coeffs_p2[0],Ncheb, 1.);
	double phi2m = D*ChebEval(&coeffs_p2[0],Ncheb, -1.);
	if(phi<=phi2m)
	{
		
		if(phi<phi1m)
		{	
			double phat = 2*pow(phi/phi1m,pPhi1)-1.;
			A = D*D*ChebEval(&coeffs_a1[0],Ncheb, phat);
		}
		else
		{
			double phat = 2.*pow((phi2m-phi)/(phi2m-phi1m),pPhi2)-1.;
			A = D*D*ChebEval(&coeffs_a2[0], Ncheb, phat);
		}
	}
	else
	{
		double duh = (phi-phi2m)/(2.*sqrt(G))+sqrt(At/Ts);
		A = (duh*duh*Ts);
	}
	return A;

}

/** 
 * Cgrav = sqrt(g*A/l) is the gravity wavespeed (l is the width of the free surface)
 * param [in] A  cross sectional area (m^2)
 * param [in] D  pipe diameter (m)
 * param [in] At transition cross-sectional area (m^2) (above this, we're in Preissman slot)
 * param [in] Ts Preissman slot width (m)
 * param [in] P  pressurization (bool; true or false)
 **/
double Cgrav(double A, double D, double At, double Ts, bool P)
{
	double c=0.;
	//check for near-zero area 
	if(A<1e-15){return 0.;}
	if(A<At)
	{
		double h = HofA(A,D,At,Ts,P);
		double l = 2.*sqrt(h/D*(1.-h/D));
		c = sqrt(G*A/l);
	}
	//in slot assign pressure wavespeed sqrt(G*Af/Ts = a) (see Sanders 2011 pg. 163)
	else  
	{
		c = sqrt(G*D*D*PI/(4.*Ts));
	} 
	return c;    
}

/**
 * \eta = g*\int_0^{h(A)}((h(A)-z)l(z)) dz  = ghA - g\bar{y}
 * where g =9.8 m/s^2 is acceleration due to gravity
 * \bar{y} is the centroid of the flow, \bar{y} = (\int_0^h(A) zl(z) dz)/A 
 * note that
 * \eta = A*\bar{p}/\rho
 *    where \bar{p} is the average hydrostatic pressure
 *    and \rho is the density 
 *
 * \param [in] A cross sectional area (m^2)
 * \param [in] D pipe diameter (m)
 * \param [in] At transition cross-sectional area (m^2) (above this, we're in Preissman slot)
 * \param [in] Ts preissman slot width (m)
 * \param [in] P pressurization (bool; true or false)   
 * */
double Eta(double A, double D, double At, double Ts, bool P)
{
	double Eta;
	if (A<At)
	{
		double y = HofA(A, D, At, Ts, P);
		Eta = G/12.*((3.*D*D-4.*D*y+4.*y*y)*sqrt(y*(D-y))
			-3.*D*D*(D-2.*y)*atan(sqrt(y)/sqrt(D-y)));
	}
	else //in slot use Sanders 2011 TPA approach
	{
		double H = (A-(PI*D*D/4.))/Ts;
		H = (A-At)/Ts;
		Eta = PI/4.*G*D*D*(H+D/2.);
	}
	return Eta;
}


/**
 * Channel is a class with dynamical variables (A, Q) in array q.
 * Has methods to update them according to de St. Venant Equations.
 * Instantiate a derived class to 
 * q is laid at as follows:
 *		q = [q(0,0), q(0,1)...,q(0,N-1), q(1,0),....q(1,N-1)], 
 * where q(0,:) is area A (m^2) and q(1,:) is discharge Q (m^3/s)
 * \param[in] Nin is the number of x gridpoints
 * \param[in] win is the width of the channel (m). Interpreted as pipe diameter for Preissman slot
 * \param[in] Lin is the length of the channel (m).
 * \param[in] Min is the number of time steps 
 * \param[in] a is the pressure wave speed (m/s) only relevant to Preissman slot
 **/
Channel::Channel(int Nin, double win, double Lin, int Min, double a, double kwin): kn(1.0),w(win),L(Lin),N(Nin),M(Min),kw(kwin) 
{
	q0 = new double[2*N];
	q = new double[2*N];
	qhat = new double[2*N];
	Cl = new double[N];
	Cl0 = new double[N];
    Clhat = new double[N];
    bfluxleft = new double[2];
	bfluxright = new double[2];
	dry = 1e-8*win;//tolerance for "dry'' pipe
   // printf("dry = %.16f\n",dry);
    //assume junctions aren't ventilated unless they're junction1s without reflection
	P.push_back(false);
	for(int i = 0; i<N+1; i++){P.push_back(false);}
	P.push_back(false);
	n = 0;
	//Array q_hist is allocated if there's less than 10 million stored variables, else complain and quit
	if(N*M<2e7){
		q_hist = new double[2*(N+2)*(M+2)]; 
		p_hist = new bool[2*(N+2)*(M+2)];
        Cl_hist = new double[(N+2)*(M+2)];
        for (int k =0; k<(N+2)*(M+2); k++)
        {
            Cl_hist[k] = 0.;
        } 
        for (int k =0; k<2*(N+2)*(M+2); k++)
        {
            p_hist[k] = 0.;
            q_hist[k] = 0.;
        }

	}
	else 
	{
		cout<<"You know, allocating "<< sizeof(w)*M*N*2<<" bytes may be a bad idea! Decrease time increment or figure out a better storage scheme.\n";
		throw("Yeah...nope.");
	}

	Mr = 0.;
	S0 = 0.;	
	int i,j;
	dx = L/((double)N);
	for(i=0; i<2*N; i++)
	{
		q0[i]=0.;
		q[i] =0.;
		qhat[i] = 0.;
	}
        for(i=0; i<N; i++)
        {
            Cl[i] = 0.;
            Cl0[i]=0.;
            Clhat[i]=0.;
        }
        bCll =0.;
        bClr=0.;
	for(j= 0; j<2;j++)
	{
		bfluxleft[j] = 0;
		bfluxright[j] = 0;
	}
    kb = 0.55/86400.;// (s^-1)=.55 days^{-1} (Rossman 1994)
}

/** Destructor */
Channel::~Channel()
{
	delete [] q0;
	delete [] q;
	delete [] q_hist;
	delete [] p_hist;
	delete [] qhat;
	delete [] Cl;
    delete [] Cl0;
    delete [] Clhat;
    delete [] Cl_hist;
    delete [] bfluxleft;
	delete [] bfluxright;
}

/** Display geometry info about uniform cross section channel */
void Cuniform::showGeom()
{
	printf("\nAbout this channel:\n"
			"Cross section is uniform\n"
			"|         |\n|<---w--->|\n|_________|\n"	
			"w = %2.1f %21s width (m)\nN = %4d %20s number of grid points\nL = %3.1f %20s length (m)\n", w,"",N,"",L,""
			"Mr = %1.3f %20s Manning roughness coefficient\nS0 = %1.2f %22s Bed slope\n", Mr, "", S0, "");
}

/** Display geometry info about Preissman slot cross section pipe */
void Cpreiss::showGeom()
{
	printf("About this channel:\n\n"
		"Cross section is Preissman Slot\n"
		"      Ts\n     <->\n     | |\n     | |      \n  /~~~~~~~\\    ^\n"
		" /         \\   |\n(<----D---->)  |  yt\n \\          /  |\n  \\________/   v\n"
		"D = %.1f %20s width (m) \nN = %d %20s number of grid points \nL = %.1f %20s length (m)\n", w,"",N,"",L,""
		"Slot width Ts= %1.5f\nTransition height yt= %1.5f%s\n",Ts,yt,""
		"Mr = %1.3f %20s Manning Roughness coeff\nS0 = %1.2f %20s bed slope\n", Mr, "", S0, "");
}


/** Display current cell values - argument is 0 for q0, 1 for q*/
void Channel::showVals(int Iwantq)
{
	int i;
	if(Iwantq)
	{
		for(i=0; i<N; i++)
		{
			printf("A0[%d] = %f and Q0[%d] = %f\n", i,q[idx(0,i)],i, q[idx(1,i)]); 
		}
	}
	else
	{
		for(i=0; i<N; i++)
		{
			printf("A[%d] = %f and Q[%d] = %f\n", i,q0[idx(0,i)],i, q0[idx(1,i)]); 
		}
	}	
}

/** Show pressure head */
void Cuniform::showp()
{

	int i;
	for(i=0; i<N; i++)
	{
		printf("h[%d] = %f and Q0[%d] = %f\n", i,q[idx(0,i)]/w,i, q[idx(1,i)]); 
	}
}

/** Show pressure head H=Eta/(g*A)*/
void Cpreiss::showp()
{

	int i;
	double h;
	for(i=0; i<N; i++)
	{
		h = fakehofA(q[idx(0,i)], P[pj(i)]);
		printf("pressure head[%d] = %f (m) = %f (psi) and Q0[%d] = %f, P = %s\n", i,h,m_to_psi*h,i, q[idx(1,i)],P[pj(i)]?"true":"false"); 
	}
}


/** Initialize q0 with constant data (A0,Q0)*/
void Channel::setq0(double A0, double Q0) 
{
	int i;
	bool p0 = (A0>=At);
	for(i=0; i<N; i++)
	{
		q0[idx(0,i)] = A0;
		q0[idx(1,i)] = Q0;
		q_hist[idx_t(0,i+1,0)] = A0;
		q_hist[idx_t(1,i+1,0)] = Q0;
		P[pj(i)] = p0;
		p_hist[pj_t(i+1,0)] = p0;
	}
	q_hist[idx_t(0,0,0)] = A0;
	q_hist[idx_t(0,N+1,0)] = A0;
	q_hist[idx_t(1,0,0)] = Q0;
	q_hist[idx_t(1,N+1,0)] = Q0;
//	P[0] = p0;
//	P[N+1] = p0;
	p_hist[pj_t(0,0)] = p0;
	p_hist[pj_t(N+1,0)] = p0;

}

/**Initialize with nonconstant array data A0, Q0*/
void Channel::setq0(double *A0, double *Q0)
{
	int i;
	for(i=0; i<N; i++)
	{
		q0[idx(0,i)] = A0[i];
		q0[idx(1,i)] = Q0[i];
		q[idx(0,i)] = A0[i];
		q[idx(1,i)] = Q0[i];
		q_hist[idx_t(0,i+1,0)] = A0[i];
		q_hist[idx_t(1,i+1,0)] = Q0[i];
		P[pj(i)] = A0[i]>=At;
		p_hist[pj_t(i+1,0)] = A0[i]>=At;
	}
	q_hist[idx_t(0,0,0)] = A0[0];
	q_hist[idx_t(0,N+1,0)] = A0[N-1];
	q_hist[idx_t(1,0,0)] = Q0[0];
	q_hist[idx_t(1,N+1,0)] = Q0[N-1];
//	P[0] = (A0[0]>=At);
//	P[N+1] = (A0[N-1]>=At);	
	p_hist[pj_t(0,0)] = (A0[0]>=At);
	p_hist[pj_t(N+1,0)] = (A0[N-1]>=At);
}

/** Initialize with nonconstant vector data A0, Q0*/
void Channel::setq(vector<double>A0, vector<double>Q0)
{
	int i;
	for(i=0; i<N; i++)
	{
		q0[idx(0,i)] = A0[i];
		q0[idx(1,i)] = Q0[i];
		q[idx(0,i)] = A0[i];
		q[idx(1,i)] = Q0[i];
		q_hist[idx_t(0,i+1,0)] = A0[i];
		q_hist[idx_t(1,i+1,0)] = Q0[i];
		P[pj(i)] = (A0[i]>=At);
		p_hist[pj_t(i+1,0)] = (A0[i]>=At);
	}
	q_hist[idx_t(0,0,0)] = A0[0];
	q_hist[idx_t(0,N+1,0)] = A0[N-1];
	q_hist[idx_t(1,0,0)] = Q0[0];
	q_hist[idx_t(1,N+1,0)] = Q0[N-1];
	P[0] = (A0[0]>=At);
	P[N+1] = (A0[N-1]>=At);	
	p_hist[pj_t(0,0)] = (A0[0]>=At);
	p_hist[pj_t(N+1,0)] = (A0[N-1]>=At);
}

/** Initialize q with constant data (A0,Q0)*/
void Channel::setq(double A0, double Q0)
{
	int i;
	bool p0 = (A0>=At);
	//cout<<"p0 = "<<p0<<" At = "<<At<<"A0 = "<<A0<<endl;
	for(i=0; i<N; i++)
	{
		q[idx(0,i)] = A0;
		q[idx(1,i)] = Q0;
		P[pj(i)] = p0;
		p_hist[pj_t(i,0)] = p0;

	}
	P[0] = p0;
	P[N+1] = p0;
	p_hist[pj_t(0,0)] = p0;
	p_hist[pj_t(N+1,0)] = p0;
}		


void Channel::setCl0(vector<double> Cl0_)
{
    for (int i = 0; i<N; i++)
    {
        double c0 = Cl0_[i];
        Cl[i] = c0;
        Cl0[i] = c0;
        Cl_hist[i+1] = c0;
       // printf(" i = %d  Cl=%f  Cl0 = %f\n",i,Cl[i], Cl0[i]);
    }
}

/** 
 * Take M Euler steps of length dt to update conservation law
 *  for left state [q1(i), q2(i)] and right state [q1p(i+1), q2(i+1)]
 *  using numflux as numerical flux
 **/
int Channel::stepEuler(double dt)
{
	cmax = 0.;  //reset CFL;
	double fplus[2] ={0};
	double fminus[2] ={0};
	double nu = dt/dx;				
	//decide when negative A forces the code throws an error over negative area (vs. just printing a warning) 
	double negtol = -dx/100.;		
	int i,k;
	//get pressurization information for leftmost cell
	Pnow = P[pj(0)]; 							              
	//update leftmost cell using externally assigned fluxes bfluxleft
	//printf("bfluxleft = [%f,%f]\n",bfluxleft[0], bfluxleft[1]);
//	printf("bfluxright = [%f,%f]\n",bfluxright[0], bfluxright[1]);
    numFlux(q0[idx(0,0)], q0[idx(0,1)], q0[idx(1,0)], q0[idx(1,1)],fplus,P[pj(0)], P[pj(1)]);
	q[idx(0,0)] = q0[idx(0,0)]-nu*(fplus[0]-bfluxleft[0]);           
	q[idx(1,0)] = q0[idx(1,0)]-nu*(fplus[1]-bfluxleft[1]);	
	P[pj(0)] = (q[idx(0,0)]>At )||( P[pj(-1)]==true && P[pj(1)]==true); 
	for(i = 1; i<N-1; i++)
	{
 		// - fluxes are the previous step's + fluxes
		fminus[0] = fplus[0];
		fminus[1] = fplus[1];                                                       
		Pnow = P[pj(i)];		
		// + fluxes need to be computed afresh
		numFlux(q0[idx(0,i)], q0[idx(0,i+1)], q0[idx(1,i)], q0[idx(1,i+1)],fplus, P[pj(i)], P[pj(i+1)]);                    
		// update conservation law
		for(k=0; k<2; k++)
		{
			q[idx(k,i)] = q0[idx(k,i)]-nu*(fplus[k]-fminus[k]);
		}
        if(q[idx(0,i)]<0)
        {
            /*
            printf("Negative area! i=%d, A= %f, Q= %f, fplus= [%f,%f] and fminus = [%f %f]\n",i,q[i], q[idx(1,i)], fplus[0],fplus[1], fminus[0],fminus[1]);
			printf("%d    %.16f   %.16f\n",-1, q_hist[idx_t(0,0,n)], q_hist[idx_t(1,0,n)]);  
            for (int i=0; i<N; i++)
            {
                printf("%d    %.16f   %.16f  %.16f   %.16f\n",i, q0[idx(0,i)], q0[idx(1,i)],q[idx(0,i)], q[idx(1,i)]);  
            }
			printf("%d    %.16f   %.16f\n",N, q_hist[idx_t(0,N+1,n)], q_hist[idx_t(1,N+1,n)]);  
            */
            return 1;
            //std::runtime_error("Oh damn. Negative area!\n");

            //
            //q[i]=0;
        }        
		//update pressurization info: if A>At-> pressurized pipe; if A<At AND a neighbor is free surface, pipe is free surface
		P[pj(i)] = (q[idx(0,i)]>At )||( P[pj(i-1)]==true && P[pj(i+1)]==true);
	}
	// update rightmost cell using externally assigned fluxes bfluxleft
	P[pj(N-1)] = (q[idx(0,N-1)]>At )||( P[pj(N-2)]==true && P[pj(N)]==true); 
	q[idx(0, N-1)] = q0[idx(0,N-1)] - nu*(bfluxright[0]-fplus[0]);                      
	q[idx(1,N-1)] = q0[idx(1,N-1)] - nu*(bfluxright[1] - fplus[1]);	
    //fold in source terms
	stepSourceTerms(dt);
        //step chlorine transport terms if tracking water quality
        if (Clplease) stepTransportTerms(dt);
        // set q0 to updated value
	for(i =0;i<N;i ++)
	{
		                                      
		//check for negative areas --do I really need this?
		if(q[idx(0,i)]<0)
		{
			/*
            printf("!!Negative area!!!\n with a[%d] = %f at time %f\n ", i, q0[i], dt*(double)n);
			printf("!!TOO MUCH Negative area!!! SHITSHITSHIT\n with a[%d] = %f at time %f\n ", i, q0[i], dt*(double)n);
			printf("%d    %.16f   %.16f\n",-1, q_hist[idx_t(0,0,n)], q_hist[idx_t(1,0,n)]);  
            for (int i=0; i<N; i++)
            {
            printf("%d    %.16f   %.16f\n",i, q0[idx(0,i)], q0[idx(1,i)]);  
            }
			printf("%d    %.16f   %.16f\n",-1, q_hist[idx_t(0,N+1,n)], q_hist[idx_t(1,N+1,n)]);  */
            return 2;
            //throw std::runtime_error("Oh damn. Negative area!\n");
		}
        for(k = 0; k<2; k++)
		{
			q0[idx(k,i)] = q[idx(k,i)];
		}     
//		if (q0[idx(0,i)]<negtol)/
//		{
//			q0[idx(0,i)] = 0.0;
//		}
	if (WTF)printf("cmax =%f and CFL=%f",cmax, dt/dx*cmax);
    }
    return 0;
}

/** Physical flux for Preissman slot*/
void physFluxBetter(double A, double Q, bool P, double D, double At, double Ts, double *flux)
{
	flux[0] = Q;
	flux[1] = (A>0? Q*Q/A:0.) +Eta(A, D, At, Ts, P); 
}

/**
 * Harten-Lax-van Leer (HLL) numerical flux
	 * (what's going on here:
	 *     sl   :     sr                 sl  sr   :
	 *      \   u*   /                    ` u*\   :
	 *       \  :   /                       `  \  :
	 *    um  \ :  / up                  um   ` \ :up
	 * ________\:/_______  or maybe    _________`\:_________ or even...
	 *          :                                 : 
	 * um, up are left and right states, u* is center state 
	 * sl, sr denote shocks )
**/ 
void Channel::numFluxHLL(double q1m, double q1p, double q2m, double q2p, double *flux, bool Pm, bool Pp)                
{
    
	double s[2]={0,0};
	double slow = 1e-5;
	//estimate speeds (Roe or fancier)
	speeds(q1m, q1p, q2m, q2p, s, Pm, Pp);
	// check for near zero speeds
	if (fabs(s[0])< slow && fabs(s[1])<slow) 
	{                                                       
		flux[0] = 0;
		flux[1] = 0;
	}
	// Update with HLL flux
	else 
	{                                                                                           
    if(s[0]>=0){physFlux(q1m,q2m,flux, Pp);}
    else if(s[0]<0 && s[1]>0)
	{
		double Fl[2];
		double Fr[2];
		physFlux(q1m,q2m,Fl, Pm);
		physFlux(q1p,q2p,Fr, Pp);
		flux[0] = (s[1]*Fl[0]-s[0]*Fr[0]+s[0]*s[1]*(q1p-q1m))/(s[1]-s[0]);
		flux[1] = (s[1]*Fl[1]-s[0]*Fr[1]+s[0]*s[1]*(q2p-q2m))/(s[1]-s[0]);
		//cout<<"fluxes are "<<flux[0]<<" "<<flux[1]<<Fl[1]<<" "<<Fr[1]<<endl;
	}
        else if(s[1]<=0)
        	physFlux(q1p, q2p, flux, Pp);
        else
          	printf("Error! Check your speeds. Something is wrong! s = [%f,%f]\n",s[0], s[1]);
	}
    return;
}

/**de St. Venant flux function (See Leon 2009)
 *This is valid for either cuniform (uniform width) or cpreiss (Preissman slot) channel classes
 */
void Channel::physFlux(double q1, double q2, double *flux, bool P)                  
{
    flux[0] = q2;
    flux[1] = (q1>0? q2*q2/q1:0.) +Eta(q1,P);  
}

/** HLL speeds (uniform channel)*/
void Cuniform::speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp)                           
{
    double um, up, hm, hp, uhat, hhat, smin, smax;
    int j;
    bool p = Pm && Pp;
    um = (q1m>dry? q2m/q1m : 0.);
    up = (q1p>dry? q2p/q1p : 0.);
    hm = HofA(q1m,Pm);
    hp = HofA(q1p,Pp);
    uhat = (hm+hp >0.)? (sqrt(hm)*um+sqrt(hp)*up)/(sqrt(hm)+sqrt(hp)) : 0. ;
    hhat = (hm+hp)/2.;
    Pnow = false;
    smin = Cgrav(hhat*w,p) +uhat;
    smax = Cgrav(hhat*w,p) +uhat; 
    for(j=0; j<2; j++){
        smin =  min3(um + pow(-1.,j)*Cgrav(hm*w,Pm), uhat + pow(-1.,j)*Cgrav(hhat*w,p),smin);
        smax =  max3(up + pow(-1.,j)*Cgrav(hp*w,Pp),uhat+ pow(-1.,j)*Cgrav(hhat*w,p),smax);
        }
    s[0] = smin;
    s[1] = smax;

    if(isnan(s[0]) || isnan(s[1]))
    {
	    printf("well shit. s is [%f,%f]\n", s[0],s[1]);
	    printf("qm = %f, qp = %f, um = %f, up = %f,\n", q1m, q1p, up, um); 
    }
}

/** Because C89 is too ghetto to have a predefined min function*/
double Channel::min3(double a, double b, double c)                       
{
	return min(min(a,b),c);
}
double Channel::max3(double a, double b, double c)
{
	return max(max(a,b),c);
}

/**
 * Step source term S in u_t+F(u)_x = S as follows:
 * q(t+dt) = q(t) +dt*S(q(t)+dt/2*S(q))
 ****/
void Channel::stepSourceTerms(double dt){
	int i;
	for (i=0; i<N; i++)
	{
		q[idx(1,i)]= q[idx(1,i)] +dt*getSourceTerms(q[idx(0,i)], q[idx(1,i)] +dt/2.*getSourceTerms(q[idx(0,i)],q[idx(1,i)]));
		q0[idx(1,i)] = q[idx(1,i)];
      //  double u= q[idx(1,i)]/q[idx(0,i)];
      //  double S1 = getSourceTerms(q[idx(0,i)],q[idx(1,i)]);
       // double q2bar = q[idx(1,i)] +dt/2.*S1;
       // double ubar = q2bar/q[idx(0,i)];
      //  double S2 =  getSourceTerms(q[idx(0,i)],q2bar);
     //   double b1= S1>0? -2.*u/S1: 0.;
    //    double b2= S2>0? -4.*ubar/(S2):0.;
    //    printf("u = %f, ubar = %f, S1 = %f, S2 == %f,dt = %f, CFL = %f",u,ubar,S1,S2, dt,min(b1,b2));
        if (q[idx(1,i)]>10) printf("whoops Q[%d] = %f\n",i,q[idx(1,i)]);
	}
	
}

/**
 * Evaluate source term S - friction slope (Manning Eqn) and actual slope S0 = -dz/dx (S0 should be a non-negative constant)
 * should I have another option for frictions slope?? Eg. Hazen Williams, Darcy-Weishback?
 * in the past this has posed some issues for very small A...be careful...
 **/
double Channel::getSourceTerms(double A, double Q){
	double Sf=0;
	double tol = At*1e-2;
	if (A>tol)
	{
		Sf = pow(Mr/kn,2)*Q*fabs(Q)/(A*A*pow(getHydRad(A),4./3.));	
	}
    //printf("Sf = %f S0 = %f A = %f, Q = %f, Rh^(4/3)=%f\n", Sf, S0, A,Q, pow(getHydRad(A),4./3.));
	return (S0-Sf)*G*A;
}

double Channel::getMassTransCoeff(double ui)//from Rossman 1994 J. Env. Eng pg 805
{
    double kf=0;
    double nu=1e-6*1.004;
    double DDiff=1.47e-7;//m^2/s (????) http://www.gsi-net.com/en/publications/gsi-chemical-database/single/107-chlorine.html
    double ReSc = ui*w/DDiff;
    double Sh = 3.65+(ui>0? (0.0668*(w/L)*ReSc)/(1+0.04*pow((w/L)*ReSc,2./3.)):0.);
    kf = Sh*DDiff/w;
    return kf;
}

double Channel::getKCl(double ai, double qi)//get Chlorine decay coefficient from Rossman 1994 J Env. Eng
{
    double Keff=0.;
    double rh = getHydRad(ai);
    double ui = ai>dry? qi/ai:0.;
    double kf = getMassTransCoeff(ui);
    if (ai<dry){Keff=kb;}
    else{
        Keff =kb+(kw*kf)/(rh*(kw+kf));
    }
    if (Keff>1){
        printf("kb = %f, kf = %f, ui= %f, ai = %f, rh = %f, Keff=%f\n",kb,kf,ui,ai,rh,Keff);}
    return Keff; 
}

void Channel::stepTransportTerms(double dt){
    //first order upwinding for Cl transport terms
    double ui =0.,ai = 0.,qi = .0,dCl=0.;
    double nu = dt/dx;
    ai = q[idx(0,0)];
    qi = q[idx(1,0)];
    ui = ai>dry? qi/ai:0.;
    for (int k =0; k<N; k++)
    if(Cl0[k]<0)
    {
        printf("Whoops! Cl[%d] = %f is negative. Rounding up to zero!\n", k,Cl0[k]);
        Cl0[k] = 0;
    }
  //  printf("n = %d, bCll =%f, bClr = %f\n",n,bCll, bClr);
    double KCl = getKCl(ai,qi);
    if (ui<0){dCl = Cl0[1]-Cl0[0];}
    else{ dCl = Cl0[0]-bCll;}
    Cl[0] = Cl0[0]-nu*ui*dCl-dt*KCl*Cl0[0];
    for(int i=1; i<N-1; i++)
    {
        ai = q[idx(0,i)];
        ui = ai>dry? q[idx(1,i)]/ai :0;//is this the right choice for u?
        KCl = getKCl(ai,qi);
        if (ui<0){dCl = (Cl0[i+1]-Cl0[i]);}
        else {dCl = Cl0[i]-Cl0[i-1];}
       // printf("t = %f, i= %d, ai = %f, ui = %f, Cl = %f, Cl0=%f, KCl =%f\n",dt*(float)n,i,ai, ui, Cl[i], Cl0[i], KCl);
        Cl[i] = Cl0[i]-nu*ui*dCl-dt*KCl*Cl0[i];
       // printf("Cl now = %f\n", Cl[i]);
        //if (Cl[i]>0){
        //}
    }
    ai = q[idx(0,N-1)];
    ui = ai>dry? q[idx(1,N-1)]/ai:0.;
    KCl = getKCl(ai,qi);
    if (ui<0){dCl = bClr-Cl0[N-1];}
    else{dCl = Cl0[N-1]-Cl0[N-2];}
    Cl[N-1] = Cl0[N-1]-nu*ui*dCl-dt*KCl*Cl0[N-1];
    //update Cl0
    for (int i=0;i<N;i++){
        Cl0[i]= Cl[i];
    }
}

/**
 * Get the volume of water in the channel
 * V = (Sum of cell values of A) x (cell width dx)
 * **/
double Channel::getTheGoddamnVolume()
{
	double evol=0.;
	for(int j = 0; j<N; j++)
	{
		evol  += dx*q[j];
	}
	return evol;
}

/**
 * Estimate kinetic energy via
 * KE = \sum_{k=0}^{N-1}u_i^2*A_i/(2g) where u_i = Q_i/A_i
 **/
double Channel::getKE(int i)
{
	double ui =0,KE=0;
	for (int j = 0; j<N; j++)
	{
		double ai = q_hist[idx_t(0,j+1,i)];
		ui = ai>dry? (q_hist[idx_t(1,j+1,i)]/ai): 0;
		KE+= ui*ui/(2*G)*ai;
	}
	return KE;
}

/*Estimate potential energy via
 *PE = Ah-Eta/(gA) 
 **/
double Channel::getPE(int i)
{
	double hi =0,PE=0;
	bool p = false;
	for (int j = 0; j<N; j++)
	{
		double ai = q_hist[idx_t(0,j+1,i)];
		PE+= (ai*HofA(ai,p)-(ai>dry?  Eta(ai,p)/(G*ai):0)); 
	}
	return PE;
}

/**
 * Estimate <dh/dx> at time t_i where h is the PRESSURE HEAD NOT WATER HEIGHT. OOPS.
 * add up |h_{k+1}-h_k| for all k (it's just |(h_{k+1}-h{k})/dx|*dx	
 **/
double Channel::getAveGradH(int i)
{
	bool p = false;
	double h1, h2,a1,a2;
	a1 = q_hist[idx_t(0,1,i)];
	a2 = q_hist[idx_t(0,2,i)];
	h1 = pbar(a1,p);
	h2 = pbar(a2,p);
	double I = 0;
	for (int k = 1; k<N+1; k++)
	{
		I+= abs(h2-h1);
		h1 = h2;
		a1 = a2;
		a2 = q_hist[idx_t(0,k+1,i)]; 
		h2 = pbar(a2,p);
	}
	return I;
}


/**
 * Write q information to (M/Mi) files 
 * These are in theory readable by smartputittogether.py or maybe something else.py...(!???)
 * */
int Channel::writeqToFile(int Mi, double dt)
{	
	
	for(int kk = 0; kk<M+1; kk+=Mi)
	{
		char fname[100];
		snprintf(fname, sizeof fname, "../movie/howdy%03d.txt", kk/Mi);
		FILE *fp = fopen(fname, "w");
		if(fp==NULL)
		{
			fputs("Error opening file!!", stderr);
			return 1;
		}
		//warning: messing with this information will screw up the python plotting script...
		fprintf(fp, " #x   h(A(x))    Q(x)  at t = %f (L = %f, dx = %f, T = %f,  dt = %f ) \n",dt*((float)kk), L,dx, ((double)M)*dt,dt); 
		for(int jj = 0; jj<N+2; jj++)
		{
			bool p = false;
		       fprintf(fp,"%d     %.10f     %.10f\n", jj, HofA(q_hist[idx_t(0,jj,kk)],p),q_hist[idx_t(1,jj,kk)]);
		}
		fclose(fp); 
	}	
	return 0;
}

/**
 * Write characteristics Q/A+sign*phi(A) to file in format readable by gnuplot
 * */
int Channel::writeRItofile(double dt, int sign) 
{
		char fname[100];
		double t,x,Q,A,RI;
		if(sign >0){snprintf(fname, sizeof fname, "../RIpics/RIplus.txt");}
		else{snprintf(fname, sizeof fname, "../RIpics/RIminus.txt");}
		FILE *fp = fopen(fname, "w");
		if(fp==NULL)
		{
			fputs("Error opening file!!", stderr);
			return 1;
		}
		fprintf(fp, "#Riemann invariants u+/-phi, T=%f, L = %f\n", dt*(double)M, L);
		for(int kk = 0; kk<M+1; kk++)
		{
			for(int jj = 0; jj<N+2; jj++)
			{
				t = kk*dt;
				x = (jj-1)*dx;
				A = q_hist[idx_t(0,jj,kk)];
				Q = q_hist[idx_t(1,jj,kk)];
				RI = Q/A +sign*PhiofA(A,false);
				if(sign==0){RI = A;}
		        	fprintf(fp,"%.10f     %.10f     %.10f\n", x, t, RI);
		        	//printf("%.10f     %.10f     %.10f\n", x, t, RI);
			}
			fprintf(fp, "\n");
		}
		fclose(fp); 
	return 0;
}

/*
 * Quickly write out some run time info
 * where is either times or places (if which[k] = 0, then where is a time; else it's a place)
 * K is the length of where
 **/
void Channel::quickWrite(double *where, int *which, int K, double T, int skip)
{
	int i, MN;
	double ds;
	for (int k=0; k<K; k++)
	{	
		//write variables at all locations at time where[k] 
		if(which[k]==0)
		{
			i = round(where[k]*(double)M/K);
			if(i>=M){i=M;}
			if(i<0){i=0;}
			printf("t = %f, index = %d\n", where[k], i);
			MN = N;
			ds = dx;
		}
		//write variables at all times at location where[k]
		else  
		{
			i = round(where[k]*(double)N/L);
			if(i>=N){i=N-1;}
			if(i<0){i=0;}
			printf("x = %f m, index = %d\n", where[k], i);
			MN = M;
			ds = T/(double)M;
		}
		printf("%s             A               h               Q \n",(which[k])?"t":"x");
		double A, h, hfake, Q;
		for(int j = 0; j<=MN; j+=skip)
		{
			if(which[k]==0)  		
			{
				A = q_hist[idx_t(0,j,i)];
				Q = q_hist[idx_t(1,j,i)];

			}
			else
			{
				A = q_hist[idx_t(0,i,j)];
				Q = q_hist[idx_t(1,i,j)];
				
			}

			h = HofA(A,false);
			hfake = fakehofA(A,true);

			printf(" %f    %f     %f    %f \n", ds*(double)j,A, h, Q);
		}
	}	
}

/**
 * Set Preissman parameters 
 * a is the pressure wave speed for this pipe -- defaults to 1200
 ***/
void Cpreiss::setGeom(double a_)

{
	a = a_;
	Af = PI*D*D/4.;
	Ts = G*Af/(a*a);
	tt = 2*(PI-asin(Ts/D));// theta such that c(A(theta)) = a
	yt = D/2.*(1-cos(tt/2.));
	At = AofH(yt,false);
	//write out a file tabulating A, g, I, c, phi, A(phi), A(h(A)) etc (use as sanity check)
	char fname1[30];
	clock_t time1 = clock();
	int duh = (int)time1;
	sprintf(fname1,"geomconfirm%d.txt",1);
	FILE *fg1 = fopen(fname1,"w");
	int Mp= 500;
	fprintf(fg1, "#D = %f, cgrav = %f, Ts = %.16f, At = %.16f, yt = %.16f \n",D, a, Ts,At,yt); 
	fprintf(fg1, "#A                  h(A)                   I(A)                  c(A)                  phi(A)                  A(phi(A))                  A(h(A))   htrue(A)\n");
	bool pp = false; 
	for(int k = 0; k<Mp; k++)
	 {
		 double tt = PI*2/Mp*(double)k;
		 double aa = D*D/8.*(tt-sin(tt));
		 double ht = 0.5*D*(1-cos(tt*.5));
		 double h = HofA(aa,pp);
		 double I = Eta(aa,pp);
		 double ah = AofH(ht,pp);
		 double c = Cgrav(aa,pp);
		 double phi = PhiofA(aa,pp);
		 double ae = AofPhi(PhiofA(aa,pp),pp);
		 fprintf(fg1,"%.16f   %.16f    %.16f   %.16f   %.16f   %.16f   %.16f    %.16f\n", aa, h, I, c, phi, ae, ah, ht); 
	 }
	 fclose(fg1);

	char fname2[30];
    sprintf(fname2, "geomconfirm%d.txt", 2);
	FILE *fg2 = fopen(fname2, "w");	
	fprintf(fg2, "A                   h(A)       Phi(A)              A(Phi(A))    A-Af\n");
	for (int i=0; i<2;i++)
	{
		for(int k = 0; k<40; k++)
		{
		 	double aa =D*D*PI/4.*(1.+pow(-1.,i+1)*pow(2.,-k-1));
		 	double h = HofA(aa,pp);
		 	fprintf(fg2, "%.16f   %.16f %.16f  %e %e   \n", aa, h, PhiofA(aa,pp),AofPhi(PhiofA(aa,pp),pp)-aa, aa-D*D*PI/4.);
		} 
	}
	fclose(fg2);
}

	
double Cpreiss::AofH(double h, bool p)
{
//	if(h<=yt && (!p))
	if(h<=yt)
	{	double t = 2.*(acos(1.-h/D*2.));
		return	 D*D*(t-sin(t))/8.;
	}
	else
		return At+Ts*(h-yt);
}

double Cpreiss::Eta(double A, bool p)
{
	return ::Eta(A,D,At,Ts,p);
}


double Cpreiss::getHydRad(double A)
	{	
		double perim,R;
		if(A<(At))
		{
			double y = HofA(A,false);
			double theta = 2*acos(1.-2*y/D);
            if (theta<1e-4)
            {
                R = D/(24.)*theta*theta; ///taylor expand for small theta
                //printf("hyrdaulic radius= %f, theta =%f, A = %f\n",R,theta,A);
            }
            else{
			    perim = D*theta/2.;
                R = A/perim;
            }
		}
		else {R = D/4.;}
        return R;
    }  

/**
 *  Roe estimates for wavespeeds
 **/
void Cpreiss::speedsRoe(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp) 
{

    double um, up, hm, hp, uhat, hhat,chat, smin, smax;
    int j;
    bool p = Pm && Pp;
    um = (q1m>dry? q2m/q1m : 0.);
    up = (q1p>dry? q2p/q1p : 0.);
    hm = HofA(q1m,Pm);
    hp = HofA(q1p,Pp);
    uhat = (hm+hp >dry)? (sqrt(hm)*um+sqrt(hp)*up)/(sqrt(hm)+sqrt(hp)) : 0. ;
    hhat = (hm+hp)/2.;
    chat = Cgrav(AofH(hhat,p),p); 
    smin = chat +uhat;
    smax = chat +uhat; 
    for(j=0; j<2; j++){
        smin =  min3(um + pow(-1.,j)*Cgrav(q1m,Pm), uhat + pow(-1.,j)*chat,smin);
        smax =  max3(up + pow(-1.,j)*Cgrav(q1p,Pp),uhat+ pow(-1.,j)*chat,smax);
        }
    s[0] = smin;
    s[1] = smax;

   // printf("well shit. s is [%f,%f]\n", s[0],s[1]);
    if(isnan(s[0]) || isnan(s[1]))
    {
	    printf("well shit. s is [%f,%f], hm=%f hp = %f\n", s[0],s[1], hm, hp);
	    printf("qm = %f, qp = %f, um = %f, up = %f,\n", q1m, q1p, up, um); 
    }
}


/**
 * HLL speed estimates from Leon 2009 
 * */
void Cpreiss::speedsHLL(double q1m, double q1p, double q2m, double q2p, double *s, bool Pm, bool Pp){
//	double dry = 1e-6*At;                                                  //pay attention to this!?
    double cbar,Astar, ym,yp,cm, cp, um =0 , up= 0;	
	double Astarp;
	//other wave speed estimates from Sanders 2011 (keep around for comparison purposes?)
	double sl1, sl2, sr1, sr2;
	ym = HofA(q1m,Pm);
	yp = HofA(q1p,Pp);
	cm = Cgrav(q1m, Pm);
	cp = Cgrav(q1p, Pp);
	//if no dry bed present
	if(fmin(ym,yp)>=dry){
		cbar = (cp+cm)/2.;
		um = q2m/q1m;
		up = q2p/q1p;
		//this verstion uses depth positivity condition
		//Astar = (q1m+q1p)/2.*(1+(um-up)/(PhiofA(q1p,Pp)+PhiofA(q1m,Pm))); 
		//this is linearized version
		Astar = (q1m+q1p)/2.*(1+( (cbar>1e-6)? (um-up)/(2.*cbar): 0));
        //if (Astar<0)
        //{
           // printf("using depth positivity condition\n");
        // Astar = (q1m+q1p)/2.*(1+(um-up)/(PhiofA(q1p,Pp)+PhiofA(q1m,Pm))); 
        //}  
		bool Ps = (Pm && Pp);
		s[0] = um - findOmega(Astar, q1m, Ps, Pm);
		s[1] = up + findOmega(Astar, q1p, Ps, Pp);
		// Sanders estimates (not very good--use them if the other estimate messes up(??!))
//		double Vs, cs, phip, phim;
//		phim = PhiofA(q1m, Pm);
//		phip = PhiofA(q1p, Pp);
//		bool Pmp = !(!Pp||!Pm);
//		Vs = 0.5*(um+up+phim-phip);
//		double phis = 0.5*(phim+phip+um-up);
//		double Astars = AofPhi(0.5*(phim+phip+um-up),Pmp);
//		if (phis>AofPhi(Af,true)){cs = a;}
//		else{cs = Cgrav(Astars,Pmp);}
//		sl1 = um - cm;
//		sl2 = Vs - cs;
//		sr1 = up + cp;
//		sr2 = Vs + cs;
//		if (WTF){
//			printf("phis = %f, cs = %f\n",phis, cs);
////			printf("Atsar = %f, Astard = %f, AstarS = %f\nwith q1m = %f and q1p = %f, um =%f, up = %f\n",Astar, Astarp, Astars, q1m, q1p, um, up);
	//		printf("sL = %f, sR = %f, Sanders speeds (vl-cl, vs-cs) = (%f, %f) and (vr+cR, vr+cr) = (%f,%f)",s[0], s[1], sl1, sl2, sr1, sr2);
	//	}
	//	s[0] = min(sl1,sl2);
   // 	s[1] = max(sr1,sr2);

			
	}
	else{
		// Both sides dry - both stay dry
		if(fmax(ym,yp)<dry)     
		{
			s[0] = 0.;
			s[1] = 0.;
		}
		//left side dry
		else if(ym<dry)  
		{
			up = q2p/q1p;
			s[0] = up - PhiofA(q1p,Pp);
			s[1] = up + cp;	
        //    printf("left side dry. s[0] = %f, s[1] = %f, up = %f, q1p = %f q2p = %f cp = %f \n", s[0], s[1], up, q1p,q2p,cp );
		}
		//right side dry
		else if(yp<dry) 
		{
			um = q2m/q1m;
			s[0] = um - cm;
			s[1] = um + PhiofA(q1m,Pm);
          //  printf("right side dry. s[0] = %f, s[1] = %f, um = %f, q1m = %f q2m = %f cm = %f \n", s[0], s[1], um, q1m,q2m,cm );
		}
		
	}
	if(isnan(s[0]) || isnan(s[1])) //check for NaNs
	{
		printf("Error!nan speeds!s[0] = %f, s[1] = %f, with y1 = %.10f, y2 = %.10f\n", s[0], s[1],ym, yp);
		printf("q1p is %f and q2p is %f, um is %f, up is %f\n", q1p, q2p, um, q2p);
		exit (1);
	}
	//check that HLL speeds have s0<s1; really, this should never bloody happen!
	if(s[0]>s[1]) 	
	{  
//	printf("The hell? s[0]>s[1], with q1m = %f, q2m = %f, q1p =%f, q2p =%f, s[0] = %f, s[1] = %f\n",q1m,q2m,q1p,q2p,s[0],s[1]);
	//use Sanders wave speed estimates instead if this happens
	double stemp = s[0];
    s[0] = s[1];
    s[1] = stemp;
    //s[0] = min(sl1,sl2);
//    s[1] = max(sr1,sr2);
	}
	cmax = max(cmax, max(fabs(s[0]),fabs(s[1]))); 		   
}

/**
 * Compute Omega (used by speedsHLL estimate)
 **/
double Cpreiss::findOmega(double Astar, double Ak, bool Ps, bool Pk)
{
	double omega;
	double eps = 1e-8;
   if (Astar>Ak+eps)
		{
			omega = sqrt((Eta(Astar, Ps)-Eta(Ak, Pk))*Astar/(Ak*(Astar-Ak)));
		} 
	else
    	{
        	omega  = Cgrav(Ak,Pk);
    	}
	return omega;
}



/**
 * The class constructor initializes parameters and allocates memory for and sets the boundary value time series.
 * \param[in] &ch0 the channel whose ghost cell is informed by this routine
 * \param[in] a_which the
 **/
Junction1::Junction1(Channel &a_ch0, int a_which, double a_bval, int a_bvaltype):ch0(a_ch0)
{
	N = ch0.N;
	bval = new double[ch0.M+1];
	Clbval = new double[ch0.M+1];
	for(int i=0;i<ch0.M+1; i++){
            bval[i] = a_bval;
            Clbval[i] = 0.;
        }
	bvaltype = a_bvaltype;
	whichend = a_which;
	w = ch0.w;
	//default is to use RI
	reflect = 0;
        //initialize boundary chlorine value
}

/**apply boundary conditions at the end of a single pipe. Several options.
 **What computation happens depends on these parameters:
 * --reflect = +1, -1, 0  (see below for details)
 * --whichend = 0, 1   (0 for left (x=0) end, 1 for right (x=L) end. Matters for signs of things.)  
 * --bvaltype = 0, 1   (0 for A, 1 for Q)  
 * Reflect: reflect all incoming waves (hard boundary), set reflect =1
 * Extrapolate: let all incoming waves through (ghost cell has same values as last interior cell), set reflect = -1
 * Specify: specify Q,A
 *			specify f(Q,A) = 0 (not implemented yet) 
 *     follows a characteristcic out of the domain to solve for unknown value, then updating boundary fluxes accordingly.
 * 
 * Summary of cases:
	  case 0: reflect everything
	  case 1: reflect nothing
	  case 2: supercritical 
				2.1 inflow: specify A and Q
				2.2 outflow: same as case 1
	  case 3: subcritical: specify A(t) or Q(t)
				3.1 can solve using RI
					3.1.0 specify Q<Qtol
					3.1.1 specify Q>=Qtol
					3.1.2 specify A
				3.2 can't solve using RI, set A = Ain
	  case 4: orifice outflow
 **/
void Junction1::boundaryFluxes()
{
	int bccase, Nq, Npi, Npe;	
	double Ain, Qin, Aext, Qext;
	double Cc=0, Cd = 0, dp0 = 0, eext=0;  //random crap for dysfunctional orifice routine
	double Qtol = 1e-12;
	bool Pin, Pext;
	double ctol = 0;		      //tolerance for trying to solve for characteristic solns
    double sign = pow(-1.,whichend+1);    //gives the sign in u (plus or minus) phi (left side: -1, right side: +1)
	//update bCl to current value if we care
    double bCl = Clbval[ch0.n];
    if (whichend)                         //if we're on the right side of pipe
	{
		Nq = N-1;
		Npe = N+1;
		Npi = N;
        if(reflect==1) bCl = ch0.Cl[N-1];
        ch0.bClr=bCl;
        ch0.Cl_hist[(N+2)*ch0.n+N+1]=bCl;
	}
	//if we're on the left side of the pipe
	else								  
	{
		Nq = 0;
		Npe = 0;
		Npi = 1;
        if (reflect ==1) bCl = ch0.Cl[0];
        ch0.bCll=bCl;
        ch0.Cl_hist[(N+2)*ch0.n]=bCl;
	}
	
	Ain = ch0.q[ch0.idx(0,Nq)];
	Qin = ch0.q[ch0.idx(1,Nq)];
	ch0.Pnow = ch0.P[Nq];
	Pin = ch0.P[Npi];
	Pext = ch0.P[Npe];
	if(reflect ==1) bccase = 0;
	else if(reflect ==-1) bccase =1;
	else
	{
		if (fabs(Qin)/Ain>ch0.Cgrav(Ain, Pin))  //first check for supercritical flow
		{
			Qext = bval[ch0.n];		
			if (fabs(Qext)<Qtol) bccase =0;      //just reflect
			else if((whichend &&(Qext>Qtol)) ||(!whichend &&Qext<-Qtol)) bccase=1;          //supercitical outflow -> extrapolate
			else bccase = 21;               //dredge up a reasonable value of Aext or Qext!
		}
		else
		{
			if (bvaltype==0)bccase = 312;		//if specifying A
			else if(bvaltype==1)				//if specifying Q
			{
				Qext = bval[ch0.n];				//sample the boundary time series at the current time step
				if(fabs(Qext)<Qtol)	bccase = 310;//if Qext != 0, problem is harder than if Q = 0 (cases treated seperately)
				else
				{
					//first use Qext to find values c1min and c1max -- range of c+/_(A, Qext) over A in [0, \infty] 
					//if c1 = c+/_(Ain, Qin) lies in  [c1min, c1max], then it is feasible to follow a characteristic out of the doman
					double c1min=0., c1max=0., c1, xs=0.;
					double uin = (Ain>ch0.dry ?Qin/Ain :0. );
					if(ch0.channeltype ==0) //uniform cross section case is easy
					{
						c1 = uin +sign*2.*sqrt(G*Ain/w);
						c1min = 3.*pow((G*fabs(Qext)/w),1./3.); //min achievable value for for c+  (outgoing on right)
						c1max = -c1min;  			//max achievable value for c_ (outgoing on left)
					}
					else  //Preissman slot requires SUBTERFUDGE (formerly, rootfinding)
					{	
						double xhat = pow(ch0.w/G*Qext*Qext,1./3.);
						c1min = Qext/xhat + ch0.PhiofA(xhat, false);//estimate bounds with uniform cross section values.
						c1max = Qext/xhat - ch0.PhiofA(xhat, false);
						c1  = uin +sign*ch0.PhiofA(Ain,Pin);
						if(WTF) printf("c1 = %f,c1 min = %f, c1max = %f, Qext = %f, Qin = %f, Ain = %f\n",c1,c1min,c1max, Qext, Qin, Ain);
					}
					if((whichend ==0 && Qext<0. && c1>c1max)||(whichend ==1 && Qext>0. && c1<c1min))//make sure left end boundary flux is enforceable 
					{
						printf("oops! Qext = %f is too %s for c1 = %f\n, setting Aext =Ain= %f.note Qin = %f\n", Qext, whichend?"large":"small",c1,Ain,Qin);
						bccase =32;
					}	
					else bccase = 311;		
				}

			}
		}
	}
	if(bvaltype==2){
		if (Ain<ch0.At) bccase =1;
		else bccase=4;
	}
	switch(bccase)
	{
		case 0://reflect everything
			Aext = Ain;
			Qext = -Qin;
			break;
		case 1://reflect nothing i.e. extrapolate
			Aext = Ain;
			Qext = Qin;
			break;
		case 21://supercritical inflow...make up correct value of Aext...?
			if (bvaltype==0)
			{
				Aext = bval[ch0.n];
				Qext = (Qin*Aext/Ain)*(ch0.Cgrav(Aext, Pin)*ch0.Cgrav(Ain, Pin));
			}
			else
			{
				Qext =bval[ch0.n];
				//Aext = Ain*Qext/Qin;
				Aext = ch0.q_hist[ch0.idx_t(0,0,max(ch0.n-1,0))];        //use previous value of Aext unless you are at the beginning...
			}
			break;
		case 310://subcritical, specify Q<Qtol	
			if(ch0.channeltype ==0) //uniform channel	
			{
				double uin = (Ain>0. ?Qin/Ain :0. );
				Aext = w/(G*4.)*pow((uin +sign*2.*sqrt(G*Ain/w)),2.);
			}
			else //Preissman slot
			{       			
				double lhs = Qin/Ain +sign*ch0.PhiofA(Ain,Pin);
				if(sign*lhs>=0)
				{
					Aext = ch0.AofPhi(sign*lhs,false);
				}
				else //I dunno, reflect?
				{
					cout<<"yeah I dunno, setting Aext =Ain and Qext = -Qin"<<endl;
					printf("Qin = %f, Ain  %f, lhs= %f, whichend = %d\n",Qin, Ain, lhs,whichend);
					Aext = Ain;
					Qext = -Qin;
					bval[ch0.n] = -Qin;
				}				
			}
			break;
		case 311://subcritical, specify Q>=Qtol
			if(ch0.channeltype ==0)				//if uniform cross section
			{	
				double c2 = 2.*sign*sqrt(G/w);
				double uin = Ain>ch0.dry? Qin/Ain :0.;
				double c1 = uin +sign*2.*sqrt(G*Ain/w);
				fRI f(c1,c2, bval[ch0.n]);
				dfRI df(c2,bval[ch0.n]);
				Aext = Newton(f,df,Ain, 1e-10, 100);
			}
			else								//Preissman slot
			{
				double uin = (Ain>ch0.dry ?Qin/Ain :0. );
				double lhs = uin +sign*ch0.PhiofA(Ain,Pin);//solve lhs = Qext/x +sign*phi(x) for x
				if(WTF)	cout<<"UHHHHHMM reflect = "<<reflect<<" Qext = "<<Qext<<"  lhs = "<<lhs<<"sign ="<<sign<<" Ain="<<Ain<<endl;
				fallpurpose fp(ch0.w, ch0.At,ch0.Ts, lhs, Qext, sign,1.,0.,false);
				if(sgn(fp(0.)*sgn(fp(ch0.Af*2.0))<0))//careful with the Ridders solver
				{
					int count;
					Aext = ::ridders(fp,0.,ch0.Af*2.0,&count, 1e-10, 1e-10);//this line crashes things quite frequently...
					double uext = (Aext>ch0.dry ?Qext/Aext :0.);
					double err = fabs(uext +sign*ch0.PhiofA(Aext,Pext)-lhs);
				}
				else Aext = Ain;
			}
			break;
		case 312://subcritical, specify A
			Aext = bval[ch0.n];
			Qext = (Qin/Ain+sign*ch0.PhiofA(Ain,Pin) - sign*ch0.PhiofA(Aext,Pext))*Aext;
			if (fabs((Qext/Aext))>ch0.Cgrav(Aext, Pext)) Qext = Qin;
			break;	
		case 32://this is the ohshitcase
			Aext = Ain;
			Qext = Qin;
			printf("warning. boundary case 3.2 (ohshit case). Extrapolating Ain=%f and Qin =%f\n",Ain, Qin);
			bval[ch0.n] = Qext;
			break;
		case 4://Orifice outflow
			Cd = 0.78; //discharge coefficient
			Cc = 0.83;//contraction coefficient (see Trajkovic 1999)
			eext = bval[ch0.n];
			dp0 = ch0.HofA(Ain, Pin)-Cc*eext;
			if (dp0>0 && Qin>0)
			{
				Aext = ch0.AofH(Cc*eext, false);
				Qext = Cd*Aext*sqrt(2.*G*dp0);
				Aext = Ain;
		//	these two lines work but are boring.
			//	Qext = Qin;
			//	Aext = Qext/(Cd*sqrt(2*G*dp0));
				cout<<"orifice eqn! Qext = "<<Qext<<"  Qin = "<<Qin<<"  Aext = "<<Aext<<"  Ain = "<<Ain<<"  eext ="<<eext<<" h(Aext)="<<ch0.HofA(Aext,Pext)<<endl;
			}
			else
			{
				Qext = Qin;
				Aext = Ain;
			}
			break;
		default:
			Aext = Ain;
			Qext = Qin;
		break;
	}	
	//compute the fluxes using numFlux
	if(whichend)//right end of pipe
	{	
		ch0.numFlux(Ain, Aext, Qin, Qext, ch0.bfluxright, ch0.P[N], ch0.P[N+1]);
		if(WTF) printf("\n in junction routine!Aext =%f, Ain = %f, Qin %f, Qext = %f, bfluxright = [%f,%f]\n",Aext, Ain, Qin, Qext,ch0.bfluxright[0],ch0.bfluxright[1]);
	}
	else
	{
		ch0.numFlux(Aext, Ain, Qext, Qin, ch0.bfluxleft, ch0.P[0], ch0.P[1]);
		if(WTF)	printf("\nin junction routine!Aext =%f, Ain = %f, Qin %f, Qext = %f, bfluxleft = [%f,%f]\n",Aext, Ain, Qin, Qext,ch0.bfluxleft[0],ch0.bfluxleft[1]);
	}
	ch0.q_hist[ch0.idx_t(0,Npe,ch0.n)] = Aext;
	ch0.q_hist[ch0.idx_t(1,Npe,ch0.n)] = Qext;
	//update boundary pressurization states
	if(reflect ==1||Qext==0){ch0.P[Npe] = true;}//for some reason this cannot be messed with, sigh.
	if(reflect ==-1){ch0.P[Npe] = ch0.P[Npi];}//for some reason this cannot be messed with, sigh.
	else if(bvaltype==0 && Aext<ch0.At){ch0.P[Npe] = false;}
	else if(Aext>ch0.At){ch0.P[Npe]= true;}
	else{ch0.P[Npe] =ch0.P[Npi];}
	ch0.p_hist[ch0.pj_t(Npe,ch0.n)] = ch0.P[Npe];
   // printf("Ain=%f Qin =%f Aext = %f Qext = %f\n",Ain, Qin, Aext, Qext);
    //	printf("Ain is %f and Qin is %f and Aext is %f and Qext is %f for end %d\n Pin is %d and Pext is %d\n", Ain, Qin, Aext, Qext, whichend, Pin, Pext);

}

double Junction1::getFlowThrough(double dt)
{

	double Qt = 0;
	int N0;
	if (whichend)
	{
		N0 = ch0.N-1;	
	}
	else
	{
		N0 = 0;
	}
	for (int i = 0; i<ch0.M+1; i++)
	{	
		Qt+=dt*ch0.q_hist[ch0.idx_t(1,N0,i)];
	}
	return Qt;
}


//overloaded functions to set boundary value time series

/** set boundary value to a constant value*/
void Junction1::setbVal(double bvalnew)
{
	for(int i =0; i<ch0.M+1; i++)
	{
		bval[i] = bvalnew;
	}
}

void Junction1::setClbVal(double Clbvalnew)
{
    for(int i =0; i<ch0.M+1; i++)
	{
		Clbval[i] = Clbvalnew;
	}
    
}

void Junction1::setClbVal(double *Clbvalnew)
{
     for(int i =0; i<ch0.M+1; i++)
	{
		Clbval[i] = Clbvalnew[i];
	}
    
}

void Junction1::setClbVal(vector<Real> x){
    for (int i=0; i<ch0.M+1;i++)
    {
        Clbval[i] = x[i];
//        printf("Clbval[%d]=%f\n",i, Clbval[i]);
    }
}

    /** set boundary value to a time series stored in various arrays*/
void Junction1::setbVal(valarray<Real> x)
{
	for(int i =0; i<ch0.M+1; i++)
	{
		bval[i] = x[i];
	}
}

void Junction1::setbVal(vector<Real> x)
{
	for(int i =0; i<ch0.M+1; i++)
	{
		bval[i] = x[i];
	}
}


void Junction1::setbVal(double*x)
{
	for(int i =0; i<ch0.M+1; i++)
	{
		bval[i] = x[i];
	}
}

Junction1::~Junction1()
{
	delete [] bval;
	delete [] Clbval;
}


Junction2::Junction2(Channel &a_ch0, Channel &a_ch1, int a_which0, int a_which1, double a_valveopen):ch0(a_ch0),ch1(a_ch1){
	whichend0 = a_which0;
	whichend1 = a_which1;
	if(whichend0){
		N0 = ch0.N-1;
		Ns0 = ch0.N+1;
	}
	else{
		N0 = 0;
		Ns0 = 0;
	}
	if(whichend1){
		N1 = ch1.N-1;
		Ns1 = ch1.N+1;
		
	}
	else{
		N1 = 0;
		Ns1= 0;
	}
	valvetimes.resize(ch0.M+1);
	valveopen = a_valveopen;
	valvetimes.resize(ch0.M+1);
	for(int i=0; i<ch0.M+1; i++)valvetimes[i]=valveopen;
	offset =0;
        dx0 = ch0.L/(double)ch0.N;
        dx1 = ch1.L/(double)ch1.N;        
}


/*fluxes across serial junction - only one reasonable way to do it:
*
*
*   ^ ~~~~~~~~            channel 0 sees h1 -offset
*   |        | 			  channel 1 sees h0 +offset			
*   h0       |
*   |        |~~~~~~~~~~  ^
*   v _______|            |  h1    
*   ^        |            | 
*  offset    |            |
* ..v........|__________  v............ground   
*            | 
*<-channel 0 | channel 1 ->
*
* */


void Junction2::setValveTimes(vector<Real>x)
{
	for (int i = 0; i<ch0.M+1; i++)
	{
		valvetimes[i]= x[i];
	}
}


void Junction2::setValveTimes(valarray<Real>x)
{
	for (int i = 0; i<ch0.M+1; i++)
	{
		valvetimes[i]= x[i];
	}
}

void Junction2::boundaryFluxes(double dt)
{	
	double q1m, q1p, q2m, q2p, q1mfake, q1pfake;
	valveopen = valvetimes[ch0.n];
	q1m = ch0.q[ch0.idx(0,N0)];
	q2m = ch0.q[ch0.idx(1,N0)];
	q1p = ch1.q[ch1.idx(0,N1)];
	q2p = ch1.q[ch1.idx(1,N1)];
	bool pm = ch0.P[ch0.pj(N0)];
	bool pp = ch1.P[ch1.pj(N1)];
    double nu0 = dt/dx0;
    double nu1 = dt/dx1;
    //attempt at incorporating valve opening coefficient
//	printf("valveopen=%f\n",valveopen);
    if(valveopen>0)
	{
		double h0f = ch1.HofA(q1p,pp)-offset;
		double h1f = ch0.HofA(q1m,pm)+offset;
	   	//what channel 0 sees
		q1pfake = ch0.AofH(h0f,pp);     
		//what channel 1 sees
		q1mfake = ch1.AofH(h1f,pm); 
        //Chlorine upwinding 
        double dCll, dClr, ul, ur;            
        ul = q1m>ch0.dry? q2m/q1m:0.;
        ur = q1p>ch1.dry? q2p/q1p:0.;
		if(whichend0)//right end of pipe 0
		{
            if(!whichend1)//left end of pipe 1
            {
                ch0.numFlux(q1m,q1pfake, q2m, q2p*valveopen, ch0.bfluxright, pm, pp);
                ch1.bfluxleft[0] = ch0.bfluxright[0];	
                ch1.bfluxleft[1] = ch0.bfluxright[1];
                if (ul<0){dCll = ch1.bCll-ch0.bClr;}
                else{dCll = ch0.bClr-ch0.Cl[N0];}
              //  dClr = dCll;
                if(ur<0){dClr = ch1.Cl[0]-ch1.bCll;}
                else{dClr = ch1.bCll-ch0.bClr;}
                ch0.bClr += -nu0*ul*dCll -dt*ch0.getKCl(q1m,ul)*ch0.bClr;
                ch1.bCll += -nu1*ur*dClr -dt*ch1.getKCl(q1p,ur)*ch1.bCll;
                ch0.Cl_hist[ch0.n*(ch0.N+2)+N0+1] = ch0.bClr;
                ch1.Cl_hist[ch1.n*(ch1.N+2)] = ch1.bCll;  
            }
            else//also right end of pipe 1
            {
                ch0.numFlux(q1m,q1pfake, q2m, q2p*valveopen, ch0.bfluxright, pm, pp);
                ch1.bfluxright[0] = ch0.bfluxright[0];	
                ch1.bfluxright[1] = ch0.bfluxright[1];
                if (ul<0){dCll = ch1.bCll-ch0.bClr;}
                else{dCll = ch0.bClr-ch0.Cl[N0];}
              //  dClr = dCll;
                if(ur<0){dClr = ch1.Cl[N1]-ch1.bClr;}
                else{dClr = ch1.bCll-ch0.bClr;}
                ch0.bClr += -nu0*ul*dCll -dt*ch0.getKCl(q1m,ul)*ch0.bClr;
                ch1.bCll += -nu1*ur*dClr -dt*ch1.getKCl(q1p,ur)*ch1.bClr;
                ch0.Cl_hist[ch0.n*(ch0.N+2)+N0+1] = ch0.bClr;
                ch1.Cl_hist[ch1.n*(ch1.N+2)+N1+1] = ch1.bClr;  
            }
        }
		else//left end of pipe0 
		{
            if(whichend1)//right end of pipe 1
            { 
                ch0.numFlux(q1pfake, q1m, q2p*valveopen, q2m, ch0.bfluxleft, pp, pm);
                ch1.bfluxright[0] = ch0.bfluxleft[0];
                ch1.bfluxright[1] = ch0.bfluxleft[1];
                if (ul<0){dCll = ch0.bCll-ch1.bClr;}
                else{dCll = ch1.bClr-ch1.Cl[N1];}
              //  dClr = dCll;
                if (ur<0){dClr = ch0.Cl[0]-ch0.bCll;}
                else{dClr = ch0.bCll-ch1.bClr;}
                ch1.bClr += -nu1*ul*dCll-dt*ch1.getKCl(q1m,ul)*ch1.bClr;
                ch0.bCll += -nu0*ur*dClr-dt*ch0.getKCl(q1p,ur)*ch0.bCll;
                ch0.Cl_hist[ch0.n*(ch0.N+2)] = ch0.bCll;
                ch1.Cl_hist[ch1.n*(ch1.N+2)+N1+1] = ch1.bClr;
            }
            else //also left end of pipe1
            {
                ch0.numFlux(q1pfake, q1m, q2p*valveopen, q2m, ch0.bfluxleft, pp, pm);
                ch1.bfluxleft[0] = ch0.bfluxleft[0];
                ch1.bfluxleft[1] = ch0.bfluxleft[1];
                if (ul<0){dCll = ch0.bCll-ch1.bCll;}
                else{dCll = ch1.bCll-ch1.Cl[0];}
              //  dClr = dCll;
                if (ur<0){dClr = ch0.Cl[0]-ch0.bCll;}
                else{dClr = ch0.bCll-ch1.bCll;}
                ch1.bCll += -nu1*ul*dCll-dt*ch1.getKCl(q1m,ul)*ch1.bCll;
                ch0.bCll += -nu0*ur*dClr-dt*ch0.getKCl(q1p,ur)*ch0.bCll;
                ch0.Cl_hist[ch0.n*(ch0.N+2)] = ch0.bCll;
                ch1.Cl_hist[ch1.n*(ch1.N+2)] = ch1.bCll;
            }
		}

		ch0.q_hist[ch0.idx_t(0,Ns0,ch0.n)] =  q1pfake;
		ch0.q_hist[ch0.idx_t(1,Ns0,ch0.n)] =  q2p*valveopen;
		ch1.q_hist[ch1.idx_t(0,Ns1,ch1.n)] =  q1mfake;
		ch1.q_hist[ch1.idx_t(1,Ns1,ch1.n)] =  q2m*valveopen;
	    if(ch0.P[ch0.pj(N0)]==false||ch1.P[ch1.pj(N1)]==false)ch1.P[Ns1] = false;
		else ch1.P[Ns1] = true;
		if(ch1.P[ch0.pj(N1)]==false||ch1.P[ch0.pj(N0)]==false) ch0.P[Ns0] =false;
		else ch0.P[Ns0] = true;
		ch0.p_hist[ch0.pj_t(Ns0,ch0.n)] = ch0.P[Ns0];
		ch1.p_hist[ch1.pj_t(Ns1,ch1.n)] = ch1.P[Ns1];
	}
	//reflect to get 0 flux
	else 
	{
		ch0.numFlux(q1m, q1m, q2m, -q2m, ch0.bfluxright, ch0.P[Ns0], ch1.P[Ns1]);
		ch1.numFlux(q1p, q1p, -q2p, q2p,  ch1.bfluxleft, ch0.P[Ns0], ch1.P[Ns1]);
	    ch0.q_hist[ch0.idx_t(0,Ns0,ch0.n)] =  q1m;
		ch0.q_hist[ch0.idx_t(1,Ns0,ch0.n)] =  0.;
		ch1.q_hist[ch1.idx_t(0,Ns1,ch1.n)] =  q1p;
		ch1.q_hist[ch1.idx_t(1,Ns1,ch1.n)] =  0.;
      //  printf("yooooooo valveopen =0\n");        
	}
}


Junction3::Junction3(Channel &ch0, Channel &ch1, Channel &ch2, int which0, int which1, int which2): 
			ch0(ch0), ch1(ch1), ch2(ch2),j2_01(ch0, ch1, which0, which1, 1), j2_12(ch1, ch2, which1, which2, 1),j2_21(ch2, ch1, which2, which1, 1), j2_20(ch0,ch2, which0, which2, 1)
{
	Ns[0] = ch0.N;
	Ns[1] = ch1.N;
	Ns[2] = ch2.N;
	whichend[0] = which0;
	whichend[1] = which1;
	whichend[2] = which2;
    printf("*******whichends = [%d,%d,%d]\n", which0, which1, which2, which2);

}
/**
 * this routine assumes you have one incoming (whichend =1) and two outgoing pipes (whichend =0) OR
 * two incoming (whichend =1) and two outging (whichend =0)
 * To do: set up machinery to have two incoming and one outgong but there's no reason to have anything else
**/
void Junction3::boundaryFluxes(double dt){
	double flux0[2], flux1[2], flux2[2];
	double Abar[3], Qbar[3];
	//pik is the percentage of flux from pipe k going into pipe i; should have sum_k pik = 1 for each i...
	double p01, p02, p10, p12, p20, p21;  
	p01 = 0.5;
	p02 = 0.5;
	p10 = 1-p01;
	p12 = 0.5;
	p20 = 1-p02;
	p21 = 1-p12;
	if(whichend[0]==1 &&whichend[1] ==0 &&whichend[2] ==0)
	{

		Abar[0] = .5*(ch1.q[ch1.idx(0,0)]+ch2.q[ch2.idx(0,0)]);
		Qbar[0] = .5*(ch1.q[ch1.idx(1,0)]+ch2.q[ch2.idx(1,0)]);
		Abar[1] = .5*(ch0.q[ch0.idx(0,ch0.N-1)]+ch2.q[ch2.idx(0,0)]);
		Qbar[1] = .5*(ch0.q[ch0.idx(1,ch0.N-1)]-ch2.q[ch2.idx(1,0)]);
		Abar[2] = .5*(ch0.q[ch0.idx(0,ch0.N-1)]+ch1.q[ch1.idx(0,0)]);	
		Qbar[2] = .5*(ch0.q[ch0.idx(1,ch0.N-1)]-ch1.q[ch1.idx(1,0)]);
	//solve fake Reimann problems between each pair
	//0-1 is easy
		j2_01.boundaryFluxes(dt);
		flux0[0] =  p01*ch0.bfluxright[0];
		flux0[1] =  p01*ch0.bfluxright[1];
		flux1[0] =  p10*ch1.bfluxleft[0];
		flux1[1] =  p10*ch1.bfluxleft[1];
	//1-2 is a pain in the neck	
		ch2.q[ch2.idx(1,0)]= -ch2.q[ch2.idx(1,0)];
		j2_12.boundaryFluxes(dt);	
		ch2.q[ch2.idx(1,0)]= -ch2.q[ch2.idx(1,0)];
		flux1[0] += p12*ch1.bfluxleft[0];
		flux1[1] += p12*ch1.bfluxleft[1];

		flux2[0] =  p21*ch2.bfluxleft[0];///formerly these two lines were bflux right. Wrong before?
		flux2[1] =  p21*ch2.bfluxleft[1];

		ch1.q[ch1.idx(1,0)]= -ch1.q[ch1.idx(1,0)];
		j2_21.boundaryFluxes(dt);
		ch1.q[ch1.idx(1,0)]= -ch1.q[ch1.idx(1,0)];
		flux2[0] = p21*ch2.bfluxleft[0];
		flux2[1] = p21*ch2.bfluxleft[1];
	

		j2_20.boundaryFluxes(dt);
		flux0[0] += p02*ch0.bfluxright[0];
		flux0[1] += p02*ch0.bfluxright[1];
		flux2[0] += p20*ch2.bfluxleft[0];
		flux2[1] += p20*ch2.bfluxleft[1];

		//set all the fluxes properly
		ch0.bfluxright[0] = flux0[0];
		ch0.bfluxright[1] = flux0[1];
		ch1.bfluxleft[0] = flux1[0];
		ch1.bfluxleft[1] = flux1[1];
		ch2.bfluxleft[0] = flux2[0];
		ch2.bfluxleft[1] = flux2[1];
       /* printf("triple junction fluxes are\n");
        printf("channel 0: [%f  %f]\n",flux0[0], flux0[1]);
        printf("channel 1: [%f  %f]\n",flux1[0], flux1[1]);
        printf("channel 2: [%f  %f]\n",ch2.bfluxleft[0], ch2.bfluxleft[1]);*/
		//store away the info
		ch0.q_hist[ch0.idx_t(0,ch0.N+1,ch0.n)] = Abar[0]; 
		ch0.q_hist[ch0.idx_t(1,ch0.N+1,ch0.n)] = Qbar[0];
		ch1.q_hist[ch1.idx_t(0,0,ch1.n)] =  Abar[1];
		ch1.q_hist[ch1.idx_t(1,0,ch1.n)] =  Qbar[1];
		ch2.q_hist[ch2.idx_t(0,0,ch2.n)] =  Abar[2];
		ch2.q_hist[ch2.idx_t(1,0,ch2.n)] =  Qbar[2];
	}
    else if( (whichend[0] ==0&& whichend[1] ==1 &&whichend[2] ==1))
	{

		Abar[0] = .5*(ch1.q[ch1.idx(0,ch1.N-1)]+ch2.q[ch2.idx(0,ch2.N-1)]);
		Qbar[0] = .5*(ch1.q[ch1.idx(1,ch1.N-1)]+ch2.q[ch2.idx(1,ch2.N-1)]);
		Abar[1] = .5*(ch0.q[ch0.idx(0,0)]+ch2.q[ch2.idx(0,ch2.N-1)]);
		Qbar[1] = .5*(ch0.q[ch0.idx(1,0)]-ch2.q[ch2.idx(1,ch2.N-1)]);
		Abar[2] = .5*(ch0.q[ch0.idx(0,0)]+ch1.q[ch1.idx(0,ch1.N-1)]);	
		Qbar[2] = .5*(ch0.q[ch0.idx(1,0)]-ch1.q[ch1.idx(1,ch1.N-1)]);
	//solve fake Reimann problems between each pair
	//0-1 is easy
		j2_01.boundaryFluxes(dt);
		flux0[0] =  p01*ch0.bfluxleft[0];
		flux0[1] =  p01*ch0.bfluxleft[1];
		flux1[0] =  p10*ch1.bfluxright[0];
		flux1[1] =  p10*ch1.bfluxright[1];
	//1-2 is a pain in the neck	
		ch2.q[ch2.idx(1,ch2.N-1)]= -ch2.q[ch2.idx(1,ch2.N-1)];
		j2_12.boundaryFluxes(dt);	
		ch2.q[ch2.idx(1,ch2.N-1)]= -ch2.q[ch2.idx(1,ch2.N-1)];
		flux1[0] += p12*ch1.bfluxright[0];
		flux1[1] += p12*ch1.bfluxright[1];
		flux2[0] =  p21*ch2.bfluxright[0];
		flux2[1] =  p21*ch2.bfluxright[1];

		ch1.q[ch1.idx(1,ch1.N-1)]= -ch1.q[ch1.idx(1,ch1.N-1)];
		j2_21.boundaryFluxes(dt);
		ch1.q[ch1.idx(1,ch1.N-1)]= -ch1.q[ch1.idx(1,ch1.N-1)];
		flux2[0] = p21*ch2.bfluxright[0];
		flux2[1] = p21*ch2.bfluxright[1];
	

		j2_20.boundaryFluxes(dt);
		flux0[0] += p02*ch0.bfluxleft[0];
		flux0[1] += p02*ch0.bfluxleft[1];
		flux2[0] += p20*ch2.bfluxright[0];
		flux2[1] += p20*ch2.bfluxright[1];

		//set all the fluxes properly
		ch0.bfluxleft[0] = flux0[0];
		ch0.bfluxleft[1] = flux0[1];
		ch1.bfluxright[0] = flux1[0];
		ch1.bfluxright[1] = flux1[1];
		ch2.bfluxright[0] = flux2[0];
		ch2.bfluxright[1] = flux2[1];
/*        printf("triple junction fluxes are\n");
        printf("channel 0: [%f  %f]\n",flux0[0], flux0[1]);
        printf("channel 1: [%f  %f]\n",flux1[0], flux1[1]);
        printf("channel 2: [%f  %f]\n",flux2[0], flux2[1]);*/
		//store away the info
		ch0.q_hist[ch0.idx_t(0,0,ch0.n)] = Abar[0]; 
		ch0.q_hist[ch0.idx_t(1,0,ch0.n)] = Qbar[0];
		ch1.q_hist[ch1.idx_t(0,ch1.N+1,ch1.n)] =  Abar[1];
		ch1.q_hist[ch1.idx_t(1,ch1.N+1,ch1.n)] =  Qbar[1];
		ch2.q_hist[ch2.idx_t(0,ch2.N+1,ch2.n)] =  Abar[2];
		ch2.q_hist[ch2.idx_t(1,ch2.N+1,ch2.n)] =  Qbar[2];

    }
	else{	
		cout<<"End 0 = "<<whichend[0]<<" End 1 = "<< whichend[1]<<" End 2 = "<<whichend[2]<<endl;
		cout<<"Have not implemented triple junctions for this configuration!Sorry!\n";
	}
}




/*
//keep around for comparison purposes, maybe?
double Cpreiss::hofAold(double A)
{
		double theta = getTheta(A,D);		
		double y = D/2.*(1.+cos(PI-theta/2.));
		return y;
}
*/
//if you want to use negative slot...?
double Cpreiss::fakehofA(double A, bool p)
{
	if (A>=At){p=true;}
	double y;
	if (p){y = yt+(A-At)/Ts;}	
	else {y = HofA(A,p);}
	return y;
}

double Cpreiss::fakeAofh(double h, bool p)
{
	if (!p){return AofH(h,p);}
	else{return (h-yt)*Ts+At;}	
}


