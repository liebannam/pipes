
#ifndef OPTDETIRTUS_H
#define OPTDETRITUS_H

#include "network.h"
#include "levmar.h"

void getTimeSeries(vector<Real> & bvals, vector<Real> &x, int m, int M, double T, int Fourier)
       	
//bvals is time series of M+1 values at t =0, T/M, ...T
//x is a length m vector. either contains Fourier modes (Fourier==1) or of Hermite spline components (Fourier ==0)
//we have m/2 values in Fourier series or m/2 spline interpolation points. Interpolation dt = 1/(m/2-1)
//
{
	if(Fourier)	//Discrete Fourier modes:
	//x[k] = a_k (k<=m/2)
	//x[k] = b_j (where j =k-m/2, for k>m/2)
	//      x =     [a_0, a_1, ...a_m/2, b_1,   b2   ...b_(m/2-1)] 
	//                |     |       |     |     |
	//index, k        0     1      m/2  m/2+1  m/2+2       m-1     
	{
		double t;
		for(int nn = 0; nn<M+1; nn++)
		{
			bvals[nn] = 0.5*x[0];
			t = (double)(nn)/(double)M;
			for (int k = 1; k<m/2; k++)
			{
				bvals[nn] += x[k]*cos(2.*PI*(double)k*t) + x[k+m/2]*sin(2.*PI*(double)k*t);
		//		cout<<k<< " "<< x[k+m/2+1]<<endl;
			}
			bvals[nn] +=0.5*x[m/2]*cos(PI*(double)m*t);
		//	printf("n = %d, t/T = %f, gah %f\n", n, t, bvals[n]);
		//	cout<<"!!!!!"<<bvals[n]<<endl;
		}
	}
	else//Hermite Spline
	//Hermite spline evaluation, adapted from C. Rycroft's code. 
	// Evaluates spline x at point t, where
	// x = [f(t0), f'(t0), f(t1), f'(t1), ...f(t_m), f'(tm)]
	//      |       |        |      |          |       |
	// i=   0       1        2      3          2*m    2*m+1
	// x0 = 0, xm = T

	{
		double t, Dt,dt;
		dt = T/M;//time series dt
		Dt = T/(m/2-1.);//spline dt. Careful! This is NOT the dt in your problem!
		for(int nn = 0; nn<M+1; nn++)
		{
			t= nn*dt;
			int j = int((t)/Dt);
			t = t/Dt-j; 
		      // cout<<Dt<<endl;	
			if(j>m+1 || j<0) printf("warning! t out of range!!");
			// Calculate square and cube, and pointer to the values to use
			double t2=t*t,t3=t2*t;
			// Calculate the value of the spline function
			if(j==m){bvals[nn] = x[2*j];}
			else{bvals[nn] = x[2*j]*(2*t3-3*t2+1)+x[2*j+1]*(t3-2*t2+t)*Dt+x[2*j+2]*(-2*t3+3*t2)+x[2*j+3]*(t3-t2)*Dt;}
		}
	}
}



//trial class to implement optimization machinery for adjusting BCs at a single node to minimize average pressure gradient (IS THIS HARDER??)
class bc_opt_dh: public levmar{
public:
	//x = [a_0, a_1,...a_m/2, b_1, ,,,,b_m/2-1] 2m values for discrete trig transform 	
        vector<Real> x0;
	vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	vector<double> bvals;//M+1 list of actual time series values
	int whichnode; //at which node we're applying this BC time series...
	Network Ntwk;
	int M;                //number of time steps
	int modetype;         //1 - Fourier   0- Hermite interpolation 
	double T;
	double dt;
	double mydelta; //for finite diff approx for J
	bc_opt_dh(int n, int M_, vector<double>x0_, Network Ntwk_i, int modetype_, double T_, int whichnode_):
		levmar(M_+1,n), x0(n),bvals(M_+1), Ntwk(Ntwk_i), modetype(modetype_), whichnode(whichnode_) 
	{

		//x = [x_0,...x_m-1 ] d.o.f in Fourier space
		//r = [r_0,...r_M] = residuals at times t=0, dt.. , M*dt
		M = Ntwk.M;
		dt = T/(double)M;
		a0.resize(Ntwk_i.channels.size());
		q0.resize(Ntwk_i.channels.size());
		for (int i = 0; i<Ntwk_i.channels.size(); i++)
		{
			int Ni =Ntwk_i.channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i.channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i.channels[i]->q0[j+Ni]);
			}
		}
		//all of this information needs to be standardized/automated!
		mydelta = 1e-6;
		for (int i=0; i<n; i++)
		{
			x0[i] = x0_[i];
			x[i] = x0[i];
		}
		
		getTimeSeries(bvals, x0, n,M,T,modetype);
		Ntwk.junction1s[whichnode]->setbVal(bvals);
	}
	void compute_r();
	void compute_J();

};

void bc_opt_dh::compute_r()
{	 
	
	for(int i=0; i<Ntwk.Nedges; i++){
		Ntwk.channels[i]->setq(a0[i],q0[i]);
		Ntwk.channels[i]->n = 0;
	}
	Ntwk.nn = 0;
	getTimeSeries(bvals, x, n,M,T,modetype);
//	for(int i =0; i<bvals.size(); i++){
//		cout<<bvals[i]<<endl;
//	}
	Ntwk.junction1s[whichnode]->setbVal(bvals);
	
	Ntwk.runForwardProblem(dt);
	for(int i=0; i<M+1; i++)
	{
		r[i] = 0;
		for(int k = 0; k<Ntwk.Nedges; k++)
		{
			r[i] += Ntwk.channels[k]->getAveGradH(i);
		}
	}
}

void bc_opt_dh::compute_J()
{
	J.reset_values(0.0);
	compute_r();
	vector<Real> r0(r);
	//Ntwk.junction1s[0]->setbVal(x);
	//printf("bvals are [%f, %f]\n",Ntwk.junction1s[0]->bval[0],Ntwk.junction1s[1]->bval[0]);
	for(int j = 0; j<n; j++)
	{
		
		x[j]+=mydelta;

		for(int i=0; i<Ntwk.Nedges; i++)
		{
			Ntwk.channels[i]->setq(a0[i],q0[i]);
			Ntwk.channels[i]->n = 0;
		}
		Ntwk.nn = 0;
		getTimeSeries(bvals, x, n,M,T,modetype);	
		Ntwk.junction1s[whichnode]->setbVal(bvals);	
		Ntwk.runForwardProblem(dt);
		for(int i=0; i<M+1; i++)
		{
			double r=0;
			for(int k = 0; k<Ntwk.Nedges; k++)
			{r += Ntwk.channels[k]->getAveGradH(i);}
			J(i,j) = (r-r0[i])/mydelta;
		}
		x[j] -= mydelta;	
	}	
}
	


#endif
