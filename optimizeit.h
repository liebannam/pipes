/**
 * \file optimizeit.h
 * \brief optimizeit.h documentation
 * contains classes for optimization via Levenberg-Marquardt
 * Class declarations AND definitions for:
 * --bc_opt_dh_c (minimize <dH/dx> by changing boundary conditions at a single node, with consrtained inflow volume
 * --mystery_bc (recover boundary conditions at a single node by minimizing deviation from a known pressure time series)
 * --bc1_opt_dc (minimize <dH/dx> by changing boundary conditions at a single node)
 * --bc_opt_dc (minimize <dH/dx> by changing boundary conditions at multiple nodes)
 * --opt_eq_outflow (attempt to equalize outflow)(not presently supported)
 *Note: the first four have many fairly redudant aspects and this could be more elegantly done with a base class for optimization at boundary nodes, with the derived classes being much simpler. But the work's already done and it's working.
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


#ifndef OPTIMIZEIT_H
#define OPTIMIZEIT_H

#include "levmar.h"
#include "setupandrun.h"
#ifdef __OMP
#include "omp.h"
#endif


/*
 * \brief Class to minimize <dH/dx> by changing boundary conditions at a single node, with constrained inflow volume
 * First component is set such that total inflow volume over simulation period = Vin
 *   --Fourier: constant mode a_0 
 *   --Hermite: first value x_0 
*/
class bc_opt_dh_c:public levmar{ 
public:
	//x = [a_0, a_1,...a_m/2, b_1, ,,,,b_m/2-1] 2m values for discrete trig transform 	
	int Ndof; //number of degrees of freedom per node
	vector<Real> x0; //size is Nn*ndof
	vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	//vector<double> bvals;//M+1 list of actual time series values
	Network Ntwk;
	int M;                //number of time steps	
	int modetype;         //1 - Fourier (must be this, or else...)   
	int whichnode;        //at which node we're applying this BC time series...
	double Vin;			  //desired total flow volutme = integral_0^T(Q(0,t))dt
	double T;
	double dt;
	double xfake0;
	double mydelta; //for finite diff approx for J
	bc_opt_dh_c(int n, int M_, vector<double>x0_, Network *Ntwk_i, double T_, int whichnode, double Vin_, int modetype_):
		levmar(M_+1,n), x0(n),Ntwk(Ntwk_i),whichnode(whichnode), Vin(Vin_)
	{

		//x = [x_1,...x_m-1 ] d.o.f in Fourier space
		//note that x0 = constant term is determined by volume constraint
		//r = [r_0,...r_M] = residuals at times t=0, dt.. , M*dt
		T = T_;
		modetype = modetype_;
		Ndof = x0.size();
		M = Ntwk.M;
		dt = T/(double)M;
		a0.resize(Ntwk_i->channels.size());
		q0.resize(Ntwk_i->channels.size());
		for (int i = 0; i<Ntwk_i->channels.size(); i++)
		{
			int Ni =Ntwk_i->channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i->channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i->channels[i]->q0[j+Ni]);
			}
		}
		mydelta = 1e-6;
		for (int i=0; i<n; i++)
		{
			x0[i] = x0_[i];
			x[i] = x0[i];
		if (WTF){printf("i = %d, x0 = %f\n ", i,x0[i]);}
		}

		vector <Real> bvals(M+1);
		setBCTimeSeries(x0,bvals);
		Ntwk.junction1s[whichnode]->setbVal(bvals);
	}
	void compute_r();
	void compute_J();
	void compute_r(vector<Real> &rr, vector <Real>&xx);
	void setBCTimeSeries(vector<Real>xx, vector <Real> &bvals);

};

void bc_opt_dh_c::compute_r()
{
	compute_r(r,x);	

}

void bc_opt_dh_c::compute_r(vector<Real> &rr, vector <Real>&xx)
{	 
	
	Network Ntwk0(Ntwk);
	for(int i=0; i<Ntwk0.Nedges; i++)
	{
		Ntwk0.channels[i]->setq(a0[i],q0[i]);
		Ntwk0.channels[i]->n = 0;
	}
	Ntwk0.nn = 0;
	vector <Real> bvals(M+1);
	setBCTimeSeries(xx,bvals);
	Ntwk0.junction1s[whichnode]->setbVal(bvals);

	Ntwk0.runForwardProblem(dt);
	for(int i=0; i<M+1; i++)
	{
		rr[i] = sqrt(Ntwk0.getAveGradH(i));
	}
}

void bc_opt_dh_c::setBCTimeSeries(vector<Real>xx, vector <Real> &bvals)
{
	vector <Real> xfake(Ndof+1);
	for (int k = 0; k<Ndof; k++)
	{
		xfake[k+1] = xx[k];
	}
	if (modetype==1){xfake0 = 2.*Vin/T;}//fourier constant mode = 2Vin/T
	else//set Q(0,0) value
	{	
		double Dt =  2.*T/((double)Ndof-1);//spline dt. 
		double sumQ = 0;
		for (int k=0;k<Ndof/2; k++)
		{
			sumQ+=x0[2*k+1];
		}
		sumQ -= x0[Ndof-2]/2.;
		xfake0 =  2*(Vin/Dt-sumQ-Dt/12.*(x[0]-x[Ndof-1]));
        //printf("sumQ = %f, x[0] = %f, Dt = %f\n", sumQ, xfake[0], Dt);
	}
	xfake[0] = xfake0;
	getTimeSeries(bvals, xfake, Ndof+1,M,T,modetype);
}

void bc_opt_dh_c::compute_J()
{
	J.reset_values(0.0);

//	compute_r();
	//Ntwk.junction1s[0]->setbVal(x);
//	printf("bvals are [%f, %f]\n",Ntwk.junction1s[0]->bval[0],Ntwk.junction1s[1]->bval[0]);
#pragma omp parallel for
	for(int j = 0; j<n; j++)
	{

	//	int NT = omp_get_num_threads();
//		printf("number of threads is %d\n",NT);
        vector<Real> rp(M+1);
		vector<Real> rm(M+1);
		vector<Real> x2(x);
		x2[j]+=mydelta;
		compute_r(rp,x2);	
		x2[j] -=2*mydelta;
		compute_r(rm, x2);
		for(int i=0; i<M+1; i++)
		{
			J(i,j) = (rp[i]-rm[i])/(2.*mydelta);
		}
	}	
}


/* \brief Class to recover boundary conditions at a single node by minimizing deviation from a known pressure time series*/
class mystery_bc: public levmar{
	public:
		vector<Real> x0;
		vector<double> hdata;
		vector<double> qfixed; //fixed boundary value time series (only 0,...delay and M+1-delay, ...M+1 matter)
		vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	   	//assume we have height data for x=xstar in pipe pj 
		int pj;
		double xstar;
		int Nstar;
		int whichnode;      //node we are varying time series at
		int nodetype;       //type of boundary conditions at that node
		Network Ntwk;
		int M;                //total number of time steps
		int delay;            //delay is # steps to let information propogate. m = M+1-2*delay is number of elements of residual
		int modetype;         //1 - Fourier   0- Hermite interpolation 
		double T;
		double dt;
		double mydelta; //for finite diff approx for J
		mystery_bc(int n, int M_, vector<double>x0_, vector<double> hdata_, Network *Ntwk_i,int modetype_ ,double T_, int pj_, double xstar_, int whichnode_, vector<double> qfixed_, int delay_=0):levmar(M_+1-2*delay_,n), hdata(hdata_), M(M_), pj(pj_), xstar(xstar_),x0(n),whichnode(whichnode_),Ntwk(Ntwk_i),T(T_)
	{
		modetype = modetype_;
		delay = delay_;
		dt= T/(double)M;
		mydelta = 1e-6;
		Nstar = (int) floor(Ntwk.channels[pj]->N*xstar/Ntwk.channels[pj]->L);
		printf("Nstar = %d and pj = %d and xstar = %f\n", Nstar, pj, xstar);
		for (int i=0; i<x0.size(); i++)
		{
			x0[i] = x0_[i];
			x[i] = x0_[i];
		}
		qfixed.resize(M+1);
		for (int i = 0; i<M+1; i++)
		{
			qfixed[i] = qfixed_[i];
		}
		a0.resize(Ntwk_i->channels.size());
		q0.resize(Ntwk_i->channels.size());
		for (int i = 0; i<Ntwk_i->channels.size(); i++)
		{
			int Ni =Ntwk_i->channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i->channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i->channels[i]->q0[j+Ni]);
			}
		}
		vector <Real> bvals(M+1);
		vector <Real> b(m);
		getTimeSeries(b, x0, n, m, T, modetype);
		for (int i = 0; i<M+1; i++)
			{
				if(i<delay || i>M+1-delay)
				{	
					bvals[i] = qfixed[i];
				}
				else
				{
					bvals[i+delay] = b[i];
				}
			}
		Ntwk.junction1s[whichnode]->setbVal(bvals);
	}
	void compute_r();
	void compute_J();
	void compute_r(vector<Real> &rr, vector <Real>&xx);
};

void mystery_bc::compute_r()
{
	compute_r(r,x);	
}

void mystery_bc::compute_r(vector<Real> &rr, vector <Real>&xx)
{	 
	Network Ntwk0(Ntwk);
	for(int i=0; i<Ntwk0.Nedges; i++){
		Ntwk0.channels[i]->setq(a0[i],q0[i]);
		Ntwk0.channels[i]->n = 0;
	}
	Ntwk0.nn = 0;
	vector <Real> bvals(M+1);
	vector <Real> b(m);
	getTimeSeries(b, xx, n, m, T, modetype);
	for (int i = 0; i<M+1; i++)
	{
		if(i<delay || i>M+1-delay)
		{	
			bvals[i] = qfixed[i];
		}
		else
		{
			bvals[i+delay] = b[i];
		}
	}

	Ntwk0.junction1s[whichnode]->setbVal(bvals);
	Ntwk0.runForwardProblem(dt);
	for(int i=0; i<m; i++)
	{
		double a = Ntwk0.channels[pj]->q_hist[Ntwk0.channels[pj]->idx_t(0,Nstar,i+delay)]; 
		double H =  Ntwk.channels[pj]->pbar(a,false);
		rr[i] =H-hdata[i+delay];
	//	printf("i= %d, hdata = %f  r =  %f  a = %f\n",i,hdata[i],rr[i], a); 
	}
}

void mystery_bc::compute_J()
{

	J.reset_values(0.0);

#pragma omp parallel for
	for(int j = 0; j<n; j++)
	{

		int NT = omp_get_num_threads();
		printf("number of threads is %d\n",NT);
        vector<Real> rp(m);
		vector<Real> rm(m);
		vector<Real> x2(x);
		x2[j]+=mydelta;
		compute_r(rp,x2);	
		x2[j] -=2*mydelta;
		compute_r(rm, x2);
		for(int i=0; i<m; i++)
		{
			J(i,j) = (rp[i]-rm[i])/(2.*mydelta);
		}
	}	

}


/*\brief Class to minimize <dH/dx> by changing boundary conditions at a single node)*/
class bc1_opt_dh: public levmar{
public:
	//x = [a_0, a_1,...a_m/2, b_1, ,,,,b_m/2-1] 2m values for discrete trig transform 	
	vector<Real> x0;
	vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	vector<double> bvals;//M+1 list of actual time series values
	Network Ntwk;
	int M;                //number of time steps
	int whichnode; //at which node we're applying this BC time series...
	int modetype;         //1 - Fourier   0- Hermite interpolation 
	double T;
	double dt;
	double mydelta; //for finite diff approx for J
	bc1_opt_dh(int n, int M_, vector<double>x0_, Network *Ntwk_i, int modetype_, double T_, int whichnode_):
		levmar(M_+1,n), x0(n),bvals(M_+1), Ntwk(Ntwk_i), whichnode(whichnode_), modetype(modetype_)
	{

		//x = [x_0,...x_m-1 ] d.o.f in Fourier space
		//r = [r_0,...r_M] = residuals at times t=0, dt.. , M*dt
		T = T_;
		M = Ntwk.M;
		dt = T/(double)M;
		a0.resize(Ntwk_i->channels.size());
		q0.resize(Ntwk_i->channels.size());
		for (int i = 0; i<Ntwk_i->channels.size(); i++)
		{
			int Ni =Ntwk_i->channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i->channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i->channels[i]->q0[j+Ni]);
			}
		}
		bvals.resize(M+1);
		//all of this information needs to be standardized/automated!
		mydelta = 1e-6;
		for (int i=0; i<n; i++)
		{
			x0[i] = x0_[i];
			x[i] = x0[i];
			printf("i = %d, x0 = %f\n ", i,x0[i]);
		}
		getTimeSeries(bvals, x0, n,M,T,modetype);
		//printf("T = %f, n = %d, M = %d, modetype = %d", T,n,M,modetype);
		Ntwk.junction1s[whichnode]->setbVal(bvals);
	}
	void compute_r();
	void compute_J();

};

void bc1_opt_dh::compute_r()
{	 
	
	for(int i=0; i<Ntwk.Nedges; i++){
		Ntwk.channels[i]->setq(a0[i],q0[i]);
		Ntwk.channels[i]->n = 0;
	}
	Ntwk.nn = 0;
	getTimeSeries(bvals, x, n,M,T,modetype);
	
	cout<<"here?\n";
	for(int i=0;i<n; i++){
		cout<<x[i]<<endl;
	}
	cout<<"OR HERE?!"<<endl;	
	for(int i =0; i<bvals.size(); i++){
		cout<<bvals[i]<<endl;
	}
	Ntwk.junction1s[whichnode]->setbVal(bvals);
	
	Ntwk.runForwardProblem(dt);
	for(int i=0; i<M+1; i++)
	{
		r[i] = dt*Ntwk.getAveGradH(i);
	//	printf("r[%d] = %f\n", i, r[i]);
	}
		//r[i] = 0;
		//for(int k = 0; k<Ntwk.Nedges; k++)
		//{r[i] += dt*Ntwk.channels[k]->getAveGradH(i);}
	
}

void bc1_opt_dh::compute_J()
{
	J.reset_values(0.0);

	vector<Real> rp(M+1);
	vector<Real> rm(M+1);
//	compute_r();
	//Ntwk.junction1s[0]->setbVal(x);
	//printf("bvals are [%f, %f]\n",Ntwk.junction1s[0]->bval[0],Ntwk.junction1s[1]->bval[0]);
	for(int j = 0; j<n; j++)
	{
		
		x[j]+=mydelta;
		compute_r();
		rp.swap(r);
		x[j] -=2*mydelta;
		compute_r();
		rm.swap(r);

	//	for(int i=0; i<Ntwk.Nedges; i++)
	//	{
	//		Ntwk.channels[i]->setq(a0[i],q0[i]);
	//		Ntwk.channels[i]->n = 0;
	//	}
	//	Ntwk.nn = 0;
	//	getTimeSeries(bvals, x, n,M,T,modetype);	
	//	Ntwk.junction1s[whichnode]->setbVal(bvals);	
	//	Ntwk.runForwardProblem(dt);
		for(int i=0; i<M+1; i++)
		{
			//double r=0;
			//for(int k = 0; k<Ntwk.Nedges; k++)
			//{r += dt*Ntwk.channels[k]->getAveGradH(i);}	
	//		double r = dt*Ntwk.getAveGradH(i);

			//printf("r[%d] = %f r0[%d] = %f\n", i, rp[i],i, rm[i]);
			J(i,j) = (rp[i]-rm[i])/(2.*mydelta);
		}
		x[j] += mydelta;	
	}	
	compute_r();
}

//same as above but consider multiple 1 junctions at once..
class bc_opt_dh:public levmar{
public:
	//x = [a_0, a_1,...a_m/2, b_1, ,,,,b_m/2-1] 2m values for discrete trig transform 	
        int Nn;  //number of nodes where bcs are being varied
	int Ndof; //number of degrees of freedom per node
	vector<Real> x0; //size is Nn*ndof
	vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	//vector<double> bvals;//M+1 list of actual time series values
	Network Ntwk;
	int M;                //number of time steps
	int skip;             //number of time steps skipped before checking r
	int modetype;         //1 - Fourier   0- Hermite interpolation 
	vector <int> whichnodes; //at which node we're applying this BC time series...
	double T;
	double dt;
	double mydelta; //for finite diff approx for J
	bc_opt_dh(int n, int M_, vector<double>x0_, Network *Ntwk_i, int modetype_, double T_, vector<int> whichnodes_, int skip_):
		levmar(M_+1,n), x0(n),Ntwk(Ntwk_i), modetype(modetype_),whichnodes(whichnodes_)	
	{

		//x = [x_0,...x_m-1 ] d.o.f in Fourier space
		//r = [r_0,...r_M] = residuals at times t=0, dt.. , M*dt
		T = T_;
		printf("did we make it in yet?");
		Nn = whichnodes.size();
		Ndof = x0.size()/Nn;  //eeek this could go SO HORRIFICALLY WRONG AHHHHHHHHHHH
		M = Ntwk.M;
		skip = skip_;
		dt = T/(double)M;
		a0.resize(Ntwk_i->channels.size());
		q0.resize(Ntwk_i->channels.size());
		for (int i = 0; i<Ntwk_i->channels.size(); i++)
		{
			int Ni =Ntwk_i->channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i->channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i->channels[i]->q0[j+Ni]);
			}
		}
		vector <Real> bvals(M+1);
		//all of this information needs to be standardized/automated!
		mydelta = 1e-6;
		for (int i=0; i<n; i++)
		{
			x0[i] = x0_[i];
			x[i] = x0[i];
			printf("i = %d, x0 = %f\n ", i,x0[i]);
		}

		vector <Real> xfake(Ndof);
		for (int i = 0; i<Nn; i++)
		{
			for (int k = 0; k<Ndof; k++)
			{
				xfake[k] = x0[i*(Ndof)+k];
			}
			getTimeSeries(bvals, xfake, Ndof,M,T,modetype);
			Ntwk.junction1s[whichnodes[i]]->setbVal(bvals);
		}
	}
	void compute_r();
	void compute_J();
	void compute_r(vector<Real> &rr, vector <Real>&xx);

};


void bc_opt_dh::compute_r()
{
	compute_r(r,x);	

}

void bc_opt_dh::compute_r(vector<Real> &rr, vector <Real>&xx)
{	 
	
	Network Ntwk0(Ntwk);
	for(int i=0; i<Ntwk0.Nedges; i++){
		Ntwk0.channels[i]->setq(a0[i],q0[i]);
		Ntwk0.channels[i]->n = 0;
	}
	Ntwk0.nn = 0;

	vector <Real> xfake(Ndof);
	vector <Real> bvals(M+1);
//	cout<<"here? (this is x)\n";	
//	for(int i=0;i<n; i++){
//			cout<<x[i]<<endl;
//	}
	for (int i = 0; i<Nn; i++)
	{
		//cout<<"Ndof ="<<Ndof<<endl;
	//	cout<<"which node? "<<whichnodes[i]<<endl;
		for (int k = 0; k<Ndof; k++)
		{
			xfake[k] = xx[i*(Ndof)+k];
		//	cout<<"k="<<k<<" xfake = "<<xfake[k]<<endl;
		}
		getTimeSeries(bvals, xfake, Ndof,M,T,modetype);
		Ntwk0.junction1s[whichnodes[i]]->setbVal(bvals);
	//	cout<<"OR HERE?! (this is bvals)"<<endl;	
	//	for(int i =0; i<bvals.size(); i++){
	//		cout<<bvals[i]<<endl;
	//	}
	}
	Ntwk0.runForwardProblem(dt);
	//for(int i=0; i<M/skip+1; i++)
	for(int i=0; i<M+1; i++)
	{
	//fuck this skip noise
	//	rr[i] = sqrt((double)skip/((float)M)*Ntwk0.getAveGradH(i*skip));
		rr[i] = sqrt(Ntwk0.getAveGradH(i));
	//	printf("r[%d] = %f\n", i, r[i]);
	}
}


void bc_opt_dh::compute_J()
{
	J.reset_values(0.0);

//	compute_r();
	//Ntwk.junction1s[0]->setbVal(x);
	printf("bvals are [%f, %f]\n",Ntwk.junction1s[0]->bval[0],Ntwk.junction1s[1]->bval[0]);
#pragma omp parallel for
	for(int j = 0; j<n; j++)
	{

	//	int NT = omp_get_num_threads();
	//	printf("number of threads is %d\n",NT);
        vector<Real> rp(M+1);
		vector<Real> rm(M+1);
		vector<Real> x2(x);
		x2[j]+=mydelta;
		compute_r(rp,x2);	
		x2[j] -=2*mydelta;
		compute_r(rm, x2);
		for(int i=0; i<M/skip+1; i++)
		{
			J(i,j) = (rp[i]-rm[i])/(2.*mydelta);
		}
//		x[j] += mydelta;	
	}	
//	compute_r();
}


/*updating this to be less shitty! (11/13/2015)*/
class opt_eq_outflow:public levmar{
public:
	//x = [x_0...x_n-1] values to describe curves for valves 	
    int Nn;  //number of nodes where bcs are being varied
	int Ndof; //number of degrees of freedom per node
	int No; //number of outflow nodes under consideration (Mm = No)
	vector<Real> x0; //size is Ndof*Nnodes = Mm
	vector< vector <double> >  a0, q0; //a0[i][j] is the jth value in edge i...I think?
	//vector<double> bvals;//M+1 list of actual time series values
	vector <int> whichnodes; //at which nodes we're applying this BC time series...
	vector <int> whichnodes_out; //at which nodes we're tracking outflow
	Network Ntwk;
	int M;                //number of time steps
	int modetype;         //1 - Fourier   0- Hermite interpolation 
	double T;
	double dt;
	double mydelta; //for finite diff approx for J
	opt_eq_outflow(int n, int No_, vector<double>x0_, Network *Ntwk_i, double T_, vector<int> whichnodes_, vector<int> whichnodes_out_):levmar(No_,n),No(No_), x0(n),whichnodes(whichnodes_),
		whichnodes_out(whichnodes_out_),Ntwk(Ntwk_i)
	{

		//x = [x_0,...x_m-1 ] d.o.f in Fourier space
		//r = [r_0,...r_M] = residuals at times t=0, dt.. , M*dt
		T = T_;
		Nn = whichnodes.size();
		Ndof = x0.size()/Nn;  //eeek this could go SO HORRIFICALLY WRONG AHHHHHHHHHHH
		M = Ntwk.M;
		dt = T/(double)M;
		a0.resize(Ntwk_i->channels.size());
		q0.resize(Ntwk_i->channels.size());
		for (int i = 0; i<Ntwk_i->channels.size(); i++)
		{
			int Ni =Ntwk_i->channels[i]->N; 
			for(int j = 0; j<Ni; j++)
			{
				a0[i].push_back(Ntwk_i->channels[i]->q0[j]);
				q0[i].push_back(Ntwk_i->channels[i]->q0[j+Ni]);
			}
		}
		vector <Real> vvals(M+1);
		//all of this information needs to be standardized/automated!
		mydelta = 1e-6;
		for (int i=0; i<n; i++)
		{
			x0[i] = x0_[i];
			x[i] = x0[i];
			printf("i = %d, x0 = %f\n ", i,x0[i]);
		}

		vector <Real> xfake(Ndof);
		for (int i = 0; i<Nn; i++)
		{
			for (int k = 0; k<Ndof; k++)
			{
				xfake[k] = x0[i*(Ndof)+k];
			}
			evaluateit1(vvals, xfake, M,T);
			Ntwk.junction2s[whichnodes[i]]->setValveTimes(vvals);
		}
	}
	void compute_r();
	void compute_r(vector<Real> &rr, vector<Real>&xx);
	void compute_J();
	void evaluateit1(vector<Real> &v, vector<Real> x, int M, double T);
};

void opt_eq_outflow::evaluateit1(vector<Real> &v, vector<Real> x, int M, double T)//evaluate time series given by sin^2(x[0]*ti+x[1]), ti  = i*T/M, i = 0,2,...M
{
	double t =0;
	double dt = T/(double)M;
	for (int i = 0; i<M+1; i++)
	{
		t = (double)i*dt;
		double s = sin(x[0]*t/T+x[1]);
		v[i] = s*s;
	//	v[i] =1.;
	//	cout<<i<<" "<<v[i]<<endl;
		void compute_r();
	void compute_J();
	void compute_r(vector<Real> &rr, vector <Real>&xx);
}
}


void opt_eq_outflow::compute_r()
{
	compute_r(r, x);
}

void opt_eq_outflow::compute_r(vector<Real> &rr, vector <Real>&xx)
{

	Network Ntwk0(Ntwk);
	for(int i=0; i<Ntwk.Nedges; i++){
		Ntwk0.channels[i]->setq(a0[i],q0[i]);
		Ntwk0.channels[i]->n = 0;
	}
	Ntwk0.nn = 0;

	vector <Real> xfake(Ndof);
	vector <Real> vvals(M+1);
	for (int i = 0; i<Nn; i++)
	{
		for (int k = 0; k<Ndof; k++)
		{
			xfake[k] = xx[i*(Ndof)+k];
		}
		evaluateit1(vvals, xfake, M,T);
	//	for(int j = 0; j<vvals.size(); j+=10)cout<<vvals[j]<<endl;
		Ntwk0.junction2s[whichnodes[i]]->setValveTimes(vvals);
	}

	Ntwk0.runForwardProblem(dt);
	double Qbar = 0;
	for (int i =0; i<No; i++)
	{
		double Qi = Ntwk0.junction1s[whichnodes_out[i]]->getFlowThrough(dt);
		Qbar += Qi;
		cout<<"i "<<i<<"  Qi= "<<Qi<<endl;
	}
	Qbar/=(double)No;
	printf("Qbar = %f\n", Qbar);
	for(int i=0; i<No; i++)
	{
		rr[i] = Ntwk0.junction1s[whichnodes_out[i]]->getFlowThrough(dt)-Qbar;
		printf("r[%d] = %f\n", i, rr[i]);
	}
}

void opt_eq_outflow::compute_J()
{
	printf("Coming soon!\n");

		J.reset_values(0.0);

//	compute_r();
	//Ntwk.junction1s[0]->setbVal(x);
	//printf("bvals are [%f, %f]\n",Ntwk.junction1s[0]->bval[0],Ntwk.junction1s[1]->bval[0]);
	#pragma omp parallel for
	for(int j = 0; j<n; j++)
	{
		vector<Real> rp(M+1);
		vector<Real> rm(M+1);
		vector<Real> x2(x);
		x2[j]+=mydelta;
		compute_r(rp,x2);	
		x2[j] -=2*mydelta;
		compute_r(rm, x2);
		for(int i=0; i<No; i++)
		{
			J(i,j) = (rp[i]-rm[i])/(2.*mydelta);
		}
//		x[j] += mydelta;	
	}	
//	compute_r();
}


#endif

