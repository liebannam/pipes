
#ifndef OPTDETIRTUS_H
#define OPTDETRITUS_H

#include "levmar.h"
#include "setupandrun.h"


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
		T = T_;
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

void bc_opt_dh::compute_r()
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

void bc_opt_dh::compute_J()
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
	


#endif
