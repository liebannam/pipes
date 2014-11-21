/////////////////////////
/*Read in input files and optimize things. THIS IS PRETTY ANNOYING...
*command line syntax is ../setupandrun example.inp example.config
*where:
*       example.inp is an EPANET-style file containing network layout information 
*       example.config is a custom configuration file with run-time parameters (e.g. number of cells, ICs, BCs)
*TO DO: make this program check everything out to ensure the .inp and .config files match!
*
*///////////////////
//#include "network.h"
#include "levmar.h"
#include "optdetritus.h"
#include "setupandrun.h"
//




int main(int argc, char *argv[] )	
{
	int writelogs = 0; //change to 1 if the height plots look dull on account of large pressure variations

	//first open .inp file and process information about network layout and components
	char * finp= argv[1], *fconfig = argv[2];
	int M, Mi;
	double T;
	int channeltype = 1;
	Network Ntwk = setupNetwork(finp, fconfig, M,Mi,T, channeltype);
	double dt = T/(double)M;
	double dx = Ntwk.channels[0]->L/(double)Ntwk.channels[0]->N;
	int Nedges = Ntwk.Nedges;
	double V0=Ntwk.getTotalVolume();
	clock_t start_t, end_t;	
	start_t = clock();
//	for(int k=0; k<Nedges; k++){
	//	Ntwk.channels[k]->showGeom();
//		Ntwk.channels[k]->showp();
//	}

//	Ntwk.runForwardProblem(dt);
	
	end_t = clock();
	printf("Elapsed Time is %f\n",(end_t-start_t)/(double)CLOCKS_PER_SEC);
	printf("Elapsed simulation time is %f\n", dt*(double)(M));
	double V = Ntwk.getTotalVolume();
	cout<<"initial volume "<<V0<< "    "<<"Final Volume " <<V<< endl;
	cout<<"dV = "<<V-V0<<endl;
	
//optimization time! at last!
	cout<<"M = "<<M<<endl;

	int ndof = 16;   // degrees of freedom (in Fourier or Hermite modes)
	int modetype = 0;
	int whichnode = 1;
	double b0 = Ntwk.junction1s[whichnode]->bval[0];
	double Dt = T/(ndof/2-1); //hermite interpolation spacing
	vector<double> h(M+1);
	vector<Real> x0(ndof+1,0);
	if (modetype)x0[0] = 2*b0;
	else{
		for(int k = 0;k<(ndof)/2+1;k++)
		{
			x0[2*k] = b0;
			x0[2*k-1] = 0.;
		}
	}
	bc_opt_dh test1(ndof, M, x0, Ntwk, modetype, T, whichnode);
	cout<<"Made it?\n";
//	cout<<T<<endl;
	double places[] = {T};
	int which[] = {0};
	test1.compute_f();
//	for (int k = 0; k<3; k++)
//	test1.Ntwk.channels[k]->quickWrite(places, which, 1, T,1);

	cout<<"f0 is!"<<test1.f<<endl;
	test1.compute_f();
	double f0 = test1.f;
//	test1.Ntwk.channels[1]->showVals(1);
	cout<<"f0 is!"<<f0<<endl;
	start_t = clock();
//	test1.chkder();
	end_t = clock();
	double chkt = (end_t-start_t)/(double)CLOCKS_PER_SEC;
	start_t = clock();
	test1.solve();
	end_t = clock();
	test1.compute_f();
	double fnew = test1.f;
	test1.dump();
	FILE *fp = fopen("test1dump.txt","w");
	test1.dump(fp);
	fclose(fp);
	double solvet = (end_t-start_t)/(double)CLOCKS_PER_SEC;
	vector <Real> bf(M+1);
	getTimeSeries(bf, test1.x, ndof,M,T,modetype);
	FILE *fb = fopen("boundaryvals.txt", "w");
	//fprintf(fb, "t        Q(Junction %d)  Q(Junction %d)\n", 1,2);
	fprintf(fb, "t        Q(Junction %d)  \n", whichnode);
	for (int j = 0; j<M+1; j++){
		fprintf(fb, "%f   %f\n", dt*(double)j, bf[j]);
	}
	fclose(fb);
	

	for(int k = 0;k<ndof;k++)cout<<test1.x[k]<<"  "<<x0[k]<<endl;
	cout<<"\n\n";
//	for (int k =0;k<M+1;k++)cout<<bf[k]<<endl;
	printf("chkder time = %f s and solve time = %f s\n",chkt, solvet);
	printf("fold = %f\n fnew = %f\n",f0, fnew);
	printf("a0 = %f, q0 = %f\n", test1.a0[0][0], test1.q0[0][0]);
	printf("dt = %f , dx = %f, CFL = %f\n",dt, dx, dt/dx*test1.Ntwk.channels[0]->a);
	
	writeOutputText(test1.Ntwk, M, Mi);
//
//////////////
//
//	
//	//Ntwk.channels[0]->showVals(1);
//	for(int j=0;j<Nedges; j++)
//	{
//		for(int k =0; k<Ns[j]; k++){cout<<lengths[j]/Ns[j]*k<<"   "<<Ntwk.channels[j]->hofA(Ntwk.channels[j]->q[k])<<endl;}
//		for(int k =0; k<Ntwk.channels[j]->N; k++){cout<<k<<"   "<<Ntwk.channels[j]->q[k]<<endl;}
//	}
//	for(int j=0;j<Nedges; j++){
//		for(int k =0; k<N; k++){
//			cout<<dx*k<<"   "<<Ntwk.channels[j]->hofA(Ntwk.channels[j]->q[k])<<" "<<Ntwk.channels[j]->q[k+N]<<endl;
//			//cout<<dx*k<<"   "<<(Ntwk.channels[j]->q[k])<<endl;
//		}
//	}
//	cout<<"Number of 1 junctions is "<<Ntwk.junction1s.size()<<endl;
//	cout<<"Number of 2 junctions is "<<Ntwk.junction2s.size()<<endl;
//	cout<<"Number of 3 junctions is "<<Ntwk.junction3s.size()<<endl;
//	cout<<"Number of edges is "<<Nedges<<endl;
//
//
//
}

