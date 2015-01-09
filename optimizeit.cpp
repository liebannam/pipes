/////////////////////////
/*Read in input files and optimize things. THIS IS PRETTY ANNOYING...
*command line syntax is ../setupandrun example.inp example.config
*where:
*       example.inp is an EPANET-style file containing network layout information 
*       example.config is a custom configuration file with run-time parameters (e.g. number of cells, ICs, BCs)
*TO DO: make this program check everything out to ensure the .inp and .config files match!
*
*///////////////////
#include "optimizeit.h"
//


int main(int argc, char *argv[] )	
{
if(0)//optimization with bc_opt_dh
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
	srand (time(NULL));
	int Nrounds = 1;//how many times to adjust the thing
	int modetype = 1;
	int whichnode = 0;
	int Nn = 2;// number of nodes varied
	vector <int> whichnodes(Nn);
	whichnodes[0] = 1;
	whichnodes[1] = 2;
	vector<Real> x0(Nn*(ndof),0);
	for (int i = 0; i<2;i++){
		double b0 = Ntwk.junction1s[whichnodes[i]]->bval[0];
		double Dt = T/(ndof/2-1); //hermite interpolation spacing
		vector<double> h(M+1);
		if (modetype)x0[i*(ndof+1)] = 2*b0;
		else{
			for(int k = 0;k<(ndof)/2+1;k++)
			{
				x0[i*(ndof)+2*k] = b0;
				x0[i*(ndof)+2*k-1] = 0.;
			}
		}
	//	for(int k = 0; k<Nn*ndof; k++){x0[k] = ((double)rand()/RAND_MAX)-.5;}
	}
	for(int dd = 0; dd<Nrounds; dd++)
	{
//	bc1_opt_dh test1(ndof, M, x0, Ntwk, modetype, T, whichnode);
	bc_opt_dh test1(ndof*Nn, M, x0, Ntwk, modetype, T, whichnodes,10);
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
	double ompstart = omp_get_wtime();
	test1.solve();
	double ompend = omp_get_wtime();
	end_t = clock();
	char file0[15];
	sprintf(file0,"test%3ddump.txt",dd); 
	FILE *fp = fopen(file0,"w");
	test1.dump(fp);
	fclose(fp);

	test1.dump();
	test1.compute_f();
	double fnew = test1.f;
	double solvet = (end_t-start_t)/(double)CLOCKS_PER_SEC;
	vector <Real> bf(M+1,0.);
	vector <Real> xfake(ndof+1,0.);
	char file1[19];
       	sprintf(file1, "boundaryvals%3d.txt",dd);
	FILE *fb = fopen(file1, "w");
	fprintf(fb,"f is %f\n",fnew);
	for (int i = 0; i<Nn; i++)
	{	
		for (int k = 0; k<ndof; k++)
		{
			xfake[k] = test1.x[i*(ndof)+k];
		}
		getTimeSeries(bf, xfake, ndof,M,T,modetype);
		fprintf (fb, "node is %d\n", whichnodes[i]);
		//fprintf(fb, "t        Q(Junction %d)  Q(Junction %d)\n", 1,2);
		fprintf(fb, "t        Q(Junction %d)  \n", whichnodes[i]);
		for (int j = 0; j<M+1; j++){
			fprintf(fb, "%f   %f\n", dt*(double)j, bf[j]);
		}
	}
	fclose(fb);
	printf("Round %d finished, f = %f\n", dd, fnew);	
	cout<<"x     x0"<<endl;
	for(int k = 0;k<x0.size();k++)cout<<test1.x[k]<<"  "<<x0[k]<<endl;
	cout<<"\n\n";
//	for (int k =0;k<M+1;k++)cout<<bf[k]<<endl;
	printf("chkder time = %f s; CPU solve time = %f s, real solve time = %f\n",chkt, solvet, ompend-ompstart);
	printf("fold = %f\n fnew = %f\n",f0, fnew);
	printf("a0 = %f, q0 = %f\n", test1.a0[0][0], test1.q0[0][0]);
	printf("T = %f, dt = %f , dx = %f, CFL = %f\n",T, dt, dx, dt/dx*test1.Ntwk.channels[0]->a);
	printf("number of nodes is %d\n", Nn);	
	for(int k = 0; k<x0.size(); k++)
	{
		x0[k] = test1.x[k]+((double)rand()/(double)RAND_MAX)*.01-.005;
	}
	//writeOutputText(test1.Ntwk, M, Mi);
	}
}

else// optimization with opt_eq_outflow (hahahahahahaha this isn't gonna work...)
{
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

	printf("T = %f, dt = %f , dx = %f, CFL = %f\n",T, dt, dx, dt/dx*Ntwk.channels[0]->a);
	Ntwk.runForwardProblem(dt);
	double times[1] = {T};
	int which[1] = {0};

	for(int i = 0; i<Ntwk.channels.size(); i++)
	{
		Ntwk.channels[i]->quickWrite(times, which, 1,T,1);
	}

	end_t = clock();
	printf("Elapsed Time is %f\n",(end_t-start_t)/(double)CLOCKS_PER_SEC);
	printf("Elapsed simulation time is %f\n", dt*(double)(M));
	double V = Ntwk.getTotalVolume();
	cout<<"initial volume "<<V0<< "    "<<"Final Volume " <<V<< endl;
	cout<<"dV = "<<V-V0<<endl;
	printf("T = %f, dt = %f , dx = %f, CFL = %f\n",T, dt, dx, dt/dx*Ntwk.channels[0]->a);
//optimization time! at last!
	cout<<"M = "<<M<<endl;	
	int Ndof = 2;   // degrees of freedom (in Fourier or Hermite modes)
	int Nn = 3;// number of nodes varied
	int No = 7;//number of outflow nodes where we're measuring
	vector <int> whichnodes(Nn);
	vector <int> whichnodes_out(No);
	whichnodes[0] = 0;
	whichnodes[1] = 1;
	whichnodes[2] = 2;
	whichnodes_out[0] = 1;
	whichnodes_out[1] = 2;
	whichnodes_out[2] = 3;
	whichnodes_out[3] = 4;
	whichnodes_out[4] = 5;
	whichnodes_out[5] = 6;
	whichnodes_out[6] = 7;
	vector<Real> x0(Nn*(Ndof),0);
	for(int i= 0; i<Nn; i++)
	{
		x0[i*Ndof] = PI;
		x0[i*Ndof+1] = 0.;//(float)i/2; 
	}
	printf("No is %d and Nn is %d and Ndof is %d and I am confused\n",No, Nn, Ndof);
	opt_eq_outflow test2(Nn*Ndof,No,x0, Ntwk, T, whichnodes, whichnodes_out);
	test2.compute_f();

	cout<<"Made it?\n";
	double fold = test2.f;
	cout<<"f = "<<test2.f<<endl;
//	test2.solve();
	cout<<"fold = "<<fold<<endl;
	cout<<"fnew ="<<test2.f<<endl;
	test2.dump();	
	char file0[15];
	sprintf(file0,"test%3ddump.txt",0); 
	FILE *fp = fopen(file0,"w");
	test2.dump(fp);
	fclose(fp);
	vector <Real> vt;
	vector <Real> xfake(Ndof);
	
	char file1[19];
       	sprintf(file1, "boundaryvals%3d.txt",0);
	FILE *fb = fopen(file1, "w");
	for(int i = 0; i<Nn; i++)
	{
	for (int k = 0; k<Ndof; k++)
		{
			xfake[k] = test2.x[i*(Ndof)+k];
		}
		test2.evaluateit1(vt, xfake, M,T);
		fprintf (fb, "node is %d\n", whichnodes[i]);
		fprintf(fb, "t        Q(Junction %d)  \n", whichnodes[i]);
		for (int j = 0; j<M+1; j++){
			fprintf(fb, "%f   %f\n", dt*(double)j, vt[j]);
		}
	}
	
	fclose(fb);


	
//

}
}

