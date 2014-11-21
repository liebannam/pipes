/////////////////////////
/*Driver for main simulation
*command line syntax is ../setupandrun example.inp example.config
*where:
*       example.inp is an EPANET-style file containing network layout information 
*       example.config is a custom configuration file with run-time parameters (e.g. number of cells, ICs, BCs)
*TO DO: make this program check everything out to ensure the .inp and .config files match!
*
*///////////////////
#include "levmar.h"
#include "setupandrun.h"
#include "basic_time_series.h"

double getTheta2(double A, double D)
{
	int count;
	ftheta th(A,D);
	double theta =::ridders(th, -.1, 2*PI+1, &count,1e-15, 1e-15);
	return theta;
}

void testallthiscrap() //print out all the crap I'm computing to see if it makes any bloody sense
{
	//also wtf is going on with cgrav???!?!?
	double D = 1;
	int M = 100;
	int N = 100;
	double L = 1000;
	Cpreiss ch0(N,D,L,M);
	double Af = PI*D*D/4.;
	double At = ch0.At;
       	double Ts = ch0.Ts;	
	int K = 500;
	double dk =1./(double)K;
	bool P = true;
	double x,hp, h,c,eta,phi, etap, phip, cp, A2;

	ch0.setGeom(1200);
	printf("#Ts = %f D = %f  a = %f \n#A               h                  eta(free)        eta(pressurized)      h(pressurized)    cgrav phi(free) phi(pressurized)\n", ch0.Ts, ch0.D, ch0.a);
	for(int j = 1; j<3; j++)
	{
	for(int k= 10; k<81; k++)
	{	
	       double p = -k/10.;
	       x = pow(-1.,j)*(pow(10.,p))+Af;
	       h = ch0.fakehofA(x,false);
	       hp = ch0.fakehofA(x,true);
	       etap = Eta(x,D,At,Ts,true);
	       eta = Eta(x,D,At,Ts,false);
	       c = Cgrav(x,D,At,Ts,false);
	       cp = Cgrav(x,D,At,Ts,true);
	       phi = PhiofA(x,D,At,Ts,false);
	       phip = PhiofA(x,D,At,Ts,true);
	       printf("%.10e    %.10f      %.10f     %.10f    %.10f    %.10f    %.10f    %.10f\n", x-Af,h,eta, etap, hp, c, phi, phip);	

	}
	}
	 bool pp = false;
	 x = Af-(1e-7)/2.;
	 etap = Eta(x,D,At,Ts,true);
	 eta = Eta(x,D,At,Ts,false);
	 c = Cgrav(x,D,At,Ts,false);
	 cp = Cgrav(x,D,At,Ts,true);
	 phi = PhiofA(x,D,At,Ts,true);
	 printf("Evaluated at Af = %.10f\n %.10f  (eta free)\n %.10f (eta press)\n %.10f (c)\n %.10f  (cpress)\n%.10f phi", Af,eta, etap, c, cp, phi);
	 printf("At-Af = %e\n", fabs(ch0.At-x));
	 printf("A                   h(A)               |A-Af|           |hnew-hold|        |A(theta) - A|\n");
	 for(int k = 0; k<40; k++)
	 {
		 double aa = D*D*PI/4.*(1.-pow(2.,-k-1));
		 double h = ch0.HofA(aa,pp);
		 double hold =ch0.hofAold(aa);
		 double th = 2*acos(1-2.*h/D);
		 printf("%.16f   %.16f   %e   %e   %e\n", aa, h,fabs(aa-D*D*PI/4.), fabs(h-hold), fabs(D*D/8.*(th-sin(th))-aa)); 
	 }
	 FILE *fg1 = fopen("geomconfirm1.txt","w");
	 int Mp= 500;
	 fprintf(fg1, "#A                  h(A)                   I(A)                  c(A)                  phi(A)                  hA(phi(A))                  A(h(A))   htrue(A)\n");
	 for(int k = 0; k<Mp; k++)
	 {
		 //double aa = D*D*PI/(4*Mp)*(double)k;
		 double tt = PI*2/Mp*(double)k;
		 double aa = D*D/8.*(tt-sin(tt));
		 double ht = 0.5*D*(1-cos(tt*.5));
		 double h = ch0.HofA(aa,pp);
		 double I = ch0.pbar(aa,pp);
		 double ah = ch0.AofH(ht,pp);
		 double c = ch0.Cgrav(aa,pp);
		 double phi = ch0.PhiofA(aa,pp);
		 double ae = ch0.AofPhi(ch0.PhiofA(aa,pp),pp);
		 fprintf(fg1,"%.16f   %.16f    %.16f   %.16f   %.16f   %.16f   %.16f    %.16f\n", aa, h, I, c, phi, ae, ah, ht); 
	 }
	 fclose(fg1);


	 //compare timings
	 clock_t t0,t1,t2,t3,t4,t5,t6;
	 t0 = clock();
	 double yy;
	 int MM = 10000;
	 for (int i = 1; i<MM; i++){yy =ch0.HofA((float)i/MM*PI*D*D/4.,false);}
	 t1 = clock();
	 for (int i = 1; i<MM; i++) {yy =ch0.hofAold((float)i/MM*PI*D*D/4.); } 
	 t2  = clock();
	 for (int i = 1; i<MM; i++){yy =ch0.pbar_old((float)i/MM*PI*D*D/4., false);}
	t3 = clock();
	for (int i = 1; i<MM; i++){yy =ch0.AofH((float)i/MM*PI*D*D/4.,false); }
	t4 = clock();
	for (int i = 1; i<MM; i++){ yy =ch0.PhiofA((float)i/MM*PI*D*D/4.,false); }
	t5 = clock();
	for (int i = 1; i<MM; i++){yy =ch0.pbar((float)i/MM*PI*D*D/4., false); }
	t6 = clock();
	printf("%d evaluations\n", MM);
	cout<<"h(A) chebyshev eval time =          "<<(t1-t0)/(double)CLOCKS_PER_SEC<<endl;
	cout<<"h(A) rootfind for theta eval time = "<<(t2-t1)/(double)CLOCKS_PER_SEC<<endl;
	cout<<"I(A) rootfind eval time =           "<<(t3-t2)/(double)CLOCKS_PER_SEC<<endl;
	cout<<"I(A) new eval time =                "<<(t6-t5)/(double)CLOCKS_PER_SEC<<endl;
	cout<<"A(h) eval time =                    "<<(t4-t3)/(double)CLOCKS_PER_SEC<<endl;
	cout<<"phi(A) eval time =                  "<<(t5-t4)/(double)CLOCKS_PER_SEC<<endl;
	
}

//namespace setup{

int main(int argc, char *argv[] )	
{
	
	clock_t start_t, end_t;
	int writelogs = 0;
	int Ntwktype = 1;
//	int stylized =1;//=1 writes out all information to text files for a funky fun plotting experience with presssure_putittogether.py
//	int output_psi =1;// write output as equivalent pressure (in psi) from height fields

	
	char *finp = argv[1];
	char *fconfig = argv[2];
	int M, Mi;
	double T;
	int channeltype = 1;
	Network Ntwk = setupNetwork(finp, fconfig, M,Mi,T, channeltype);
	double dt = T/(double)M;
	double dx = Ntwk.channels[0]->L/(double)Ntwk.channels[0]->N;
	int Nedges = Ntwk.Nedges;
	double V0=Ntwk.getTotalVolume();
	start_t = clock();
//	for(int k=0; k<Nedges; k++){
	//	Ntwk.channels[k]->showGeom();
//		Ntwk.channels[k]->showp();
//	}
	printf("h = .014, A = %.10f\n",Ntwk.channels[1]->HofA(0.4,false));
	printf("h = .008, A = %.10f\n",Ntwk.channels[0]->HofA(0.45,false));
	Ntwk.runForwardProblem(dt);
	end_t = clock();
//	for(int k=0; k<Nedges; k++){
	//	Ntwk.channels[k]->showGeom();
//		Ntwk.channels[k]->showp();
//	}
//	for(int k=0; k<Ntwk.channels[0]->N; k++)
//	{ 
//		printf("%d   %f\n", k, Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q[(Ntwk.channels[0]->idx(0,k))], true));
//	}
	printf("Elapsed real time is %f\n", (end_t-start_t)/(double)CLOCKS_PER_SEC);	
	printf("Elapsed simulation time is %f\n", dt*(double)(M));
	double f = 0;
	for (int i=0; i<M+1; i++)f+=pow(dt*Ntwk.getAveGradH(i),2)/2.;
	double V = Ntwk.getTotalVolume();
	cout<<"initial volume "<<V0<< "    "<<"Final Volume " <<V<< endl;
	cout<<"dV = "<<V-V0<<endl;
	cout<<"f = "<<f<<endl;
	cout<<"maximum wave speed is "<<Ntwk.channels[0]->Cgrav(PI*.25/4.,false)<<endl;
	//
//printf("t    H(valve)  h(valve-phys)  Q(reservoir)\n");   
//for (int k = 0; k<M/Mi; k++)
//{
//	printf("%f     %f   \n", dt*float(k*Mi), Ntwk.channels[1]->hofA(Ntwk.channels[1]->q_hist[Ntwk.channels[1]->idx_t(0,N/2,k*Mi)]));}
//negprestest3.config
//	printf("%f     %.10f   %10f   %f\n", dt*float(k*Mi), Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(0,199,k*Mi)],true), Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(0,199,k*Mi)],false), Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(1,1,k*Mi)]);
//trajkovic.config
//	printf("%f     %.10f   %10f   %f\n", dt*float(k*Mi), Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(0,140,k*Mi)],true), Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(0,140,k*Mi)],false), Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(1,140,k*Mi)]);
//}


printf("h = 150, A = %.10f\n",Ntwk.channels[0]->AofH(150,false));
printf("h = 1, A = %.10f\n",Ntwk.channels[0]->AofH(1,false));
printf("Af is %f\n", Ntwk.channels[0]->At);
printf("dt = %f , dx = %f, CFL = %f\n",dt, dx, dt/dx*Ntwk.channels[0]->a);
/*for (int i=0;i<Ntwk.channels[0]->N; i++){
	//printf("%f   %f   \n",i*Ntwk.channels[0]->dx, Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q[Ntwk.channels[0]->idx(0,i)],Ntwk.channels[0]->P[i+1]));
	printf("%f   %f   \n",i*Ntwk.channels[0]->dx, Ntwk.channels[0]->fakehofA(Ntwk.channels[0]->q[Ntwk.channels[0]->idx(0,i)],Ntwk.channels[0]->P[i+1]));
}*/


double places0[1] = {0};
double places1[1] = {2.75};
double times[1] = {T};
int which[1] = {1};

//}
writeOutputTarga(Ntwk, M, Mi,T, 0);

Ntwk.channels[0]->quickWrite(places1, which, 1,T,100); 
Ntwk.channels[2]->quickWrite(places1, which, 1,T,100); 
Ntwk.channels[2]->quickWrite(times, which, 1,T,100); 
writeOutputText(Ntwk, M, Mi);
//testallthiscrap();
//printf("Coefficients!!\n");

//for(int i = 0; i<Ntwk.channels[0]->Ncheb+1; i++)
//printf("%d   %.15f    %.15f    %.15f    %.15f    %.15f   \n", i, coeffs_h[i],coeffs_p1[i],coeffs_p2[i], Ntwk.channels[0]->coeffs_a1[i],coeffs_a2[i]);
}//}


//optimization crap	
//	int ndof = 16;   // degrees of freedom (in Fourier or Hermite modes)
////	int modetype;
////	double Dt = T/(ndof/2-1); //hermite interpolation spacing
////	vector<double> h(M+1);
////	vector<Real> x0(ndof, 0.);
////	vector<Real> bvals(M+1,0);
//////////////
//
//	
//	//Ntwk.channels[0]->showVals(1);
//	for(int j=0;j<Nedges; j++)
//	{
//		for(int k =0; k<Ns[j]; k++){cout<<lengths[j]/Ns[j]*k<<"   "<<Ntwk.channels[j]->hofA(Ntwk.channels[j]->q[k])<<endl;}
//		for(int k =0; k<N; k++){cout<<dx*k<<"   "<<Ntwk.channels[j]->q[k]<<endl;}
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
////	cout<<"Final volume distribution ="<<Ntwk.channels[0]->getTheGoddamnVolume()<<" "<<Ntwk.channels[1]->getTheGoddamnVolume()<<" "<<Ntwk.channels[2]->getTheGoddamnVolume()<<endl; 
//	printf("triple juncton values are A0 = %f, A1 = %f, A2 = %f\n", Ntwk.channels[0]->q[Ns[0]-1], Ntwk.channels[1]->q[0], Ntwk.channels[2]->q[0]);
//
//



Network setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_)
{
	//first open .inp file and process information about network layout and components
	ifstream file1(finp);
	char BC_filename[50];
	string stuff;
	vector<int> jIDs, pIDs, conns;
	vector<double> lengths, diams, Mrs, S0s, xcoords, ycoords, elevs;
	int jflag =0, pflag = 0, cflag =0;
	int first;
	while (getline(file1, stuff, '\n')) 
	 {
		 if (stuff[0]==';')
			 continue;
		 else{ 
			if(strncmp(stuff.c_str(), "[JUNCTIONS]", 11 )==0)
			{	
				jflag =1;
				first = 1;
				pflag = 0;
				cflag = 0;
			}
			if(strncmp(stuff.c_str(), "[PIPES]", 7 )==0)
			{
				pflag =1;
				first =1;
				jflag =0;
				cflag =0;
			}
			if(strncmp(stuff.c_str(), "[COORDINATES]", 13 )==0)
			{
				cflag =1;
				first =1;
				jflag =0;
				pflag =0;
			}
			if(jflag)
			{      	
				stringstream ss(stuff); // Insert the string into a stream
				ss>>appendTo(jIDs)>>appendTo(elevs);
				if(stuff[0] =='[' && first==0){jflag =0;}
				first = 0;
			}
			if(pflag)
			{      	
				stringstream ss(stuff); 
				ss>>appendTo(pIDs)>>appendTo(conns)>>appendTo(conns)>>appendTo(lengths)>>appendTo(diams)>>appendTo(Mrs);
				if(stuff[0] =='[' && first==0){pflag =0;}
				first = 0;
			}
			if (cflag)
			{
				int tmp;
				stringstream ss(stuff);
				ss>>tmp>>appendTo(xcoords)>>appendTo(ycoords);
				if(stuff[0] =='[' && first==0){cflag =0;}
				first = 0;
			}
		 }
	 }
	file1.close();
	////
	//figure out slopes from elevation info
	for(int k = 0; k<pIDs.size(); k++)
	{	
		int lk = find(jIDs.begin(), jIDs.end(), conns[2*k]   ) - jIDs.begin();
		int rk = find(jIDs.begin(), jIDs.end(), conns[2*k+1] ) - jIDs.begin();
		S0s.push_back((elevs[lk] - elevs[rk])/lengths[k]);
	}
	//convert diamters to m
//	for(int k = 0; k<diams.size(); k++){
//		diams[k];
//		cout<<diams[k]<<endl;
//	}
	//convert labels to 0 based indexing (if not already)
	if(jIDs[0]>0)
	{ 
		int shift = jIDs[0];
		cout<<"shifting indexing by "<<shift<<endl;
		for(int k=0; k<jIDs.size(); k++) jIDs[k]-=shift;
		for(int k = 0; k<pIDs.size(); k++){
			pIDs[k]-=shift;
			conns[2*k]-=shift;
			conns[2*k+1]-=shift;
		}
	}

	//now read  .config file to set run parameters
	ifstream file2(fconfig);
	string morestuff;

	jflag =0;
	pflag = 0;
	first = 0;
	int tflag = 0;
	vector<int> bvaltypes;
	vector<double>bvals;
	vector<double>reflects;
	vector<double> offsets;
	vector<double> valveopens;
	vector<double> offset01s;
	vector<double> offset02s;
	vector<double> offset12s;
	vector<int> Ns;
	vector<double> h0s;
	vector<double> q0s;
	vector <Real > xbval;
	int whichnode;
	int modetype;
	while (getline(file2, morestuff, '\n')) 
	{
		 if (morestuff[0]==';')
			 continue;
		 else{ 
			if(strncmp(morestuff.c_str(), "[JUNCTION_INFO]", 15 )==0)
			{	
				jflag =1;
				first = 1;
				pflag = 0;
				tflag = 0;
			}
			if(strncmp(morestuff.c_str(), "[PIPE_INFO]", 11 )==0)
			{
				pflag =1;
				first =1;
				jflag =0;
				tflag = 0;
			}
			if(strncmp(morestuff.c_str(), "[TIME_INFO]", 11 )==0)
			{
				tflag =1;
				first =1;
				jflag =0;
				pflag = 0;
			}

			if(jflag)
			{      	
				int tmp,type;
				stringstream ss(morestuff); // Insert the string into a stream
				ss>>tmp>>type;
				if(type==1)//next three columns matter
				{
					ss>>appendTo(bvaltypes)>>appendTo(bvals)>>appendTo(reflects);
				}
				else if(type ==2)//skip next three columns; next 2 matter
				{	
					ss>>tmp>>tmp>>tmp>>appendTo(offsets)>>appendTo(valveopens);
				}
				else if(type ==3)//skip next 5, next 2 matter
				{
					ss>>tmp>>tmp>>tmp>>tmp>>tmp>>appendTo(offset01s)>>appendTo(offset02s)>>appendTo(offset12s);
				}

				if(morestuff[0] =='[' && first==0){jflag =0;}
				first = 0;
			}
			if(pflag)
			{      	
				int tmp;
				stringstream ss(morestuff); 
				ss>>tmp>>appendTo(Ns)>>appendTo(h0s)>>appendTo(q0s);
				if(morestuff[0] =='[' && first==0){pflag =0;}
				first = 0;
			}
			if(tflag)
			{
				stringstream ss(morestuff);	
				ss>>T>>M>>Mi;
				cout<<"T = "<<T<<endl;
				if(morestuff[0] =='[' && first==0){tflag =0;}
				first = 0;

			}
		 }
		 if(strncmp(morestuff.c_str(),"[BC_FILENAME]",13)==0)
		 {
			string tmp(morestuff.c_str());
			int L = tmp.length()/sizeof(tmp[0]);
			cout<<"L="<<L<<endl;
			int go = 0;
			int count  =0;
			for (int i = 0; i<L; i++){
				if (go>0){
					if(tmp[i]==' ')
					{
						count ++;
					}
					else{
						BC_filename[i-count] = tmp[i];
					}
				}
				else{
					count++;
//					char tmp2 = tmp[i];
					if (tmp[i]==']'){
						cout<<"yay\n";
						go=1;
					}
				}
			}
			cout<<"BC Filename is:"<<BC_filename<<endl;
			string evenmorestuff;
			ifstream fbc(BC_filename);
			string trash;
			//ifstream fbc("../indata/bcs2.txt");
			while (getline(fbc, evenmorestuff, '\n'))
				{
					stringstream ss(evenmorestuff);
					if(strncmp(evenmorestuff.c_str(),"modetype",8)==0)
					{
						ss>>trash>>modetype;
						cout<<"modetype = "<<modetype<<endl;
					}
					if(strncmp(evenmorestuff.c_str(),"whichnode",9)==0)
					{
						ss>>trash>>whichnode;
						cout<<"whichnode = "<<whichnode<<endl;
					}
					else{
						ss>>appendTo(xbval);
					}
					//cout<<evenmorestuff;
				}	
			fbc.close();
		}

	 }
	cout<<"Junction Info:\njunction ID     elevation\n";
	for(int k= 0; k<jIDs.size(); k++)
	{    cout<< jIDs[k]<<"                "<<elevs[k]<<"   "<<xcoords[k]<<" "<<ycoords[k]<<endl;}
	cout<<"Pipe Info:\npipe ID   left_end   right_end    length(m)   diam(mm)    manning coeff  slope\n";
	for(int k= 0; k<pIDs.size(); k++)
	{ 
		printf("%d          %d          %d            %.1f      %.3f        %.4f          %.4f     %.2f    %.2f\n", pIDs[k], conns[2*k], conns[2*k+1],lengths[k], diams[k], Mrs[k], S0s[k], h0s[k], q0s[k]); }


	int channeltype =channeltype_;  //0 for uniform cross section, 1 for Preissman slot

	int Nnodes = jIDs.size();
	int Nedges = pIDs.size();

	double dt = T/(double)M;

	//make the damn return network
	Network Ntwk(Nnodes, conns, Nedges, Ns, diams, lengths, S0s, Mrs, h0s, q0s, M, channeltype);
	
	
	cout<<"Number of 1 junctions is "<<Ntwk.junction1s.size()<<endl;
	cout<<"Number of 2 junctions is "<<Ntwk.junction2s.size()<<endl;
	cout<<"Number of 3 junctions is "<<Ntwk.junction3s.size()<<endl;
	cout<<"Number of edges is "<<Nedges<<endl;
	for(int k = 0; k<Ntwk.channels.size(); k++)
	{
		double a0 = Ntwk.channels[k]->AofH(h0s[k],false);
		printf("h0 = %f, a0 = %.10f, h(a0) = %f\n", h0s[k],a0, Ntwk.channels[k]->HofA(a0,false) );
		Ntwk.channels[k]->setq(a0, q0s[k]);
		Ntwk.channels[k]->setq0(a0, q0s[k]);
	}	
	for(int k = 0; k<Ntwk.junction1s.size(); k++)
	{
		Ntwk.junction1s[k]->bvaltype = bvaltypes[k];
		Ntwk.junction1s[k]->setbVal(bvals[k]);
		Ntwk.junction1s[k]->reflect = reflects[k];
	}
	for(int k = 0; k<Ntwk.junction2s.size(); k++)
	{
		Ntwk.junction2s[k]->offset = offsets[k];
		Ntwk.junction2s[k]->valveopen = valveopens[k];
	}
	for(int k = 0; k<Ntwk.junction3s.size(); k++)
	{
		Ntwk.junction3s[k]->j2_01.offset = offset01s[k];
		Ntwk.junction3s[k]->j2_20.offset = offset02s[k];
		Ntwk.junction3s[k]->j2_12.offset = offset12s[k];
		Ntwk.junction3s[k]->j2_21.offset = offset12s[k];
	}
	if (xbval.size()>0){
		printf("Setting boundary values for node %d using %s modes:\n",whichnode, (modetype?"Fourier":"Hermite"));
		for (int ii = 0; ii<xbval.size();ii++)cout<<"i = "<<ii<<" m[i]= "<<xbval[ii]<<endl;
		vector<Real> bvals(M+1);
		getTimeSeries(bvals, xbval, xbval.size(), M, T, modetype);
	//	for(int ii= 0; ii<bvals.size(); ii++)
	//	{
	//		printf("%f  %f\n", ii*T/M,bvals[ii]);
	//	}
		Ntwk.junction1s[whichnode]->setbVal(bvals);	

	}

	
	char mdata[] = "../output_data/mapdata.txt";
	FILE *fm = fopen(mdata, "w");
	if(fm==NULL)
	{
		fputs("Error opening file!!", stderr);
		throw;
	}

	////warning mess with the following line at your peril because the plot script will epically *&%^ up...
	fprintf(fm, "%f\n", T);
	for(int k =0; k<Nedges; k++){fprintf(fm, "%d   %d   ", Ntwk.conns[2*k], Ntwk.conns[2*k+1]);}
	fprintf(fm, "\n");
	for(int k = 0; k<Nedges; k++){fprintf(fm, "%f  ", (Ntwk.channels[k]->w)/2.);}	
	fprintf(fm, "\n");
	for(int k =0; k<Nnodes; k++){
		fprintf(fm, "%d   %f    %f    %lf\n  ",jIDs[k], xcoords[k], ycoords[k],elevs[k]);
	}
	fclose(fm);
	printf("network data written to %s\n", mdata);
	return Ntwk;
}


////output heightfields and a textfile "runinfo.txt" of maxvalues to accompany them.
void writeOutputTarga(Network &Ntwk, int M, int Mi, double T, int writelogs)
{
	char sdata[] = "../output_data/scalings.txt";
	int Nedges = Ntwk.Nedges;
	int Nnodes = Ntwk.Nnodes;
	FILE *fp = fopen(sdata, "w");
		if(fp==NULL)
		{
			fputs("Error opening file!!", stderr);
			throw;
		}

	for(int ii=0;ii<M; ii+=Mi)
	{
	//	printf("i=%d\n", ii);
		fprintf(fp, "%d    ", ii/Mi);
		for(int kk = 0; kk<Nedges;kk++)
		{
			char filename[100];
			sprintf(filename,"../output_data/out%d_%03d.tga", kk,ii/Mi);
			int mm = Ntwk.channels[kk]->N;                        //x direction
			int nn = mm/2;	                        //y direction
			double myfld[mm*nn];
			double val;	
			double zmin =0, zmax = 0;
			for(int j = 0;j<nn; j++)
			{
				for (int i= 0; i<mm; i++)
				{
				//	if(output_psi==1)
				//	{
				//		val = m_to_psi*Ntwk.channels[kk]->fakehofA(Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)], Ntwk.channels[kk]->P[ii+1]);
				//	}
					if(writelogs)
					{
						bool p =  false;//Ntwk.channels[kk]->P_hist[Ntwk.channels[kk]->pidx_t(i+1, ii)]
						double a = Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)];
						val = log(Ntwk.channels[kk]->HofA(a,p)+1);
					}
					else{
						bool p = false;
						double a = Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)];
						val = Ntwk.channels[kk]->HofA(a,p);
					}	
					myfld[i+mm*j] = val;
				       //cout<<val<<"   ";	
					zmin = fmin(zmin, val);
					zmax = fmax(zmax, val);
				}
			}
			
			fprintf(fp, "%f   %f   ", zmax, zmin);
			zmin = 0;
		//	zmax = Ntwk.channels[kk]->w;
		//	printf("zmax =  %f, zmin =%f  \n", zmax, zmin);
			w3d_targa_output_surface(filename, myfld,mm,nn,zmin,zmax);
			}
	  // printf("\n");
	   fprintf(fp, "\n");
	}

	fclose(fp);
	printf("Writing %s to targa files for plotting\n", writelogs?"log(height fields)":"height fields");
	printf("scalings data writen to %s\n",sdata);
	cout<<"number of writes is  "<<M/Mi<<endl;

}

void writeOutputText(Network &Ntwk, int M, int Mi)	
{
	
	int Nedges= Ntwk.Nedges;
	cout<<"number of writes is  "<<M/Mi<<endl;
	for(int ii=0;ii<M; ii+=Mi)
	{
		double val;
		char fname[100];
		sprintf(fname,"../output_data/txt_out%03d.txt",ii/Mi);		
		FILE *fd = fopen(fname, "w");
		if(fd==NULL)
		{
			fputs("Error opening file!!", stderr);
			throw;	
		}
		for(int kk=0; kk<Nedges; kk++)
		{
			int NN = Ntwk.channels[kk]->N;
			for(int j = 0;j<NN; j++)
			{
				val = m_to_psi*Ntwk.channels[kk]->fakehofA(Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,j+1, ii)], false);
				fprintf(fd, "%.10f     ", val);

			}
		fprintf(fd, "\n");	
		}
		fclose(fd);
	}
	printf("Writing psi data to, e.g., ../output_data/txt_out000.txt\n number of writes is %d\n", M/Mi); 
/*	double Al = Ntwk.channels[0]->q[49];
	double Ar =  Ntwk.channels[0]->q[49];
	printf("Al = %f, Ar = %f\n", Al, Ar);
	testcrappy_f(Al, Ar, Ntwk.channels[0]->w, Ntwk.channels[0]->Ts,Ntwk.channels[0]->At);*/
}


