#include "setupandrun.h"

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



