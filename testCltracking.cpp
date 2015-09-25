#include "setupandrun.h"


int main(int argc, char* argv[])
{

    clock_t start_t, end_t;
    double A = 0.5;
    double B = 0.5;
    double KCl = 0.2;
    char *finp = argv[1];
	char *fconfig = argv[2];
	int M, Mi;
	double T;
	int channeltype = 1;
	
    Network * Ntwk = setupNetwork(finp, fconfig, M,Mi,T, channeltype);
	int N = Ntwk->channels[0]->N;
    if (N>1) printf("whoops, check your .inp file. this network has %d pipes\n", N);
    double L =Ntwk->channels[0]->L;
	double dx = L/(double)N;
    double dt = T/(double)M;
    double Clvals0[M+1];
    double Clvals1[M+1];
    double pi = 4.*atan(1.);
    double u = Ntwk->channels[0]->q[N]/Ntwk->channels[0]->q[0];
    double t;
    for (int i = 0; i<M+1; i++){
        t = dt*(float)i;
        Clvals0[i] = exp(-KCl*t)*(A+B*cos(-2*pi*u/L*t));
        Clvals1[i] = exp(-KCl*t)*(A+B*cos(2*pi/L*(L-u*t)));
    }
    double c0[N];
    for (int i = 0; i<N; i++){
        c0[i] = A+B*cos(2*pi*dx/L*(double)i);
        Ntwk->channels[0]->Cl[i] = c0[i];
        Ntwk->channels[0]->Cl0[i] = c0[i];
    }
    Ntwk->junction1s[0]->setClbval(Clvals1);
    Ntwk->junction1s[1]->setClbval(Clvals0);
    Ntwk->channels[0]->setKCl(KCl);
    int Nedges = Ntwk->Nedges;
	double V0=Ntwk->getTotalVolume();
	start_t = clock();
	for(int k=0; k<Nedges; k++)
	{
		printf("Channel %d info:\n",k);
		Ntwk->channels[k]->showGeom();
	}
    for(int k = 0; k<N;k++)
    printf("%f   %f   %f\n",(float)k*dx, Ntwk->channels[0]->Cl[k], exp(-KCl*0.)*(A+B*cos(2*pi/L*(dx*(float)k))));
   	Ntwk->runForwardProblem(dt); 
    printf("Elapsed simulation time is %f\n", dt*(double)(M));
	double f = 0;
	for (int i=0; i<M+1; i++)f+=pow(dt*Ntwk->getAveGradH(i),2.)/2.;
	double V = Ntwk->getTotalVolume();
	cout<<"initial volume "<<V0<< "    "<<"Final Volume " <<V<< endl;
	cout<<"dV = "<<V-V0<<endl;
	cout<<"f = "<<f<<endl;
	cout<<"maximum wave speed is "<<Ntwk->channels[0]->Cgrav(PI*.25/4.,false)<<endl;
	printf("Af is %f\n", Ntwk->channels[0]->At);
	printf("dt = %f , dx = %f, CFL = %f\n",dt, dx, dt/dx*Ntwk->channels[0]->a);
    double times[1] = {0};
	int which[1] = {0};
    double Cest, Ctrue, err, e2=0.;
   // Ntwk->channels[0]->quickWrite(times, which, 1,T,1);
    for(int k = 0; k<N;k++){
        Cest = Ntwk->channels[0]->Cl[k];
        Ctrue =exp(-KCl*T)*(A+B*cos(2*pi/L*(dx*(float)k - u*T)));
        err = Cest-Ctrue;
        e2+=dx*err*err;
        printf("%f  %f  %f  %f  %f\n",(float)k*dx, c0[k], Cest, Ctrue, err);
    }
    printf("Chlorine decay const = %f and u = %f\n",KCl,u);
    printf("||c-ctrue||_2 = %f\n",e2);
}
