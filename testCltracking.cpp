#include "setupandrun.h"


int main(int argc, char* argv[])
{

    clock_t start_t, end_t;
    
    double A = 0.5;
    double B = 0.5;
    double KCl = 0.55;
    double kw = 0.0;
    char *finp = argv[1];
	char *fconfig = argv[2];
	int M, Mi;
	
    double T,L,u,e2=0;
	int channeltype = 1;
	
    Network * Ntwk = setupNetwork(finp, fconfig, M,Mi,T, channeltype);
    int Nedges = Ntwk->Nedges;
    int N = Ntwk->channels[0]->N;
    double dt = T/(double)M;
    double Clvals0[M+1];
    double Clvals1[M+1];
    double pi = 4.*atan(1.);
    u = Ntwk->channels[0]->q[Ntwk->channels[0]->N]/Ntwk->channels[0]->q[0];
    double Lloc = Ntwk->channels[0]->L;
    if (Nedges<2){ L =Lloc*Nedges;}
    else {L = Lloc*2;}
    double t,dx;
    for (int i = 0; i<M+1; i++)
    {
        t = dt*(float)i;
        Clvals0[i] = exp(-KCl*t)*(A+B*cos(-2*pi*u/L*t));
        Clvals1[i] = exp(-KCl*t)*(A+B*cos(2*pi/L*(L-u*t)));
    }


    if (Nedges==3)
    {
        int Nf = N*2;
        double c0[Nf];
        dx = Lloc/(double)N;
        for (int i = 0; i<Nf; i++)
        {
            c0[i] = A+B*cos(2*pi*dx/Lloc*(double)i);
        }
        for (int i = 0; i<M+1;i++)
        {
            Clvals0[i] = exp(-KCl*t)*(A+B*cos(-2*pi*u/Lloc*t));
            Clvals1[i] = 0.;
        }
        for (int i =0; i<N; i++)
        {
            Ntwk->channels[0]->Cl[i] =  c0[i];
            Ntwk->channels[0]->Cl0[i] = c0[i];
            Ntwk->channels[1]->Cl[i] =  c0[i+N];
            Ntwk->channels[1]->Cl0[i] = c0[i+N];
            Ntwk->channels[2]->Cl[i] =  c0[i+N];
            Ntwk->channels[2]->Cl0[i] = c0[i+N];
        }
        Ntwk->junction1s[0]->setClbVal(Clvals0);
        Ntwk->junction1s[1]->setClbVal(Clvals1);
        Ntwk->junction1s[2]->setClbVal(Clvals1);
        Ntwk->channels[0]->setClkw(kw);
        Ntwk->channels[1]->setClkw(kw);
        Ntwk->channels[2]->setClkw(kw);
        Ntwk->runForwardProblem(dt); 
        double times[1] = {0};
        int which[1] = {0};
        for (int j=0; j<3;j++)      
        {printf("j = %d\n\n",j);
            Ntwk->channels[j]->quickWrite(times, which, 1,T,1);
        }
        printf("Elapsed simulation time is %f\n", dt*(double)(M));
        double err, Cest, Ctrue,x;
        for(int j =0; j<3; j++)
        {
            int jfake = min(j,1);
            //printf("jfake=%d\n\n!\n",jfake);
            for(int k = 0; k<Ntwk->channels[0]->N;k++)
            {
                Cest = Ntwk->channels[j]->Cl[k];
                x = dx*(float)k+jfake*Lloc; 
                Ctrue =exp(-KCl*T)*(A+B*cos(2*pi/Lloc*(x- u*T)));
                err = Cest-Ctrue;
                e2+=dx*err*err;
                printf("%f  %f  %f  %f  %f\n",(float)k*dx, c0[k], Cest, Ctrue, err);
            }
        }

    }
    else
    {
        int Nf =N*Nedges;
        double c0[Nf];
        dx = Lloc/(double)N;
        if (Nedges==2) 
        {
            printf("N = %d Nf = %d \n",N,Nf);
    //      printf("whoops, check your .inp file. this network has %d pipes\n", N);
            for (int i = 0; i<Nf; i++)
            {
                c0[i] = A+B*cos(2*pi*dx/L*(double)i);
            }
            for (int i =0; i<N; i++)
            {
                Ntwk->channels[0]->Cl[i] =  c0[i];
                Ntwk->channels[0]->Cl0[i] = c0[i];
                Ntwk->channels[1]->Cl[i] =  c0[i+N];
                Ntwk->channels[1]->Cl0[i] = c0[i+N];
            }
            Ntwk->junction1s[0]->setClbVal(Clvals0);
            Ntwk->junction1s[1]->setClbVal(Clvals1);
            Ntwk->channels[0]->setClkw(kw);
            Ntwk->channels[1]->setClkw(kw);

        }
        else
        {
            for (int i = 0; i<N; i++)
            {
                c0[i] = A+B*cos(2*pi*dx/L*(double)i);
                Ntwk->channels[0]->Cl[i] = c0[i];
                Ntwk->channels[0]->Cl0[i] = c0[i];
            }
            Ntwk->junction1s[0]->setClbVal(Clvals0);
            Ntwk->junction1s[1]->setClbVal(Clvals1);
            Ntwk->channels[0]->setClkw(kw);
        }
            double V0=Ntwk->getTotalVolume();
            start_t = clock();
            for(int k=0; k<Nedges; k++)
            {
                printf("Channel %d info:\n",k);
                Ntwk->channels[k]->showGeom();
            }
        for(int k = 0; k<N;k++)printf("%f   %f   %f\n",(float)k*dx, Ntwk->channels[0]->Cl[k], exp(-KCl*0.)*(A+B*cos(2*pi/L*(dx*(float)k))));
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
        double Cest, Ctrue, err=0., x;
       // Ntwk->channels[0]->quickWrite(times, which, 1,T,1);
        for(int j =0; j<Nedges; j++)
        {
            for(int k = 0; k<Ntwk->channels[0]->N;k++)
            {
                Cest = Ntwk->channels[j]->Cl[k];
                x = dx*(float)k+j*Lloc; 
                Ctrue =exp(-KCl*T)*(A+B*cos(2*pi/L*(x- u*T)));
                err = Cest-Ctrue;
                e2+=dx*err*err;
                printf("%f  %f  %f  %f  %f\n",(float)k*dx+j*Lloc, c0[k+j*N], Cest, Ctrue, err);
            }
        }
    }
    printf("Npipes = %d\n",Nedges);
    printf("Chlorine decay const = %f and u = %f and dx = %f\n",KCl,u,dx);
    printf("||c-ctrue||_2 = %.16f\n",e2);
}

