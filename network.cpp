/******
 * implementation of network.h
 * for details, see comment at top of network.h
 ****/

#include "network.h"

int find_nth(std::vector<int> values,  int desired, int n_ordinal, int N)
{
	int retval = -1;
	int count =0;
	for(int i = 0; i<N; i++)
	{
		if(values[i]== desired)
		{
			count+=1;
			retval = i;
		}
		if(count == n_ordinal){return retval;}
	}		
}
/*
Network_params::Network_params(std::vector<int> Ns_, std::vector<double> ws_, std::vector<double> Ls_, 
					std::vector<double> S0s_, std::vector<double> Mrs_, std::vector<double> 
					a0s_, std::vector<double> q0s_,double a_):a(a_)
{ 
	Ns = Ns_;
	ws = ws_;
	Ls = Ls_;
	S0s = S0s_;
	Mrs = Mrs_;
	a0s = a0s_;
	q0s = q0s_;
			
}*/
/*Network::Network(int Nnodes_, std::vector<int> conns_, int Nedges_, int M_,  int channeltype_, Network_params p):
	        Nnodes(Nnodes_), Nedges(Nedges_), M(M_), channeltype(channeltype_)*/

Network::Network(int Nnodes_, std::vector<int> conns_, int Nedges_, std::vector<int> Ns, std::vector<double> ws, std::vector<double> Ls, 
		std::vector<double> S0s, std::vector<double> Mrs, std::vector<double> a0s, std::vector<double> q0s, int M_,  int channeltype_, double a):
	        Nnodes(Nnodes_), Nedges(Nedges_), M(M_), channeltype(channeltype_)
{
	
	nn = 0; //start with time =0	
	//initialize internal array of connectivity data 
	//conns = new int[Nedges*2];
	for(int k = 0; k<Nedges*2; k++)
	{
		conns.push_back(conns_[k]);
	}
	//figure out the type of each node (is it connected to 1, 2 or 3 pipes)
//	nodeTypes = new int[2*Nedges];
	for (int k=0; k<2*Nedges; k++){nodeTypes.push_back(0.);}
	for(int k=0; k<2*Nedges; k++)
	{
		nodeTypes[conns[k]]+=1;
	}
	//fill channels with data from Ns, ws, Ls, etc
	for(int i = 0; i<Nedges; i++)
	{	
		if(channeltype==0){channels.push_back(new Cuniform(Ns[i], ws[i], Ls[i],M,a));}
		else{channels.push_back(new Cpreiss(Ns[i], ws[i], Ls[i],M,a));}
		channels[i]->setq(a0s[i],q0s[i]);
		channels[i]->setq0(a0s[i], q0s[i]);
		channels[i]->Mr = Mrs[i];
		channels[i]->S0 = S0s[i];
	//	channels[i]->showGeom();
	}

	//now can assign channels to the correct junctions
	for(int j=0; j<Nnodes; j++)
	{
	//	printf("node %d gets a junction%d\n", j,nodeTypes[j]);
		if(nodeTypes[j] ==1)
		{
			int idx =find_nth(conns, j, 1, 2*Nedges);
	//		printf(" index is %d, row is %d, column is %d\n", idx, idx/2, idx%2);
			printf("\njunction1!!!\n edge is %d and whichend is %d\n ", idx/2, idx%2); 
		        junction1s.push_back(new Junction1(*channels[idx/2], idx%2 ,0.,1)); //row in conns corresponds to edge; column corresponds to left or right end.
		}
		else if(nodeTypes[j] ==2)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			printf("\njunction2!!! index1 is %d, edge0 is %d, which0 is %d\n", idx1, idx1/2, idx1%2);
			printf(" index2 is %d, edge1 is %d, which1 is %d\n", idx2, idx2/2, idx2%2);
			junction2s.push_back(new Junction2(*channels[idx1/2], *channels[idx2/2], idx1%2, idx2%2, 1.0));
		}

		else if(nodeTypes[j] ==3)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			int idx3 =find_nth(conns, j, 3, 2*Nedges);
			printf("\njunction3!!!\n index1 is %d, row is %d, column is %d\n", idx1, idx1/2, idx1%2);
			printf(" index2 is %d, row is %d, column is %d\n", idx2, idx2/2, idx2%2);	
			printf(" index3 is %d, row is %d, column is %d\n", idx3, idx3/2, idx3%2);
			//set it up so that either you have [end0, end1, end2] = [1,0,0] or = [0,1,1];
			//[0,1,1] case
			if(idx1%2+idx2%2+idx3%2 ==2){
				if (idx1%2 == 0){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );}
				else if(idx2%2 == 0){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );}

			}
			//[1,0,0] case
			else if(idx1%2+idx2%2+idx3%2 ==1){
				if(idx1%2 ==1){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );}
				else if(idx2%2 ==1){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
				printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );
				}
			}
			else {
				cout<<"Ruh-roh, boundary conditions not implemented for this configuration yet. Your simulation is fairly certainily going to suck.\n";
			}
			
			
		}			
		else {printf("Well, then I have no idea. --Keanu Reeves  (nodetypes[j] = %d)\n", nodeTypes[j]);
			for (int k = 0; k<Nedges; k++)cout<<conns[2*k]<<" "<<conns[2*k+1]<<endl;
		}

	}
}

///WELCOME TO THE COPY CONSTRUCTOR! IT'S A...work in progress 
//
Network::Network(const Network &N_old):Nnodes(N_old.Nnodes), Nedges(N_old.Nedges), M(N_old.M)
	      {
	nn = 0;
	double a = N_old.channels[0]->a;
	channeltype = N_old.channeltype;
	//copy nodes and connectivity info
	for( int i =0; i<Nedges*2; i++){conns.push_back(N_old.conns[i]);}
	for(int i =0; i<Nedges*2; i++){nodeTypes.push_back(N_old.nodeTypes[i]);}
     	//copy channels, their values and their parameters N, w, L, M, a, M, S0
	for(int i = 0; i<Nedges; i++)
	{	
		int Ni = N_old.channels[i]->N;
		if(N_old.channeltype ==0){channels.push_back(new Cuniform(Ni, N_old.channels[i]->w, N_old.channels[i]->L,M,a));}
		else{channels.push_back(new Cpreiss(Ni, N_old.channels[i]->w, N_old.channels[i]->L,M,a));}
		for(int k=0; k<2*Ni; k++){
			channels[i]->q[k] = N_old.channels[i]->q[k];
			channels[i]->q0[k] = N_old.channels[i]->q0[k];
		}
		channels[i]->Mr = N_old.channels[i]->Mr;
		channels[i]->S0 = N_old.channels[i]->S0;
	}
	//now can assign channels to the correct junctions
	for(int j=0; j<Nnodes; j++)
	{
	//	printf("node %d gets a junction%d\n", j,nodeTypes[j]);
		if(nodeTypes[j] ==1)
		{
			int idx =find_nth(conns, j, 1, 2*Nedges);
	//		printf(" index is %d, row is %d, column is %d\n", idx, idx/2, idx%2);
	//		printf("\njunction1!!!\n edge is %d and whichend is %d\n ", idx/2, idx%2); 
		        junction1s.push_back(new Junction1(*channels[idx/2], idx%2 ,0.,1)); //row in conns corresponds to edge; column corresponds to left or right end.
		}
		else if(nodeTypes[j] ==2)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			//printf("\njunction2!!! index1 is %d, edge0 is %d, which0 is %d\n", idx1, idx1/2, idx1%2);
			//printf(" index2 is %d, edge1 is %d, which1 is %d\n", idx2, idx2/2, idx2%2);
			junction2s.push_back(new Junction2(*channels[idx1/2], *channels[idx2/2], idx1%2, idx2%2, 1.0));
		}

		else if(nodeTypes[j] ==3)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			int idx3 =find_nth(conns, j, 3, 2*Nedges);
	///		printf("\njunction3!!!\n index1 is %d, row is %d, column is %d\n", idx1, idx1/2, idx1%2);
	//		printf(" index2 is %d, row is %d, column is %d\n", idx2, idx2/2, idx2%2);	
	//		printf(" index3 is %d, row is %d, column is %d\n", idx3, idx3/2, idx3%2);
			//set it up so that either you have [end0, end1, end2] = [1,0,0] or = [0,1,1];
			//[0,1,1] case
			if(idx1%2+idx2%2+idx3%2 ==2){
				if (idx1%2 == 0){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );
				}
				else if(idx2%2 == 0){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );
				}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );
				}

			}
			//[1,0,0] case
			else if(idx1%2+idx2%2+idx3%2 ==1){
				if(idx1%2 ==1){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );
				}
				else if(idx2%2 ==1){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );
				}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );
				}
			}
			else {
				cout<<"Ruh-roh, boundary conditions not implemented for this configuration yet. Your simulation is fairly certainily going to suck.\n";
			}
			
			
		}			
		else {printf("Well, then I have no idea. --Keanu Reeves  (nodetypes[j] = %d)\n", nodeTypes[j]);
			for (int k = 0; k<Nedges; k++)cout<<conns[2*k]<<" "<<conns[2*k+1]<<endl;
		}

	}
	//copy other junction info
	for(int k = 0; k<junction1s.size(); k++)
	{
		junction1s[k]->bvaltype = N_old.junction1s[k]->bvaltype;
		junction1s[k]->setbVal(N_old.junction1s[k]->bval);
		//cout<<junction1s[k]->bval[0]<<"yupgoddamnitfuckyou\n";
		junction1s[k]->reflect = N_old.junction1s[k]->reflect;
	}
	for(int k = 0; k<junction2s.size(); k++)
	{
		junction2s[k]->offset = N_old.junction2s[k]->offset;
		junction2s[k]->valveopen = N_old.junction2s[k]->valveopen;
		junction2s[k]->setValveTimes(N_old.junction2s[k]->valvetimes);
	}
	
	for(int k = 0; k<junction3s.size(); k++)
	{
		junction3s[k]->j2_01.offset = N_old.junction3s[k]->j2_01.offset; 
		junction3s[k]->j2_20.offset = N_old.junction3s[k]->j2_20.offset;
		junction3s[k]->j2_12.offset = N_old.junction3s[k]->j2_12.offset;
		junction3s[k]->j2_21.offset = N_old.junction3s[k]->j2_21.offset;
	}


}

//copy constructor from pointer to an old network...
Network::Network(Network *N_old):Nnodes(N_old->Nnodes), Nedges(N_old->Nedges), M(N_old->M)
	      {
	nn = 0;
	double a = N_old->channels[0]->a;
	channeltype = N_old->channeltype;
	//copy nodes and connectivity info
	for( int i =0; i<Nedges*2; i++){conns.push_back(N_old->conns[i]);}
	for(int i =0; i<Nedges*2; i++){nodeTypes.push_back(N_old->nodeTypes[i]);}
     	//copy channels, their values and their parameters N, w, L, M, a, M, S0
	for(int i = 0; i<Nedges; i++)
	{	
		int Ni = N_old->channels[i]->N;
		if(N_old->channeltype ==0){channels.push_back(new Cuniform(Ni, N_old->channels[i]->w, N_old->channels[i]->L,M,a));}
		else{channels.push_back(new Cpreiss(Ni, N_old->channels[i]->w, N_old->channels[i]->L,M,a));}
		for(int k=0; k<2*Ni; k++){
			channels[i]->q[k] = N_old->channels[i]->q[k];
			channels[i]->q0[k] = N_old->channels[i]->q0[k];
		}
		channels[i]->Mr = N_old->channels[i]->Mr;
		channels[i]->S0 = N_old->channels[i]->S0;
	}
	//now can assign channels to the correct junctions
	for(int j=0; j<Nnodes; j++)
	{
	//	printf("node %d gets a junction%d\n", j,nodeTypes[j]);
		if(nodeTypes[j] ==1)
		{
			int idx =find_nth(conns, j, 1, 2*Nedges);
	//		printf(" index is %d, row is %d, column is %d\n", idx, idx/2, idx%2);
	//		printf("\njunction1!!!\n edge is %d and whichend is %d\n ", idx/2, idx%2); 
		        junction1s.push_back(new Junction1(*channels[idx/2], idx%2 ,0.,1)); //row in conns corresponds to edge; column corresponds to left or right end.
		}
		else if(nodeTypes[j] ==2)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			//printf("\njunction2!!! index1 is %d, edge0 is %d, which0 is %d\n", idx1, idx1/2, idx1%2);
			//printf(" index2 is %d, edge1 is %d, which1 is %d\n", idx2, idx2/2, idx2%2);
			junction2s.push_back(new Junction2(*channels[idx1/2], *channels[idx2/2], idx1%2, idx2%2, 1.0));
		}

		else if(nodeTypes[j] ==3)
		{
			int idx1 =find_nth(conns, j, 1, 2*Nedges);
			int idx2 =find_nth(conns, j, 2, 2*Nedges);
			int idx3 =find_nth(conns, j, 3, 2*Nedges);
	///		printf("\njunction3!!!\n index1 is %d, row is %d, column is %d\n", idx1, idx1/2, idx1%2);
	//		printf(" index2 is %d, row is %d, column is %d\n", idx2, idx2/2, idx2%2);	
	//		printf(" index3 is %d, row is %d, column is %d\n", idx3, idx3/2, idx3%2);
			//set it up so that either you have [end0, end1, end2] = [1,0,0] or = [0,1,1];
			//[0,1,1] case
			if(idx1%2+idx2%2+idx3%2 ==2){
				if (idx1%2 == 0){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );
				}
				else if(idx2%2 == 0){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );
				}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );
				}

			}
			//[1,0,0] case
			else if(idx1%2+idx2%2+idx3%2 ==1){
				if(idx1%2 ==1){junction3s.push_back(new Junction3(*channels[idx1/2], *channels[idx2/2], *channels[idx3/2], idx1%2, idx2%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx1/2, idx2/2, idx3/2 );
				}
				else if(idx2%2 ==1){junction3s.push_back(new Junction3(*channels[idx2/2], *channels[idx1/2], *channels[idx3/2], idx2%2, idx1%2, idx3%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx2/2, idx1/2, idx3/2 );
				}
				else {junction3s.push_back(new Junction3(*channels[idx3/2], *channels[idx2/2], *channels[idx1/2], idx3%2, idx2%2, idx1%2));
			//	printf("ch0 is %d ch1 is %d ch2 is %d", idx3/2, idx2/2, idx1/2 );
				}
			}
			else {
				cout<<"Ruh-roh, boundary conditions not implemented for this configuration yet. Your simulation is fairly certainily going to suck.\n";
			}
			
			
		}			
		else {printf("Well, then I have no idea. --Keanu Reeves  (nodetypes[j] = %d)\n", nodeTypes[j]);
			for (int k = 0; k<Nedges; k++)cout<<conns[2*k]<<" "<<conns[2*k+1]<<endl;
		}

	}
	//copy other junction info
	for(int k = 0; k<junction1s.size(); k++)
	{
		junction1s[k]->bvaltype = N_old->junction1s[k]->bvaltype;
		junction1s[k]->setbVal(N_old->junction1s[k]->bval);
		//cout<<junction1s[k]->bval[0]<<"yupgoddamnitfuckyou\n";
		junction1s[k]->reflect = N_old->junction1s[k]->reflect;
	}
	for(int k = 0; k<junction2s.size(); k++)
	{
		junction2s[k]->offset = N_old->junction2s[k]->offset;
		junction2s[k]->valveopen = N_old->junction2s[k]->valveopen;
		junction2s[k]->setValveTimes(N_old->junction2s[k]->valvetimes);
	}
	
	for(int k = 0; k<junction3s.size(); k++)
	{
		junction3s[k]->j2_01.offset = N_old->junction3s[k]->j2_01.offset; 
		junction3s[k]->j2_20.offset = N_old->junction3s[k]->j2_20.offset;
		junction3s[k]->j2_12.offset = N_old->junction3s[k]->j2_12.offset;
		junction3s[k]->j2_21.offset = N_old->junction3s[k]->j2_21.offset;
	}


}



Network:: ~Network()
{
//#pragma omp parallel for
	for(unsigned int k=0; k<channels.size(); k++){delete channels[k];}	
	for(unsigned int k=0; k<junction1s.size(); k++){delete junction1s[k];}
	for(unsigned int k=0; k<junction2s.size(); k++){delete junction2s[k];}
	for(unsigned int k=0; k<junction3s.size(); k++){delete junction3s[k];}
//	delete [] conns;
//	delete [] nodeTypes;	
}


void Network::EulerStep(double dt)
{
	//old way that was first order
	//
	for (unsigned int j = 0; j<junction3s.size(); j++){junction3s[j]->boundaryFluxes();}
	for (unsigned int j = 0; j<junction1s.size(); j++){junction1s[j]->boundaryFluxes();}
	for (unsigned int j = 0; j<junction2s.size(); j++){junction2s[j]->boundaryFluxes();}
	for (int k = 0;k<Nedges; k++){
		//cout<<"k= "<<k<<endl;
		channels[k]->stepEuler(dt);
	}
}

void Network::stepRK3_SSP(double dt)
{
	for(int j=0; j<Nedges; j++)
	{
		for(int i=0; i<channels[j]->N*2; i++)
		{ 
			channels[j]->qhat[i] = channels[j]->q0[i];
		}
	}
	EulerStep(dt);
	EulerStep(dt);
	for(int j=0; j<Nedges; j++)
	{
		channels[j]->n++;
		//cout<<"n thinks it is "<<channels[j]->n<<endl;
		for(int i = 0; i<channels[j]->N; i++)
		{
			for (int k=0; k<2; k++)
			{
				channels[j]->q[channels[j]->idx(k,i)] = 0.75*channels[j]->qhat[channels[j]->idx(k,i)] 
					+.25*channels[j]->q[channels[j]->idx(k,i)];
				channels[j]->q0[channels[j]->idx(k,i)] = channels[j]->q[channels[j]->idx(k,i)];
			}
		}
	}
	EulerStep(dt);
	for(int j=0; j<Nedges; j++)
	{
		for(int i = 0; i<channels[j]->N; i++)
		{
			for (int k=0; k<2; k++)
			{
				channels[j]->q[channels[j]->idx(k,i)] = 1./3.*channels[j]->qhat[channels[j]->idx(k,i)]
				       	+2./3.*channels[j]->q[channels[j]->idx(k,i)];
				channels[j]->q0[channels[j]->idx(k,i)] = channels[j]->q[channels[j]->idx(k,i)];
				channels[j]->q_hist[channels[j]->idx_t(k,i+1,nn)] = channels[j]->q[channels[j]->idx(k,i)];
			}
			channels[j]->p_hist[channels[j]->pj_t(i+1,nn)] = channels[j]->P[i+1];
		}
	}

}


void Network::runForwardProblem(double dt)
{
	nn = 0;
	//int Mi = M<500?1:M/500;
	for(int i=0; i<M; i++)
	{
			
	//	if(i%Mi==0)
		{
	//		printf("current time is = %f s ", (double)nn*dt);
	//		printf("Average Gradient is %f \n", getAveGradH(i));	
	if (WTF)
			printf("%f %d %f \n",0.,i ,getAveGradH(i));	
		}
		nn ++;
		stepRK3_SSP(dt);
		//EulerStep(dt);
	
	}
}
double Network::getTotalVolume()
{
	double v = 0;
	for(unsigned int j = 0; j<channels.size(); j++)
	{
		v += channels[j]->getTheGoddamnVolume();
	}
	return v;
}

double Network::getAveGradH(int i)
{
	double dhdx= 0.;
	for(unsigned int j = 0; j <channels.size(); j++)
	{
		dhdx += channels[j]->getAveGradH(i); 
	}
	return dhdx;
}


