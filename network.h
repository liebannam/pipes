/****************
* \file network.h
* \brief network.h documentation
* Network creates a network of channels
* constructor:
	Network(int Nnodes_, int *conns_, int Nedges_, int * Ns, double *ws, double *Ls, double *S0s, double *Mrs, double *a0s, double *q0s);

*/

//Things I'd like to add:
//draw() --- draw an ascii pic of the network, with labels
//timeStep(dt) -- possibly refining/taking multiple steps on any pieces that are violating CFL)

#ifndef NETWORK_H
#define NETWORK_H
#include <vector>
#include "channel.h"
#include "omp.h"

/**\brief Class for setting up networks of connected channels.*/
class Network{
	public:
	//number of nodes
	int Nnodes;     
	//number of edges	
	int Nedges;     
	//labels each node as type 1, 2 or 3 (is it associated with a junction1, junction2, etc..)
	std::vector<int> nodeTypes; 
   	//array of connectivity data- form is  [leftend0, rightend0, leftend1, rightend1,...] 
	std::vector<int> conns;    
	//vector of channel class instances 
	std::vector<Channel*> channels;	
	//vector of "1 junctions" (dead ends or outlets where external B.C. is specified)
	std::vector<Junction1*> junction1s; 
	//vector of "2 junctions" (intersection between two edges)
	std::vector<Junction2*> junction2s; 
	//vector of "3 junctions" (intersection between three edges)	
	std::vector<Junction3*> junction3s; 
	//number of timesteps (wow is this sloppy!)
	int M; 
	//which timestep we're currently at
	int nn;  
	//0 for uniform, 1 for Preissman slot
	int channeltype; 
	//constructor
	Network(int Nnodes_, std::vector<int> conns_, int Nedges_, std::vector<int> Ns, std::vector<double> ws, std::vector<double> Ls, std::vector<double> S0s, std::vector<double> Mrs, std::vector<double> a0s, std::vector<double> q0s, int M_,  int channeltype_, double a = 1200.);
	//copy constructor
	Network(const Network &N_old);
	Network(Network *N_old);
	//destructor
	~Network();
	//other handy functions 
	void EulerStep(double dt);
	void stepRK3_SSP(double dt);
	//step everybody in the network from t =0 to t = T
	void runForwardProblem(double dt);
	//total volume (because this is an FV code, V= sum(A(i))*dx))
	double getTotalVolume();
	//average gradient of H(x) over lengths
	double getAveGradH(int i); 
	//get total kinetic energy
	double getKE(int i);
	//get total potential energy
	double getPE(int i);
};
#endif


