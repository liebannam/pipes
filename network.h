/****************
* class for creating networks of channels
* constructor:

	Network(int Nnodes_, int *conns_, int Nedges_, int * Ns, double *ws, double *Ls, double *S0s, double *Mrs, double *a0s, double *q0s);

* Nnodes is the number of nodes, which are presumed to be numbered [0,1, ...Nnodes]. 
* conns = [leftend0, rightend0, leftend1, rightend1,...] (an array of length 2*Nedges)
* so if an element of conns has an (odd/even) index, it corresponds to a (left/right) end. 
* Nedges = length(conns)/2. 
* the edge corresponding to the ith row of conns is stored in the ith element of channels.
* Ns, ws, Ls, etc are arrays of length Nedges containing parameters for each edge:
* Ns - has number of grid points
* ws - has channel widths
* Ls - has channel lengths
* S0s - has channel slopes
* Mrs - has Manning roughness coeffs
* a0s- has initial area values
* q0s - has ininital q0 values
* M = number of time steps
* destructor is kind of involved because the class does shady things with vectors of pointers to classes, and many creatures born with "new"
 *******/

//Things I'd like to add:
//draw() --- draw an ascii pic of the network, with labels
//timeStep(dt) -- possibly refining/taking multiple steps on any pieces that are violating CFL)
//write(filename) -- write all data to a text file.  Should have a corresponding python script to turn this into properly labeled figures
//eventually this class needs to be templated-- i.e. you should be able to specify an edge object which could be a channel or something totally different!!!


#ifndef NETWORK_H
#define NETWORK_H

#include <vector>
#include "channel.h"
class Network{


	public:
	int Nnodes;     //number of nodes
	int Nedges;     //number of edges	
	std::vector<int> nodeTypes; //labels each node as type 1, 2 or 3 (is it associated with a junction1, junction2, etc..)
	std::vector<int> conns;     //array of connectivity data- form is  [leftend0, rightend0, leftend1, rightend1,...] 
	std::vector<Channel*> channels;	//vector of channel class instances 
	std::vector<Junction1*> junction1s; //vector of "1 junctions" (dead ends)
	std::vector<Junction2*> junction2s; //vector of "2 junctions" (intersection between two edges)
	std::vector<Junction3*> junction3s; //vector of "3 junctions" (intersection between three edges)	
	int M; //number of timesteps (wow is this sloppy!)
	int nn;  //which timestep we're currently at...
	int channeltype; //0 for uniform, 1 for Preissman slot
	//constructor
	Network(int Nnodes_, std::vector<int> conns_, int Nedges_, std::vector<int> Ns, std::vector<double> ws, std::vector<double> Ls, 
			std::vector<double> S0s, std::vector<double> Mrs, std::vector<double> a0s, std::vector<double> q0s, int M_,  int channeltype_);
	//copy constructor goddamnit, I have to do this now, and do it perfectly :(
	Network(const Network &N_old);
	//destructor
	~Network();
	//other handy functions coming soon...
	void EulerStep(double dt);
	void stepRK3_SSP(double dt);
	void runForwardProblem(double dt);//step dynamic variables from t =0 to t = T
	double getTotalVolume();//total volume (because this is an FV code, V= sum(A(i))*dx))
	double getAveGradH(double dt); //average gradient of h(x) over lengths
};



#endif


