#ifndef SETUPANDRUN_H
#define SETUPANDRUN_H



#include "network.h"
/*
class Coords{
	public: 
		vector<double>xcoords;
		vector<double>ycoords;
		vector<double>elevs;
	Coords(){
		xcoords(xin);
		ycoords(yin);
		elevs(ein);
	}
	void update(vector<double> xin, vector<double> yin, vector<double> ein){
		for(int i=1; i<xin_.size(); i++){
			xcoords[i]= xin[i];
			ycoords[i] = yin[i];
			elevs[i]= ein[i];
		}
	}


};
*/
//for writing output to labeled .tga files for use in POVRAY
void w3d_targa_output_surface(const char* filename,double *fld,int m,int n,double zmin,double zmax);
Network setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_);
void writeOutputTarga(Network &Ntwk, int M, int Mi,vector <int> jIDs, vector<double>xcoords, vector<double>ycoords, vector<double> elevs,double T, int writelogs);
void writeOutputText(Network &Ntwk, int M, int Mi);

//using namespace std;

//Templating magic for streaming into vectors (THANKS ROB!!!!!)
template<typename T>
struct AppendToVector
{
    std::vector<T>& vec;
    AppendToVector(std::vector<T>& vec) : vec(vec) {}
};

template<typename T>
std::istream& operator>>(std::istream& s, const AppendToVector<T>& app)
{
    T val;
    if (s >> val)
        app.vec.push_back(val);
    return s;
}

template<typename T>
AppendToVector<T> appendTo(std::vector<T>& vec)
{
    return AppendToVector<T>(vec);
}

#endif

