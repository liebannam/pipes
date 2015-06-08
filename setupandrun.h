#ifndef SETUPANDRUN_H
#define SETUPANDRUN_H


#include "network.h"
#include "string.h"

//for dealing with Hermite/Fourier representations
void getTimeSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier);
void getCoeffSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier);  //in an ideal world this is the inverse of getTimeSeries...

//for writing output to labeled .tga files for use in POVRAY
void w3d_targa_output_surface(const char* filename,double *fld,int m,int n,double zmin,double zmax);
Network* setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_);
void writeOutputTarga(Network *Ntwk, int M, int Mi,double T, int writelogs);
void writeOutputText(Network *Ntwk, int M, int Mi);


#endif

