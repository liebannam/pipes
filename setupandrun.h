/*\file setupandrun.h
 *\brief setupandrun.h file documentation
 * Contains various useful functions for setting up networks and outputting data:
 *  	--read in network data from .inp file and .config files
 *		--get time series from Hermite or Fourier coefficients
 *		--write output to targa files
 *		--write output to text files
 **/

/*This file is part of Pipes.

    Pipes is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Pipes is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Pipes.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef SETUPANDRUN_H
#define SETUPANDRUN_H


#include "network.h"
#include "string.h"

/**\Recover time series from a Hermite or Fourier representation*/
void getTimeSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier);
/**\brief (Not implemented yet) In an ideal world this is the inverse of getTimeSeries...*/
void getCoeffSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier);  
/** Output targa files (thanks to Chris Rycroft for this function) */
void w3d_targa_output_surface(const char* filename,double *fld,int m,int n,double zmin,double zmax);
/** \brief Set up a pointer to a network class based on input files*/
Network* setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_);
/** \brief Write pressure data to labeled .tga files for use in POVRAY. use smarterputittogether.py to parse the output.**/
void writeOutputTarga(Network *Ntwk, int M, int Mi,double T, int writelogs);
/**\brief Write pressure data to labeled text files for perusal or use in other plotting software.*/
void writeOutputText(Network *Ntwk, int M, int Mi);

#endif

