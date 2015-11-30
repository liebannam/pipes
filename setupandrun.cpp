/**\file setupandrun.cpp
 * \brief function definitions for setupandrun.h
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


#include "setupandrun.h"

/*recover time series
 * \param[in] bvals: time series of M+1 values at t =0, T/M, ...T
 * \param[in] x: length m vector. 
 *				either contains Fourier modes (Fourier==1) 
 *				or Hermite spline components (Fourier ==0)
 *\param[in] m: interpolation points. 
 *				have m/2 values in Fourier series or m/2 spline interpolation points. 
 *				Interpolation dt = 1/(m/2-1)
 *\param[in] M: output time series with M+1 points
 *\param[in] T: ending time. i.e. time series is assumed to be [0,T]
 *\param[in] Fourier: flag for interpolation type.
				Fourier =1 means Fourier interpolation
				Fourier =0 means Hermite interpoloation
 */
void getTimeSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier)     	
{
	/**Discrete Fourier modes:
	x[k] = a_k (k<=m/2)
	x[k] = b_j (where j =k-m/2, for k>m/2)
	      x =     [a_0, a_1, ...a_m/2, b_1,   b2   ...b_(m/2-1)] 
	                |     |       |     |     |
	index, k        0     1      m/2  m/2+1  m/2+2       m-1 */    
	if(Fourier)	
	{
		double t;
		for(int nn = 0; nn<M+1; nn++)
		{
			bvals[nn] = 0.5*x[0];
			t = (double)(nn)/(double)M;
			for (int k = 1; k<m/2; k++)
			{
				bvals[nn] += x[k]*cos(2.*PI*(double)k*t) + x[k+m/2]*sin(2.*PI*(double)k*t);
			}
			bvals[nn] +=0.5*x[m/2]*cos(PI*(double)m*t);
		}
	}
	/**Hermite Spline
	 Hermite spline evaluation, adapted from C. Rycroft's code. 
	 Evaluates spline x at point t, where
	 x = [f(t0), f'(t0), f(t1), f'(t1), ...f(t_m), f'(tm)]
	      |       |        |      |          |       |
	 i=   0       1        2      3          2*m    2*m+1
	 x0 = 0, xm = T*/
	else
	{
		double t, Dt,dt;
		dt = T/M;//time series dt
		Dt = T/((double)m/2.-1.);//spline dt. Careful! This is NOT the dt in your problem!
		for(int nn = 0; nn<M+1; nn++)
		{
			t= nn*dt;
			int j = int((t)/Dt);
			t = t/Dt-j; 
			if(j>m+1 || j<0) printf("warning! t=%f out of range!!, Dt = %f, T = %f, m = %d\n",t,Dt,T,m);
			// Calculate square and cube, and pointer to the values to use
			double t2=t*t,t3=t2*t;
			// Calculate the value of the spline function
			if(j==m){bvals[nn] = x[2*j];}
			else{
				bvals[nn] = x[2*j]*(2*t3-3*t2+1)+x[2*j+1]*(t3-2*t2+t)*Dt+x[2*j+2]*(-2*t3+3*t2)+x[2*j+3]*(t3-t2)*Dt;
			}
		}
	}
}

//turn time series into hermite or Fourier interpolateion (not implemented yet)
void getCoeffSeries(vector< vector <Real> > & bvals, vector< vector <Real> > &x, const int m, const int M, double T, int Fourier)
// in an ideal world this is the inverse of getTimeSeries...
{
	if (Fourier){printf("Whoops! not implemented yet");}
	else//Hermite Spline Coeffs (holy shit how do I do this!?)
	{
		//yikes.
	}
}

/**set up a network based on contents of files finp, fconfig.
 * \param[in] finp: path to input file with extension .inp
 * \param[in] fconfig: path to configuration file with extension .config
 * \param[in] M: number of simulaiton time steps
 * \param[in] Mi: number of time steps between writes to output files (1<=Mi<=M)
 * \param[in] T: total simulation time (s)
 * \param[in] T: which channel cross section to use. 0 for uniform, 1 for Preissman slot
 * */
Network* setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_)
{
	//first open .inp file and process information about network layout and components
	double default_a = 12.; //pressure wave speed if not specified in INP file
	ifstream file1(finp);
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
	double a = 0;
	ifstream file2(fconfig);
	string morestuff;
	int Nmodes = 0;
	jflag =0;
	pflag = 0;
	first = 0;
	int tflag = 0;
	int count2 = 0;
	vector<int> bvaltypes;
	vector<double> bvals;
	vector< vector<Real> > bvalsfancy;
	vector<double> reflects;
	vector<double> offsets;
	vector<double> valveopens;
	vector<double> offset01s;
	vector<double> offset02s;
	vector<double> offset12s;
	vector<int> Ns;
	vector<double> h0s;
	vector<double> q0s;
	vector < vector <Real > > xbval(jIDs.size(),vector<Real> (0,0));
	vector <int> whichnode;
	vector <int> modetype;
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
				int tmp=0,btype=0;
				stringstream ss(morestuff); // Insert the string into a stream
				ss>>tmp>>btype;
				printf("btype = %d\n",btype);
				if(btype==1)//next three columns matter
				{
					ss>>appendTo(bvaltypes)>>appendTo(bvals)>>appendTo(reflects);
				}
				else if(btype ==2)//skip next three columns; next 2 matter
				{	
					ss>>tmp>>tmp>>tmp>>appendTo(offsets)>>appendTo(valveopens);
				}
				else if(btype ==3)//skip next 5, next 2 matter
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
				ss>>T>>M>>Mi>>a;
				if (a ==0) a = default_a;
				//cout<<morestuff<<endl;
				cout<<"T = "<<T<<" tflag ="<<tflag<<endl;
				if(first==0){tflag =0;}
				first = 0;

			}
		 }
		 if(strncmp(morestuff.c_str(),"[BC_FILENAME]",13)==0)
		 {
			string tmp(morestuff.c_str());
			char BC_filename[200];
			strcpy(BC_filename,morestuff.c_str()+13);
			cout<<"BC Filename is:"<<BC_filename<<endl;
			//cout<<strncmp(BC_filename,"../indata/bcs2.txt",18)<<"FUCK YOU"<<endl;
			//cout<<strncmp(BC_filename,"~/Dropbox/Research/Network7.0/indata/bcs2.txt",18)<<"FUCK YOU"<<endl;
			//cout<<sizeof(BC_filename)/sizeof(BC_filename[0])<<endl;
			string evenmorestuff;
			string trash;
			vector <Real> tmp_x;
			count2++;
			//ifstream fbc("../indata/bcs2.txt");	
			ifstream fbc(BC_filename);
			if (!fbc){cout<<"Never mind! I don't feel like opening"<<BC_filename<<endl;}
			else{
			while (getline(fbc, evenmorestuff, '\n'))
				{
					//cout<<evenmorestuff<<endl;
					stringstream ss(evenmorestuff);
					if(strncmp(evenmorestuff.c_str(),"modetype",8)==0)
					{
						ss>>trash>>appendTo(modetype);
						cout<<"modetype = "<<modetype[0]<<endl;
					}
					else if(strncmp(evenmorestuff.c_str(),"whichnode",9)==0)
					{
						ss>>trash>>appendTo(whichnode);
						cout<<"whichnode = "<<whichnode[0]<<endl;
											}
					else{
						ss>>appendTo(tmp_x);
					}
				}	
			fbc.close();
			Nmodes = tmp_x.size();
			for (int jj = 0; jj<tmp_x.size();jj++){
				xbval[whichnode[count2-1]].push_back(tmp_x[jj]);
			}
			//cout<<"hmmm"<<count2-1<<endl;
			//cout<<"tmp_x.size = "<<tmp_x.size()<<"  xbval.size() ="<<xbval.size()<<endl;
			}
		}

	 }
	cout<<"Junction Info:\njunction ID   elev.  xcoord     ycoord \n";
	for(int k= 0; k<jIDs.size(); k++)
	{    cout<< jIDs[k]<<"                "<<elevs[k]<<"      "<<xcoords[k]<<"     "<<ycoords[k]<<endl;}
	cout<<"Pipe Info:\npipe ID   left_end   right_end    length(m)   diam(m)    manning coeff    slope    h0        q0\n";
	for(int k= 0; k<pIDs.size(); k++)
	{ 
		printf("%d          %d          %d            %.1f      %.3f        %.4f          %.4f     %.2f    %.2f\n", pIDs[k], conns[2*k], conns[2*k+1],lengths[k], diams[k], Mrs[k], S0s[k], h0s[k], q0s[k]); }


	int channeltype =channeltype_;  //0 for uniform cross section, 1 for Preissman slot

	int Nnodes = jIDs.size();
	int Nedges = pIDs.size();

//	double dt = T/(double)M;
	cout<<"T = "<<T<<" M = "<<M<<endl;
	cout<<"a = "<<a<<endl;
	//make the damn return network
	Network *Ntwk = new Network(Nnodes, conns, Nedges, Ns, diams, lengths, S0s, Mrs, h0s, q0s, M, channeltype,a);
	//proclaim some stuff about the damn return network
	cout<<"Number of 1 junctions is "<<Ntwk->junction1s.size()<<endl;
	cout<<"Number of 2 junctions is "<<Ntwk->junction2s.size()<<endl;
	cout<<"Number of 3 junctions is "<<Ntwk->junction3s.size()<<endl;
	cout<<"Number of edges is "<<Nedges<<endl;
	for(int k = 0; k<Ntwk->channels.size(); k++)
	{
		double a0 = Ntwk->channels[k]->AofH(h0s[k],false);
		//printf("h0 = %f, a0 = %.10f, h(a0) = %f\n", h0s[k],a0, Ntwk->channels[k]->HofA(a0,false) );
		Ntwk->channels[k]->setq(a0, q0s[k]);
		Ntwk->channels[k]->setq0(a0, q0s[k]);
	}
	for(int k = 0; k<Ntwk->junction1s.size(); k++)
	{
		Ntwk->junction1s[k]->bvaltype = bvaltypes[k];
		Ntwk->junction1s[k]->setbVal(bvals[k]);
		Ntwk->junction1s[k]->reflect = reflects[k];
		cout<<"k = "<<k<<" reflect = "<<reflects[k]<<endl;
	}
	for(int k = 0; k<Ntwk->junction2s.size(); k++)
	{
		Ntwk->junction2s[k]->offset = offsets[k];
		Ntwk->junction2s[k]->valveopen = valveopens[k];
	}
	for(int k = 0; k<Ntwk->junction3s.size(); k++)
	{
		Ntwk->junction3s[k]->j2_01.offset = offset01s[k];
		Ntwk->junction3s[k]->j2_20.offset = offset02s[k];
		Ntwk->junction3s[k]->j2_12.offset = offset12s[k];
		Ntwk->junction3s[k]->j2_21.offset = offset12s[k];
	}
	cout<<"number of specified boundaries is "<<whichnode.size()<<endl;
	for(int ii= 0; ii<whichnode.size(); ii++)
		{

			cout<<"m = "<<Nmodes<<"T= "<<T<<endl;
			vector<Real> bvalsfancy(M+1,0.);
			printf("Setting boundary values for using %s modes:\n", (modetype[0]?"Fourier":"Hermite"));
			int DAMN = xbval[whichnode[ii]].size();
			vector <Real>xfake(DAMN);
			if(DAMN>0){
			for (int damn =0; damn<DAMN;damn++){
			        xfake[damn] = xbval[whichnode[ii]][damn];
			}
			getTimeSeries(bvalsfancy, xfake, Nmodes, M, T, modetype[ii]);
			Ntwk->junction1s[whichnode[ii]]->setbVal(bvalsfancy);
			cout<<"Setting bvals for node "<<whichnode[ii]<<endl;
		}
		}
	
	#ifdef commandline
	char mdata[] = "../output_data/mapdata.txt";
	#else
	char mdata[] = "output_data/mapdata.txt";
	#endif
	cout<<mdata<<endl;
	FILE *fm = fopen(mdata, "w");
	if(fm==NULL)
	{
		fputs("Error opening file!!", stderr);
		throw;
	}

	//printf("are we there yet!?!?\n");
	////warning mess with the following line at your peril because the plot script will epically *&%^ up...
	fprintf(fm, "%f\n", T);
	for(int k =0; k<Nedges; k++){fprintf(fm, "%d   %d   ", Ntwk->conns[2*k], Ntwk->conns[2*k+1]);}
	fprintf(fm, "\n");
	for(int k = 0; k<Nedges; k++){fprintf(fm, "%f  ", (Ntwk->channels[k]->w)/2.);}	
	fprintf(fm, "\n");
	for(int k =0; k<Nnodes; k++){
		fprintf(fm, "%d   %f    %f    %lf\n  ",jIDs[k], xcoords[k], ycoords[k],elevs[k]);
	}
	fclose(fm);
	printf("network data written to %s\n", mdata);
	return Ntwk;
}


/**output heightfields and a textfile "runinfo.txt" of maxvalues to accompany them.
 * \param[in] Ntwk a network class whose data we want to write
 * \param[in] M: number of simulation time steps
 * \param[in] Mi: number of time steps between writes to output files (1<=Mi<=M)
 * \param[in] T: total simulation time (s)
 * \param[in] writelogs: whether to write log(1+p) or just p, where p is pressure (to get 
 *			the output data consists of either 
 *			writelogs=0: pbar (the average pressure head, reported in m) 
 *		 or writelogs=1: log(pbar+1) to make colorbars more informative if there's large variation 
 * */
void writeOutputTarga(Network *Ntwk, int M, int Mi, double T, int writelogs)
{
//	char sdata[] = "~/Dropbox/Research/Network7.0/output_data/scalings.txt";
	#ifdef commandline
	char sdata[] = "../output_data/scalings.txt";
	#else
	char sdata[] = "output_data/scalings.txt";
	#endif
	int Nedges = Ntwk->Nedges;
	//int Nnodes = Ntwk.Nnodes;
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
			#ifdef commandline
			sprintf(filename,"../output_data/out%d_%03d.tga", kk,ii/Mi);
			#else
			sprintf(filename,"output_data/out%d_%03d.tga", kk,ii/Mi);
			#endif
			int mm = Ntwk->channels[kk]->N;                        //x direction
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
						double a = Ntwk->channels[kk]->q_hist[Ntwk->channels[kk]->idx_t(0,i+1, ii)];
						val = log(Ntwk->channels[kk]->pbar(a,p)+1);
					}
					else{
						bool p = false;
						double a = Ntwk->channels[kk]->q_hist[Ntwk->channels[kk]->idx_t(0,i+1, ii)];
						val = Ntwk->channels[kk]->pbar(a,p);
					//	printf("pipe %d, i = %d a = %f pbar = %f ",kk,i,a,val);	
					}	
					myfld[i+mm*j] = val;
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

/**Output pressure head (m) to text files
 * \param[in] Ntwk a network class whose data we want to write
 * \param[in] M: number of simulation time steps
 * \param[in] Mi: number of time steps between writes to output files (1<=Mi<=M)
  */
void writeOutputText(Network *Ntwk, int M, int Mi)	
{
	
	int Nedges= Ntwk->Nedges;
	cout<<"number of writes is  "<<M/Mi<<endl;
	for(int ii=0;ii<M; ii+=Mi)
	{
		double val;
		char fname[100];
		#ifdef commandline
		sprintf(fname,"../output_data/txt_out%03d.txt",ii/Mi);		
		#else
		sprintf(fname,"output_data/txt_out%03d.txt",ii/Mi);		
		#endif
		FILE *fd = fopen(fname, "w");
		if(fd==NULL)
		{
			fputs("Error opening file!!", stderr);
			throw;	
		}
		for(int kk=0; kk<Nedges; kk++)
		{
			int NN = Ntwk->channels[kk]->N;
			for(int j = 0;j<NN; j++)
			{
				val = Ntwk->channels[kk]->pbar(Ntwk->channels[kk]->q_hist[Ntwk->channels[kk]->idx_t(0,j+1, ii)], false);
				fprintf(fd, "%.10f     ", val);
			}
		fprintf(fd, "\n");	
		}
		fclose(fd);
	}
	printf("Writing pressure head data to, e.g., ../output_data/txt_out000.txt\n number of writes is %d\n", M/Mi); 
}


