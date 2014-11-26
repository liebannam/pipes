/////////////////////////
/*Driver for main simulation
*command line syntax is ../setupandrun example.inp example.config
*where:
*       example.inp is an EPANET-style file containing network layout information 
*       example.config is a custom configuration file with run-time parameters (e.g. number of cells, ICs, BCs)
*TO DO: make this program check everything out to ensure the .inp and .config files match!
*
*///////////////////
#include "setupandrun.h"

//void getTimeSeries(vector< vector <Real> > & bvals, vector< vector <Real> > &x, const int m, const int M, double T, int Fourier)
//same thing as getTimeSeries, only now does it for multiple nodes.
//
//{
//	vector<Real> tmp_b(M+1);
//	int K = xbvals.size();
//	for (int k = 0; k<K; k++)
//	{
//		vector<Real> tmp_x(x[k]);
//		getTimeSeries(tmp_b, tmp_x, m, M, T, Fourier);
//	}
//}
void getTimeSeries(vector<Real> & bvals, vector<Real> &x, const int m, const int M, double T, int Fourier)     	
//bvals is time series of M+1 values at t =0, T/M, ...T
//x is a length m vector. either contains Fourier modes (Fourier==1) or of Hermite spline components (Fourier ==0)
//we have m/2 values in Fourier series or m/2 spline interpolation points. Interpolation dt = 1/(m/2-1)
//
{
	if(Fourier)	//Discrete Fourier modes:
	//x[k] = a_k (k<=m/2)
	//x[k] = b_j (where j =k-m/2, for k>m/2)
	//      x =     [a_0, a_1, ...a_m/2, b_1,   b2   ...b_(m/2-1)] 
	//                |     |       |     |     |
	//index, k        0     1      m/2  m/2+1  m/2+2       m-1     
	{
		double t;
		for(int nn = 0; nn<M+1; nn++)
		{
			bvals[nn] = 0.5*x[0];
			t = (double)(nn)/(double)M;
			for (int k = 1; k<m/2; k++)
			{
				bvals[nn] += x[k]*cos(2.*PI*(double)k*t) + x[k+m/2]*sin(2.*PI*(double)k*t);
		//		cout<<k<< " "<< x[k+m/2+1]<<endl;
			}
			bvals[nn] +=0.5*x[m/2]*cos(PI*(double)m*t);
		//	printf("n = %d, t/T = %f, gah %f\n", n, t, bvals[n]);
		//	cout<<"!!!!!"<<bvals[n]<<endl;
		}
	}
	else//Hermite Spline
	//Hermite spline evaluation, adapted from C. Rycroft's code. 
	// Evaluates spline x at point t, where
	// x = [f(t0), f'(t0), f(t1), f'(t1), ...f(t_m), f'(tm)]
	//      |       |        |      |          |       |
	// i=   0       1        2      3          2*m    2*m+1
	// x0 = 0, xm = T

	{
		double t, Dt,dt;
		dt = T/M;//time series dt
		Dt = T/((double)m/2.-1.);//spline dt. Careful! This is NOT the dt in your problem!
		for(int nn = 0; nn<M+1; nn++)
		{
			t= nn*dt;
			int j = int((t)/Dt);
			t = t/Dt-j; 
		        cout<<Dt<<endl;	
			if(j>m+1 || j<0) printf("warning! t=%f out of range!!, Dt = %f, T = %d, m = %d\n",t,Dt,T,m);
			// Calculate square and cube, and pointer to the values to use
			double t2=t*t,t3=t2*t;
			// Calculate the value of the spline function
			if(j==m){bvals[nn] = x[2*j];}
			else{bvals[nn] = x[2*j]*(2*t3-3*t2+1)+x[2*j+1]*(t3-2*t2+t)*Dt+x[2*j+2]*(-2*t3+3*t2)+x[2*j+3]*(t3-t2)*Dt;}
		}
	}
}




Network setupNetwork(char *finp, char *fconfig, int &M, int &Mi, double &T, int channeltype_)
{
	//first open .inp file and process information about network layout and components
	ifstream file1(finp);
	char BC_filename[50];
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
	ifstream file2(fconfig);
	string morestuff;

	jflag =0;
	pflag = 0;
	first = 0;
	int tflag = 0;
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
	vector < vector <Real > > xbval;
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
				int tmp,type;
				stringstream ss(morestuff); // Insert the string into a stream
				ss>>tmp>>type;
				if(type==1)//next three columns matter
				{
					ss>>appendTo(bvaltypes)>>appendTo(bvals)>>appendTo(reflects);
				}
				else if(type ==2)//skip next three columns; next 2 matter
				{	
					ss>>tmp>>tmp>>tmp>>appendTo(offsets)>>appendTo(valveopens);
				}
				else if(type ==3)//skip next 5, next 2 matter
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
				ss>>T>>M>>Mi;
				cout<<morestuff<<endl;
				cout<<"T = "<<T<<" tflag ="<<tflag<<endl;
				if(first==0){tflag =0;}
				first = 0;

			}
		 }
		 if(strncmp(morestuff.c_str(),"[BC_FILENAME]",13)==0)
		 {
			string tmp(morestuff.c_str());
			int L = tmp.length()/sizeof(tmp[0]);
			cout<<"L="<<L<<endl;
			int go = 0;
			int count  =0;
			for (int i = 0; i<L; i++){
				if (go>0){
					if(tmp[i]==' ')
					{
						count ++;
					}
					else{
						BC_filename[i-count] = tmp[i];
					}
				}
				else{
					count++;
//					char tmp2 = tmp[i];
					if (tmp[i]==']'){
						cout<<"yay\n";
						go=1;
					}
				}
			}
			cout<<"BC Filename is:"<<BC_filename<<endl;
			string evenmorestuff;
			ifstream fbc(BC_filename);
			string trash;
			vector <Real> tmp_x;
			//ifstream fbc("../indata/bcs2.txt");
			while (getline(fbc, evenmorestuff, '\n'))
				{
					stringstream ss(evenmorestuff);
					if(strncmp(evenmorestuff.c_str(),"modetype",8)==0)
					{
						ss>>trash>>appendTo(modetype);
						cout<<"modetype = "<<modetype[0]<<endl;
					}
					if(strncmp(evenmorestuff.c_str(),"whichnode",9)==0)
					{
						ss>>trash>>appendTo(whichnode);
						cout<<"whichnode = "<<whichnode[0]<<endl;
					}
					else{
						ss>>appendTo(tmp_x);
					}
					//cout<<evenmorestuff;
				}	
			fbc.close();
			xbval.push_back(tmp_x);
		}

	 }
	cout<<"Junction Info:\njunction ID     elevation\n";
	for(int k= 0; k<jIDs.size(); k++)
	{    cout<< jIDs[k]<<"                "<<elevs[k]<<"   "<<xcoords[k]<<" "<<ycoords[k]<<endl;}
	cout<<"Pipe Info:\npipe ID   left_end   right_end    length(m)   diam(mm)    manning coeff  slope\n";
	for(int k= 0; k<pIDs.size(); k++)
	{ 
		printf("%d          %d          %d            %.1f      %.3f        %.4f          %.4f     %.2f    %.2f\n", pIDs[k], conns[2*k], conns[2*k+1],lengths[k], diams[k], Mrs[k], S0s[k], h0s[k], q0s[k]); }


	int channeltype =channeltype_;  //0 for uniform cross section, 1 for Preissman slot

	int Nnodes = jIDs.size();
	int Nedges = pIDs.size();

	double dt = T/(double)M;
	cout<<"T = "<<T<<" M = "<<M<<endl;
	//make the damn return network
	Network Ntwk(Nnodes, conns, Nedges, Ns, diams, lengths, S0s, Mrs, h0s, q0s, M, channeltype);
	
	
	cout<<"Number of 1 junctions is "<<Ntwk.junction1s.size()<<endl;
	cout<<"Number of 2 junctions is "<<Ntwk.junction2s.size()<<endl;
	cout<<"Number of 3 junctions is "<<Ntwk.junction3s.size()<<endl;
	cout<<"Number of edges is "<<Nedges<<endl;
	for(int k = 0; k<Ntwk.channels.size(); k++)
	{
		double a0 = Ntwk.channels[k]->AofH(h0s[k],false);
		printf("h0 = %f, a0 = %.10f, h(a0) = %f\n", h0s[k],a0, Ntwk.channels[k]->HofA(a0,false) );
		Ntwk.channels[k]->setq(a0, q0s[k]);
		Ntwk.channels[k]->setq0(a0, q0s[k]);
	}	
	for(int k = 0; k<Ntwk.junction1s.size(); k++)
	{
		Ntwk.junction1s[k]->bvaltype = bvaltypes[k];
		Ntwk.junction1s[k]->setbVal(bvals[k]);
		Ntwk.junction1s[k]->reflect = reflects[k];
	}
	for(int k = 0; k<Ntwk.junction2s.size(); k++)
	{
		Ntwk.junction2s[k]->offset = offsets[k];
		Ntwk.junction2s[k]->valveopen = valveopens[k];
	}
	for(int k = 0; k<Ntwk.junction3s.size(); k++)
	{
		Ntwk.junction3s[k]->j2_01.offset = offset01s[k];
		Ntwk.junction3s[k]->j2_20.offset = offset02s[k];
		Ntwk.junction3s[k]->j2_12.offset = offset12s[k];
		Ntwk.junction3s[k]->j2_21.offset = offset12s[k];
	}
	if (xbval.size()>0){
		printf("Setting boundary values for using %s modes:\n", (modetype[0]?"Fourier":"Hermite"));
	//	for (int ii = 0; ii<xbval.size();ii++)cout<<"i = "<<ii<<" m[i]= "<<xbval[ii]<<endl;
		int mm = xbval[0].size();
		vector<Real> bvalsfancy(M+1,0.);
		cout<<"m = "<<mm<<"T= "<<T<<endl;;
		for(int ii= 0; ii<xbval.size(); ii++)
		{
			getTimeSeries(bvalsfancy, xbval[ii], mm, M, T, modetype[0]);
			Ntwk.junction1s[whichnode[ii]]->setbVal(bvalsfancy);
			cout<<"Setting bvals for node"<<whichnode[ii]<<endl;
			for (int kk = 0; kk<bvalsfancy.size(); kk++)
				cout<<bvalsfancy[ii]<<endl;	
		}

	}

	
	char mdata[] = "../output_data/mapdata.txt";
	FILE *fm = fopen(mdata, "w");
	if(fm==NULL)
	{
		fputs("Error opening file!!", stderr);
		throw;
	}

	////warning mess with the following line at your peril because the plot script will epically *&%^ up...
	fprintf(fm, "%f\n", T);
	for(int k =0; k<Nedges; k++){fprintf(fm, "%d   %d   ", Ntwk.conns[2*k], Ntwk.conns[2*k+1]);}
	fprintf(fm, "\n");
	for(int k = 0; k<Nedges; k++){fprintf(fm, "%f  ", (Ntwk.channels[k]->w)/2.);}	
	fprintf(fm, "\n");
	for(int k =0; k<Nnodes; k++){
		fprintf(fm, "%d   %f    %f    %lf\n  ",jIDs[k], xcoords[k], ycoords[k],elevs[k]);
	}
	fclose(fm);
	printf("network data written to %s\n", mdata);
	return Ntwk;
}


////output heightfields and a textfile "runinfo.txt" of maxvalues to accompany them.
void writeOutputTarga(Network &Ntwk, int M, int Mi, double T, int writelogs)
{
	char sdata[] = "../output_data/scalings.txt";
	int Nedges = Ntwk.Nedges;
	int Nnodes = Ntwk.Nnodes;
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
			sprintf(filename,"../output_data/out%d_%03d.tga", kk,ii/Mi);
			int mm = Ntwk.channels[kk]->N;                        //x direction
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
						double a = Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)];
						val = log(Ntwk.channels[kk]->HofA(a,p)+1);
					}
					else{
						bool p = false;
						double a = Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)];
						val = Ntwk.channels[kk]->HofA(a,p);
					}	
					myfld[i+mm*j] = val;
				       //cout<<val<<"   ";	
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

void writeOutputText(Network &Ntwk, int M, int Mi)	
{
	
	int Nedges= Ntwk.Nedges;
	cout<<"number of writes is  "<<M/Mi<<endl;
	for(int ii=0;ii<M; ii+=Mi)
	{
		double val;
		char fname[100];
		sprintf(fname,"../output_data/txt_out%03d.txt",ii/Mi);		
		FILE *fd = fopen(fname, "w");
		if(fd==NULL)
		{
			fputs("Error opening file!!", stderr);
			throw;	
		}
		for(int kk=0; kk<Nedges; kk++)
		{
			int NN = Ntwk.channels[kk]->N;
			for(int j = 0;j<NN; j++)
			{
				val = m_to_psi*Ntwk.channels[kk]->fakehofA(Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,j+1, ii)], false);
				fprintf(fd, "%.10f     ", val);

			}
		fprintf(fd, "\n");	
		}
		fclose(fd);
	}
	printf("Writing psi data to, e.g., ../output_data/txt_out000.txt\n number of writes is %d\n", M/Mi); 
/*	double Al = Ntwk.channels[0]->q[49];
	double Ar =  Ntwk.channels[0]->q[49];
	printf("Al = %f, Ar = %f\n", Al, Ar);
	testcrappy_f(Al, Ar, Ntwk.channels[0]->w, Ntwk.channels[0]->Ts,Ntwk.channels[0]->At);*/
}


