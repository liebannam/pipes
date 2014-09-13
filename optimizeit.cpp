/////////////////////////
/*Read in input files and optimize things. THIS IS PRETTY ANNOYING...
*command line syntax is ../setupandrun example.inp example.config
*where:
*       example.inp is an EPANET-style file containing network layout information 
*       example.config is a custom configuration file with run-time parameters (e.g. number of cells, ICs, BCs)
*TO DO: make this program check everything out to ensure the .inp and .config files match!
*
*///////////////////
//#include "network.h"
//#include "levmar.h"
#include "optdetritus.h"
#include "setupandrun.h"
//




int main(int argc, char *argv[] )	
	{
	int writelogs = 0; //change to 1 if the height plots look dull on account of large pressure variations

	//first open .inp file and process information about network layout and components
	char * f1= argv[1], *f2 = argv[2];
	ifstream file1(f1);
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
	//convert diamters to meters....should check for this, maybe?
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
	ifstream file2(f2);
	string morestuff;

	jflag =0;
	pflag = 0;
	first = 0;
	int tflag = 0;
	vector<int> bvaltypes;
	vector<double>bvals;
	vector<double>reflects;
	vector<double> offsets;
	vector<double> valveopens;
	vector<double> offset01s;
	vector<double> offset02s;
	vector<double> offset12s;
	vector<int> Ns;
	vector<double> h0s;
	vector<double> q0s;
	int M, Mi;
	double T;
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
				cout<<"T = "<<T<<endl;
				if(morestuff[0] =='[' && first==0){tflag =0;}
				first = 0;

			}
		 }
	 }
	cout<<"Junction Info:\njunction ID     elevation\n";
	for(int k= 0; k<jIDs.size(); k++)
	{    cout<< jIDs[k]<<"                "<<elevs[k]<<"   "<<xcoords[k]<<" "<<ycoords[k]<<endl;}
	cout<<"Pipe Info:\npipe ID   left_end   right_end    length(m)   diam(mm)    manning coeff  slope\n";
	for(int k= 0; k<pIDs.size(); k++)
	{ printf("%d          %d          %d            %.1f      %.3f        %.4f          %.4f     %.2f    %.2f\n", pIDs[k], conns[2*k], conns[2*k+1],lengths[k], diams[k], Mrs[k], S0s[k], h0s[k], q0s[k]); }


		clock_t start_t, end_t;
		int channeltype =1;  //0 for uniform cross section, 1 for Preissman slot

		int Nnodes = jIDs.size();
		int Nedges = pIDs.size();

		double dt = T/(double)M;

		//make the damn network
		Network Ntwk(Nnodes, conns, Nedges, Ns, diams, lengths, S0s, Mrs, h0s, q0s, M, channeltype);
		
		
		cout<<"Number of 1 junctions is "<<Ntwk.junction1s.size()<<endl;
		cout<<"Number of 2 junctions is "<<Ntwk.junction2s.size()<<endl;
		cout<<"Number of 3 junctions is "<<Ntwk.junction3s.size()<<endl;
		cout<<"Number of edges is "<<Nedges<<endl;
		for(int k = 0; k<Ntwk.channels.size(); k++)
		{
			double a0 = Ntwk.channels[k]->Aofh(h0s[k]);
			Ntwk.channels[k]->setq(a0, q0s[k]);
			Ntwk.channels[k]->setq0(a0, q0s[k]);
		}	
		for(int k = 0; k<Ntwk.junction1s.size(); k++)
		{
			Ntwk.junction1s[k]->bvaltype = bvaltypes[k];
			Ntwk.junction1s[k]->setbVal(bvals[k]);
			cout<<"k = "<<k<<"bval= "<<bvals[k]<<endl;
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
		
		cout<<Ntwk.junction1s[0]->bval;	
		double V0=Ntwk.getTotalVolume();
		start_t = clock();
			
//		Network Ntwk2(Ntwk);
		Ntwk.runForwardProblem(dt);
		end_t = clock();
/*		for(int k=0; k<Nedges; k++){
		//	Ntwk.channels[k]->showGeom();
			Ntwk.channels[k]->showVals(0);
		}
*/
		printf("Elapsed real time is %f\n", (end_t-start_t)/(double)CLOCKS_PER_SEC);	
		printf("Elapsed simulation time is %f\n", dt*(double)(M));
		double V = Ntwk.getTotalVolume();

		cout<<"Network 1\n";
		cout<<"initial volume "<<V0<< "    "<<"Final Volume " <<V<< endl;
		cout<<"dV = "<<V-V0<<endl;
//		double V20 = Ntwk2.getTotalVolume();
//		Ntwk2.runForwardProblem(dt);		
//		double V2 = Ntwk2.getTotalVolume();
//		cout<<"Network 2\n";
//		cout<<"initial volume "<<V20<< "    "<<"Final Volume " <<V2<< endl;
//		cout<<"dV2 = "<<V2-V20<<endl;

//optimization time! at last!
	cout<<"M = "<<M<<endl;
	int ndof = 16;   // degrees of freedom (in Fourier or Hermite modes)
	int modetype = 0;
	int whichnode = 0;
	double Dt = T/(ndof/2-1); //hermite interpolation spacing
	vector<double> h(M+1);
	vector<Real> x0(ndof, 0.);

	bc_opt_dh test1(ndof, M, x0, Ntwk, modetype, T, whichnode);
	cout<<"Made it?\n";
	test1.compute_f();
	double f0 = test1.f;
	cout<<"f0 is!"<<f0<<endl;
//	test1.chkder();
//	test1.solve();
//	test1.dump();
//////////////
//
//	
//	//Ntwk.channels[0]->showVals(1);
	for(int j=0;j<Nedges; j++)
	{
//		for(int k =0; k<Ns[j]; k++){cout<<lengths[j]/Ns[j]*k<<"   "<<Ntwk.channels[j]->hofA(Ntwk.channels[j]->q[k])<<endl;}
		for(int k =0; k<Ntwk.channels[j]->N; k++){cout<<k<<"   "<<Ntwk.channels[j]->q[k]<<endl;}
	}
//	for(int j=0;j<Nedges; j++){
//		for(int k =0; k<N; k++){
//			cout<<dx*k<<"   "<<Ntwk.channels[j]->hofA(Ntwk.channels[j]->q[k])<<" "<<Ntwk.channels[j]->q[k+N]<<endl;
//			//cout<<dx*k<<"   "<<(Ntwk.channels[j]->q[k])<<endl;
//		}
//	}
//	cout<<"Number of 1 junctions is "<<Ntwk.junction1s.size()<<endl;
//	cout<<"Number of 2 junctions is "<<Ntwk.junction2s.size()<<endl;
//	cout<<"Number of 3 junctions is "<<Ntwk.junction3s.size()<<endl;
//	cout<<"Number of edges is "<<Nedges<<endl;
//
//
//
////	cout<<"Final volume distribution ="<<Ntwk.channels[0]->getTheGoddamnVolume()<<" "<<Ntwk.channels[1]->getTheGoddamnVolume()<<" "<<Ntwk.channels[2]->getTheGoddamnVolume()<<endl; 
	printf("triple juncton values are A0 = %f, A1 = %f, A2 = %f\n", Ntwk.channels[0]->q[Ns[0]-1], Ntwk.channels[1]->q[0], Ntwk.channels[2]->q[0]);
//
//
////output heightfields and a textfile "runinfo.txt" of maxvalues to accompany them.
char sdata[] = "../output_data/scalings.txt";
char mdata[] = "../output_data/mapdata.txt";

FILE *fp = fopen(sdata, "w");
	if(fp==NULL)
	{
		fputs("Error opening file!!", stderr);
		return 1;
	}
FILE *fm = fopen(mdata, "w");
	if(fp==NULL)
	{
		fputs("Error opening file!!", stderr);
		return 1;
	}
////warning mess with the following line at your peril because the plot script will epically *&%^ up...
fprintf(fm, "%f\n", T);
for(int k =0; k<Nedges; k++){fprintf(fm, "%d   %d   ", conns[2*k], conns[2*k+1]);}
fprintf(fm, "\n");
for(int k = 0; k<Nedges; k++){fprintf(fm, "%f  ", diams[k]/2.);cout<<diams[k]/3.<<endl;}	
fprintf(fm, "\n");
for(int k =0; k<Nnodes; k++){
	fprintf(fm, "%d   %f    %f    %lf\n  ",jIDs[k], xcoords[k], ycoords[k],elevs[k]);
}

cout<<"number of writes is  "<<M/Mi<<endl;
for(int ii=0;ii<M; ii+=Mi)
{
//	printf("i=%d\n", ii);
	fprintf(fp, "%d    ", ii/Mi);
	for(int kk = 0; kk<Nedges;kk++)
	{
		char filename[100];
		sprintf(filename,"../output_data/out%d_%03d.tga", kk,ii/Mi);
		int mm = Ns[kk];                        //x direction
		int nn = mm/2;	                        //y direction
		double myfld[mm*nn];
		double val;	
		double zmin =0, zmax = 0;
		for(int j = 0;j<nn; j++)
		{
			for (int i= 0; i<mm; i++)
			{
				if(writelogs)
				{
					val = Ntwk.channels[kk]->hofA(Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)]);
				}
				else
				{
					val = log(Ntwk.channels[kk]->hofA(Ntwk.channels[kk]->q_hist[Ntwk.channels[kk]->idx_t(0,i+1, ii)])+1);
				}
				myfld[i+mm*j] = val;
			       //cout<<val<<"   ";	
				zmin = fmin(zmin, val);
				zmax = fmax(zmax, val);
			}
		}

		fprintf(fp, "%f      ", zmax);
		zmin = 0;
	//	zmax = Ntwk.channels[kk]->w;
	//	printf("zmax =  %f  \n", zmax);
		w3d_targa_output_surface(filename, myfld,mm,nn,zmin,zmax);
		}
  // printf("\n");
   fprintf(fp, "\n");
}
fclose(fp);
fclose(fm);

//for (int k = 0; k<M/Mi; k++)
//{
//	printf("%f     %f   \n", dt*float(k*Mi), Ntwk.channels[1]->hofA(Ntwk.channels[1]->q_hist[Ntwk.channels[1]->idx_t(0,N/2,k*Mi)]));}
//	printf("%f     %f   \n", dt*float(k*Mi), Ntwk.channels[0]->hofA(Ntwk.channels[0]->q_hist[Ntwk.channels[0]->idx_t(0,N/2,k*Mi)]));}

}

