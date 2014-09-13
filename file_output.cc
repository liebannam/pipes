#include <cmath>

#include "file_output.hh"

/** Outputs a two dimensional array as a Targa image that can be read as a
 * POV-Ray height field. The routine treats the array as periodic and repeats
 * the first row and column in the output image.
 * \param[in] filename the name of the file to write to.
 * \param[in] fld the array to use.
 * \param[in] m the horizontal grid size.
 * \param[in] n the vertical grid size.
 * \param[in] zmin the field value to encode to 0.
 * \param[in] zmax the field value to encode to 16777215. */ 
void w3d_targa_output_surface(const char* filename,double *fld,int m,int n,double zmin,double zmax) {
	char *buf(new char[m<5?18:3*(m+1)]),*bp(buf);
	int i,j;double *fldp;
	double zfac(1/(zmax-zmin));

	// Try opening the output file. If it is unsuccessful, then just return
	// without doing anything.
	FILE *outf(fopen(filename,"wb"));
	if(outf==NULL) return;

	// Output the Targa header
	w3d_targa_create_header(buf,m,n);
	fwrite(buf,1,18,outf);

	// Loop over the lines of the surface
	for(j=0;j<n;j++) {
		bp=buf;
		fldp=fld+j*m;
		for(i=0;i<m;i++) w3d_targa_encode_value(bp,(fldp[i]-zmin)*zfac);
		//w3d_targa_encode_value(bp,(*fldp-zmin)*zfac);
		fwrite(buf,3,m,outf);
	}

	// Write first line again
//	bp=buf;
//	for(i=0;i<m;i++) w3d_targa_encode_value(bp,(fld[i]-zmin)*zfac);
//	w3d_targa_encode_value(bp,(*fld-zmin)*zfac);
//	fwrite(buf,3,m+1,outf);

	// Close the output file and free the temporary buffer
	fclose(outf);
	delete [] buf;
}

/** Computes the minimum and maximum values of an array.
 * \param[in] fld a pointer to the array to consider.
 * \param[in] mn the number of elements in the array.
 * \param[out] zmin a reference in which to return the minimum height.
 * \param[out] zmax a reference in which to return the maximum height. */
void w3d_compute_min_max(double *fld,int mn,double &zmin,double &zmax) {
	double *flde=fld+mn;zmin=zmax=*(fld++);
	while(fld<flde) {
		if(*fld<zmin) zmin=*fld;
		if(*fld>zmax) zmax=*fld;
		fld++;
	}
}

/** Outputs a two dimensional array in the Gnuplot binary format. The routine
 * treats the array as periodic and repeats the first row and column.
 * \param[in] filename the name of the file to write to.
 * \param[in] fld the array to use.
 * \param[in] m the horizontal grid size.
 * \param[in] n the vertical grid size. */
void w3d_output(const char* filename,double *fld,int m,int n) {
	const double tpi=8*atan(1.0);
	double dx=tpi/m,dy=tpi/n;
	int i,j;

	// Open file and write header line 
	FILE *outf;
	outf=fopen(filename,"wb");
	float *fbuf(new float[m+2]),*fp(fbuf);
	double *pp;
	*(fp++)=m+1;for(i=0;i<m;i++) *(fp++)=i*dx;*fp=tpi;
	sfwrite(fbuf,sizeof(float),m+2,outf);
	
	// Write field entries line-by-line
	for(j=0;j<n;j++) {

		// Write header entry
		fp=fbuf;*(fp++)=j*dy;

		// Write a horizontal line to the buffer
		pp=fld+j*m;for(i=0;i<m;i++) *(fp++)=*(pp++);
		*fp=fld[j*m];
		sfwrite(fbuf,sizeof(float),m+2,outf);
	}

	// Since we are working in periodic coordinates, output the first line
	// again at y=2*pi
	fp=fbuf;*(fp++)=tpi;
	pp=fld;for(i=0;i<m;i++) *(fp++)=*(pp++);
	*fp=*fld;
	sfwrite(fbuf,sizeof(float),m+2,outf);

	// Remove temporary memory and close file
	delete [] fbuf;
	fclose(outf);
}
