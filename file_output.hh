#include <cstdio>
#include <cstdlib>
using namespace std;

#ifndef WATER3D_FILE_OUTPUT_HH
#define WATER3D_FILE_OUTPUT_HH

inline void sfwrite(const void *bu,size_t size,size_t count,FILE *f) {
	if(fwrite(bu,size,count,f)!=count) {
		fprintf(stderr,"Error writing data to file");
		exit(1);
	}
}

/** Creates the Truevision Targa header.
 * \param[in] bp a temporary buffer of 18 bytes in which to write the header. 
 * \param[in] tm the horizontal size of the image.
 * \param[in] tn the vertical size of the image. */
inline void w3d_targa_create_header(char *bp,int tm,int tn) {
	*(bp++)=0;*(bp++)=0;*(bp++)=2;*(bp++)=0;
	*(bp++)=0;*(bp++)=0;*(bp++)=0;*(bp++)=0;
	*(bp++)=0;*(bp++)=0;*(bp++)=0;*(bp++)=0;
	*(bp++)=tm&255;*(bp++)=(tm>>8)&255;
	*(bp++)=tn&255;*(bp++)=(tn>>8)&255;
	*(bp++)=24;*(bp++)=0;
}

/** Codes a floating point number as a 24-bit integer within the RGB channels
 * of a Targa file.
 * \param[in,out] bp a pointer to the buffer to write to.
 * \param[in] val the value to code. */
inline void w3d_targa_encode_value(char *&bp,double val) {
	int ht(int(16777215.0*val+0.5));
	if(ht<0) ht=0;if(ht>16777215) ht=16777215;

	// Code the height in the red, green, and blue
	// channels
	*(bp++)=ht&255;
	*(bp++)=(ht>>8)&255;
	*(bp++)=(ht>>16)&255;
}

void w3d_targa_output_surface(const char* filename,double *fld,int m,int n,double zmin,double zmax);
void w3d_compute_min_max(double *fld,int mn,double &zmin,double &zmax);
void w3d_output(const char* filename,double *fld,int m,int n);

#endif
