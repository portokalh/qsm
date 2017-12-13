/* calcphase.c

This is a C implementation of part of Chunlei Lu's Matlab MRI phase imaging program MapFreq.m.
Some options have been removed (multiple coils) because the users told me they aren't needed.
It is expected but not required that this will be called from a Matlab program that reads the
p-file header and finds the appropriate size parameters for the arrays. The expected command-line
arguments are clearly identified in the code below.

This program makes use of the the FFTW library, which must be installed before this program
is compiled. In particular, the single-precision library (fftw3f) must be present. Also, for
big arrays, it is necesary to specify that a 64 bit architecture is to be used. So the configure
call might look like this:
./configure --enable-float --with-gcc-arch=x86_64 CFLAGS="-arch x86_64" CC="gcc" 
See http://www.fftw.org/fftw3_doc/Installation-on-Unix.html#Installation-on-Unix for details.
You will need root access to install.

To compile with gcc, try:
gcc -o calcphase calcphase.c -lfftw3f -lm -arch x86_64

Sam Johnston
110519
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fftw3.h>
#include <time.h>

float PI = M_PI;
float PI2 = 2*M_PI;

int main(int argc, char *argv[])
{

	unsigned long i, j, k, ii, jj, kk, ind, ind1, ind2, nx, ny, nz, nx2, ny2, nz2, nxy, nxyz, offset, baseline_spacing, baseline_skip, pval_type, pval_size, time;
	long x, y, z; 
	int intval;
	short shortval;
	float fermir, fermiw, val, ang, scale, te;
	float *phase;
	char *filenamein, *filenameout;
	FILE *fp;
	fftwf_complex *img, *img_low;
	fftwf_plan fft3d_reconstruct, fft3d_forward, fft3d_backward;
	//clock_t ticks;
	
	// read in arguments
	filenamein = argv[1]; // name of pfile 
	filenameout = argv[2]; // name of file to write phase array
	offset = atoi(argv[3]); // number of bytes to skip in beginning of pfile
	baseline_spacing = atoi(argv[4]); // number of elements between each baseline
	pval_type = atoi(argv[5]); // data type in pfile: (0) short or (1) long
	nx = atoi(argv[6]); // number of elements in first dimension of array (columns)
	ny = atoi(argv[7]); // number of elements in second dimension of array (rows)
	nz = atoi(argv[8]); // number of elements in third dimension of array (slices)
	fermir = atof(argv[9]); // radius of Fermi window
	fermiw = atof(argv[10]); // width of Fermi window
	te = atof(argv[11]); // echo time
	//time = atoi(argv[12]); // should we report time? (0) no or (1) yes

	nx2 = nx/2;
	ny2 = ny/2;
	nz2 = nz/2;
	
	nxy = nx*ny;
	nxyz = nxy*nz;
	
	scale = 1.0/(te*PI2);
	
	if (pval_type) {
		pval_size = sizeof(int);
	} else {
		pval_size = sizeof(short);
	}
	
	// allocate memory
	phase = (float *) malloc(sizeof(float)*nxyz);
	img = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nxyz);
	img_low = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex)*nxyz);

	// Fourier transform plans
	fft3d_reconstruct = fftwf_plan_dft_3d(nz, ny, nx, img, img, FFTW_BACKWARD, FFTW_ESTIMATE);
	fft3d_forward = fftwf_plan_dft_3d(nz, ny, nx, img_low, img_low, FFTW_FORWARD, FFTW_ESTIMATE);
	fft3d_backward = fftwf_plan_dft_3d(nz, ny, nx, img_low, img_low, FFTW_BACKWARD, FFTW_ESTIMATE);
	
	// open pfile
	fp = fopen(filenamein,"r");
	
	//if (time) printf("Time before read: %d\n",clock()/CLOCKS_PER_SEC);

	// read in raw data from file
	baseline_skip = 2*nx*pval_size;
	fseek(fp,offset,SEEK_CUR);
	fseek(fp,baseline_skip,SEEK_CUR);
	for (i=0; i<nz; i++) { // for each slice
		for (j=0; j<ny; j++) { // for each row
			for (k=0; k<nx; k++) { // for each column
				ind = k + j*nx + i*nxy;
				if (pval_type) {
					fread(&intval,pval_size,1,fp); // imaginary
					img[ind][1] = 1.0*intval;
					fread(&intval,pval_size,1,fp); // real
					img[ind][0] = 1.0*intval;
				} else {
					fread(&shortval,pval_size,1,fp); // imaginary
					img[ind][1] = 1.0*shortval;
					fread(&shortval,pval_size,1,fp); // real
					img[ind][0] = 1.0*shortval;	
				}	
			}
			if ((((i*ny)+j+1)%baseline_spacing)==0) {
				fseek(fp,baseline_skip,SEEK_CUR);
			}
		}
	}	
	
	// close pfile
	fclose(fp);
	
	//if (time) printf("Time before shift 1: %d\n",clock()/CLOCKS_PER_SEC);
	
	// fft shift
	for (i=0; i<nz2; i++) { // for each slice (doing half of all slices because contents are swapped with other half)
		ii = (i+nz2)%nz;
		for (j=0; j<ny; j++) { // for each row
			jj = (j+ny2)%ny;
			for (k=0; k<nx; k++) { // for each column
				kk = (k+nx2)%nx;
				ind1 = k + j*nx + i*nxy;
				ind2 = kk + jj*nx + ii*nxy;
				// swap real
				val = img[ind1][0];
				img[ind1][0] = img[ind2][0];
				img[ind2][0] = val;
				// swap imaginary
				val = img[ind1][1];
				img[ind1][1] = img[ind2][1];
				img[ind2][1] = val;
			}
		}
	}
	
	//if (time) printf("Time before ifft: %d\n",clock()/CLOCKS_PER_SEC);
	
	fftwf_execute(fft3d_reconstruct); // out of k space
	
	//if (time) printf("Time before shift 2: %d\n",clock()/CLOCKS_PER_SEC);

	// fft shift
	for (i=0; i<nz2; i++) { // for each slice (doing half of all slices because contents are swapped with other half)
		ii = (i+nz2)%nz;
		for (j=0; j<ny; j++) { // for each row
			jj = (j+ny2)%ny;
			for (k=0; k<nx; k++) { // for each column
				kk = (k+nx2)%nx;
				ind1 = k + j*nx + i*nxy;
				ind2 = kk + jj*nx + ii*nxy;
				// swap real
				val = img[ind1][0];
				img[ind1][0] = img[ind2][0];
				img[ind2][0] = val;
				// swap imaginary
				val = img[ind1][1];
				img[ind1][1] = img[ind2][1];
				img[ind2][1] = val;
			}
		}
	}
	
	
	//
	// Construct a low resolution image by filtering with a Fermi window in the Fourier domain.
	//
	
	//if (time) printf("Time before copy: %s\n",clock()/CLOCKS_PER_SEC);
	
	// copy data to low res array
	for (i=0; i<nxyz; i++) {
		img_low[i][0] = img[i][0];
		img_low[i][1] = img[i][1];
	}
	
	//if (time) printf("Time before fft: %s\n",clock()/CLOCKS_PER_SEC);

	fftwf_execute(fft3d_forward); // into Fourier domain
	
	//if (time) printf("Time before filter: %s\n",clock()/CLOCKS_PER_SEC);
	
	// multiply by Fermi window
	// using modulo indices because of fft shift
	// scaling filter by number of elements because fftw doesn't scale
	for (i=0; i<nz; i++) { // for each slice
		z = ((i+nz2)%nz)-nz2;
		for (j=0; j<ny; j++) { // for each row
			y = ((j+ny2)%ny)-ny2;
			for (k=0; k<nx; k++) { // for each column
				x = ((k+nx2)%nx)-nx2;
				ind = k + j*nx + i*nxy;
				val = 1.0/(nxyz*(1.0 + exp((sqrt((1.0*x*x) + (1.0*y*y) + (1.0*z*z))-fermir)/fermiw)));
				img_low[ind][0] *= val; // real
				img_low[ind][1] *= val; // imaginary
			}
		}
	}
	
	//if (time) printf("Time before ifft: %s\n",clock()/CLOCKS_PER_SEC);
	
	fftwf_execute(fft3d_backward); // out of Fourier domain
	
	//
	// Compute the phase angle at each element, subtract the corresponding angle from 
	// the low resolution image.
	//
	
	//if (time) printf("Time before phase: %s\n",clock()/CLOCKS_PER_SEC);
	
	// compute phase
	for (i=0; i<nxyz; i++) { // for each element	
		ang = atan2(img[i][1],img[i][0]) - atan2(img_low[i][1],img_low[i][0]);
		// force between -pi and pi
		ang = ang>PI ? ang-PI2 : ang; 
		ang = ang<-PI ? ang+PI2 : ang;
		phase[i] = ang*scale;
	}
	
	//if (time) printf("Time before write: %s\n",clock()/CLOCKS_PER_SEC);
	
	// write phase file
	fp = fopen(filenameout,"w");
	fwrite(phase,sizeof(float),nxyz,fp);
	fclose(fp);
	
	// deallocate memory
	fftwf_destroy_plan(fft3d_forward);
	fftwf_destroy_plan(fft3d_backward);
	fftwf_free(img);
	fftwf_free(img_low);
	free(phase);

	return 0;
}

/* junk

	FILE *fp2;
	fp2 = fopen("low2.bin","w");
	for (i=0; i<nxyz; i++) {
		fwrite(&img_low[i][0],sizeof(float),1,fp2); 
	}
	for (i=0; i<nxyz; i++) {
		fwrite(&img_low[i][1],sizeof(float),1,fp2); 
	}
	fclose(fp2);
	printf("%f %f\n",img[0][0],img[0][1]);
	printf("%f %f\n",img_low[0][0],img_low[0][1]);
	printf("%f %f %f\n",phase[3094],atan2(img[3094][1],img[3094][0]),atan2(img_low[3094][1],img_low[3094][0]));
	
	*/
	

