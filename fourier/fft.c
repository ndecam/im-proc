#include <float.h>
#include <stdlib.h>
#include <math.h>

#include <fft.h>

fftw_complex *
forward(int rows, int cols, unsigned short* g_img)
{
  /*(void)rows;
  (void)cols;
  (void)g_img;
  return NULL;*/

  fftw_complex* img_complex = (fftw_complex*) malloc(rows*cols*sizeof(fftw_complex));
  for(int cmp=0; cmp< rows*cols; cmp++ ){
    img_complex[cmp] = (double)g_img[cmp];
  }

  fftw_complex* struc_complex = (fftw_complex*) malloc(rows*cols*sizeof(fftw_complex));

  fftw_plan plan = fftw_plan_dft_2d(rows,cols,img_complex,struc_complex,FFTW_FORWARD,FFTW_ESTIMATE );
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  free(img_complex);
  return struc_complex;
}

unsigned short *
backward(int rows, int cols, fftw_complex* freq_repr)
{
  /*
  unsigned short* g_img = (unsigned short) malloc(rows*cols*sizeoff(unsigned short));
  (void)rows;
  (void)cols;
  (void)freq_repr;
  return NULL;*/
  fftw_complex* struc_complex = (fftw_complex*) malloc(rows*cols*sizeof(fftw_complex));

  fftw_plan plan = fftw_plan_dft_2d(rows,cols,freq_repr,struc_complex,FFTW_BACKWARD,FFTW_ESTIMATE );
  fftw_execute(plan);
  fftw_destroy_plan(plan);

  unsigned short* g_img = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));

  for(int cmp=0; cmp< rows*cols; cmp++ ){
    //printf("%f\n",creal(struc_complex[cmp]) / 65353);
    g_img[cmp] = creal(struc_complex[cmp]) / (65535*4) ;
    //printf("%d\n", g_img[cmp] );
  }

  free(struc_complex);

  return g_img;
}

void
freq2spectra(int rows, int cols, fftw_complex* freq_repr, float* as, float* ps)
{
  /*(void)rows;
  (void)cols;
  (void)freq_repr;
  (void)as;
  (void)ps;*/
  for(int i=0;i<rows*cols;i++){
    as[i] = cabs(freq_repr[i]);
    ps[i] = carg(freq_repr[i]);
  }
}

void
spectra2freq(int rows, int cols, float* as, float* ps, fftw_complex* freq_repr)
{
  /*(void)rows;
  (void)cols;
  (void)as;
  (void)ps;
  (void)freq_repr;*/
  for(int i = 0;i<rows*cols;i++)
    freq_repr[i] = as[i] + ps[i]*_Complex_I;
}

