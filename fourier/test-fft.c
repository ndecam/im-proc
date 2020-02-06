/**
 * @file test-fft.c
 * @brief test the behaviors of functions in fft module
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include <fft.h>
#include <math.h>
#define M_PI       3.14159265358979323846
/**
 * @brief test the forward and backward functions
 * @param char* name, the input image file name
 */


void
test_forward_backward(char* name)
{
  fprintf(stderr, "test_forward_backward:\n");
  pnm ims = pnm_load(name);
  
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);
  for(int channel = 0;channel <3;channel++){
    unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
    pnm_get_channel(ims,data,channel);
    
    fftw_complex* img_complex = forward(rows,cols,data);
    data = backward(rows,cols,img_complex);
    
    pnm_set_channel(ims,data,channel);
    free(data);
  }
  char* newname = (char*) malloc(50* sizeof(char));
  sprintf(newname, "FB-%s",name);
  pnm_save(ims,PnmRawPpm,newname);


  (void)name;
  fprintf(stderr, "OK\n");
  free(newname);
}

/**
 * @brief test image reconstruction from of magnitude and phase spectrum
 * @param char *name: the input image file name
 */
void
test_reconstruction(char* name)
{
  fprintf(stderr, "test_reconstruction: ");
  pnm ims = pnm_load(name);

  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);
  for(int channel = 0;channel <3;channel++){
    unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
    pnm_get_channel(ims,data,channel);

    fftw_complex* img_complex = forward(rows,cols,data);

    float as[rows*cols],ps[rows*cols];

    freq2spectra(rows,cols,img_complex,as,ps);
    spectra2freq(rows,cols,as,ps,img_complex);

    data = backward(rows,cols,img_complex);
    pnm_set_channel(ims,data,channel);
    free(data);
  }
  char* newname = (char*) malloc(50* sizeof(char));
  
  sprintf(newname, "FB-ASPM-%s",name);
  pnm_save(ims,PnmRawPpm,newname);


  (void)name;
  fprintf(stderr, "OK\n");
  free(newname);

}

void
re_center(int rows, int cols,unsigned short* data){
  unsigned short output[rows*cols];
  int middle = (rows*cols)/2;
  for(int i=0;i<rows*cols;i++){
    output[(i+middle)%(rows*cols)] = data[i];
  }

  for(int i=0;i<rows*cols;i++){
    data[i]=output[(i-1+cols/2)%(rows*cols)];
  }
}

/**
 * @brief test construction of magnitude and phase images in ppm files
 * @param char* name, the input image file name
 */
void
test_display(char* name)
{
  fprintf(stderr, "test_display: ");
  pnm ims = pnm_load(name);
  pnm ims2 = pnm_load(name);
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);
  float k = 0.2;
  for(int channel = 0;channel <3;channel++){
    unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
    pnm_get_channel(ims,data,channel);

    fftw_complex* img_complex = forward(rows,cols,data);
    float as[rows*cols],ps[rows*cols];
    unsigned short norm_as[rows*cols], norm_ps[rows*cols];
    freq2spectra(rows,cols,img_complex,as,ps);
    float amax = 0;
    for(int i=0;i<rows*cols;i++)
      if(as[i] > amax) amax = as[i];
    for(int i=0;i<rows*cols;i++){
      as[i] = pow( as[i] / amax, k) * 255;
      norm_as[i] = as[i];
      norm_ps[i] = ps[i];
    }

    re_center(rows,cols,norm_as);
    re_center(rows,cols,norm_ps);

    pnm_set_channel(ims,norm_as,channel);
    pnm_set_channel(ims2,norm_ps,channel);
    free(data);
  }
  char* newname = (char*) malloc(50* sizeof(char));
  
  sprintf(newname, "AS-%s",name);
  pnm_save(ims,PnmRawPpm,newname);

  sprintf(newname, "PS-%s",name);
  pnm_save(ims2,PnmRawPpm,newname);

  (void)name;
  fprintf(stderr, "OK\n");
  free(newname);
}
/**
 * @brief test the modification of magnitude by adding a periodic functions
          on both vertical and horizontal axis, and
 *        construct output images
 * @param char* name, the input image file name
 */
void

test_add_frequencies(char* name)
{
  fprintf(stderr, "test_add_frequencies: ");
  pnm ims = pnm_load(name);
  pnm ims2 = pnm_load(name);
  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);
  float k = 0.2;
  for(int channel = 0;channel <3;channel++){
    unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
    pnm_get_channel(ims,data,channel);

    fftw_complex* img_complex = forward(rows,cols,data);

    float as[rows*cols],ps[rows*cols];
    unsigned short norm_as[rows*cols];

    freq2spectra(rows,cols,img_complex,as,ps);

    float amax = 0;
    for(int i=0;i<rows*cols;i++)
      if(as[i] > amax) amax = as[i];
    
    //i est vertical
    //j est horizontal

    for(int i=0;i<rows;i++)
      for(int j=0;j<cols;j++){
        float sinx = 0.25 * amax * sinf(2*M_PI*i);
        float siny = 0.25 * amax * sinf(2*M_PI*j);
        as[(i*cols)+j] = as[(i*cols)+j] + sinx + siny;
        norm_as[(i*cols)+j] = pow( as[i] / amax, k) * 255;
      }
    
    re_center(rows,cols,norm_as);
    pnm_set_channel(ims2,norm_as,channel);
    

    spectra2freq(rows,cols,as,ps,img_complex);
    data = backward(rows,cols,img_complex);

    pnm_set_channel(ims,data,channel);
    free(data);
  }

  char* newname = (char*) malloc(50* sizeof(char));
  sprintf(newname, "FREQ-%s",name);
  pnm_save(ims,PnmRawPpm,newname);

  sprintf(newname, "FAS-%s",name);
  pnm_save(ims2,PnmRawPpm,newname);

  (void)name;
  fprintf(stderr, "OK\n");
}

void
run(char* name)
{
  test_forward_backward(name);
  test_reconstruction(name);
  test_display(name);
  test_add_frequencies(name);
}



void
usage(const char *s)
{
  fprintf(stderr, "Usage: %s <ims> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 1
int
main(int argc, char *argv[])
{
  if (argc != PARAM+1) usage(argv[0]);
  run(argv[1]);
  /*test_forward_backward(argv[1]);
  test_reconstruction(argv[1]);
  test_display(argv[1]);*/
  return EXIT_SUCCESS;
}
