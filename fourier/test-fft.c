/**
 * @file test-fft.c
 * @brief test the behaviors of functions in fft module
 *
 */
#include <stdlib.h>
#include <stdio.h>

#include <fft.h>

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

  unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
  
  pnm_get_channel(ims,data,0);
  printf("%d %d %d\n",data[0],data[1],data[2]);
  
  
  fftw_complex* img_complex = forward(rows,cols,data);
  
  data = backward(rows,cols,img_complex);
  char* newname = (char*) malloc(50* sizeof(char));
  
  pnm_set_channel(ims,data,0);
  pnm_set_channel(ims,data,1);
  pnm_set_channel(ims,data,2);
  printf("%d %d %d\n",data[0],data[1],data[2]);
  sprintf(newname, "FB-%s.ppm",name);
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
  unsigned short *data = pnm_get_image(ims);

  int cols = pnm_get_width(ims);
  int rows = pnm_get_height(ims);
  fftw_complex* img_complex = forward(rows,cols,data);

  float as[rows*cols],ps[rows*cols];

  freq2spectra(rows,cols,img_complex,as,ps);
  spectra2freq(rows,cols,as,ps,img_complex);

  data = backward(rows,cols,img_complex);
  char* newname = (char*) malloc(50* sizeof(char));
  sprintf(newname, "FB-%s.ppm",name);
  pnm_save(ims,PnmRawPpm,newname);


  (void)name;
  fprintf(stderr, "OK\n");
  free(newname);

}

/**
 * @brief test construction of magnitude and phase images in ppm files
 * @param char* name, the input image file name
 */
void
test_display(char* name)
{
  fprintf(stderr, "test_display: ");
  (void)name;
  fprintf(stderr, "OK\n");
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
  //run(argv[1]);
  test_forward_backward(argv[1]);
  //test_reconstruction(argv[1]);
  return EXIT_SUCCESS;
}
