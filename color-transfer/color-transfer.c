/**
 * @file color-transfert
 * @brief transfert color from source image to target image.
 *        Method from Reinhard et al. : 
 *        Erik Reinhard, Michael Ashikhmin, Bruce Gooch and Peter Shirley, 
 *        'Color Transfer between Images', IEEE CGA special issue on 
 *        Applied Perception, Vol 21, No 5, pp 34-41, September - October 2001
 */

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <bcl.h>

#define D 3

float RGB2LMS[D][D] = {
  {0.3811, 0.5783, 0.0402}, 
  {0.1967, 0.7244, 0.0782},  
  {0.0241, 0.1288, 0.8444}
};


float LMS2RGB[D][D] = {
  {4.4679,-3.5873,0.1193},
  {-1.2186,2.3809,-0.1624},
  {0.0497,-0.2439,1.2045}
};

void LMS_from_RGB(float R,float G,float B,float* LMS){
  LMS[0] = log(R*RGB2LMS[0][0] + G*RGB2LMS[0][1] + B*RGB2LMS[0][2]);
  LMS[1] = log(R*RGB2LMS[1][0] + G*RGB2LMS[1][1] + B*RGB2LMS[1][2]);
  LMS[2] = log(R*RGB2LMS[2][0] + G*RGB2LMS[2][1] + B*RGB2LMS[2][2]);
}

//matrixes multiplication might be wrong, recheck that later if it doesn't work
void Lalphabeta_from_LMS(float L,float M,float S,float* lalphabeta){
  lalphabeta[0] = L*0.57735026919	+ M * 0.40824829046 + S * 0.70710678118;
  lalphabeta[1] = L*0.57735026919 + M * 0.40824829046 + S * -1.41421356236;
  lalphabeta[2] = L*0.57735026919 + M * -0.40824829046	+ S * 0;
}

void LMS_from_lalphabeta(float l,float alpha,float beta,float* LMS){
  LMS[0] = l * 0.57735026919 + alpha * 0.40824829046 + beta * 0.70710678118;
  LMS[1] = l * 0.57735026919 + alpha * 0.40824829046 + beta * -0.70710678118;
  LMS[2] = l * 0.57735026919 + alpha * -0.81649658092	+ beta*0;
}

void RGB_from_LMS(float L,float M,float S, float* RGB){
  RGB[0] = log(L*LMS2RGB[0][0] + M*LMS2RGB[0][1] + S*LMS2RGB[0][2]);
  RGB[1] = log(L*LMS2RGB[1][0] + M*LMS2RGB[1][1] + S*LMS2RGB[1][2]);
  RGB[2] = log(L*LMS2RGB[2][0] + M*LMS2RGB[2][1] + S*LMS2RGB[2][2]);
}

float* lalphabeta_means(float** lalphabetas, int pixels,float* lalphabeta_means_result){
  for(int i;i<pixels;i++)
    for(int y=0;y<3;y++)
      lalphabeta_means_result[y]+=lalphabetas[i][y]/pixels;
}

void lalphabeta_star(float* l,float* alpha, float* beta,float* means){
  *l-=-means[0];
  *alpha-=means[1];
  *beta-=means[2];
}

void lalphabeta_prime(float* l_star,float* alpha_star, float* beta_star,float* means)




void
process(char *ims, char *imt, char* imd){
  (void) ims;
  (void) imt;
  (void) imd;

}

void
usage (char *s){
  fprintf(stderr, "Usage: %s <ims> <imt> <imd> \n", s);
  exit(EXIT_FAILURE);
}

#define PARAM 3
int
main(int argc, char *argv[]){
  if (argc != PARAM+1) usage(argv[0]);
  process(argv[1], argv[2], argv[3]);
  return EXIT_SUCCESS;
}
