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
  LMS[0] = log10(R*RGB2LMS[0][0] + G*RGB2LMS[0][1] + B*RGB2LMS[0][2]);
  LMS[1] = log10(R*RGB2LMS[1][0] + G*RGB2LMS[1][1] + B*RGB2LMS[1][2]);
  LMS[2] = log10(R*RGB2LMS[2][0] + G*RGB2LMS[2][1] + B*RGB2LMS[2][2]);
}

//matrixes multiplication might be wrong, recheck that later if it doesn't work
void Lalphabeta_from_LMS(float L,float M,float S,float* lalphabeta){
  lalphabeta[0] = L*0.57735026919	+ M *0.57735026919 + S *0.57735026919;
  lalphabeta[1] = L*0.40824829046 + M *0.40824829046 + S *-0.81649658092;
  lalphabeta[2] = L*0.70710678118 + M *-0.70710678118	+ S *0;
}

void LMS_from_lalphabeta(float l,float alpha,float beta,float* LMS){
  LMS[0] = pow(10,l * 0.57735026919 + alpha * 0.40824829046 + beta * 0.70710678118);
  LMS[1] = pow(10,l * 0.57735026919 + alpha * 0.40824829046 + beta * -0.70710678118);
  LMS[2] = pow(10,l * 0.57735026919 + alpha * -0.81649658092	+ beta*0);
}

void RGB_from_LMS(float L,float M,float S, unsigned short* RGB){
  RGB[0] = L*LMS2RGB[0][0] + M*LMS2RGB[0][1] + S*LMS2RGB[0][2];
  RGB[1] = L*LMS2RGB[1][0] + M*LMS2RGB[1][1] + S*LMS2RGB[1][2];
  RGB[2] = L*LMS2RGB[2][0] + M*LMS2RGB[2][1] + S*LMS2RGB[2][2];
}

void lalphabeta_means(float** lalphabetas, int pixels,float* lalphabeta_means_result){
  for(int i;i<pixels;i++)
    for(int y=0;y<3;y++)
      lalphabeta_means_result[y]+=lalphabetas[i][y]/pixels;
}

void lalphabeta_star(float* l,float* alpha, float* beta,float* means){
  *l-=-means[0];
  *alpha-=means[1];
  *beta-=means[2];
}

void lalphabeta_prime(float lstar,float alphastar,float betastar,float* sd_source,float* sd_target,float * lalphabetaprime){
  lalphabetaprime[0] = (sd_target[0]/sd_source[0])*lstar;
  lalphabetaprime[1] = (sd_target[1]/sd_source[1])*alphastar;
  lalphabetaprime[2] = (sd_target[2]/sd_source[2])*betastar;
}

//void lalphabeta_prime(float* l_star,float* alpha_star, float* beta_star,float* means)
void image_to_lalphabeta_data_points(pnm image,float lalphabeta_data_points[pnm_get_width(image)][pnm_get_height(image)][3]){
  int cols = pnm_get_width(image);
  int rows = pnm_get_height(image);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      unsigned short R = pnm_get_component(image,i,j,0);
      unsigned short G = pnm_get_component(image,i,j,1);
      unsigned short B = pnm_get_component(image,i,j,2);
      float LMS[3];
      LMS_from_RGB(R,G,B,LMS);
      float Lalphabeta[3];
      Lalphabeta_from_LMS(LMS[0],LMS[1],LMS[2],Lalphabeta);
      lalphabeta_data_points[i][j][0] = Lalphabeta[0];
      lalphabeta_data_points[i][j][1] = Lalphabeta[1];
      lalphabeta_data_points[i][j][2] = Lalphabeta[2];
    }
}

void lalphabeta_data_points_to_image(pnm imageimd,float lalphabeta_data_points[pnm_get_width(imageimd)][pnm_get_height(imageimd)][3]){
  int cols = pnm_get_width(imageimd);
  int rows = pnm_get_height(imageimd);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      float LMS[3];
      LMS_from_lalphabeta(lalphabeta_data_points[i][j][0],lalphabeta_data_points[i][j][1],lalphabeta_data_points[i][j][2],LMS);
      unsigned short RGB[3];
      RGB_from_LMS(LMS[0],LMS[1],LMS[2],RGB);
      pnm_set_component(imageimd,i,j,0,RGB[0]);
      pnm_set_component(imageimd,i,j,1,RGB[1]);
      pnm_set_component(imageimd,i,j,2,RGB[2]);
    }
}

void lalphabeta_image_mean_standard_deviation(int rows,int cols,float lalphabeta_data_points[rows][cols][3],float* mean,float* standard_deviation){
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int channel = 0;channel<3;channel++)
        mean[channel] += lalphabeta_data_points[i][j][channel]/(rows*cols);
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int channel = 0;channel<3;channel++)
        standard_deviation[channel] += (lalphabeta_data_points[i][j][channel]-mean[channel])/(rows*cols);
    
}



void
process(char *ims, char *imt, char* imd){
  pnm imageims = pnm_load(ims);
  int cols = pnm_get_width(imageims);
  int rows = pnm_get_height(imageims);
  pnm imageimd = pnm_new(cols, rows, PnmRawPpm);

  pnm imageimt = pnm_load(imt);
  int cols_imt = pnm_get_width(imageimt);
  int rows_imt = pnm_get_height(imageimt);
  

  float lalphabeta_source[rows][cols][3];
  float lalphabeta_target[rows_imt][cols_imt][3];
  image_to_lalphabeta_data_points(imageims,lalphabeta_source);
  image_to_lalphabeta_data_points(imageimt,lalphabeta_target);

  float sd_source[3];
  float mean_source[3];
  float sd_target[3];
  float mean_target[3];

  lalphabeta_image_mean_standard_deviation(rows,cols,lalphabeta_source,mean_source,sd_source);
  lalphabeta_image_mean_standard_deviation(rows_imt,cols_imt,lalphabeta_target,mean_target,sd_target);

  lalphabeta_data_points_to_image(imageimd,lalphabeta_source);
  pnm_save(imageimd,PnmRawPpm,imd);
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
