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

double RGB2LMS[D][D] = {
  {0.3811, 0.5783, 0.0402}, 
  {0.1967, 0.7244, 0.0782},  
  {0.0241, 0.1288, 0.8444}
};


double LMS2RGB[D][D] = {
  {4.4679,-3.5873,0.1193},
  {-1.2186,2.3809,-0.1624},
  {0.0497,-0.2439,1.2045}
};

void LMS_from_RGB(double R,double G,double B,double* LMS){
  if(R ==0.0 && G==0.0 && B==0.0){
    R=0.0001;
    G=0.0001;
    B=0.0001;
  }
  LMS[0] = log10(R*RGB2LMS[0][0] + G*RGB2LMS[0][1] + B*RGB2LMS[0][2]);
  LMS[1] = log10(R*RGB2LMS[1][0] + G*RGB2LMS[1][1] + B*RGB2LMS[1][2]);
  LMS[2] = log10(R*RGB2LMS[2][0] + G*RGB2LMS[2][1] + B*RGB2LMS[2][2]);
}

void Lalphabeta_from_LMS(double L,double M,double S,double* lalphabeta){
  lalphabeta[0] = L*0.57735026919	+ M *0.57735026919 + S *0.57735026919;
  lalphabeta[1] = L*0.40824829046 + M *0.40824829046 + S *-0.81649658092;
  lalphabeta[2] = L*0.70710678118 + M *-0.70710678118	+ S *0;
}

void LMS_from_lalphabeta(double l,double alpha,double beta,double* LMS){
  LMS[0] = pow(10,l * 0.57735026919 + alpha * 0.40824829046 + beta * 0.70710678118);
  LMS[1] = pow(10,l * 0.57735026919 + alpha * 0.40824829046 + beta * -0.70710678118);
  LMS[2] = pow(10,l * 0.57735026919 + alpha * -0.81649658092	+ beta*0);
}

void RGB_from_LMS(double L,double M,double S, unsigned short* RGB){
  RGB[0] = L*LMS2RGB[0][0] + M*LMS2RGB[0][1] + S*LMS2RGB[0][2];
  RGB[1] = L*LMS2RGB[1][0] + M*LMS2RGB[1][1] + S*LMS2RGB[1][2];
  RGB[2] = L*LMS2RGB[2][0] + M*LMS2RGB[2][1] + S*LMS2RGB[2][2];
}



//void lalphabeta_prime(double* l_star,double* alpha_star, double* beta_star,double* means)
void image_to_lalphabeta_data_points(pnm image,double*** lalphabeta_data_points){
  printf("Entered!\n");
  int cols = pnm_get_width(image);
  int rows = pnm_get_height(image);
  printf("Rows:%d,cols:%d\n",rows,cols);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      unsigned short R = pnm_get_component(image,i,j,0);
      unsigned short G = pnm_get_component(image,i,j,1);
      unsigned short B = pnm_get_component(image,i,j,2);
      double LMS[3];
      LMS_from_RGB(R,G,B,LMS);
      double Lalphabeta[3];
      Lalphabeta_from_LMS(LMS[0],LMS[1],LMS[2],Lalphabeta);
      lalphabeta_data_points[i][j][0] = Lalphabeta[0];
      lalphabeta_data_points[i][j][1] = Lalphabeta[1];
      lalphabeta_data_points[i][j][2] = Lalphabeta[2];
    }
}

void lalphabeta_data_points_to_image(pnm imageimd,double*** lalphabeta_data_points){
  int cols = pnm_get_width(imageimd);
  int rows = pnm_get_height(imageimd);
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++){
      double LMS[3];
      LMS_from_lalphabeta(lalphabeta_data_points[i][j][0],lalphabeta_data_points[i][j][1],lalphabeta_data_points[i][j][2],LMS);
      unsigned short RGB[3];
      RGB_from_LMS(LMS[0],LMS[1],LMS[2],RGB);
      pnm_set_component(imageimd,i,j,0,RGB[0]);
      pnm_set_component(imageimd,i,j,1,RGB[1]);
      pnm_set_component(imageimd,i,j,2,RGB[2]);
    }
}

void lalphabeta_image_mean_standard_deviation(int rows,int cols,double*** lalphabeta_data_points,double* mean,double* standard_deviation){
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int channel = 0;channel<3;channel++){
        mean[channel] += lalphabeta_data_points[i][j][channel];
      }
  for(int channel = 0;channel<3;channel++)
    mean[channel]/=(rows*cols);
  
  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int channel = 0;channel<3;channel++)
        standard_deviation[channel] += pow(lalphabeta_data_points[i][j][channel] - mean[channel], 2);
  
  for(int channel = 0;channel<3;channel++){
    standard_deviation[channel]/=(double)(rows*cols);
    standard_deviation[channel] = sqrt(standard_deviation[channel]);
  }
    
}
//applies colors from the target image in lalphabeta form on source.
void apply_color_transfer(double*** lalphabeta_data_points_source,int rows,int cols,double*** lalphabeta_data_points_target,int rows_imt,int cols_imt){
  double sd_source[3] = {0.0,0.0,0.0};
  double mean_source[3] =  {0.0,0.0,0.0};
  double sd_target[3] =  {0.0,0.0,0.0};
  double mean_target[3] =  {0.0,0.0,0.0};

  lalphabeta_image_mean_standard_deviation(rows,cols,lalphabeta_data_points_source,mean_source,sd_source);
  lalphabeta_image_mean_standard_deviation(rows_imt,cols_imt,lalphabeta_data_points_target,mean_target,sd_target);

  for(int i=0;i<rows;i++)
    for(int j=0;j<cols;j++)
      for(int channel=0;channel<3;channel++){
        lalphabeta_data_points_source[i][j][channel]-=mean_source[channel];
        lalphabeta_data_points_source[i][j][channel]*=(sd_target[channel]/sd_source[channel]);
        lalphabeta_data_points_source[i][j][channel]+=mean_target[channel];
      }

}

double*** malloc_3dim_array(int x,int y,int z){
  double *** array = (double ***)malloc(x*sizeof(double**));

  for (int i = 0; i< x; i++) {
    array[i] = (double **) malloc(y*sizeof(double *));
    for (int j = 0; j < y; j++) {
        array[i][j] = (double *)malloc(z*sizeof(double));
    }
  }
  return array;
}


void
process(char *ims, char *imt, char* imd){
  pnm imageims = pnm_load(ims);
  int cols = pnm_get_width(imageims);
  int rows = pnm_get_height(imageims);
  pnm imageimt = pnm_load(imt);
  int cols_imt = pnm_get_width(imageimt);
  int rows_imt = pnm_get_height(imageimt);
  pnm imageimd = pnm_new(cols_imt, rows_imt, PnmRawPpm);
  
  
  
  printf("Rows:%d,cols:%d\n",rows,cols);
  double*** lalphabeta_source=malloc_3dim_array(rows,cols,3);
  double*** lalphabeta_target=malloc_3dim_array(rows_imt,cols_imt,3);
  
  image_to_lalphabeta_data_points(imageims,lalphabeta_source);
  image_to_lalphabeta_data_points(imageimt,lalphabeta_target);

  
  apply_color_transfer(lalphabeta_target,rows_imt,cols_imt,lalphabeta_source,rows,cols);

  lalphabeta_data_points_to_image(imageimd,lalphabeta_target);
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
