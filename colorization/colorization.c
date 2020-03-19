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

int** malloc_2dim_array(int x,int y){
  int ** array = (int **)malloc(x*sizeof(int*));

  for (int i = 0; i< x; i++) {
    array[i] = (int *) malloc(y*sizeof(int *));
  }
  return array;
}

void neighborhood_preprocessing(double*** lalphabeta_data_points, double*** luminance_mean_sd, int rows, int cols){
    int patch_size = 5;
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            for(int patched_i=i-(patch_size/2);patched_i<i+(patch_size/2+1);patched_i++)
                for(int patched_j=j-(patch_size/2);patched_j<j+(patch_size/2+1);patched_j++)
                    if(patched_i>0 && patched_i<rows && patched_j>0 && patched_j<cols)
                        luminance_mean_sd[i][j][0]+=(lalphabeta_data_points[patched_i][patched_j][0]/(patch_size*patch_size));
                

    
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++){
            for(int patched_i=i-(patch_size/2);patched_i<i+(patch_size/2+1);patched_i++)
                for(int patched_j=j-(patch_size/2);patched_j<j+(patch_size/2+1);patched_j++)
                    if(patched_i>0 && patched_i<rows && patched_j>0 && patched_j<cols)
                        luminance_mean_sd[i][j][1] += pow(lalphabeta_data_points[patched_i][patched_j][0] - luminance_mean_sd[i][j][0], 2);
        luminance_mean_sd[i][j][1]/=(double)(patch_size*patch_size);
        luminance_mean_sd[i][j][1] = sqrt(luminance_mean_sd[i][j][1]);
        }
    
}
void select_best_swatch(int** random_swatches,double*** neighborhood_source, double* neighborhood_target,int* swatch_row,int* swatch_col){
    
    double best_result = (neighborhood_target[0]-neighborhood_source[0][0][0])+(neighborhood_target[1]-neighborhood_source[0][0][1]);
    for(int swatch=0;swatch<200;swatch++){
        int swatch_row_t = random_swatches[swatch][0];
        int swatch_col_t = random_swatches[swatch][1];
        int neighborhood_comparison = (neighborhood_target[0]-neighborhood_source[swatch_row_t][swatch_col_t][0])+(neighborhood_target[1]-neighborhood_source[swatch_row_t][swatch_col_t][1]);
        if(   neighborhood_comparison < best_result){
            *swatch_row=swatch_row_t;
            *swatch_col=swatch_col_t;
            best_result = neighborhood_comparison;
        }
    }
}

void 
apply_colorization(double*** lalphabeta_source,double*** neighborhood_source,int rows,int cols, double*** lalphabeta_target,double*** neighborhood_target,int rows_imt,int cols_imt){
    int** random_swatches = malloc_2dim_array(200,2) ;
    for(int swatch = 0;swatch<200;swatch++){
        int row = rand()%rows;
        int col = rand()%cols;
        random_swatches[swatch][0] = row;
        random_swatches[swatch][1] = col;
        printf("%d %d\n",row,col);
    }

    for(int i=0;i<rows_imt;i++)
        for(int j=0;j<cols_imt;j++){
            int swatch_row;
            int swatch_col;
            select_best_swatch(random_swatches,neighborhood_source,neighborhood_target[i][j],&swatch_row,&swatch_col);
            for(int channel=1;channel<3;channel++){
                lalphabeta_target[i][j][channel]=lalphabeta_source[swatch_row][swatch_col][channel];
            }
        }
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
  
  
  
  double*** lalphabeta_source=malloc_3dim_array(rows,cols,3);
  double*** lalphabeta_target=malloc_3dim_array(rows_imt,cols_imt,3);
  double*** neighborhood_source = malloc_3dim_array(rows,cols,2);
  double*** neighborhood_target = malloc_3dim_array(rows_imt,cols_imt,2);

  image_to_lalphabeta_data_points(imageims,lalphabeta_source);
  image_to_lalphabeta_data_points(imageimt,lalphabeta_target);
  printf("Preprocessing neighborhood\n");
  neighborhood_preprocessing(lalphabeta_source,neighborhood_source,rows,cols);
  neighborhood_preprocessing(lalphabeta_target,neighborhood_target,rows_imt,cols_imt);
  printf("Applying colorization\n");
  apply_colorization(lalphabeta_source,neighborhood_source,rows,cols, lalphabeta_target,neighborhood_target,rows_imt,cols_imt);
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