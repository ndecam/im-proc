#include <stdlib.h>
#include <stdio.h>
#include <bcl.h>
#include<string.h>
#include <math.h>



double box(double pixel){
  if(pixel > (-1/2) && pixel<(1/2)) return 1;
  return 0;
}

double tent(double x){
  if(x>-1 && x<1) return 1-abs(x);
  return 0;
}

double bell(double x){
  if(x<=0.5) return (x*x)+0.75;
  if(abs(x)>0.5 && abs(x)<1.5) return 0.5*(abs(x) -1.5)*(abs(x) -1.5);
  return 0;
}

double mitchell_netravali(double x){
  if(x>-1 && x<1) return (7/6)*pow(abs(x),3) - 2*x*x + (8/9);
  if( (x>-2 && x<-1) || (x>1 && x<2) ) return (-(7/18))*pow(abs(x),3)+2*x*x - (10/3)*abs(x) + (16/9);
  return 0.f;
}

double h(double x,char* filter){
  printf("Filter is %s|\n",filter);
  if(strcmp(filter, "box")==0){
    printf("Box returned %f\n",box(x));
    return box(x);
  }
  if(strcmp(filter, "tent")==0){
    return tent(x);
  }
  if(strcmp(filter, "bell")==0){
    return bell(x);
  }
  if(strcmp(filter, "mitch")==0){
    return mitchell_netravali(x);
  }
  return 0;
}
void process(int factor, char* ims, char* imd,char* filter){
    pnm imageims = pnm_load(ims);
    int cols = pnm_get_width(imageims);
    int rows = pnm_get_height(imageims);
    int domain;
    if(strcmp(filter, "box")==0){
      domain=1;
    }
    if(strcmp(filter, "tent")==0){
      domain = 2;
    }
    if(strcmp(filter, "bell")==0){
      domain = 3;
    }
    if(strcmp(filter, "mitch")==0){
      domain = 4;
    }
    pnm imageimd = pnm_new(cols*factor, rows, PnmRawPpm);
    printf("Doing rotation on cols = %d , rows = %d\n",cols,rows);
    for(int channel=0;channel <3;channel ++){
      for(int i =0;i<rows;i++)
        for(int j_prime=0;j_prime<cols*factor;j_prime++){
          printf("J_prime = %d\n",j_prime);
          double j = j_prime/factor;
          float left = j-(domain/2);
          float right = j+(domain/2);
          if(left<0) left = 0;
          if(right>cols) right = cols;
          unsigned short S;
          printf("Left is %f,right is %f\n",left,right);
          for(int k=left;k<right;k++){
            S = S + (pnm_get_component(imageims,i,k,channel)*h(k-j,filter));
          }
          printf("New component in %d\n",S);
          pnm_set_component(imageimd,i,j_prime,channel,S);
        }
    }
    pnm_save(imageimd,PnmRawPpm,imd);

}

void rotate(char* ims, char* imd){
  pnm imageims = pnm_load(ims);
  int cols = pnm_get_width(imageims);
  int rows = pnm_get_height(imageims);
  pnm imageimd = pnm_new(rows, cols, PnmRawPpm);
  for(int channel = 0; channel <3; channel++){
    for(int i=0; i< rows; i++){
      for(int j=0; j< cols; j++ ){
        int newi = cols - i-1;
        int newj = i;
        unsigned short value = pnm_get_component(imageims, i, j, channel);
        pnm_set_component(imageimd,newi,newj,channel,value);


      }
    }
  }

  pnm_save(imageimd,PnmRawPpm,imd);
}




void usage(){
    printf("filter <factor> <filter-name> <ims> <imd>");
    exit(EXIT_FAILURE);
}

#define PARAM 4
int
main(int argc, char *argv[])
{
  if (argc != PARAM+1) usage(argv[0]);
  int factor =  atoi(argv[1]);
  char* filter = argv[2];
  char* ims = argv[3];
  char* imd = argv[4];

  process(factor,ims,"intermediare.ppm",filter);
  printf("first process done\n");

  rotate("intermediare.ppm","rotate.ppm");
  printf("rotate done\n");

  process(factor,"rotate.ppm",imd,filter);
  printf("second process done\n");

  rotate(imd,imd);
  printf("rotate done\n");
  rotate(imd,imd);
  printf("rotate done\n");
  rotate(imd,imd);
  printf("rotate done\n");

  return EXIT_SUCCESS;
}
