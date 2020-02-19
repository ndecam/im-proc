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
  if(strcmp(filter, "box")){
    return box(x);
  }
  if(strcmp(filter, "tent")){
    return tent(x);
  }
  if(strcmp(filter, "bell")){
    return bell(x);
  }
  if(strcmp(filter, "mitch")){
    return mitchell_netravali(x);
  }
  return 0;
}
void process(int factor, char* ims, char* imd,char* filter){
    pnm imageims = pnm_load(ims);
    int cols = pnm_get_width(imageims);
    int rows = pnm_get_height(imageims);
    int domain;
    if(strcmp(filter, "box")){
      domain=1;
    }
    if(strcmp(filter, "tent")){
      domain = 2;
    }
    if(strcmp(filter, "bell")){
      domain = 3;
    }
    if(strcmp(filter, "mitch")){
      domain = 4;
    }
    pnm imageimd = pnm_new(cols*factor, rows*factor, PnmRawPpm);

    //TO DO:
    for(int channel;channel <3;channel ++){
      for(int i =0;i<rows;i++)
        for(int j_prime=0;j_prime<cols;j_prime++){
          double j = j_prime/factor;
          int left = j-(domain/2);
          int right = j+(domain/2);
          unsigned short S;
          for(int k=left;k<right;k++){
            S = S + pnm_get_component(imageims,i,k,channel)*h(k-j,filter);
          }
          pnm_set_component(imageimd,i,j_prime,channel,S);
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
  process(factor,ims,imd,filter);

  return EXIT_SUCCESS;
}
