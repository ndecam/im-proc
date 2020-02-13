#include <stdlib.h>
#include <stdio.h>
#include <bcl.h>
#include "fft.h"

void 
padding(int rows,int cols,int factor,fftw_complex* freq_repr,fftw_complex* freq_repr_with_padding){
    printf("Image is of size %d*%d\n",rows,cols);
    int topleftcorner_row = ((rows*factor)-rows)/2;
    int topleftcorner_col = ((cols*factor)-cols)/2;
    printf("Bigger image is of size %d*%d\n",rows*factor,cols*factor);
    printf("Top-left corner is (%d,%d)\n",topleftcorner_row,topleftcorner_col);
    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++){
            printf("%d %d modifies index %d\n",i,j, (i+topleftcorner_row) * (cols*factor) + topleftcorner_col + j);
            printf("It is supposed to be %d*%d\n",(i+topleftcorner_row),topleftcorner_col + j);
            freq_repr_with_padding[ (i+topleftcorner_row) * (cols*factor) + topleftcorner_col + j ] = freq_repr[(i*cols)+j];
        }
    
    /*for(int i=0;i<cols*factor;i++){
        printf("[");
        for(int j=0;j<rows*factor;j++)
            printf("%f",creal(freq_repr_with_padding[(i*rows)+j]));
        printf("]\n");
    }*/
}

void process(int factor, char* ims, char* imd){
    pnm imageims = pnm_load(ims);
    int cols = pnm_get_width(imageims);
    int rows = pnm_get_height(imageims);
    pnm imageimd = pnm_new(cols*factor,rows*factor,PnmRawPpm);

    unsigned short *data = (unsigned short*) malloc(rows*cols*sizeof(unsigned short));
    unsigned short *returned_data = (unsigned short*) malloc(factor*rows*cols*factor*sizeof(unsigned short));
    
    for(int channel = 0;channel <3;channel++){
        fftw_complex* padding_complex = (fftw_complex*)calloc(rows*cols*factor*factor,sizeof(fftw_complex));

        pnm_get_channel(imageims,data,channel);
        fftw_complex* img_complex = forward(rows,cols,data);

        padding(rows,cols,factor,img_complex,padding_complex);
        printf("Padding done\n");

        returned_data = backward(rows*factor,cols*factor,padding_complex);
        /*for(int i =0;i<rows*cols;i++)
            printf("%d\n",returned_data[i]);*/
        pnm_set_channel(imageimd,returned_data,channel);
        free(padding_complex);
    }
    free(data);
    free(returned_data);
    pnm_save(imageimd,PnmRawPpm,imd);


            
}



void usage(){
    printf("padding <factor> <source> <target>");
    exit(EXIT_FAILURE);
}





#define PARAM 3
int
main(int argc, char *argv[])
{
  if (argc != PARAM+1) usage(argv[0]);
  int factor =  atoi(argv[1]);
  char* ims = argv[2];
  char* imd = argv[3];
  process(factor, ims, imd );


  return EXIT_SUCCESS;
}
