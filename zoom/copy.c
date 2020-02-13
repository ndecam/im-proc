#include <stdlib.h>
#include <stdio.h>
#include <bcl.h>

void process(int factor, char* ims, char* imd){
    pnm imageims = pnm_load(ims);
    int cols = pnm_get_width(imageims);
    int rows = pnm_get_height(imageims);

    pnm imageimd = pnm_new(cols*factor, rows*factor, PnmRawPpm);

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            for(int fi =0;fi<factor;fi++)
                for(int fi2=0;fi2<factor;fi2++)
                    for(int channel = 0; channel<3; channel++){
                        pnm_set_component(imageimd,(i*factor)+fi,(j*factor)+fi2,channel,pnm_get_component(imageims,i,j,channel));
                    }

    pnm_save(imageimd,PnmRawPpm,imd);
            
}




void usage(){
    printf("copy <factor> <source> <target>");
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
