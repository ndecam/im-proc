#include <stdlib.h>
#include <stdio.h>

#include <bcl.h>

void
usage (char *s)
{
  fprintf(stderr,"Usage: %s %s", s, "<i> <j> <rows> <cols> <ims> <imd>\n");
  exit(EXIT_FAILURE);
}



void
process(const size_t rows, const size_t cols, const size_t i, const size_t j, char* ims, char* imd)
{
  pnm image_imd = pnm_new(cols, rows, PnmRawPpm);
  pnm image_ims = pnm_load(ims);

  //int ims_rows = pnm_get_height(image_ims);
  //int ims_cols = pnm_get_width(image_ims);

  for(size_t cmpcol = 0; cmpcol < cols; cmpcol ++ ){
    for(size_t cmprow = 0; cmprow < rows; cmprow ++ ){
       short redValue = pnm_get_component(image_ims,i+ cmpcol,j + cmprow, PnmRed);
       short greenValue = pnm_get_component(image_ims,i+ cmpcol,j + cmprow, PnmGreen);
       short blueValue = pnm_get_component(image_ims,i+ cmpcol,j + cmprow, PnmBlue);

       pnm_set_component(image_imd,cmpcol, cmprow, PnmRed, redValue);
       pnm_set_component(image_imd,cmpcol, cmprow, PnmGreen, greenValue);
       pnm_set_component(image_imd,cmpcol, cmprow, PnmBlue, blueValue);
    }
  }
  pnm_save(image_imd, PnmRawPpm, imd);

  pnm_free(image_imd);
  pnm_free(image_ims);

}



#define PARAM 6
int
main(int argc, char *argv[])
{
  if (argc != PARAM+1) usage(argv[0]);
  size_t i =  atoi(argv[1]);
  size_t j =  atoi(argv[2]);
  size_t rows = atoi(argv[3]);
  size_t cols = atoi(argv[4]);
  char* ims = argv[5];
  char* imd = argv[6];
  process(rows, cols ,i , j, ims, imd );


  return EXIT_SUCCESS;
}
