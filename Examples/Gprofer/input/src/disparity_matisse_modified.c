/* Implementation file for disparity_matisse_modified */

#include "disparity_matisse_modified.h"
#include "getDisparity.h"
#include "lib/matisse.h"
#include "lib/matisse_string.h"
#include "lib/tensor.h"
#include "lib/tensor_struct.h"
#include <stdlib.h>


/**
 */
tensor_d* disparity_matisse_modified_tctdtd_row_2d_2d_1(tensor_c* type, tensor_d* imleft, tensor_d* imright, tensor_d** restrict imDispOwn)
{
   int WIN_SZ;
   int SHIFT;
   tensor_d* getDisparity_arg1 = NULL;
   int numel_result;
   int iter_1;
   tensor_d* getDisparity_arg2 = NULL;
   int numel_result_1;
   int iter;

   // path(path,common);
   WIN_SZ = 8;
   SHIFT = 64;
   if(strcmp_partially_checked_1tc(type, "test")){
      WIN_SZ = 2;
      SHIFT = 1;
   }else if(strcmp_partially_checked_1tc(type, "sim_fast")){
      WIN_SZ = 4;
      SHIFT = 4;
   }else if(strcmp_partially_checked_1tc(type, "sim")){
      WIN_SZ = 4;
      SHIFT = 8;
   }
   
   // outFile = [resultDir, '/', 'out', '.bmp'];
   // file = [dataDir, '/1.bmp'];
   // imleft = readImage(file);
   // imright = readImage([dataDir, '/2.bmp']);
   // fprintf(1,'Input size\t\t- (%dx%d)\n', rows, cols);
   // start = photonStartTiming;
   new_array_helper_d(imleft->shape, imleft->dims, &getDisparity_arg1);
   numel_result = imleft->length;
   for(iter_1 = 0; iter_1 < numel_result; ++iter_1){
      getDisparity_arg1->data[iter_1] = (double) imleft->data[iter_1];
   }
   
   new_array_helper_d(imright->shape, imright->dims, &getDisparity_arg2);
   numel_result_1 = imright->length;
   for(iter = 0; iter < numel_result_1; ++iter){
      getDisparity_arg2->data[iter] = (double) imright->data[iter];
   }
   
   getDisparity_tdtdii_2d_2d_1(getDisparity_arg1, getDisparity_arg2, WIN_SZ, SHIFT, imDispOwn);
   // stop = photonEndTiming;
   // elapsed = photonReportTiming(start, stop);
   // writeMatrix(imDispOwn, dataDir);
   // imwrite(uint8(minSAD), outFile, 'bmp');
   // photonPrintTiming(elapsed);
   tensor_free_d(&getDisparity_arg1);
   tensor_free_d(&getDisparity_arg2);
   
   return *imDispOwn;
}
