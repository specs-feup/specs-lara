/* Implementation file for getDisparity */

#include "getDisparity.h"
#include "lib/array_creators_alloc.h"
#include "lib/general_matrix.h"
#include "lib/matisse.h"
#include "lib/matrix.h"
#include "lib/tensor.h"
#include "lib/tensor_struct.h"
#include "padarray_specialized_pre_post_tdd1x2_2.h"
#include "padarray_specialized_pre_tdi1x2s_2.h"
#include <stdlib.h>


/**
 */
tensor_d* getDisparity_tdtdii_2d_2d_1(tensor_d* Ileft, tensor_d* Iright, int win_sz, int max_shift, tensor_d** restrict retDisparity)
{
   tensor_d* Ileft_1 = NULL;
   int numel_result;
   int iter_1;
   tensor_d* Iright_1 = NULL;
   int numel_result_1;
   int iter;
   int nr;
   int nc;
   int nb;
   double half_win_sz;
   double matrix_1[2];
   tensor_d* IleftPadded = NULL;
   double matrix[2];
   tensor_d* IrightPadded = NULL;
   tensor_d* minSAD = NULL;
   int mtimes;
   int minSAD_numel;
   int iter_3;
   int retDisparity_numel;
   int iter_2;
   int k;
   tensor_d* retSAD = NULL;
   int j;
   int i;
   double a;
   double b;

   new_array_helper_d(Ileft->shape, Ileft->dims, &Ileft_1);
   numel_result = Ileft->length;
   for(iter_1 = 0; iter_1 < numel_result; ++iter_1){
      Ileft_1->data[iter_1] = (double) Ileft->data[iter_1];
   }
   
   new_array_helper_d(Iright->shape, Iright->dims, &Iright_1);
   numel_result_1 = Iright->length;
   for(iter = 0; iter < numel_result_1; ++iter){
      Iright_1->data[iter] = (double) Iright->data[iter];
   }
   
   size_multiargs_td_3_generic(Ileft_1, &nr, &nc, &nb);
   if(win_sz > 1){
      half_win_sz = (double) win_sz / 2.0;
      // Removed no-op function call: matisse_new_array_from_dims
      matrix_1[0] = half_win_sz;
      matrix_1[1] = half_win_sz;
      padarray_specialized_pre_post_tdd1x2_2_tdd1x2_2d_1(Ileft_1, matrix_1, &IleftPadded);
      // Removed no-op function call: matisse_new_array_from_dims
      matrix[0] = half_win_sz;
      matrix[1] = half_win_sz;
      padarray_specialized_pre_post_tdd1x2_2_tdd1x2_2d_1(Iright_1, matrix, &IrightPadded);
   }else{
      if(Ileft_1 != NULL){
         new_array_helper_d(Ileft_1->shape, Ileft_1->dims, &IleftPadded);
         copy_td_ptd(Ileft_1, &IleftPadded);
      }
      
      if(Iright_1 != NULL){
         new_array_helper_d(Iright_1->shape, Iright_1->dims, &IrightPadded);
         copy_td_ptd(Iright_1, &IrightPadded);
      }
      
   }
   
   new_array_d_2(nr, nc, &minSAD);
   mtimes = 255 * 255 * 255;
   minSAD_numel = minSAD->length;
   for(iter_3 = 0; iter_3 < minSAD_numel; ++iter_3){
      minSAD->data[iter_3] = (double) mtimes;
   }
   
   new_array_d_2(nr, nc, retDisparity);
   retDisparity_numel = (*retDisparity)->length;
   for(iter_2 = 0; iter_2 < retDisparity_numel; ++iter_2){
      (*retDisparity)->data[iter_2] = (double) max_shift;
   }
   
   for(k = 1; k <= max_shift; ++k){
      correlateSAD_tdtdii_undef_undef_1(IleftPadded, IrightPadded, win_sz, k - 1, &retSAD);
      //  Find Disparity
      //  findDisparity(retSAD, minSAD, retDisp, k, nr, nc)
      for(j = 0; j < nc; ++j){
         for(i = 0; i < nr; ++i){
            a = retSAD->data[(i) + (j) * retSAD->shape[0]];
            b = minSAD->data[(i) + (j) * minSAD->shape[0]];
            if(a < b){
               minSAD->data[(i) + (j) * minSAD->shape[0]] = a;
               (*retDisparity)->data[(i) + (j) * (*retDisparity)->shape[0]] = (double) k;
            }
            
         }
         
      }
      
   }
   
   tensor_free_d(&IleftPadded);
   tensor_free_d(&Ileft_1);
   tensor_free_d(&IrightPadded);
   tensor_free_d(&Iright_1);
   tensor_free_d(&minSAD);
   tensor_free_d(&retSAD);
   
   return *retDisparity;
}

/**
 */
tensor_d* correlateSAD_tdtdii_undef_undef_1(tensor_d* Ileft, tensor_d* Iright, int win_sz, int disparity, tensor_d** restrict retSAD)
{
   int matrix[2];
   tensor_d* Iright_moved = NULL;
   int rows;
   int cols;
   tensor_d* SAD = NULL;
   int end;
   int i;
   double diff;
   tensor_d* integralImg = NULL;
   int colon_arg1_6;
   int colon_arg1_7;
   int range_size;
   int range_size_1;
   tensor_d* plus_arg1 = NULL;
   int colon_arg2_1;
   tensor_i* size_1 = NULL;
   tensor_d* minus_arg1 = NULL;
   int colon_arg1_3;
   tensor_i* size_2 = NULL;
   tensor_d* minus_arg1_6 = NULL;
   int colon_arg1;
   tensor_i* size = NULL;
   int iter_1;
   int integralImg_index_4;
   int integralImg_index_1;
   int integralImg_index_3;
   int integralImg_index_6;
   int iter;

   // Removed no-op function call: matisse_new_array_from_dims
   matrix[0] = 0;
   matrix[1] = disparity;
   padarray_specialized_pre_tdi1x2s_2_tdi1x2_undef_1(Iright, matrix, &Iright_moved);
   size_multiargs_td_2_generic(Ileft, &rows, &cols);
   zeros_d2(rows, cols, &SAD);
   end = Ileft->length;
   for(i = 0; i < end; ++i){
      diff = Ileft->data[i] - Iright_moved->data[i];
      SAD->data[i] = diff * diff;
   }
   
   // 2D scan.
   integralImage2D_td_2d_1(SAD, &integralImg);
   colon_arg1_6 = win_sz + 1;
   colon_arg1_7 = win_sz + 1;
   range_size = integralImg->shape[0] - colon_arg1_6 + 1;
   range_size_1 = integralImg->shape[1] - colon_arg1_7 + 1;
   new_array_d_3(range_size, range_size_1, 1, &plus_arg1);
   colon_arg2_1 = integralImg->shape[1] - win_sz + 1;
   if(!(integralImg->shape[0] - win_sz + 1 <= integralImg->shape[0])){
      abort();
   }
   
   if(!(colon_arg2_1 <= integralImg->shape[1])){
      abort();
   }
   
   size_d(plus_arg1, &size_1);
   new_array_ti_d(size_1, &minus_arg1);
   colon_arg1_3 = win_sz + 1;
   if(!(integralImg->shape[0] - win_sz + 1 <= integralImg->shape[0])){
      abort();
   }
   
   size_d(minus_arg1, &size_2);
   new_array_ti_d(size_2, &minus_arg1_6);
   colon_arg1 = win_sz + 1;
   if(!(integralImg->shape[1] - win_sz + 1 <= integralImg->shape[1])){
      abort();
   }
   
   size_d(minus_arg1_6, &size);
   new_array_ti_d(size, retSAD);
   for(iter_1 = 1; iter_1 <= range_size_1; ++iter_1){
      integralImg_index_4 = iter_1 + colon_arg1_7 - 1;
      integralImg_index_1 = iter_1 + 2 - 1;
      integralImg_index_3 = iter_1 + colon_arg1_3 - 1;
      integralImg_index_6 = iter_1 + 2 - 1;
      for(iter = 1; iter <= range_size; ++iter){
         (*retSAD)->data[(iter - 1) + (iter_1 - 1) * (*retSAD)->shape[0]] = integralImg->data[(iter + colon_arg1_6 - 1 - 1) + (integralImg_index_4 - 1) * integralImg->shape[0]] + integralImg->data[(iter + 2 - 1 - 1) + (integralImg_index_1 - 1) * integralImg->shape[0]] - integralImg->data[(iter + 2 - 1 - 1) + (integralImg_index_3 - 1) * integralImg->shape[0]] - integralImg->data[(iter + colon_arg1 - 1 - 1) + (integralImg_index_6 - 1) * integralImg->shape[0]];
      }
      
   }
   
   tensor_free_d(&Iright_moved);
   tensor_free_d(&SAD);
   tensor_free_d(&integralImg);
   tensor_free_d(&minus_arg1);
   tensor_free_d(&minus_arg1_6);
   tensor_free_d(&plus_arg1);
   tensor_free_i(&size);
   tensor_free_i(&size_1);
   tensor_free_i(&size_2);
   
   return *retSAD;
}

/**
 */
tensor_d* integralImage2D_td_2d_1(tensor_d* I, tensor_d** restrict retImg)
{
   int nr;
   int nc;
   int nb;
   int I_end;
   int iter_2;
   int i;
   int retImg_arg1;
   int retImg_end;
   int iter_1;
   int j;
   int retImg_arg2;
   int retImg_end_1;
   int iter;

   size_multiargs_td_3_generic(I, &nr, &nc, &nb);
   zeros_d3(nr, nc, nb, retImg);
   if(!(1 <= I->shape[0])){
      abort();
   }
   
   I_end = I->shape[1];
   for(iter_2 = 0; iter_2 < I_end; ++iter_2){
      (*retImg)->data[iter_2] = I->data[iter_2];
   }
   
   for(i = 2; i <= nr; ++i){
      retImg_arg1 = i - 1;
      retImg_end = (*retImg)->shape[1];
      for(iter_1 = 0; iter_1 < retImg_end; ++iter_1){
         (*retImg)->data[(i - 1) + (iter_1) * (*retImg)->shape[0]] = (*retImg)->data[(retImg_arg1 - 1) + (iter_1) * (*retImg)->shape[0]] + I->data[(i - 1) + (iter_1) * I->shape[0]];
      }
      
   }
   
   // vtuneResumeMex;
   for(j = 2; j <= nc; ++j){
      retImg_arg2 = j - 1;
      retImg_end_1 = (*retImg)->shape[0];
      for(iter = 0; iter < retImg_end_1; ++iter){
         (*retImg)->data[(iter) + (j - 1) * (*retImg)->shape[0]] = (*retImg)->data[(iter) + (retImg_arg2 - 1) * (*retImg)->shape[0]] + (*retImg)->data[(iter) + (j - 1) * (*retImg)->shape[0]];
      }
      
   }
   
   // vtunePauseMex;
   
   return *retImg;
}
