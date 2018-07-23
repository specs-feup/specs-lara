/* Implementation file for padarray_specialized_pre_tdi1x2s_2 */

#include "lib/general_matrix.h"
#include "lib/matisse.h"
#include "lib/matlab_general.h"
#include "lib/tensor.h"
#include "lib/tensor_struct.h"
#include "new_from_matrix_double_0_d.h"
#include "padarray_specialized_pre_tdi1x2s_2.h"
#include <stdlib.h>


/**
 */
tensor_d* padarray_specialized_pre_tdi1x2s_2_tdi1x2_undef_1(tensor_d* A, int pad[2], tensor_d** restrict y)
{
   int A_ndims;
   int pad_numel;
   int y_ndims;
   tensor_d* y_size = NULL;
   int i;
   int size_at_i;
   int y_numel;
   int iter;
   int s1;
   int s2;
   int i2;
   int i1;

   //  TODO: Validate pad is vector.
   A_ndims = ndims_alloc_d(A);
   pad_numel = 2;
   y_ndims = max_scalars_dec_ii(A_ndims, pad_numel);
   new_array_d_2(1, y_ndims, &y_size);
   for(i = 1; i <= y_ndims; ++i){
      size_at_i = A->dims < i ? 1 : A->shape[i - 1];
      if(i <= pad_numel){
         size_at_i += (pad[i - 1]) * (1 + 0);
      }
      
      y_size->data[i - 1] = (double) size_at_i;
   }
   
   new_from_matrix_double_0_d_td_row_1(y_size, y);
   y_numel = (*y)->length;
   for(iter = 0; iter < y_numel; ++iter){
      (*y)->data[iter] = 0.0;
   }
   
   size_multiargs_td_2_generic(A, &s1, &s2);
   for(i2 = 1; i2 <= s2; ++i2){
      for(i1 = 1; i1 <= s1; ++i1){
         (*y)->data[(i1 + (pad[0]) - 1) + (i2 + (pad[1]) - 1) * (*y)->shape[0]] = A->data[(i1 - 1) + (i2 - 1) * A->shape[0]];
      }
      
   }
   
   tensor_free_d(&y_size);
   
   return *y;
}
