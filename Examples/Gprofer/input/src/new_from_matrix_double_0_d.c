/* Implementation file for new_from_matrix_double_0_d */

#include "lib/matisse.h"
#include "lib/tensor.h"
#include "lib/tensor_struct.h"
#include "new_from_matrix_double_0_d.h"
#include <stdlib.h>


/**
 */
tensor_d* new_from_matrix_double_0_d_td_row_1(tensor_d* shape, tensor_d** restrict new_matrix)
{
   int end;
   int i;

   // function new_matrix = new_from_matrix_double_<VALUE_STRING>(shape)
   //  Create matrix
   new_array_td_d(shape, new_matrix);
   //  Initialize matrix
   end = (*new_matrix)->length;
   for(i = 0; i < end; ++i){
      (*new_matrix)->data[i] = 0.0;
   }
   
   
   return *new_matrix;
}
