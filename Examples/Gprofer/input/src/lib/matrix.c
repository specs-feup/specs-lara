/* Implementation file for lib/matrix */

#include "matisse.h"
#include "matrix.h"
#include <stdlib.h>
#include "tensor_struct.h"


/**
 */
tensor_d* copy_td_ptd(tensor_d* source_matrix, tensor_d** restrict target_matrix)
{
   int i;

   for(i = 0; i < source_matrix->length; ++i){
      (*target_matrix)->data[i] = source_matrix->data[i];
   }
   
   
   return *target_matrix;
}
