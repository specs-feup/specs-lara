/* Implementation file for lib/general_matrix */

#include "general_matrix.h"
#include "matisse.h"
#include <stdlib.h>
#include "tensor.h"
#include "tensor_struct.h"


/**
 */
void size_multiargs_td_3_generic(tensor_d* in, int* restrict out_1, int* restrict out_2, int* restrict out_3)
{
   int result;
   int i;

   *out_1 = in->shape[0];
   *out_2 = in->shape[1];
   result = 1;
   for(i = 2; i < (ndims_alloc_d(in)); ++i){
      result *= in->shape[i];
   }
   
   *out_3 = result;
}

/**
 */
void size_multiargs_td_2_of_2(tensor_d* in, int* restrict out_1, int* restrict out_2)
{

   *out_1 = in->shape[0];
   *out_2 = in->shape[1];
}

/**
 */
void size_multiargs_td_2_generic(tensor_d* in, int* restrict out_1, int* restrict out_2)
{
   int result;
   int i;

   *out_1 = in->shape[0];
   result = 1;
   for(i = 1; i < (ndims_alloc_d(in)); ++i){
      result *= in->shape[i];
   }
   
   *out_2 = result;
}
