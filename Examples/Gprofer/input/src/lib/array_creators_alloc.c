/* Implementation file for lib/array_creators_alloc */

#include "array_creators_alloc.h"
#include "matisse.h"
#include <stdlib.h>
#include "tensor.h"
#include "tensor_struct.h"


/**
 */
tensor_d* zeros_d2(int dim_1, int dim_2, tensor_d** restrict t)
{

   
   return new_const_array_d2(dim_1, dim_2, 0.0, t);
}

/**
 */
tensor_d* zeros_d3(int dim_1, int dim_2, int dim_3, tensor_d** restrict t)
{

   
   return new_const_array_d3(dim_1, dim_2, dim_3, 0.0, t);
}

/**
 */
tensor_c* zeros_c2(int dim_1, int dim_2, tensor_c** restrict t)
{

   
   return new_const_array_c2(dim_1, dim_2, 0, t);
}
