/* Header file for lib/array_creators_alloc */

#ifndef LIB_ARRAY_CREATORS_ALLOC_H
#define LIB_ARRAY_CREATORS_ALLOC_H

#include "matisse.h"
#include <stdlib.h>
#include "tensor_struct.h"

/**
 */
tensor_d* zeros_d2(int dim_1, int dim_2, tensor_d** restrict t);

/**
 */
tensor_d* zeros_d3(int dim_1, int dim_2, int dim_3, tensor_d** restrict t);

/**
 */
tensor_c* zeros_c2(int dim_1, int dim_2, tensor_c** restrict t);

#endif
