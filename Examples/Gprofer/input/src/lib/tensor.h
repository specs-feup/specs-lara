/* Header file for lib/tensor */

#ifndef LIB_TENSOR_H
#define LIB_TENSOR_H

#include "matisse.h"
#include <stdlib.h>
#include "tensor_struct.h"

/**
 */
tensor_d* new_array_helper_d(int* shape, int dims, tensor_d** restrict t);

/**
 */
int is_same_shape_alloc_d(tensor_d* t, int* restrict shape, int dims);

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_d(tensor_d** restrict t);

/**
 */
tensor_d* new_array_d_2(int dim_1, int dim_2, tensor_d** restrict t);

/**
 */
int ndims_alloc_d(tensor_d* t);

/**
 */
tensor_d* new_array_td_d(tensor_d* shape, tensor_d** restrict t);

/**
 */
tensor_d* new_array_ti_d(tensor_i* shape, tensor_d** restrict t);

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_i(tensor_i** restrict t);

/**
 */
tensor_i* size_d(tensor_d* t, tensor_i** restrict size_array);

/**
 */
tensor_i* new_array_i_2(int dim_1, int dim_2, tensor_i** restrict t);

/**
 */
tensor_i* new_array_helper_i(int* shape, int dims, tensor_i** restrict t);

/**
 */
int is_same_shape_alloc_i(tensor_i* t, int* restrict shape, int dims);

/**
 */
tensor_d* new_array_d_3(int dim_1, int dim_2, int dim_3, tensor_d** restrict t);

/**
 */
tensor_d* new_const_array_d2(int dim_1, int dim_2, double value, tensor_d** restrict t);

/**
 */
void set_matrix_values_alloc_d(tensor_d* t, double value);

/**
 */
tensor_d* new_const_array_d3(int dim_1, int dim_2, int dim_3, double value, tensor_d** restrict t);

/**
 */
tensor_c* new_const_array_c2(int dim_1, int dim_2, char value, tensor_c** restrict t);

/**
 */
void set_matrix_values_alloc_c(tensor_c* t, char value);

/**
 */
tensor_c* new_array_c_2(int dim_1, int dim_2, tensor_c** restrict t);

/**
 */
tensor_c* new_array_helper_c(int* shape, int dims, tensor_c** restrict t);

/**
 */
int is_same_shape_alloc_c(tensor_c* t, int* restrict shape, int dims);

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_c(tensor_c** restrict t);

#endif
