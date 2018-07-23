/* Header file for lib/general_matrix */

#ifndef LIB_GENERAL_MATRIX_H
#define LIB_GENERAL_MATRIX_H

#include "matisse.h"
#include <stdlib.h>
#include "tensor_struct.h"

/**
 */
void size_multiargs_td_3_generic(tensor_d* in, int* restrict out_1, int* restrict out_2, int* restrict out_3);

/**
 */
void size_multiargs_td_2_of_2(tensor_d* in, int* restrict out_1, int* restrict out_2);

/**
 */
void size_multiargs_td_2_generic(tensor_d* in, int* restrict out_1, int* restrict out_2);

#endif
