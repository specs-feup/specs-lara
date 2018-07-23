/* Header file for lib/matrix */

#ifndef LIB_MATRIX_H
#define LIB_MATRIX_H

#include "matisse.h"
#include <stdlib.h>
#include "tensor_struct.h"

/**
 */
tensor_d* copy_td_ptd(tensor_d* source_matrix, tensor_d** restrict target_matrix);

#endif
