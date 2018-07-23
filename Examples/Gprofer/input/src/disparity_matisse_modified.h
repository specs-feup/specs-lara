/* Header file for disparity_matisse_modified */

#ifndef DISPARITY_MATISSE_MODIFIED_H
#define DISPARITY_MATISSE_MODIFIED_H

#include "lib/matisse.h"
#include "lib/tensor_struct.h"
#include <stdlib.h>

/**
 */
tensor_d* disparity_matisse_modified_tctdtd_row_2d_2d_1(tensor_c* type, tensor_d* imleft, tensor_d* imright, tensor_d** restrict imDispOwn);

#endif
