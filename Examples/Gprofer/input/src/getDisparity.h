/* Header file for getDisparity */

#ifndef GETDISPARITY_H
#define GETDISPARITY_H

#include "lib/matisse.h"
#include "lib/tensor_struct.h"
#include <stdlib.h>

/**
 */
tensor_d* getDisparity_tdtdii_2d_2d_1(tensor_d* Ileft, tensor_d* Iright, int win_sz, int max_shift, tensor_d** restrict retDisparity);

/**
 */
tensor_d* correlateSAD_tdtdii_undef_undef_1(tensor_d* Ileft, tensor_d* Iright, int win_sz, int disparity, tensor_d** restrict retSAD);

/**
 */
tensor_d* integralImage2D_td_2d_1(tensor_d* I, tensor_d** restrict retImg);

#endif
