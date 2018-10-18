/* Header file for lib/load */

#ifndef LIB_LOAD_H
#define LIB_LOAD_H

#include <inttypes.h>
#include "matisse.h"
#include <stdlib.h>
#include "tensor_struct.h"

/**
 */
tensor_d* load_matrix_variable_2_td(const char * varname, int dim_0, int dim_1, tensor_d** restrict out);

/**
 */
void load_raw_data_from_file_td(const char * varname, uint32_t num_elements, tensor_d* out);

/**
 */
char * get_absolute_filename(char * filename);

#ifdef _WIN32
#include <windows.h>
#endif


#endif
