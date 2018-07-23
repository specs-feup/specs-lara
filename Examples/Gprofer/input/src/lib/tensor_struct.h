/* Header file for lib/tensor_struct */

#ifndef LIB_TENSOR_STRUCT_H
#define LIB_TENSOR_STRUCT_H

/**
 * @struct tensor
 *
 * Represents a tensor.
 * Has information about the shape of the tensor and the length of its linearized version.
 *
 */
typedef struct tensor_struct_d {

	double* data;
	int length;

	int* shape;
	int dims;
	int owns_data;

} tensor_d;

/**
 * @struct tensor
 *
 * Represents a tensor.
 * Has information about the shape of the tensor and the length of its linearized version.
 *
 */
typedef struct tensor_struct_c {

	char* data;
	int length;

	int* shape;
	int dims;
	int owns_data;

} tensor_c;

/**
 * @struct tensor
 *
 * Represents a tensor.
 * Has information about the shape of the tensor and the length of its linearized version.
 *
 */
typedef struct tensor_struct_i {

	int* data;
	int length;

	int* shape;
	int dims;
	int owns_data;

} tensor_i;

#endif
