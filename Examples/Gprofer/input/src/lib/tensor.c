/* Implementation file for lib/tensor */

#include "matisse.h"
#include <stdio.h>
#include <stdlib.h>
#include "tensor.h"
#include "tensor_struct.h"


/**
 */
tensor_d* new_array_helper_d(int* shape, int dims, tensor_d** restrict t)
{
/* new_row(int* shape, int dims, tensor** t) */

	int TRIES = 2;
	int i;
	int length;
	int sameShape;
	int* previous_shape = NULL;
	double* previous_data = NULL;
	
	/* Check if matrix is already allocated */
	if(*t != NULL) {
	
		/* If shape is the same, return current matrix, */
		/* even if it does not own the data (so it can implement window-writing) */
		sameShape = is_same_shape_alloc_d(*t, shape, dims);	
		//if(sameShape && (*t)->owns_data) {
		if(sameShape) {
			return *t;
		} 		
		
		/* Save pointers to previous shape and data */
		/* Only free data if tensor owns it */
		if((*t)->owns_data) {
			previous_data = (*t)->data;
		}
		
		previous_shape = (*t)->shape;		
	}

	/* Create the tensor and return it. */
	free(*t);
	
	for(i=0; i<TRIES; i++) {
		*t = (tensor_d*) malloc(sizeof(tensor_d));
		if (*t != NULL) {
			break;
		}
	}
	/**t = (tensor_d*) malloc(sizeof(tensor_d)); */

	if (*t == NULL) {
       printf("ERROR: Could not allocate memory for the matrix structure\n");
	   exit(EXIT_FAILURE);
	}

	/* Calculate the length of the linearized version of the tensor. */
	length = 1;
	for (i = 0; i < dims; i++) {
		length *= shape[i];
	}

	free(previous_data);
	if(length == 0) {
		(*t)->data = NULL;
	} else {
		(*t)->data = (double*) malloc(sizeof(double) * length);
		if((*t)->data == NULL) {
			printf("ERROR: Could not allocate memory for the matrix elements (%d elements)\n", length);
			exit(EXIT_FAILURE);
        }
	}
	

	(*t)->length = length;

	free(previous_shape);
	(*t)->shape = (int*) malloc(sizeof(int) * dims);
	if((*t)->shape == NULL) {
		printf("ERROR: Could not allocated memory for the matrix shape\n");
		exit(EXIT_FAILURE);
    }	
	
	for (i = 0; i < dims; ++i) {
		(*t)->shape[i] = shape[i];
	}

	while (dims > 2 && shape[dims - 1] == 1) {
		--dims;
	}

	(*t)->dims = dims;

	/* Data is owned by the tensor, since it allocated it */
	(*t)->owns_data = 1;
	
	return *t;
}

/**
 */
int is_same_shape_alloc_d(tensor_d* t, int* restrict shape, int dims)
{
	int i;
	int result;
	
	/* Check if it has the same number of dimensions */
	if(t->dims != dims) {
		return 0;
	}

	/* As default, result is one */
	result = 1;
	
	/* Check if all dimensions of the tensor are the same */
	for(i=0; i<dims; i++) {
		if(t->shape[i] != shape[i]) {
			result = 0;
		}
	}

	/* Return result */
	
	return result;}

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_d(tensor_d** restrict t)
{

	/* If already null, return */
	if(*t == NULL) {
		return;
	}

	/* Free the values, if tensor owns the data */
	if((*t)->owns_data) {
		free((*t)->data);
		(*t)->data = NULL;
	}
	
	/* Free the shape data */
	free((*t)->shape);
	(*t)->shape = NULL;
	
	/* Free the tensor itself */
	free(*t);
	
	/* Set the pointer to null */
	*t = NULL;
}

/**
 */
tensor_d* new_array_d_2(int dim_1, int dim_2, tensor_d** restrict t)
{
/* new_row(int dim1, int dim2..., tensor** t) */
	/* int* shape; */
	int shape[2];
	int dims;

	shape[0] = dim_1 > 0 ? dim_1: 0;
	shape[1] = dim_2 > 0 ? dim_2: 0;

/*	shape = (int[2]){<DIM_NAMES>}; */
	dims = 2;
	
	return new_array_helper_d(shape,  dims, t);	
}

/**
 */
int ndims_alloc_d(tensor_d* t)
{
   return t->dims;
}

/**
 */
tensor_d* new_array_td_d(tensor_d* shape, tensor_d** restrict t)
{
/*tensor_d** new_array_d(tensor_i* shape, tensor_d** t) */
	int dims = shape->length;
	int* newShape = malloc(sizeof(int) * dims);
	int i;

	for(i=0;i<dims;i++) {
		newShape[i] = shape->data[i];

	}

	new_array_helper_d(newShape,  dims, t);
	free(newShape);
	return *t;
}

/**
 */
tensor_d* new_array_ti_d(tensor_i* shape, tensor_d** restrict t)
{
/*tensor_d** new_array_d(tensor_i* shape, tensor_d** t) */
	int dims = shape->length;
	int* newShape = malloc(sizeof(int) * dims);
	int i;

	for(i=0;i<dims;i++) {
		newShape[i] = shape->data[i];

	}

	new_array_helper_d(newShape,  dims, t);
	free(newShape);
	return *t;
}

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_i(tensor_i** restrict t)
{

	/* If already null, return */
	if(*t == NULL) {
		return;
	}

	/* Free the values, if tensor owns the data */
	if((*t)->owns_data) {
		free((*t)->data);
		(*t)->data = NULL;
	}
	
	/* Free the shape data */
	free((*t)->shape);
	(*t)->shape = NULL;
	
	/* Free the tensor itself */
	free(*t);
	
	/* Set the pointer to null */
	*t = NULL;
}

/**
 */
tensor_i* size_d(tensor_d* t, tensor_i** restrict size_array)
{
   int i;

   new_array_i_2(1, ndims_alloc_d(t), size_array);
   for(i = 0; i < (ndims_alloc_d(t)); ++i){
      (*size_array)->data[i] = t->shape[i];
   }
   
   
   return *size_array;
}

/**
 */
tensor_i* new_array_i_2(int dim_1, int dim_2, tensor_i** restrict t)
{
/* new_row(int dim1, int dim2..., tensor** t) */
	/* int* shape; */
	int shape[2];
	int dims;

	shape[0] = dim_1 > 0 ? dim_1: 0;
	shape[1] = dim_2 > 0 ? dim_2: 0;

/*	shape = (int[2]){<DIM_NAMES>}; */
	dims = 2;
	
	return new_array_helper_i(shape,  dims, t);	
}

/**
 */
tensor_i* new_array_helper_i(int* shape, int dims, tensor_i** restrict t)
{
/* new_row(int* shape, int dims, tensor** t) */

	int TRIES = 2;
	int i;
	int length;
	int sameShape;
	int* previous_shape = NULL;
	int* previous_data = NULL;
	
	/* Check if matrix is already allocated */
	if(*t != NULL) {
	
		/* If shape is the same, return current matrix, */
		/* even if it does not own the data (so it can implement window-writing) */
		sameShape = is_same_shape_alloc_i(*t, shape, dims);	
		//if(sameShape && (*t)->owns_data) {
		if(sameShape) {
			return *t;
		} 		
		
		/* Save pointers to previous shape and data */
		/* Only free data if tensor owns it */
		if((*t)->owns_data) {
			previous_data = (*t)->data;
		}
		
		previous_shape = (*t)->shape;		
	}

	/* Create the tensor and return it. */
	free(*t);
	
	for(i=0; i<TRIES; i++) {
		*t = (tensor_i*) malloc(sizeof(tensor_i));
		if (*t != NULL) {
			break;
		}
	}
	/**t = (tensor_i*) malloc(sizeof(tensor_i)); */

	if (*t == NULL) {
       printf("ERROR: Could not allocate memory for the matrix structure\n");
	   exit(EXIT_FAILURE);
	}

	/* Calculate the length of the linearized version of the tensor. */
	length = 1;
	for (i = 0; i < dims; i++) {
		length *= shape[i];
	}

	free(previous_data);
	if(length == 0) {
		(*t)->data = NULL;
	} else {
		(*t)->data = (int*) malloc(sizeof(int) * length);
		if((*t)->data == NULL) {
			printf("ERROR: Could not allocate memory for the matrix elements (%d elements)\n", length);
			exit(EXIT_FAILURE);
        }
	}
	

	(*t)->length = length;

	free(previous_shape);
	(*t)->shape = (int*) malloc(sizeof(int) * dims);
	if((*t)->shape == NULL) {
		printf("ERROR: Could not allocated memory for the matrix shape\n");
		exit(EXIT_FAILURE);
    }	
	
	for (i = 0; i < dims; ++i) {
		(*t)->shape[i] = shape[i];
	}

	while (dims > 2 && shape[dims - 1] == 1) {
		--dims;
	}

	(*t)->dims = dims;

	/* Data is owned by the tensor, since it allocated it */
	(*t)->owns_data = 1;
	
	return *t;
}

/**
 */
int is_same_shape_alloc_i(tensor_i* t, int* restrict shape, int dims)
{
	int i;
	int result;
	
	/* Check if it has the same number of dimensions */
	if(t->dims != dims) {
		return 0;
	}

	/* As default, result is one */
	result = 1;
	
	/* Check if all dimensions of the tensor are the same */
	for(i=0; i<dims; i++) {
		if(t->shape[i] != shape[i]) {
			result = 0;
		}
	}

	/* Return result */
	
	return result;}

/**
 */
tensor_d* new_array_d_3(int dim_1, int dim_2, int dim_3, tensor_d** restrict t)
{
/* new_row(int dim1, int dim2..., tensor** t) */
	/* int* shape; */
	int shape[3];
	int dims;

	shape[0] = dim_1 > 0 ? dim_1: 0;
	shape[1] = dim_2 > 0 ? dim_2: 0;
	shape[2] = dim_3 > 0 ? dim_3: 0;

/*	shape = (int[3]){<DIM_NAMES>}; */
	dims = 3;
	
	return new_array_helper_d(shape,  dims, t);	
}

/**
 */
tensor_d* new_const_array_d2(int dim_1, int dim_2, double value, tensor_d** restrict t)
{

   new_array_d_2(dim_1, dim_2, t);
   if((*t)->owns_data) {
	   set_matrix_values_alloc_d(*t, value);
   } else {
      printf("WARNING (new_const_array): Call to zeros for an array that does not own its data.\n");
   }
   
   return *t;
}

/**
 */
void set_matrix_values_alloc_d(tensor_d* t, double value)
{
/* set_matrix_values(tensor* t, elementType value) */

	int i;
	
	/* Set the values inside the tensor */
	for(i = 0; i<t->length; i = i+1)
   {
      t->data[i] = value;
   }
}

/**
 */
tensor_d* new_const_array_d3(int dim_1, int dim_2, int dim_3, double value, tensor_d** restrict t)
{

   new_array_d_3(dim_1, dim_2, dim_3, t);
   if((*t)->owns_data) {
	   set_matrix_values_alloc_d(*t, value);
   } else {
      printf("WARNING (new_const_array): Call to zeros for an array that does not own its data.\n");
   }
   
   return *t;
}

/**
 */
tensor_c* new_const_array_c2(int dim_1, int dim_2, char value, tensor_c** restrict t)
{

   new_array_c_2(dim_1, dim_2, t);
   if((*t)->owns_data) {
	   set_matrix_values_alloc_c(*t, value);
   } else {
      printf("WARNING (new_const_array): Call to zeros for an array that does not own its data.\n");
   }
   
   return *t;
}

/**
 */
void set_matrix_values_alloc_c(tensor_c* t, char value)
{
/* set_matrix_values(tensor* t, elementType value) */

	int i;
	
	/* Set the values inside the tensor */
	for(i = 0; i<t->length; i = i+1)
   {
      t->data[i] = value;
   }
}

/**
 */
tensor_c* new_array_c_2(int dim_1, int dim_2, tensor_c** restrict t)
{
/* new_row(int dim1, int dim2..., tensor** t) */
	/* int* shape; */
	int shape[2];
	int dims;

	shape[0] = dim_1 > 0 ? dim_1: 0;
	shape[1] = dim_2 > 0 ? dim_2: 0;

/*	shape = (int[2]){<DIM_NAMES>}; */
	dims = 2;
	
	return new_array_helper_c(shape,  dims, t);	
}

/**
 */
tensor_c* new_array_helper_c(int* shape, int dims, tensor_c** restrict t)
{
/* new_row(int* shape, int dims, tensor** t) */

	int TRIES = 2;
	int i;
	int length;
	int sameShape;
	int* previous_shape = NULL;
	char* previous_data = NULL;
	
	/* Check if matrix is already allocated */
	if(*t != NULL) {
	
		/* If shape is the same, return current matrix, */
		/* even if it does not own the data (so it can implement window-writing) */
		sameShape = is_same_shape_alloc_c(*t, shape, dims);	
		//if(sameShape && (*t)->owns_data) {
		if(sameShape) {
			return *t;
		} 		
		
		/* Save pointers to previous shape and data */
		/* Only free data if tensor owns it */
		if((*t)->owns_data) {
			previous_data = (*t)->data;
		}
		
		previous_shape = (*t)->shape;		
	}

	/* Create the tensor and return it. */
	free(*t);
	
	for(i=0; i<TRIES; i++) {
		*t = (tensor_c*) malloc(sizeof(tensor_c));
		if (*t != NULL) {
			break;
		}
	}
	/**t = (tensor_c*) malloc(sizeof(tensor_c)); */

	if (*t == NULL) {
       printf("ERROR: Could not allocate memory for the matrix structure\n");
	   exit(EXIT_FAILURE);
	}

	/* Calculate the length of the linearized version of the tensor. */
	length = 1;
	for (i = 0; i < dims; i++) {
		length *= shape[i];
	}

	free(previous_data);
	if(length == 0) {
		(*t)->data = NULL;
	} else {
		(*t)->data = (char*) malloc(sizeof(char) * length);
		if((*t)->data == NULL) {
			printf("ERROR: Could not allocate memory for the matrix elements (%d elements)\n", length);
			exit(EXIT_FAILURE);
        }
	}
	

	(*t)->length = length;

	free(previous_shape);
	(*t)->shape = (int*) malloc(sizeof(int) * dims);
	if((*t)->shape == NULL) {
		printf("ERROR: Could not allocated memory for the matrix shape\n");
		exit(EXIT_FAILURE);
    }	
	
	for (i = 0; i < dims; ++i) {
		(*t)->shape[i] = shape[i];
	}

	while (dims > 2 && shape[dims - 1] == 1) {
		--dims;
	}

	(*t)->dims = dims;

	/* Data is owned by the tensor, since it allocated it */
	(*t)->owns_data = 1;
	
	return *t;
}

/**
 */
int is_same_shape_alloc_c(tensor_c* t, int* restrict shape, int dims)
{
	int i;
	int result;
	
	/* Check if it has the same number of dimensions */
	if(t->dims != dims) {
		return 0;
	}

	/* As default, result is one */
	result = 1;
	
	/* Check if all dimensions of the tensor are the same */
	for(i=0; i<dims; i++) {
		if(t->shape[i] != shape[i]) {
			result = 0;
		}
	}

	/* Return result */
	
	return result;}

/**
 *  Deallocates the space used by the tensor;
 * 
 *  @param t
 *  		the tensor to free
 */
void tensor_free_c(tensor_c** restrict t)
{

	/* If already null, return */
	if(*t == NULL) {
		return;
	}

	/* Free the values, if tensor owns the data */
	if((*t)->owns_data) {
		free((*t)->data);
		(*t)->data = NULL;
	}
	
	/* Free the shape data */
	free((*t)->shape);
	(*t)->shape = NULL;
	
	/* Free the tensor itself */
	free(*t);
	
	/* Set the pointer to null */
	*t = NULL;
}
