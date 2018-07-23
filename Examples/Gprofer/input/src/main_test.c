/* Implementation file for main_test */

#include "disparity_matisse_modified.h"
#include "lib/array_creators_alloc.h"
#include "lib/load.h"
#include "lib/tensor_struct.h"
#include <stdlib.h>


#include <stdio.h>
#include <time.h>

#define MULTI_EXEC 1

/**
 */
int main(int argc, char** argv)
{
   tensor_d* imleft = NULL;
   tensor_d* imright = NULL;
   tensor_c* type = NULL;
   char type_static[6] = {102,117,108,108,104,100};
   tensor_d* imDispOwn = NULL;

   load_matrix_variable_2_td("imleft", 1920, 1080, &imleft);
   load_matrix_variable_2_td("imright", 1920, 1080, &imright);
   zeros_c2(1, 6, &type);
   free(type->data);
   type->data = type_static;
   
   // Initialize random seed, in case rand is used
   srand((unsigned) time(NULL));
   
   double timeElapsed;
   double start, end;
   #if MULTI_EXEC
   int iterations_idx;   
   int iterations = 1;
   #endif
   #if MULTI_EXEC
   if(argc>1)
   {
      iterations = atoi(argv[1]);
   }
   #endif
   
   #if MULTI_EXEC
   for(iterations_idx = 0; iterations_idx<iterations; iterations_idx++)
   {
   #endif
   disparity_matisse_modified_tctdtd_row_2d_2d_1(type, imleft, imright, &imDispOwn);
   #if MULTI_EXEC
   }
   #endif
   
   return (int) imDispOwn->data[0];
}
