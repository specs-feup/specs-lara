/* Implementation file for lib/matisse_string */

#include "matisse_string.h"
#include <stdlib.h>
#include "tensor_struct.h"


/**
 */
int strcmp_partially_checked_1tc(tensor_c* in_m, const char * in_s)
{
   int matrix_size;
   int i;

   matrix_size = in_m->length;
   if(matrix_size != 4){
      
      return 0;
   }
   
   for(i = 0; i < 4; ++i){
      if(in_m->data[i] != in_s[i]){
         
         return 0;
      }
      
   }
   
   
   return 1;
}
