#include "lib.h"
void range_update( unsigned int id, double value );
extern double range_min[ 2 ]; extern double range_max[ 2 ];

int lib_call(int a) {
   int b;
   /* monitoring b */
   b = a + 1;
   range_update(1, b);
   
   return b;
}
