#include "lib.h"
#include "../range_values_monitor1.h"

int lib_call(int a) {
   int b;
   /* monitoring b */
   b = a + 1;
   monitor1_range_update(1, b);
   
   return b;
}
