#include "lib/lib.h"
#include "range_values_monitor1.h"

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      /* monitoring a */
      a += bar();
      monitor1_range_update(0, a);
   }
   /* monitoring a */
   a = lib_call(a);
   monitor1_range_update(0, a);
   
   return a;
}


int main() {
   monitor1_range_init();
   foo();
   monitor1_print_ranges();
}
