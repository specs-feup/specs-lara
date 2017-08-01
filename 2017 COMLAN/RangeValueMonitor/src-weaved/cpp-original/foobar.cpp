#include "lib/lib.h"
double range_min[ 2 ] = {0};
double range_max[ 2 ] = {0};

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      /* monitoring a */
      a += bar();
      range_update(0, a);
   }
   /* monitoring a */
   a = lib_call(a);
   range_update(0, a);
   
   return a;
}



void range_init() {

	unsigned int i;
	for(i=0; i < 2; i++) {

		range_min[i] = 1.0/0.0;
		range_max[i] = -1.0/0.0;
	}
}

void range_update(unsigned int id, double value) {

	if(value < range_min[id]) range_min[id] = value;
	if(value > range_max[id]) range_max[id] = value;
}

void print_ranges() {
	printf("foo\n");
printf("\ta: {%f, %f}\n", range_min[0], range_max[0]);
printf("lib_call\n");
printf("\tb: {%f, %f}\n", range_min[1], range_max[1]);

}
		

int main() {
   atexit(print_ranges);
   range_init();
   foo();
}
