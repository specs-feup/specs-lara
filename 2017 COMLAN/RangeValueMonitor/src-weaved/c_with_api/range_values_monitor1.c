#include <stdio.h>
#include <stdlib.h>
double monitor1_range_min[2] = {0};
double monitor1_range_max[2] = {0};
void monitor1_range_init() {

	unsigned int i;
	for(i=0; i < 2; i++) {

		monitor1_range_min[i] = 1.0/0.0;
		monitor1_range_max[i] = -1.0/0.0;
	}
}void monitor1_range_update(unsigned int id, double value) {
	if(value < monitor1_range_min[id]) monitor1_range_min[id] = value;
	if(value > monitor1_range_max[id]) monitor1_range_max[id] = value;
}

void monitor1_print_ranges() {
   FILE *log_file_0 = fopen("ranges_foobar.c.txt", "a+");
   if (log_file_0 == NULL)
   {
       printf("Error opening file ranges_foobar.c.txt\n");
       exit(1);
   } 
   fprintf(log_file_0, "foo\n");
   fprintf(log_file_0, "\ta: {%f, %f}\n", monitor1_range_min[0], monitor1_range_max[0]);
   fprintf(log_file_0, "lib_call\n");
   fprintf(log_file_0, "\tb: {%f, %f}\n", monitor1_range_min[1], monitor1_range_max[1]);
   // End of function
   fclose(log_file_0);
}
