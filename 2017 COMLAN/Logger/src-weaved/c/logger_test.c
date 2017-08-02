#include <stdio.h>
#include <stdlib.h>

double bar() {
   
   return 1.0;
}


double foo() {
   FILE * log_file_0 = fopen("log.txt", "a+");
   if (log_file_0 == NULL){
   printf("Error opening file log.txt\n");
   exit(1);
   }
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      fprintf(log_file_0, "Logging to a file\n");
      a += bar();
      printf("Printing again\n");
      fprintf(log_file_0, "Logging again to a file\n");
      printf("Print double %f after bar\n", 2.0);
   }
   fclose(log_file_0);
   
   return a;
}


int main() {
   FILE * log_file_1 = fopen("log.txt", "a+");
   if (log_file_1 == NULL){
   printf("Error opening file log.txt\n");
   exit(1);
   }
   fprintf(log_file_1, "Logging to a file\n");
   foo();
   printf("Printing again\n");
   fprintf(log_file_1, "Logging again to a file\n");
   printf("Print double %f after bar\n", 2.0);
   fclose(log_file_1);
}
