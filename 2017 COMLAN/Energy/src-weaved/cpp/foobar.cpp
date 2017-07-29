#include "rapl.h"
#include <iostream>

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      rapl_monitor_start();
      a += bar();
      double clava_energy_output_0 = rapl_monitor_report();
      std::cout << "Energy:" << clava_energy_output_0;
   }
   
   return a;
}


int main() {
   rapl_monitor_start();
   foo();
   double clava_energy_output_1 = rapl_monitor_report();
   std::cout << "Energy:" << clava_energy_output_1;
}
