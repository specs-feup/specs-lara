#include <chrono>
#include <iostream>

double bar() {
   
   return 1.0;
}


double foo() {
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      std::chrono::high_resolution_clock::time_point clava_timing_start_0 = std::chrono::high_resolution_clock::now();
      a += bar();
      std::chrono::high_resolution_clock::time_point clava_timing_end_0 = std::chrono::high_resolution_clock::now();
      auto clava_timing_duration_0 = std::chrono::duration_cast<std::chrono::microseconds>(clava_timing_end_0 - clava_timing_start_0).count();
      std::cout << "Time:" << clava_timing_duration_0 << "microseconds" << "\n";
   }
   
   return a;
}


int main() {
   std::chrono::high_resolution_clock::time_point clava_timing_start_1 = std::chrono::high_resolution_clock::now();
   foo();
   std::chrono::high_resolution_clock::time_point clava_timing_end_1 = std::chrono::high_resolution_clock::now();
   auto clava_timing_duration_1 = std::chrono::duration_cast<std::chrono::microseconds>(clava_timing_end_1 - clava_timing_start_1).count();
   std::cout << "Time:" << clava_timing_duration_1 << "microseconds" << "\n";
}
