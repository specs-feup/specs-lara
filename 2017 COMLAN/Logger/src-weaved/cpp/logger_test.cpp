#include <iostream>
#include <fstream>

double bar() {
   
   return 1.0;
}


double foo() {
   std::ofstream log_file_0;
   log_file_0.open("log.txt", std::ios_base::app);
   double a = 0;
   for(int i = 0; i < 1000; i++) {
      std::cout << "Print double " << 2 << " after bar" << "\n";
      log_file_0 << "Logging to a file" << "\n";
      a += bar();
      log_file_0 << "Logging again to a file" << "\n";
      std::cout << "Printing again" << "\n";
   }
   
   return a;
}


int main() {
   std::ofstream log_file_1;
   log_file_1.open("log.txt", std::ios_base::app);
   std::cout << "Print double " << 2 << " after foo" << "\n";
   log_file_1 << "Logging to a file" << "\n";
   foo();
   log_file_1 << "Logging again to a file" << "\n";
   std::cout << "Printing again" << "\n";
}
