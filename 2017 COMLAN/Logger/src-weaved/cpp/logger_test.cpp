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
      log_file_0 << "Logging to a file" << "\n";
      std::cout << "Print double " << 2 << " after bar" << "\n";
      a += bar();
      std::cout << "Printing again";
      log_file_0 << "Logging again to a file" << "\n";
      std::cout << "Print double " << 2 << " after bar" << "\n";
   }
   
   return a;
}


int main() {
   std::ofstream log_file_1;
   log_file_1.open("log.txt", std::ios_base::app);
   log_file_1 << "Logging to a file" << "\n";
   std::cout << "Print double " << 2 << " after bar" << "\n";
   foo();
   std::cout << "Printing again";
   log_file_1 << "Logging again to a file" << "\n";
   std::cout << "Print double " << 2 << " after foo" << "\n";
}
