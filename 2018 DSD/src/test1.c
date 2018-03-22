double bar(double a) {
   
   return a + 1.0;
}


double foo() {
   double a = 0;
   #pragma omp parallel for
   for(int i = 0; i < 1000; i++) {
      a += bar(a);
      #pragma omp parallel for
      for(int j = 0; j < 1000; j++) {
         a += bar(a);
      }
   }
   
   return a;
}


int main() {
   foo();
}