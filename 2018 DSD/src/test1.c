#define NUM_ATOMS 1000

double bar(double a) {
   double b = 1.0;
   return a + b;
}


double foo(double* atoms, int numAtoms) {
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
	
	double atoms[NUM_ATOMS];		
	
	foo(atoms, NUM_ATOMS);
}