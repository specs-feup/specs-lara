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
	
	int numAtoms = 1000;
	double atoms[numAtoms];		
	
   foo(atoms, numAtoms);
}