#include <stdio.h>
#include <stdlib.h>

void foo(double* a, double* b, double* c, int numElems) {
	
		// Loop to parallelize
		for(int i=0; i<numElems; i++) {
			c[i] = a[i] + b[i];
		}
	
}

int main() {
		
	int numElems = 1011;	
		
	double* a;
	double* b;
	double* c;
	
	a = (double*) malloc(sizeof(double)*numElems);
	b = (double*) malloc(sizeof(double)*numElems);
	c = (double*) malloc(sizeof(double)*numElems);
	
	for(int i=0; i<numElems; i++) {
		a[i] = i;
		b[i] = i + 1;
	}
	
	foo(a, b, c, numElems);
	
	// test output
	double acc = 0;
	for(int i=0; i<numElems; i++) {
		acc += c[i];
	}
	
	printf("Result: %f\n", acc);
	
	return 0;
}