#include <stdio.h>
#define N 1011


void foo(double a[N], double b[N], double c[N]) {
	
	{
		// Loop to parallelize
		#pragma clava mpi_scather_gather <additional parameters...>
		for(int i=0; i<N; i++) {
			c[i] = a[i] + b[i];
		}
	}
	
}

int main() {
		
	double a[N], b[N], c[N];
	
	for(int i=0; i<N; i++) {
		a[i] = i;
		b[i] = i + 1;
	}
	
	foo(a, b, c);
	
	// test output
	double acc = 0;
	for(int i=0; i<N; i++) {
		acc += c[i];
	}
	
	printf("Result: %f\n", acc);
	
	return 0;
}