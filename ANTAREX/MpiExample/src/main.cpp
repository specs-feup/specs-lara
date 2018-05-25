#include <stdio.h>
#define N 1011

int main(int argc, char **argv) {
		
	double a[N], b[N], c[N];
	
	for(int i=0; i<N; i++) {
		a[i] = i;
		b[i] = i + 1;
	}
	
	
	// Loop to parallelize
	#pragma clava parallel	
	for(int i=0; i<N; i++) {
		c[i] = a[i] + b[i];
	}

	// test output
	double acc = 0;
	for(int i=0; i<N; i++) {
		acc += c[i];
	}
	
	printf("Result: %f\n", acc);
	
	return 0;
}