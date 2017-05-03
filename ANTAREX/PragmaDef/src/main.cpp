#define N 1000

int main() {
	
	int a[N];
	
	#pragma clava $loop isParallel, iterations = 1000
	for(int i=0; i<N; i++) {
		a[i] = i + 10;
	}
	
	#pragma clava iterations = 1000
	for(int i=0; i<N; i++) {
		a[i] = i + 10;
	}
	
	#pragma clava $loop isParallel
	for(int i=0; i<N; i++) {
		a[i] = i + 10;
	}
}