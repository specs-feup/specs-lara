__kernel void mat_mul_kernel(global double* A, global double* B, global double* C, uint N,  uint M, uint K)
{
	size_t i = get_global_id(0U);
	size_t j = get_global_id(1U);

	if(i < N && j < K){
		
		double acc;
		int pos;
		
		acc = 0.f;
		for(pos = 0; pos < M; ++pos){
			acc += B[j + (pos) * K] * A[pos + (i) * M];
		}

		C[j + (i) * K] = acc;
	}
}
