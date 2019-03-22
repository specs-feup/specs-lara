void kernel_gemm(int ni, int nj, int nk,
   float alpha,
   float beta,
   float C[ ni + 0][nj + 0],
   float A[ ni + 0][nk + 0],
   float B[ nk + 0][nj + 0])
{
  int i, j, k;

#pragma scop
  for (i = 0; i < ni; i++) {
    for (j = 0; j < nj; j++)
 C[i][j] *= beta;
    for (k = 0; k < nk; k++) {
       for (j = 0; j < nj; j++)
   C[i][j] += alpha * A[i][k] * B[k][j];
    }
  }
#pragma endscop

}
