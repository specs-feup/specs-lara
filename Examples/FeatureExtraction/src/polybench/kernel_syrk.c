void kernel_syrk(int n, int m,
   float alpha,
   float beta,
   float C[ n + 0][n + 0],
   float A[ n + 0][m + 0])
{
  int i, j, k;







#pragma scop
  for (i = 0; i < n; i++) {
    for (j = 0; j <= i; j++)
      C[i][j] *= beta;
    for (k = 0; k < m; k++) {
      for (j = 0; j <= i; j++)
        C[i][j] += alpha * A[i][k] * A[j][k];
    }
  }
#pragma endscop

}
