
void kernel_syr2k(int n, int m,
    float alpha,
    float beta,
    float C[ n + 0][n + 0],
    float A[ n + 0][m + 0],
    float B[ n + 0][m + 0])
{
  int i, j, k;







#pragma scop
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      C[i][j] *= beta;
  for (i = 0; i < n; i++)
    for (k = 0; k < m; k++) {
      for (j = 0; j < n; j++)
 {
   C[i][j] += A[j][k] * alpha*B[i][k] + B[j][k] * alpha*A[i][k];
 }
     }
#pragma endscop

}

