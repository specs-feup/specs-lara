void kernel_symm(int m, int n,
   float alpha,
   float beta,
   float C[ m + 0][n + 0],
   float A[ m + 0][m + 0],
   float B[ m + 0][n + 0])
{
  int i, j, k;
  float temp2;
# 75 "symm.c"
#pragma scop
   for (i = 0; i < m; i++)
      for (j = 0; j < n; j++ )
      {
        temp2 = 0;
        for (k = 0; k < i; k++) {
           C[k][j] += alpha*B[i][j] * A[i][k];
           temp2 += B[k][j] * A[i][k];
        }
        C[i][j] = beta * C[i][j] + alpha*B[i][j] * A[i][i] + alpha * temp2;
     }
#pragma endscop

}
