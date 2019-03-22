void kernel_trmm(int m, int n,
   float alpha,
   float A[ m + 0][m + 0],
   float B[ m + 0][n + 0])
{
  int i, j, k;
  float temp;
# 69 "trmm.c"
#pragma scop
  for (i = 0; i < m; i++)
     for (j = 0; j < n; j++) {
        for (k = i+1; k < m; k++)
           B[i][j] += A[k][i] * B[k][j];
        B[i][j] = alpha * B[i][j];
     }
#pragma endscop

}

