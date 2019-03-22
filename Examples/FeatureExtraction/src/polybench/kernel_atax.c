void kernel_atax(int m, int n,
   float A[ m + 0][n + 0],
   float x[ n + 0],
   float y[ n + 0],
   float tmp[ m + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i++)
    y[i] = 0;
  for (i = 0; i < m; i++)
    {
      tmp[i] = 0.0f;
      for (j = 0; j < n; j++)
 tmp[i] = tmp[i] + A[i][j] * x[j];
      for (j = 0; j < n; j++)
 y[j] = y[j] + A[i][j] * tmp[i];
    }
#pragma endscop

}
