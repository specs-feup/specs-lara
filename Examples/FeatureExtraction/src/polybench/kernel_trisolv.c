void kernel_trisolv(int n,
      float L[ n + 0][n + 0],
      float x[ n + 0],
      float b[ n + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i++)
    {
      x[i] = b[i];
      for (j = 0; j <i; j++)
        x[i] -= L[i][j] * x[j];
      x[i] = x[i] / L[i][i];
    }
#pragma endscop

}
