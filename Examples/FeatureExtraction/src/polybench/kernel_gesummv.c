void kernel_gesummv(int n,
      float alpha,
      float beta,
      float A[ n + 0][n + 0],
      float B[ n + 0][n + 0],
      float tmp[ n + 0],
      float x[ n + 0],
      float y[ n + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i++)
    {
      tmp[i] = 0.0f;
      y[i] = 0.0f;
      for (j = 0; j < n; j++)
 {
   tmp[i] = A[i][j] * x[j] + tmp[i];
   y[i] = B[i][j] * x[j] + y[i];
 }
      y[i] = alpha * tmp[i] + beta * y[i];
    }
#pragma endscop

}
