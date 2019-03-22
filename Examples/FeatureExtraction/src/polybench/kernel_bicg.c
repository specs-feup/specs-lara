void kernel_bicg(int m, int n,
   float A[ n + 0][m + 0],
   float s[ m + 0],
   float q[ n + 0],
   float p[ m + 0],
   float r[ n + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < m; i++)
    s[i] = 0;
  for (i = 0; i < n; i++)
    {
      q[i] = 0.0f;
      for (j = 0; j < m; j++)
 {
   s[j] = s[j] + r[i] * A[i][j];
   q[i] = q[i] + A[i][j] * p[j];
 }
    }
#pragma endscop

}
