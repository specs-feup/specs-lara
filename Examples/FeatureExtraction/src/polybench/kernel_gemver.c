void kernel_gemver(int n,
     float alpha,
     float beta,
     float A[ n + 0][n + 0],
     float u1[ n + 0],
     float v1[ n + 0],
     float u2[ n + 0],
     float v2[ n + 0],
     float w[ n + 0],
     float x[ n + 0],
     float y[ n + 0],
     float z[ n + 0])
{
  int i, j;

#pragma scop

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      A[i][j] = A[i][j] + u1[i] * v1[j] + u2[i] * v2[j];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      x[i] = x[i] + beta * A[j][i] * y[j];

  for (i = 0; i < n; i++)
    x[i] = x[i] + z[i];

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      w[i] = w[i] + alpha * A[i][j] * x[j];

#pragma endscop
}
