void kernel_mvt(int n,
  float x1[ n + 0],
  float x2[ n + 0],
  float y_1[ n + 0],
  float y_2[ n + 0],
  float A[ n + 0][n + 0])
{
  int i, j;

#pragma scop
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      x1[i] = x1[i] + A[i][j] * y_1[j];
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      x2[i] = x2[i] + A[j][i] * y_2[j];
#pragma endscop

}
