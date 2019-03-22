void kernel_covariance(int m, int n,
         float float_n,
         float data[ n + 0][m + 0],
         float cov[ m + 0][m + 0],
         float mean[ m + 0])
{
  int i, j, k;

#pragma scop
  for (j = 0; j < m; j++)
    {
      mean[j] = 0.0f;
      for (i = 0; i < n; i++)
        mean[j] += data[i][j];
      mean[j] /= float_n;
    }

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      data[i][j] -= mean[j];

  for (i = 0; i < m; i++)
    for (j = i; j < m; j++)
      {
        cov[i][j] = 0.0f;
        for (k = 0; k < n; k++)
   cov[i][j] += data[k][i] * data[k][j];
        cov[i][j] /= (float_n - 1.0f);
        cov[j][i] = cov[i][j];
      }
#pragma endscop

}
