void kernel_2mm(int ni, int nj, int nk, int nl,
  float alpha,
  float beta,
  float tmp[ ni + 0][nj + 0],
  float A[ ni + 0][nk + 0],
  float B[ nk + 0][nj + 0],
  float C[ nj + 0][nl + 0],
  float D[ ni + 0][nl + 0])
{
  int i, j, k;

#pragma scop

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      {
 tmp[i][j] = 0.0f;
 for (k = 0; k < nk; ++k)
   tmp[i][j] += alpha * A[i][k] * B[k][j];
      }
  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      {
 D[i][j] *= beta;
 for (k = 0; k < nj; ++k)
   D[i][j] += tmp[i][k] * C[k][j];
      }
#pragma endscop

}

