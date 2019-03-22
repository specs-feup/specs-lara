void kernel_3mm(int ni, int nj, int nk, int nl, int nm,
  float E[ ni + 0][nj + 0],
  float A[ ni + 0][nk + 0],
  float B[ nk + 0][nj + 0],
  float F[ nj + 0][nl + 0],
  float C[ nj + 0][nm + 0],
  float D[ nm + 0][nl + 0],
  float G[ ni + 0][nl + 0])
{
  int i, j, k;

#pragma scop

  for (i = 0; i < ni; i++)
    for (j = 0; j < nj; j++)
      {
 E[i][j] = 0.0f;
 for (k = 0; k < nk; ++k)
   E[i][j] += A[i][k] * B[k][j];
      }

  for (i = 0; i < nj; i++)
    for (j = 0; j < nl; j++)
      {
 F[i][j] = 0.0f;
 for (k = 0; k < nm; ++k)
   F[i][j] += C[i][k] * D[k][j];
      }

  for (i = 0; i < ni; i++)
    for (j = 0; j < nl; j++)
      {
 G[i][j] = 0.0f;
 for (k = 0; k < nj; ++k)
   G[i][j] += E[i][k] * F[k][j];
      }
#pragma endscop

}
