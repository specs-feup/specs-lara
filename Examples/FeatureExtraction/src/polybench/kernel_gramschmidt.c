#include <math.h>

void kernel_gramschmidt(int m, int n,
   float A[ m + 0][n + 0],
   float R[ n + 0][n + 0],
   float Q[ m + 0][n + 0])
{
  int i, j, k;

  float nrm;

#pragma scop
  for (k = 0; k < n; k++)
    {
      nrm = 0.0f;
      for (i = 0; i < m; i++)
        nrm += A[i][k] * A[i][k];
      R[k][k] = sqrtf(nrm);
      for (i = 0; i < m; i++)
        Q[i][k] = A[i][k] / R[k][k];
      for (j = k + 1; j < n; j++)
 {
   R[k][j] = 0.0f;
   for (i = 0; i < m; i++)
     R[k][j] += Q[i][k] * A[i][j];
   for (i = 0; i < m; i++)
     A[i][j] = A[i][j] - Q[i][k] * R[k][j];
 }
    }
#pragma endscop

}
