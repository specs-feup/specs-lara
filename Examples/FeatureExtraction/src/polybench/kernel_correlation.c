#include <math.h>

void kernel_correlation(int m, int n,
   float float_n,
   float data[ n + 0][m + 0],
   float corr[ m + 0][m + 0],
   float mean[ m + 0],
   float stddev[ m + 0])
{
  int i, j, k;

  float eps = 0.1f;


#pragma scop
  for (j = 0; j < m; j++)
    {
      mean[j] = 0.0f;
      for (i = 0; i < n; i++)
 mean[j] += data[i][j];
      mean[j] /= float_n;
    }


   for (j = 0; j < m; j++)
    {
      stddev[j] = 0.0f;
      for (i = 0; i < n; i++)
        stddev[j] += (data[i][j] - mean[j]) * (data[i][j] - mean[j]);
      stddev[j] /= float_n;
      stddev[j] = sqrtf(stddev[j]);



      stddev[j] = stddev[j] <= eps ? 1.0f : stddev[j];
    }


  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++)
      {
        data[i][j] -= mean[j];
        data[i][j] /= sqrtf(float_n) * stddev[j];
      }


  for (i = 0; i < m-1; i++)
    {
      corr[i][i] = 1.0f;
      for (j = i+1; j < m; j++)
        {
          corr[i][j] = 0.0f;
          for (k = 0; k < n; k++)
            corr[i][j] += (data[k][i] * data[k][j]);
          corr[j][i] = corr[i][j];
        }
    }
  corr[m-1][m-1] = 1.0f;
#pragma endscop

}

