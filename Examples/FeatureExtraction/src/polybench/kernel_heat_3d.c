void kernel_heat_3d(int tsteps,
        int n,
        float A[ n + 0][n + 0][n + 0],
        float B[ n + 0][n + 0][n + 0])
{
  int t, i, j, k;

#pragma scop
    for (t = 1; t <= 20; t++) {
        for (i = 1; i < n-1; i++) {
            for (j = 1; j < n-1; j++) {
                for (k = 1; k < n-1; k++) {
                    B[i][j][k] = 0.125f * (A[i+1][j][k] - 2.0f * A[i][j][k] + A[i-1][j][k])
                                 + 0.125f * (A[i][j+1][k] - 2.0f * A[i][j][k] + A[i][j-1][k])
                                 + 0.125f * (A[i][j][k+1] - 2.0f * A[i][j][k] + A[i][j][k-1])
                                 + A[i][j][k];
                }
            }
        }
        for (i = 1; i < n-1; i++) {
           for (j = 1; j < n-1; j++) {
               for (k = 1; k < n-1; k++) {
                   A[i][j][k] = 0.125f * (B[i+1][j][k] - 2.0f * B[i][j][k] + B[i-1][j][k])
                                + 0.125f * (B[i][j+1][k] - 2.0f * B[i][j][k] + B[i][j-1][k])
                                + 0.125f * (B[i][j][k+1] - 2.0f * B[i][j][k] + B[i][j][k-1])
                                + B[i][j][k];
               }
           }
       }
    }
#pragma endscop

}
