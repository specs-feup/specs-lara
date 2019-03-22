void kernel_doitgen(int nr, int nq, int np,
      float A[ nr + 0][nq + 0][np + 0],
      float C4[ np + 0][np + 0],
      float sum[ np + 0])
{
  int r, q, p, s;

#pragma scop
  for (r = 0; r < nr; r++)
    for (q = 0; q < nq; q++) {
      for (p = 0; p < np; p++) {
 sum[p] = 0.0f;
 for (s = 0; s < np; s++)
   sum[p] += A[r][q][s] * C4[s][p];
      }
      for (p = 0; p < np; p++)
 A[r][q][p] = sum[p];
    }
#pragma endscop

}
