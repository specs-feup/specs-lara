void kernel_adi(int tsteps, int n,
  float u[ n + 0][n + 0],
  float v[ n + 0][n + 0],
  float p[ n + 0][n + 0],
  float q[ n + 0][n + 0])
{
  int t, i, j;
  float DX, DY, DT;
  float B1, B2;
  float mul1, mul2;
  float a, b, c, d, e, f;

#pragma scop

  DX = 1.0f/(float)n;
  DY = 1.0f/(float)n;
  DT = 1.0f/(float)tsteps;
  B1 = 2.0f;
  B2 = 1.0f;
  mul1 = B1 * DT / (DX * DX);
  mul2 = B2 * DT / (DY * DY);

  a = -mul1 / 2.0f;
  b = 1.0f +mul1;
  c = a;
  d = -mul2 / 2.0f;
  e = 1.0f +mul2;
  f = d;

 for (t=1; t<=tsteps; t++) {

    for (i=1; i<n-1; i++) {
      v[0][i] = 1.0f;
      p[i][0] = 0.0f;
      q[i][0] = v[0][i];
      for (j=1; j<n-1; j++) {
        p[i][j] = -c / (a*p[i][j-1]+b);
        q[i][j] = (-d*u[j][i-1]+(1.0f +2.0f*d)*u[j][i] - f*u[j][i+1]-a*q[i][j-1])/(a*p[i][j-1]+b);
      }

      v[n-1][i] = 1.0f;
      for (j=n-2; j>=1; j--) {
        v[j][i] = p[i][j] * v[j+1][i] + q[i][j];
      }
    }

    for (i=1; i<n-1; i++) {
      u[i][0] = 1.0f;
      p[i][0] = 0.0f;
      q[i][0] = u[i][0];
      for (j=1; j<n-1; j++) {
        p[i][j] = -f / (d*p[i][j-1]+e);
        q[i][j] = (-a*v[i-1][j]+(1.0f +2.0f*a)*v[i][j] - c*v[i+1][j]-d*q[i][j-1])/(d*p[i][j-1]+e);
      }
      u[i][n-1] = 1.0f;
      for (j=n-2; j>=1; j--) {
        u[i][j] = p[i][j] * u[i][j+1] + q[i][j];
      }
    }
  }
#pragma endscop
}
