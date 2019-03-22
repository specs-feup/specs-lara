void kernel_durbin(int n,
     float r[ n + 0],
     float y[ n + 0])
{
 float z[40];
 float alpha;
 float beta;
 float sum;

 int i,k;

#pragma scop
 y[0] = -r[0];
 beta = 1.0f;
 alpha = -r[0];

 for (k = 1; k < n; k++) {
   beta = (1-alpha*alpha)*beta;
   sum = 0.0f;
   for (i=0; i<k; i++) {
      sum += r[k-i-1]*y[i];
   }
   alpha = - (r[k] + sum)/beta;

   for (i=0; i<k; i++) {
      z[i] = y[i] + alpha*y[k-i-1];
   }
   for (i=0; i<k; i++) {
     y[i] = z[i];
   }
   y[k] = alpha;
 }
#pragma endscop

}
