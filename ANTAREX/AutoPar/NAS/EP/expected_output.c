#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
//---------------------------------------------------------------------
// program EP
//---------------------------------------------------------------------
//----------
//  Class S:
//----------
//----------
//  Class W:
//----------
//----------
//  Class A:
//----------
//----------
//  Class B:
//----------
//----------
//  Class C:
//----------
//----------
//  Class D:
//----------
//----------
//  Class E:
//----------

struct anon_NAS_EP_c_69 {
   double real;
   double imag;
};

typedef struct anon_NAS_EP_c_69 dcomplex;
double randlc(double *x, double a);
void vranlc(int n, double *x, double a, double y[]);
void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified);
double start[64];
double elapsed[64];
double elapsed_time();
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void wtime(double *t);
int main() {
   double Mops, t1, t2, t3, t4, x1, x2;
   double sx, sy, tm, an, tt, gc;
   double sx_verify_value, sy_verify_value, sx_err, sy_err;
   int np;
   int i, ik, kk, l, k, nit;
   int k_offset, j;
   int verified;
   double dum[3] = {1.0, 1.0, 1.0};
   char size[16];
   double x[131072];
   double q[10];
   FILE *fp;
   //--------------------------------------------------------------------
   //  Because the size of the problem is too large to store in a 32-bit
   //  integer for some classes, we put it into a string (for printing).
   //  Have to strip off the decimal point put in there by the floating
   //  point print statement (internal file)
   //--------------------------------------------------------------------
   sprintf(size, "%15.0lf", pow(2.0, 25 + 1));
   j = 14;
   if(size[j] == '.') j--;
   size[j + 1] = '\0';
   printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER-C) - EP Benchmark\n");
   printf("\n Number of random numbers generated: %15s\n", size);
   verified = 0;
   //--------------------------------------------------------------------
   //  Compute the number of "batches" of random number pairs generated
   //  per processor. Adjust if the number of processors does not evenly
   //  divide the total number
   //--------------------------------------------------------------------
   np = (1 << (25 - 16));
   //--------------------------------------------------------------------
   //  Call the random number generator functions and initialize
   //  the x-array to reduce the effects of paging on the timings.
   //  Also, call all mathematical functions that are used. Make
   //  sure these initializations cannot be eliminated as dead code.
   //--------------------------------------------------------------------
   vranlc(0, &dum[0], dum[1], &dum[2]);
   dum[0] = randlc(&dum[1], dum[2]);
   #pragma omp parallel for default(shared) private(i)
   for(i = 0; i < 2 * (1 << 16); i++) {
      x[i] = -1.0e99;
   }
   Mops = log(sqrt(fabs((((1.0) > (1.0)) ? (1.0) : (1.0)))));
   timer_clear(0);
   timer_clear(1);
   timer_clear(2);
   timer_start(0);
   t1 = 1220703125.0;
   vranlc(0, &t1, 1220703125.0, x);
   //--------------------------------------------------------------------
   //  Compute AN = A ^ (2 * NK) (mod 2^46).
   //--------------------------------------------------------------------
   t1 = 1220703125.0;
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 16 + 1; i++) {
      t2 = randlc(&t1, t1);
   }
   an = t1;
   tt = 271828183.0;
   gc = 0.0;
   sx = 0.0;
   sy = 0.0;
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 10; i++) {
      q[i] = 0.0;
   }
   //--------------------------------------------------------------------
   //  Each instance of this loop may be performed independently. We compute
   //  the k offsets separately to take into account the fact that some nodes
   //  have more numbers to generate than others
   //--------------------------------------------------------------------
   k_offset = -1;
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess q	 use : RW
   ****************************************/
   for(k = 1; k <= np; k++) {
      kk = k_offset + k;
      t1 = 271828183.0;
      t2 = an;
      // Find starting seed t1 for this kk.
      /*************** Clava msgError **************
      Loop contains Invalid Statement -> BreakStmt#189
      ****************************************/
      for(i = 1; i <= 100; i++) {
         ik = kk / 2;
         if((2 * ik) != kk) t3 = randlc(&t1, t2);
         if(ik == 0) break;
         t3 = randlc(&t2, t2);
         kk = ik;
      }
      //--------------------------------------------------------------------
      //  Compute uniform pseudorandom numbers.
      //--------------------------------------------------------------------
      vranlc(2 * (1 << 16), &t1, 1220703125.0, x);
      //--------------------------------------------------------------------
      //  Compute Gaussian deviates by acceptance-rejection method and
      //  tally counts in concentri//square annuli.  This loop is not
      //  vectorizable.
      //--------------------------------------------------------------------
      /*************** Clava msgError **************
      unsolved dependency for arrayAccess q	 use : RW
      ****************************************/
      for(i = 0; i < (1 << 16); i++) {
         x1 = 2.0 * x[2 * i] - 1.0;
         x2 = 2.0 * x[2 * i + 1] - 1.0;
         t1 = x1 * x1 + x2 * x2;
         if(t1 <= 1.0) {
            t2 = sqrt(-2.0 * log(t1) / t1);
            t3 = (x1 * t2);
            t4 = (x2 * t2);
            l = (((fabs(t3)) > (fabs(t4))) ? (fabs(t3)) : (fabs(t4)));
            q[l] = q[l] + 1.0;
            sx = sx + t3;
            sy = sy + t4;
         }
      }
   }
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 10; i++) {
      gc = gc + q[i];
   }
   timer_stop(0);
   tm = timer_read(0);
   nit = 0;
   verified = 1;
   if(25 == 24) {
      sx_verify_value = -3.247834652034740e+3;
      sy_verify_value = -6.958407078382297e+3;
   }
   else if(25 == 25) {
      sx_verify_value = -2.863319731645753e+3;
      sy_verify_value = -6.320053679109499e+3;
   }
   else if(25 == 28) {
      sx_verify_value = -4.295875165629892e+3;
      sy_verify_value = -1.580732573678431e+4;
   }
   else if(25 == 30) {
      sx_verify_value = 4.033815542441498e+4;
      sy_verify_value = -2.660669192809235e+4;
   }
   else if(25 == 32) {
      sx_verify_value = 4.764367927995374e+4;
      sy_verify_value = -8.084072988043731e+4;
   }
   else if(25 == 36) {
      sx_verify_value = 1.982481200946593e+5;
      sy_verify_value = -1.020596636361769e+5;
   }
   else if(25 == 40) {
      sx_verify_value = -5.319717441530e+05;
      sy_verify_value = -3.688834557731e+05;
   }
   else {
      verified = 0;
   }
   if(verified) {
      sx_err = fabs((sx - sx_verify_value) / sx_verify_value);
      sy_err = fabs((sy - sy_verify_value) / sy_verify_value);
      verified = ((sx_err <= 1.0e-8) && (sy_err <= 1.0e-8));
   }
   Mops = pow(2.0, 25 + 1) / tm / 1000000.0;
   printf("\nEP Benchmark Results:\n\n");
   printf("CPU Time =%10.4lf\n", tm);
   printf("N = 2^%5d\n", 25);
   printf("No. Gaussian Pairs = %15.0lf\n", gc);
   printf("Sums = %25.15lE %25.15lE\n", sx, sy);
   printf("Counts: \n");
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 10; i++) {
      printf("%3d%15.0lf\n", i, q[i]);
   }
   print_results("EP", 'W', 25 + 1, 0, 0, nit, tm, Mops, "Random numbers generated", verified);
   int exitValue = verified ? 0 : 1;
   
   return exitValue;
}

double randlc(double *x, double a) {
   //--------------------------------------------------------------------
   //
   //  This routine returns a uniform pseudorandom double precision number in the
   //  range (0, 1) by using the linear congruential generator
   //
   //  x_{k+1} = a x_k  (mod 2^46)
   //
   //  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   //  before repeating.  The argument A is the same as 'a' in the above formula,
   //  and X is the same as x_0.  A and X must be odd double precision integers
   //  in the range (1, 2^46).  The returned value RANDLC is normalized to be
   //  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
   //  the new seed x_1, so that subsequent calls to RANDLC using the same
   //  arguments will generate a continuous sequence.
   //
   //  This routine should produce the same results on any computer with at least
   //  48 mantissa bits in double precision floating point data.  On 64 bit
   //  systems, double precision should be disabled.
   //
   //  David H. Bailey     October 26, 1990
   //
   //--------------------------------------------------------------------
   // r23 = pow(0.5, 23.0);
   ////  pow(0.5, 23.0) = 1.1920928955078125e-07
   // r46 = r23 * r23;
   // t23 = pow(2.0, 23.0);
   ////  pow(2.0, 23.0) = 8.388608e+06
   // t46 = t23 * t23;
   double const r23 = 1.1920928955078125e-07;
   double const r46 = r23 * r23;
   double const t23 = 8.388608e+06;
   double const t46 = t23 * t23;
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   double r;
   //--------------------------------------------------------------------
   //  Break A into two parts such that A = 2^23 * A1 + A2.
   //--------------------------------------------------------------------
   t1 = r23 * a;
   a1 = (int) t1;
   a2 = a - t23 * a1;
   //--------------------------------------------------------------------
   //  Break X into two parts such that X = 2^23 * X1 + X2, compute
   //  Z = A1 * X2 + A2 * X1  (mod 2^23), and then
   //  X = 2^23 * Z + A2 * X2  (mod 2^46).
   //--------------------------------------------------------------------
   t1 = r23 * (*x);
   x1 = (int) t1;
   x2 = *x - t23 * x1;
   t1 = a1 * x2 + a2 * x1;
   t2 = (int) (r23 * t1);
   z = t1 - t23 * t2;
   t3 = t23 * z + a2 * x2;
   t4 = (int) (r46 * t3);
   *x = t3 - t46 * t4;
   r = r46 * (*x);
   
   return r;
}

void vranlc(int n, double *x, double a, double y[]) {
   //--------------------------------------------------------------------
   //
   //  This routine generates N uniform pseudorandom double precision numbers in
   //  the range (0, 1) by using the linear congruential generator
   //
   //  x_{k+1} = a x_k  (mod 2^46)
   //
   //  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
   //  before repeating.  The argument A is the same as 'a' in the above formula,
   //  and X is the same as x_0.  A and X must be odd double precision integers
   //  in the range (1, 2^46).  The N results are placed in Y and are normalized
   //  to be between 0 and 1.  X is updated to contain the new seed, so that
   //  subsequent calls to VRANLC using the same arguments will generate a
   //  continuous sequence.  If N is zero, only initialization is performed, and
   //  the variables X, A and Y are ignored.
   //
   //  This routine is the standard version designed for scalar or RISC systems.
   //  However, it should produce the same results on any single processor
   //  computer with at least 48 mantissa bits in double precision floating point
   //  data.  On 64 bit systems, double precision should be disabled.
   //
   //--------------------------------------------------------------------
   // r23 = pow(0.5, 23.0);
   ////  pow(0.5, 23.0) = 1.1920928955078125e-07
   // r46 = r23 * r23;
   // t23 = pow(2.0, 23.0);
   ////  pow(2.0, 23.0) = 8.388608e+06
   // t46 = t23 * t23;
   double const r23 = 1.1920928955078125e-07;
   double const r46 = r23 * r23;
   double const t23 = 8.388608e+06;
   double const t46 = t23 * t23;
   double t1, t2, t3, t4, a1, a2, x1, x2, z;
   int i;
   //--------------------------------------------------------------------
   //  Break A into two parts such that A = 2^23 * A1 + A2.
   //--------------------------------------------------------------------
   t1 = r23 * a;
   a1 = (int) t1;
   a2 = a - t23 * a1;
   //--------------------------------------------------------------------
   //  Generate N results.   This loop is not vectorizable.
   //--------------------------------------------------------------------
   /*************** Clava msgError **************
   Variable x could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(i = 0; i < n; i++) {
      //--------------------------------------------------------------------
      //  Break X into two parts such that X = 2^23 * X1 + X2, compute
      //  Z = A1 * X2 + A2 * X1  (mod 2^23), and then
      //  X = 2^23 * Z + A2 * X2  (mod 2^46).
      //--------------------------------------------------------------------
      t1 = r23 * (*x);
      x1 = (int) t1;
      x2 = *x - t23 * x1;
      t1 = a1 * x2 + a2 * x1;
      t2 = (int) (r23 * t1);
      z = t1 - t23 * t2;
      t3 = t23 * z + a2 * x2;
      t4 = (int) (r46 * t3);
      *x = t3 - t46 * t4;
      y[i] = r46 * (*x);
   }
   
   return;
}

void print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int verified) {
   char size[16];
   int j;
   printf("\n\n %s Benchmark Completed.\n", name);
   printf(" Class           =             %12c\n", class);
   // If this is not a grid-based problem (EP, FT, CG), then
   // we only print n1, which contains some measure of the
   // problem size. In that case, n2 and n3 are both zero.
   // Otherwise, we print the grid size n1xn2xn3
   if((n2 == 0) && (n3 == 0)) {
      if((name[0] == 'E') && (name[1] == 'P')) {
         sprintf(size, "%15.0lf", pow(2.0, n1));
         j = 14;
         if(size[j] == '.') {
            size[j] = ' ';
            j--;
         }
         size[j + 1] = '\0';
         printf(" Size            =          %15s\n", size);
      }
      else {
         printf(" Size            =             %12d\n", n1);
      }
   }
   else {
      printf(" Size            =           %4dx%4dx%4d\n", n1, n2, n3);
   }
   printf(" Iterations      =             %12d\n", niter);
   printf(" Time in seconds =             %12.2lf\n", t);
   printf(" Mop/s total     =          %15.2lf\n", mops);
   printf(" Operation type  = %24s\n", optype);
   if(verified) printf(" Verification    =             %12s\n", "SUCCESSFUL");
   else printf(" Verification    =             %12s\n", "UNSUCCESSFUL");
}

void wtime(double *t) {
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (void *) 0);
   if(sec < 0) sec = tv.tv_sec;
   *t = (tv.tv_sec - sec) + 1.0e-6 * tv.tv_usec;
}

/*****************************************************************/
/******         E  L  A  P  S  E  D  _  T  I  M  E          ******/
/*****************************************************************/
double elapsed_time() {
   double t;
   wtime(&t);
   
   return (t);
}

/*****************************************************************/
/******            T  I  M  E  R  _  C  L  E  A  R          ******/
/*****************************************************************/
void timer_clear(int n) {
   elapsed[n] = 0.0;
}

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  A  R  T          ******/
/*****************************************************************/
void timer_start(int n) {
   start[n] = elapsed_time();
}

/*****************************************************************/
/******            T  I  M  E  R  _  S  T  O  P             ******/
/*****************************************************************/
void timer_stop(int n) {
   double t, now;
   now = elapsed_time();
   t = now - start[n];
   elapsed[n] += t;
}

/*****************************************************************/
/******            T  I  M  E  R  _  R  E  A  D             ******/
/*****************************************************************/
double timer_read(int n) {
   
   return (elapsed[n]);
}
