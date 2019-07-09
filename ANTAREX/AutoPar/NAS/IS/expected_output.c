#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
//---------------------------------------------------------------------

// program IS

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

/*************************************/

/*Typedef: if necessary, change the*/

/*size of int here by changing the*/

/*int type to, say, long*/

/*************************************/

typedef int INT_TYPE;

struct anon_NAS_IS_c_97 {
   double real;
   double imag;
};

typedef struct anon_NAS_IS_c_97 dcomplex;
/********************/

/*Some global info*/

/********************/

INT_TYPE *key_buff_ptr_global;
/*used by full_verify to get*/

/*copies of rank info*/

int passed_verification;
/************************************/

/*These are the three main arrays.*/

/*See SIZE_OF_BUFFERS def above*/

/************************************/

INT_TYPE key_array[1048576];
INT_TYPE key_buff1[65536];
INT_TYPE key_buff2[1048576];
INT_TYPE partial_verify_vals[5];
/**********************/

/*Partial verif info*/

/**********************/

INT_TYPE test_index_array[5];
INT_TYPE test_rank_array[5];
INT_TYPE S_test_index_array[5] = {48427, 17148, 23627, 62548, 4431};
INT_TYPE S_test_rank_array[5] = {0, 18, 346, 64917, 65463};
INT_TYPE W_test_index_array[5] = {357773, 934767, 875723, 898999, 404505};
INT_TYPE W_test_rank_array[5] = {1249, 11698, 1039987, 1043896, 1048018};
INT_TYPE A_test_index_array[5] = {2112377, 662041, 5336171, 3642833, 4250760};
INT_TYPE A_test_rank_array[5] = {104, 17523, 123928, 8288932, 8388264};
INT_TYPE B_test_index_array[5] = {41869, 812306, 5102857, 18232239, 26860214};
INT_TYPE B_test_rank_array[5] = {33422937, 10244, 59149, 33135281, 99};
INT_TYPE C_test_index_array[5] = {44172927, 72999161, 74326391, 129606274, 21736814};
INT_TYPE C_test_rank_array[5] = {61147, 882988, 266290, 133997595, 133525895};
INT_TYPE D_test_index_array[5] = {1317351170, 995930646, 1157283250, 1503301535, 1453734525};
INT_TYPE D_test_rank_array[5] = {1, 36538729, 1978098519, 2145192618, 2147425337};
/***********************/
/*function prototypes*/
/***********************/
double randlc(double *X, double *A);
void full_verify();
void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int passed_verification);
double start[64];
double elapsed[64];
double elapsed_time();
void timer_clear(int n);
void timer_start(int n);
void timer_stop(int n);
double timer_read(int n);
void wtime(double *t);
/*****************************************************************/
/*************           R  A  N  D  L  C             ************/
/*************                                        ************/
/*************    portable random number generator    ************/
/*****************************************************************/
double randlc(double *X, double *A) {
   int KS = 0;
   double R23, R46, T23, T46;
   double T1, T2, T3, T4;
   double A1;
   double A2;
   double X1;
   double X2;
   double Z;
   int i, j;
   if(KS == 0) {
      R23 = 1.0;
      R46 = 1.0;
      T23 = 1.0;
      T46 = 1.0;
      /*************** Clava msgError **************
      Loop Iteration number is too low
      ****************************************/
      for(i = 1; i <= 23; i++) {
         R23 = 0.50 * R23;
         T23 = 2.0 * T23;
      }
      /*************** Clava msgError **************
      Loop Iteration number is too low
      ****************************************/
      for(i = 1; i <= 46; i++) {
         R46 = 0.50 * R46;
         T46 = 2.0 * T46;
      }
      KS = 1;
   }
   /*Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.*/
   T1 = R23 * *A;
   j = T1;
   A1 = j;
   A2 = *A - T23 * A1;
   /*Break X into two parts such that X = 2^23 * X1 + X2, compute
   Z = A1 * X2 + A2 * X1  (mod 2^23), and then
   X = 2^23 * Z + A2 * X2  (mod 2^46).*/
   T1 = R23 * *X;
   j = T1;
   X1 = j;
   X2 = *X - T23 * X1;
   T1 = A1 * X2 + A2 * X1;
   j = R23 * T1;
   T2 = j;
   Z = T1 - T23 * T2;
   T3 = T23 * Z + A2 * X2;
   j = R46 * T3;
   T4 = j;
   *X = T3 - T46 * T4;
   
   return (R46 * *X);
}

/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/
void create_seq(double seed, double a) {
   double x;
   INT_TYPE i, k;
   k = (1 << 16) / 4;
   /*************** Clava msgError **************
   Variable seed could not be categorized into any OpenMP Variable Scopeuse : RWR
   ****************************************/
   for(i = 0; i < (1 << 20); i++) {
      x = randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);
      x += randlc(&seed, &a);
      key_array[i] = k * x;
   }
}

/*****************************************************************/
/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
/*****************************************************************/
void full_verify() {
   INT_TYPE i, j;
   /*Now, finally, sort the keys:*/
   /*key_buff2[] already has the proper information, so do nothing*/
   /*Copy keys into work array; keys in key_array will be reassigned.*/
   #pragma omp parallel for default(shared) private(i) firstprivate(key_array)
   for(i = 0; i < (1 << 20); i++) key_buff2[i] = key_array[i];
   /*************** Clava msgError **************
   Array access key_array[--key_buff_ptr_global[key_buff2[i]]] which is used for writing has subscript of arrayType --key_buff_ptr_global[key_buff2[i]]
   ****************************************/
   for(i = 0; i < (1 << 20); i++) key_array[--key_buff_ptr_global[key_buff2[i]]] = key_buff2[i];
   /*Confirm keys correctly sorted: count incorrectly sorted keys, if any*/
   j = 0;
   #pragma omp parallel for default(shared) private(i) firstprivate(key_array) reduction(+ : j)
   for(i = 1; i < (1 << 20); i++) if(key_array[i - 1] > key_array[i]) j++;
   if(j != 0) {
      printf("Full_verify: number of keys out of sort: %ld\n", (long) j);
   }
   else passed_verification++;
}

/*****************************************************************/
/*************             R  A  N  K             ****************/
/*****************************************************************/
void rank(int iteration) {
   INT_TYPE i, k;
   INT_TYPE *key_buff_ptr, *key_buff_ptr2;
   key_array[iteration] = iteration;
   key_array[iteration + 10] = (1 << 16) - iteration;
   /*Determine where the partial verify test keys are, load into*/
   /*top of array bucket_size*/
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(i = 0; i < 5; i++) partial_verify_vals[i] = key_array[test_index_array[i]];
   /*Initialize*/
   /*Determine the number of keys in each bucket*/
   /*Accumulative bucket sizes are the bucket pointers*/
   /*Sort into appropriate bucket*/
   key_buff_ptr2 = key_array;
   /*Clear the work array*/
   #pragma omp parallel for default(shared) private(i)
   for(i = 0; i < (1 << 16); i++) key_buff1[i] = 0;
   /*Ranking of all keys occurs in this section:*/
   key_buff_ptr = key_buff1;
   /*In this section, the keys themselves are used as their
   own indexes to determine how many of each there are: their
   individual population*/
   /*************** Clava msgError **************
   Array access key_buff_ptr[key_buff_ptr2[i]] which is used for writing has subscript of arrayType key_buff_ptr2[i]
   ****************************************/
   for(i = 0; i < (1 << 20); i++) key_buff_ptr[key_buff_ptr2[i]]++;
   /*Now they have individual key*/
   /*population*/
   /*To obtain ranks of each key, successively add the individual key
   population*/
   /*************** Clava msgError **************
   unsolved dependency for arrayAccess key_buff_ptr	 use : RWR
   ****************************************/
   for(i = 0; i < (1 << 16) - 1; i++) key_buff_ptr[i + 1] += key_buff_ptr[i];
   /*This is the partial verify test section*/
   /*Observe that test_rank_array vals are*/
   /*shifted differently for different cases*/
   /*************** Clava msgError **************
   Loop contains Invalid Statement -> BreakStmt#486
   ****************************************/
   for(i = 0; i < 5; i++) {
      k = partial_verify_vals[i];
      /*test vals were put here*/
      if(0 < k && k <= (1 << 20) - 1) {
         INT_TYPE key_rank = key_buff_ptr[k - 1];
         int failed = 0;
         switch ('W') {
            case 'S':
            if(i <= 2) {
               if(key_rank != test_rank_array[i] + iteration) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - iteration) failed = 1;
               else passed_verification++;
            }
            break;
            case 'W':
            if(i < 2) {
               if(key_rank != test_rank_array[i] + (iteration - 2)) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - iteration) failed = 1;
               else passed_verification++;
            }
            break;
            case 'A':
            if(i <= 2) {
               if(key_rank != test_rank_array[i] + (iteration - 1)) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - (iteration - 1)) failed = 1;
               else passed_verification++;
            }
            break;
            case 'B':
            if(i == 1 || i == 2 || i == 4) {
               if(key_rank != test_rank_array[i] + iteration) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - iteration) failed = 1;
               else passed_verification++;
            }
            break;
            case 'C':
            if(i <= 2) {
               if(key_rank != test_rank_array[i] + iteration) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - iteration) failed = 1;
               else passed_verification++;
            }
            break;
            case 'D':
            if(i < 2) {
               if(key_rank != test_rank_array[i] + iteration) failed = 1;
               else passed_verification++;
            }
            else {
               if(key_rank != test_rank_array[i] - iteration) failed = 1;
               else passed_verification++;
            }
            break;
         }
         if(failed == 1) printf("Failed partial verification: iteration %d, test key %d\n", iteration, (int) i);
      }
   }
   /*Make copies of rank info for use by full_verify: these variables
   in rank are local; making them global slows down the code, probably
   since they cannot be made register by compiler*/
   if(iteration == 10) key_buff_ptr_global = key_buff_ptr;
}

/*****************************************************************/
/*************             M  A  I  N             ****************/
/*****************************************************************/
int main(int argc, char **argv) {
   int i, iteration;
   double timecounter;
   FILE *fp;
   /*Initialize timers*/
   timer_clear(0);
   /*Initialize the verification arrays if a valid class*/
   /*************** Clava msgError **************
   Loop contains Invalid Statement -> BreakStmt#608
   ****************************************/
   for(i = 0; i < 5; i++) switch ('W') {
      case 'S':
      test_index_array[i] = S_test_index_array[i];
      test_rank_array[i] = S_test_rank_array[i];
      break;
      case 'A':
      test_index_array[i] = A_test_index_array[i];
      test_rank_array[i] = A_test_rank_array[i];
      break;
      case 'W':
      test_index_array[i] = W_test_index_array[i];
      test_rank_array[i] = W_test_rank_array[i];
      break;
      case 'B':
      test_index_array[i] = B_test_index_array[i];
      test_rank_array[i] = B_test_rank_array[i];
      break;
      case 'C':
      test_index_array[i] = C_test_index_array[i];
      test_rank_array[i] = C_test_rank_array[i];
      break;
      case 'D':
      test_index_array[i] = D_test_index_array[i];
      test_rank_array[i] = D_test_rank_array[i];
      break;
   }
   ;
   /*Printout initial NPB info*/
   printf("\n\n NAS Parallel Benchmarks (NPB3.3-SER) - IS Benchmark\n\n");
   printf(" Size:  %ld  (class %c)\n", (long) (1 << 20), 'W');
   printf(" Iterations:   %d\n", 10);
   /*Generate random number sequence and subsequent keys on all procs*/
   create_seq(314159265.00, 1220703125.00);
   /*Random number gen seed*/
   /*Random number gen mult*/
   /*Do one interation for free (i.e., untimed) to guarantee initialization of
   all data and code pages and respective tables*/
   rank(1);
   /*Start verification counter*/
   passed_verification = 0;
   if('W' != 'S') printf("\n   iteration\n");
   /*Start timer*/
   timer_start(0);
   /*This is the main iteration*/
   /*************** Clava msgError **************
   Loop Iteration number is too low
   ****************************************/
   for(iteration = 1; iteration <= 10; iteration++) {
      if('W' != 'S') printf("        %d\n", iteration);
      rank(iteration);
   }
   /*End of timing, obtain maximum time of all processors*/
   timer_stop(0);
   timecounter = timer_read(0);
   /*This tests that keys are in sequence: sorting of last ranked key seq
   occurs here, but is an untimed operation*/
   full_verify();
   /*The final printout*/
   if(passed_verification != 5 * 10 + 1) passed_verification = 0;
   c_print_results("IS", 'W', (int) ((1 << 20) / 64), 64, 0, 10, timecounter, ((double) (10 * (1 << 20))) / timecounter / 1000000., "keys ranked", passed_verification);
   int exitValue = passed_verification ? 0 : 1;
   
   return exitValue;
}

/**************************/
/*E N D  P R O G R A M*/
/**************************/
void c_print_results(char *name, char class, int n1, int n2, int n3, int niter, double t, double mops, char *optype, int passed_verification) {
   printf("\n\n %s Benchmark Completed\n", name);
   printf(" Class           =                        %c\n", class);
   if(n3 == 0) {
      long nn = n1;
      if(n2 != 0) nn *= n2;
      printf(" Size            =             %12ld\n", nn);
   }
   else printf(" Size            =             %4dx%4dx%4d\n", n1, n2, n3);
   printf(" Iterations      =             %12d\n", niter);
   printf(" Time in seconds =             %12.2f\n", t);
   printf(" Mop/s total     =             %12.2f\n", mops);
   printf(" Operation type  = %24s\n", optype);
   if(passed_verification < 0) printf(" Verification    =            NOT PERFORMED\n");
   else if(passed_verification) printf(" Verification    =               SUCCESSFUL\n");
   else printf(" Verification    =             UNSUCCESSFUL\n");
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
