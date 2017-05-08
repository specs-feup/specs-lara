#include <omp.h>
#include <math.h>
#define ceild(n,d)  ceil(((double)(n))/((double)(d)))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define max(x,y)    ((x) > (y)? (x) : (y))
#define min(x,y)    ((x) < (y)? (x) : (y))

void kernel_2mm( int ni, int nj, int nk, int nl,
                 double alpha,
                 double beta,
                 double tmp[ ni + 0][nj + 0],
                 double A[ ni + 0][nk + 0],
                 double B[ nk + 0][nj + 0],
                 double C[ nj + 0][nl + 0],
                 double D[ ni + 0][nl + 0] )
{
    int i, j, k;

    int t1, t2, t3, t4, t5;
    int lb, ub, lbp, ubp, lb2, ub2;
    register int lbv, ubv;

    if ( ni >= 1 ) {
        if ( ( nj >= nl + 1 ) && ( nl >= 1 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                lbv = 0;
                ubv = nl - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    D[t2][t3] *= beta;;
                    tmp[t2][t3] = 0.0;;
                }
                lbv = nl;
                ubv = nj - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    tmp[t2][t3] = 0.0;;
                }
            }
        }
        if ( ( nj >= 1 ) && ( nj <= nl - 1 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                lbv = 0;
                ubv = nj - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    D[t2][t3] *= beta;;
                    tmp[t2][t3] = 0.0;;
                }
                lbv = nj;
                ubv = nl - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    D[t2][t3] *= beta;;
                }
            }
        }
        if ( ( nj >= 1 ) && ( nj == nl ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                lbv = 0;
                ubv = nj - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    D[t2][t3] *= beta;;
                    tmp[t2][t3] = 0.0;;
                }
            }
        }
        if ( ( nj >= 1 ) && ( nl <= 0 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                lbv = 0;
                ubv = nj - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    tmp[t2][t3] = 0.0;;
                }
            }
        }
        if ( ( nj <= 0 ) && ( nl >= 1 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                lbv = 0;
                ubv = nl - 1;
#pragma ivdep
#pragma vector always
                for ( t3 = lbv; t3 <= ubv; t3++ ) {
                    D[t2][t3] *= beta;;
                }
            }
        }
        if ( ( nj >= 1 ) && ( nk >= 1 ) && ( nl >= 1 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                for ( t3 = 0; t3 <= nj - 1; t3++ ) {
                    for ( t5 = 0; t5 <= nk - 1; t5++ ) {
                        tmp[t2][t3] += alpha * A[t2][t5] * B[t5][t3];;
                    }
                    lbv = 0;
                    ubv = nl - 1;
#pragma ivdep
#pragma vector always
                    for ( t5 = lbv; t5 <= ubv; t5++ ) {
                        D[t2][t5] += tmp[t2][t3] * C[t3][t5];;
                    }
                }
            }
        }
        if ( ( nj >= 1 ) && ( nk >= 1 ) && ( nl <= 0 ) ) {
            lbp = 0;
            ubp = ni - 1;
            
            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                for ( t3 = 0; t3 <= nj - 1; t3++ ) {
                    for ( t5 = 0; t5 <= nk - 1; t5++ ) {
                        tmp[t2][t3] += alpha * A[t2][t5] * B[t5][t3];;
                    }
                }
            }
        }
        if ( ( nj >= 1 ) && ( nk <= 0 ) && ( nl >= 1 ) ) {
            lbp = 0;
            ubp = ni - 1;

            for ( t2 = lbp; t2 <= ubp; t2++ ) {
                for ( t3 = 0; t3 <= nj - 1; t3++ ) {
                    lbv = 0;
                    ubv = nl - 1;
#pragma ivdep
#pragma vector always
                    for ( t5 = lbv; t5 <= ubv; t5++ ) {
                        D[t2][t5] += tmp[t2][t3] * C[t3][t5];;
                    }
                }
            }
        }
    }
}

