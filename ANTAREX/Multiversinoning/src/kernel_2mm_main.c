#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#define N_EXEC 1

void kernel_2mm( int ni, int nj, int nk, int nl,
                 double alpha,
                 double beta,
                 double tmp[ ni + 0][nj + 0],
                 double A[ ni + 0][nk + 0],
                 double B[ nk + 0][nj + 0],
                 double C[ nj + 0][nl + 0],
                 double D[ ni + 0][nl + 0] );


static
void *
xmalloc ( size_t num )
{
    void* cur = NULL;
    int ret = posix_memalign ( &cur, 32, num );
    if ( ! cur || ret ) {
        fprintf ( stderr, "[PolyBench] posix_memalign: cannot allocate memory" );
        exit ( 1 );
    }
    return cur;
}


void* polybench_alloc_data( unsigned long long int n, int elt_size )
{
    /// FIXME: detect overflow!
    size_t val = n;
    void* ret;
    val *= elt_size;
    ret = xmalloc ( val );

    return ret;
}


static
void init_array( int ni, int nj, int nk, int nl,
                 double *alpha,
                 double *beta,
                 double A[ ni + 0][nk + 0],
                 double B[ nk + 0][nj + 0],
                 double C[ nj + 0][nl + 0],
                 double D[ ni + 0][nl + 0] )
{
    int i, j;

    *alpha = 1.5;
    *beta = 1.2;
    for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nk; j++ ) {
            A[i][j] = ( double ) ( i * j % ni ) / ni;
        }
    for ( i = 0; i < nk; i++ )
        for ( j = 0; j < nj; j++ ) {
            B[i][j] = ( double ) ( i * ( j + 1 ) % nj ) / nj;
        }
    for ( i = 0; i < nj; i++ )
        for ( j = 0; j < nl; j++ ) {
            C[i][j] = ( double ) ( i * ( j + 3 ) % nl ) / nl;
        }
    for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nl; j++ ) {
            D[i][j] = ( double ) ( i * ( j + 2 ) % nk ) / nk;
        }
}


static
void print_array( int ni, int nl,
                  double D[ ni + 0][nl + 0] )
{
    int i, j;

    fprintf( stderr, "==BEGIN DUMP_ARRAYS==\n" );
    fprintf( stderr, "begin dump: %s", "D" );
    for ( i = 0; i < ni; i++ )
        for ( j = 0; j < nl; j++ ) {
            if ( ( i * ni + j ) % 20 == 0 ) {
                fprintf ( stderr, "\n" );
            }
            fprintf ( stderr, "%0.6lf ", D[i][j] );
        }
    fprintf( stderr, "\nend   dump: %s\n", "D" );
    fprintf( stderr, "==END   DUMP_ARRAYS==\n" );
}


void main( int argc, char** argv )
{

    int ni = 400;
    int nj = 450;
    int nk = 550;
    int nl = 600;


    double alpha;
    double beta;
    int n_exec_i;
    double ( *tmp )[ni + 0][nj + 0];
    double ( *A )[ni + 0][nk + 0];
    double ( *B )[nk + 0][nj + 0];
    double ( *C )[nj + 0][nl + 0];
    double ( *D )[ni + 0][nl + 0];

    tmp = ( double( * )[ni + 0][nj + 0] )polybench_alloc_data ( ( ni + 0 ) * ( nj + 0 ), sizeof( double ) );
    A = ( double( * )[ni + 0][nk + 0] )polybench_alloc_data ( ( ni + 0 ) * ( nk + 0 ), sizeof( double ) );
    B = ( double( * )[nk + 0][nj + 0] )polybench_alloc_data ( ( nk + 0 ) * ( nj + 0 ), sizeof( double ) );
    C = ( double( * )[nj + 0][nl + 0] )polybench_alloc_data ( ( nj + 0 ) * ( nl + 0 ), sizeof( double ) );
    D = ( double( * )[ni + 0][nl + 0] )polybench_alloc_data ( ( ni + 0 ) * ( nl + 0 ), sizeof( double ) );

    init_array ( ni, nj, nk, nl, &alpha, &beta,
                 *A,
                 *B,
                 *C,
                 *D );

    for ( n_exec_i = 0; n_exec_i < N_EXEC; n_exec_i++ ) {
        kernel_2mm ( ni, nj, nk, nl,
                     alpha, beta,
                     *tmp,
                     *A,
                     *B,
                     *C,
                     *D );
    }

    if ( argc > 42 && ! strcmp( argv[0], "" ) ) {
        print_array( ni, nl, *D );
    }

    free( ( void* )tmp );
    free( ( void* )A );
    free( ( void* )B );
    free( ( void* )C );
    free( ( void* )D );
}
