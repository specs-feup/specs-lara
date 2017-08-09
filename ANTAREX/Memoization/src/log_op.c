#include <math.h>

#define SIZE 100
#define C 45.0f

/* to remove */
int log_wrapper(int n) {

	int result;

	result = log(n);
	
	return result;
}
/* to remove */

unsigned char logarithm( const unsigned char value, float c )
{

    return c * log( 1 + value );
}

int main()
{

    unsigned char const input[SIZE] = {
        230, 122, 63, 217, 230, 56, 204, 91, 167, 137, 53, 132, 165, 29, 251, 13, 181, 195, 147, 245, 242, 8, 173, 29, 141, 160, 191, 28, 140, 149, 177, 143, 166, 76, 207, 81, 0, 170, 91, 193, 156, 5, 24, 251, 188, 207, 53, 163, 157, 129, 69, 90, 0, 107, 251, 136, 7, 118, 173, 126, 42, 96, 19, 157, 45, 0, 111, 26, 64, 77, 119, 226, 196, 18, 112, 193, 176, 80, 164, 124, 55, 8, 78, 38, 151, 249, 246, 79, 122, 253, 141, 146, 88, 13, 164, 250, 25, 65, 202, 178
    };

    unsigned char output[SIZE];

    unsigned int i;
    for( i = 0; i < SIZE; i++ ) {

        output[i] = logarithm( input[i], C );
    }

    return 0;
}
