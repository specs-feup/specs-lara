#include <math.h>
#include <stdio.h>

#define SIZE 100
#define C 45.0

double logarithm( const double value, double c )
{

    return c * log( 1 + value );
}

int main()
{

    double const input[SIZE] = {
        230.0, 122.0, 63.0, 217.0, 230.0, 56.0, 204.0, 91.0, 167.0, 137.0, 53.0, 132.0, 165.0, 29.0, 251.0, 13.0, 181.0, 195.0, 147.0, 245.0, 242.0, 8.0, 173.0, 29.0, 141.0, 160.0, 191.0, 28.0, 140.0, 149.0, 177.0, 143.0, 166.0, 76.0, 207.0, 81.0, 0.0, 170.0, 91.0, 193.0, 156.0, 5.0, 24.0, 251.0, 188.0, 207.0, 53.0, 163.0, 157.0, 129.0, 69.0, 90.0, 0.0, 107.0, 251.0, 136.0, 7.0, 118.0, 173.0, 126.0, 42.0, 96.0, 19.0, 157.0, 45.0, 0.0, 111.0, 26.0, 64.0, 77.0, 119.0, 226.0, 196.0, 18.0, 112.0, 193.0, 176.0, 80.0, 164.0, 124.0, 55.0, 8.0, 78.0, 38.0, 151.0, 249.0, 246.0, 79.0, 122.0, 253.0, 141.0, 146.0, 88.0, 13.0, 164.0, 250.0, 25.0, 65.0, 202.0, 178.0
    };

    double output[SIZE];
    double sum;

    unsigned int i;
    for( i = 0; i < SIZE; i++ ) {

        output[i] = logarithm( input[i], C );
        sum += output[i];
    }

	printf("%f\n", sum);

    return 0;
}
