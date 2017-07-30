#include "lib/lib.h"

double bar() {
    return 1.0;
}

double foo() {
    double a = 0;
    
    for(int i=0; i<1000; i++) {
        a += bar();
    }

	a = lib_call(a);
    
    return a;
}

int main() {
    
    foo();

}