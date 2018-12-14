#include "lib.h"

double bar() {
    return 1.0;
}

double foo() {
    double a = 0;
    
    for(int i=0; i<1000; i++) {
        a += bar();
    }
    
    return a;
}

int main() {
    
    foo();
    int a = 0;
    if(a == 0) {
        a++;
    } else {
        a--;
    }
    

}