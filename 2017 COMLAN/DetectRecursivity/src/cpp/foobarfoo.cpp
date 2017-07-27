double foo(int seed);

double bar(int a) {
    if(a > 0) {
        return a;
    }
    
    return foo(a+1);
}

double foo(int seed) {
    double a = 0;
    
    for(int i=seed; i<1000; i++) {
        a += bar(i);
    }


    return a;
}

int main() {
    
    foo(0);

}