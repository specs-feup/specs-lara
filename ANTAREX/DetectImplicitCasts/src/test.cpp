double foo() {
    return 1.0;
}

int main() {
    
    float a = 0;
    
    for(int i=0; i<1000; i++) {
        a += foo();
    }
}