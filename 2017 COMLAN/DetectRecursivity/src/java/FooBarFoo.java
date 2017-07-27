public class FooBarFoo {

    public double bar(int a) {
        if(a > 0) {
			return a;
		}
		
		return foo(a+1);
    }

    public double foo(int seed) {
        double a = 0;

        for (int i = seed; i < 5; i++) {
            a += bar(i);
        }

        return a;
    }

    public static void main(String[] args) {
        new FooBarFoo().foo(0);
    }
}
