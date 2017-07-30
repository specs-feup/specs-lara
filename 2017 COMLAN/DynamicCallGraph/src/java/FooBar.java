public class FooBar {

    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;

        for (int i = 0; i < 5; i++) {
            a += bar();
        }

        return a;
    }

    public static void main(String[] args) {
       double foo = new FooBar().foo();

       System.out.println(foo);
    }
}
