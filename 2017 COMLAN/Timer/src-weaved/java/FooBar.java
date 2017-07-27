


public class FooBar {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            long kadabra_timing_start_0 = System.nanoTime();
            a += bar();
            double kadabra_timing_interval_0 = (double)(System.nanoTime() - kadabra_timing_start_0) /  (double)1000;
            System.out.printf("Time:%fus\n", kadabra_timing_interval_0);
        }
        return a;
    }

    public static void main(String[] args) {
        long kadabra_timing_start_1 = System.nanoTime();
        new FooBar().foo();
        double kadabra_timing_interval_1 = (double)(System.nanoTime() - kadabra_timing_start_1) /  (double)1000;
        System.out.printf("Time:%fus\n", kadabra_timing_interval_1);
    }
}

