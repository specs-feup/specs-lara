


public class FooBar {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        kadabra.utils.RangeMonitors.monitor_range_update(0, a);
        for (int i = 0; i < 5; i++) {
            a += bar();
            kadabra.utils.RangeMonitors.monitor_range_update(0, a);
        }
        return a;
    }

    public static void main(String[] args) {
        new FooBar().foo();
        kadabra.utils.RangeMonitors.monitor_range_print();
    }
}

