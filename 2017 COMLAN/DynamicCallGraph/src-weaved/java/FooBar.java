


public class FooBar {
    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            a += bar();
            kadabra.utils.Graphs.dcg_CallGraph[0]++;
        }
        return a;
    }

    public static void main(String[] args) {
        double foo = new FooBar().foo();
        kadabra.utils.Graphs.dcg_CallGraph[1]++;
        System.out.println(foo);
        kadabra.utils.Graphs.dcg_CallGraph[2]++;
        kadabra.utils.Graphs.print_dcg_CallGraph();
    }
}

