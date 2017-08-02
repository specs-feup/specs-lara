

import pt.up.fe.specs.lara.kadabra.utils.CountingMonitorList;

public class FooBar {
    private static CountingMonitorList kadabraCountingMonitorList0 = new CountingMonitorList(3);

    public double bar() {
        return 1.0;
    }

    public double foo() {
        double a = 0;
        for (int i = 0; i < 5; i++) {
            kadabraCountingMonitorList0.increment(0);
            a += bar();
        }
        return a;
    }

    public static void main(String[] args) {
        kadabraCountingMonitorList0.increment(1);
        double foo = new FooBar().foo();
        kadabraCountingMonitorList0.increment(2);
        System.out.println(foo);
        FooBar.printCallGraph()
    }

    public static void printCallGraph() {
        int value;
        System.out.println("digraph DCG {\n");
            System.out.println("\tsubgraph cluster_FooBar {\n");
                System.out.println("\t\tlabel=\"FooBar\";\n");
                System.out.println("\t\tcolor=blue;\n");
                System.out.println("\t\tFooBar_bar [label=\"bar\"]\n");
                System.out.println("\t\tFooBar_foo [label=\"foo\"]\n");
                System.out.println("\t\tFooBar_main [label=\"main\"]\n");
                System.out.println("\t}\n");
            System.out.println("\tsubgraph cluster_PrintStream {\n");
                System.out.println("\t\tlabel=\"PrintStream\";\n");
                System.out.println("\t\tcolor=gray;\n");
                System.out.println("\t\tPrintStream_println [label=\"println\"]\n");
                System.out.println("\t}\n");
            value = kadabraCountingMonitorList0.get(0);
            System.out.println("FooBar_foo->FooBar_bar [label=\""+value+"\"]\n");
            value = kadabraCountingMonitorList0.get(1);
            System.out.println("FooBar_main->FooBar_foo [label=\""+value+"\"]\n");
            value = kadabraCountingMonitorList0.get(2);
            System.out.println("FooBar_main->PrintStream_println [label=\""+value+"\"]\n");
            System.out.println("}");
        }
    }

    