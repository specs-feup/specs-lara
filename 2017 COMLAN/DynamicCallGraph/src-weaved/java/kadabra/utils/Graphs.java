package kadabra.utils;


public class Graphs {
    public static int[] dcg_CallGraph = new int[ 3 ];

    public static void print_dcg_CallGraph(){
        System.out.println("digraph DCG {");
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
            if (kadabra.utils.Graphs.dcg_CallGraph[ 0 ] != 0) {
                System.out.printf("\tFooBar_foo -> FooBar_bar [label=\"%d\"];\n", kadabra.utils.Graphs.dcg_CallGraph[0]);
            }	if (kadabra.utils.Graphs.dcg_CallGraph[ 1 ] != 0) {
                System.out.printf("\tFooBar_main -> FooBar_foo [label=\"%d\"];\n", kadabra.utils.Graphs.dcg_CallGraph[1]);
            }	if (kadabra.utils.Graphs.dcg_CallGraph[ 2 ] != 0) {
                System.out.printf("\tFooBar_main -> PrintStream_println [label=\"%d\"];\n", kadabra.utils.Graphs.dcg_CallGraph[2]);
            }
            System.out.println("\tSource [shape=rectangle,color=blue];");
            System.out.println("\tApi [shape=rectangle,color=gray];");
            System.out.println("}");
        
    }
}

