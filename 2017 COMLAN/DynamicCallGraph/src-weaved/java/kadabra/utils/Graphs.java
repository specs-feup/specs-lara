package kadabra.utils;


public class Graphs {
    public static int[] dcg_CallGraph = new int[ 2 ];

    public static void print_dcg_CallGraph(){
        System.out.println("digraph CallGraph {");
            if (dcg_CallGraph[ 0 ] != 0) {
                System.out.printf("\tFooBar_foo -> FooBar_bar [label=\"%d\"];\n", dcg_CallGraph[0]);
            }	if (dcg_CallGraph[ 1 ] != 0) {
                System.out.printf("\tFooBar_main -> FooBar_foo [label=\"%d\"];\n", dcg_CallGraph[1]);
            }
            System.out.println("}");
    }
}

