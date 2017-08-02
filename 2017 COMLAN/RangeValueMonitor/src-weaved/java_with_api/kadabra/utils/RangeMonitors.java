package kadabra.utils;


public class RangeMonitors {
    public static double[] monitor1Max = new double[ 1 ];

    public static double[] monitor1Min = new double[ 1 ];

    static {
        kadabra.utils.RangeMonitors.monitor1_range_init();
    }

    public static void monitor1_range_init() {
        
        for(int i=0; i < 1; i++) {
            
            kadabra.utils.RangeMonitors.monitor1Min[i] = Double.POSITIVE_INFINITY;
            kadabra.utils.RangeMonitors.monitor1Max[i] = Double.NEGATIVE_INFINITY;
        }
    }

    public static void monitor1_range_update(int id, double value) {
        if(value < monitor1Min[id]) monitor1Min[id] = value;
        if(value > monitor1Max[id]) monitor1Max[id] = value;
    }

    public static void monitor1_range_print() {
        System.out.printf("foo\n");
        System.out.printf("a: {%f, %f}\n", kadabra.utils.RangeMonitors.monitor1Min[0], kadabra.utils.RangeMonitors.monitor1Max[0]);
    }
}

