package kadabra.utils;


public class RangeMonitors {
    public static double[] monitorMax = new double[ 1 ];

    public static double[] monitorMin = new double[ 1 ];

    static {
        kadabra.utils.RangeMonitors.monitor_range_init();
    }

    public static void monitor_range_init() {
        
        for(int i=0; i < 1; i++) {
            
            kadabra.utils.RangeMonitors.monitorMin[i] = Double.POSITIVE_INFINITY;
            kadabra.utils.RangeMonitors.monitorMax[i] = Double.NEGATIVE_INFINITY;
        }
    }

    public static void monitor_range_update(int id, double value) {
        if(value < monitorMin[id]) monitorMin[id] = value;
        if(value > monitorMax[id]) monitorMax[id] = value;
    }

    public static void monitor_range_print() {
        System.out.printf("foo\n");
        System.out.printf("a: {%f, %f}\n", kadabra.utils.RangeMonitors.monitorMin[0], kadabra.utils.RangeMonitors.monitorMax[0]);
    }
}

