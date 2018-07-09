package algorithms;

public class CountingSort {

    public static void sort(int[] values) {

        // this code is assuming positive values
        int dimension = Integer.MIN_VALUE;
        int minValue = Integer.MAX_VALUE;
        /**
         * Get maximum and minimum value
         */
        for (int value : values) {
            if (value > dimension) {
                dimension = value;
            }
            if (value < minValue) {
                minValue = value;
            }
        }
        dimension -= minValue - 1;
        // System.out.println(minValue + "->" + dimension);
        int[] sortingTable = new int[dimension];
        for (int i = 0; i < values.length; i++) {

            sortingTable[values[i] - minValue]++;
        }

        int pos = 0;
        for (int i = 0; i < dimension; i++) {
            int repeat = sortingTable[i];
            while (repeat > 0) {
                values[pos++] = i + minValue;
                repeat--;
            }
        }
    }

}
