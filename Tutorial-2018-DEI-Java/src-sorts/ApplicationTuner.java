import java.util.Arrays;
import java.util.Random;

import algorithms.Quicksort;

public class ApplicationTuner {
    private static final int[] SIZES = { 5, 20, 50, 100, 200, 500 };

    private static final int NUM_EXECS = 200;

    private static final int MAX_INT = 512;

    public static void main(String[] args) {

        Random random = new Random();
        for (int size : SIZES) {
            int[] values = new int[size];
            for (int it = 0; it < NUM_EXECS; it++) {
                for (int i = 0; i < size; i++) {
                    values[i] = random.nextInt(MAX_INT);
                }
                System.out.println("Sorting Array #" + (it + 1) + ": " + Arrays.toString(values));
                Quicksort.sort(values); // sorting unsorted array
            }
        }
    }
}