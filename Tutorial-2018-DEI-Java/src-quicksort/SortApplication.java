import java.util.Random;

import sort.Quicksort;

public class SortApplication {
    private static final int NUM = 1000000;

    public static void main(String[] args) {
        Random random = new Random();
        int[] values = new int[NUM];

        for (int i = 0; i < NUM; i++) {
            values[i] = random.nextInt();
        }

        Quicksort sorter = new Quicksort();
        sorter.sort(values); // sorting unsorted array
        sorter.sort(values); // sorting already sorted array
    }
}