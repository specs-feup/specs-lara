package algorithms;

public class Quicksort {

    public static void sort(int[] values) {
        if (values == null || values.length == 0) {
            return;
        }
        quicksort(values, 0, values.length - 1);
    }

    private static void quicksort(int[] numbers, int low, int high) {
        int i = low, j = high;
        int pivot = numbers[low + (high - low) / 2];

        while (i <= j) {
            while (numbers[i] < pivot) {
                i++;
            }
            while (numbers[j] > pivot) {
                j--;
            }

            if (i <= j) {
                exchange(numbers, i, j);
                i++;
                j--;
            }
        }
        if (low < j) {
            quicksort(numbers, low, j);
        }
        if (i < high) {
            quicksort(numbers, i, high);
        }
    }

    private static void exchange(int[] numbers, int i, int j) {
        int temp = numbers[i];
        numbers[i] = numbers[j];
        numbers[j] = temp;
    }
}
