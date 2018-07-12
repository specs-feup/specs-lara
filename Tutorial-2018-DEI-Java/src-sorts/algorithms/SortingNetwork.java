package algorithms;

public class SortingNetwork {

    public static void sort(int[] valuesArray) {

        int i;
        int length = valuesArray.length;
        int halfLength = length / 2;
        for (i = 0; i < halfLength; i++) {

            /* odd */
            compareLevel(valuesArray, 0);

            /* even */
            compareLevel(valuesArray, 1);
        }

        if (length % 2 != 0) {

            compareLevel(valuesArray, 0);
        }
    }

    static void swapMax(int[] valuesArray, int index) {

        if (valuesArray[index] > valuesArray[index + 1]) {

            int temp = valuesArray[index];
            valuesArray[index] = valuesArray[index + 1];
            valuesArray[index + 1] = temp;
        }
    }

    static void compareLevel(int[] values, int start) {

        int i;
        int length = values.length;
        for (i = start; i + 1 < length; i += 2) {

            swapMax(values, i);
        }
    }

    static void compareAndSwap(int[] valuesArray, boolean[] mask) {

        int i;
        int maskLength = mask.length;
        for (i = 0; i < maskLength; i++) {

            if (mask[i]) {
                swapMax(valuesArray, i);
            }
        }
    }

}
