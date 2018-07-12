import java.util.Random;

import matrix.MMult;

public class ApplicationMxM {
    private static final int NUM_EXECS = 100;
    private static final int LENGTH = 512;

    public static void main(String[] args) {
        for (int i = 0; i < NUM_EXECS; i++) {
            System.out.println("Multiplication #" + (i + 1));
            int[][] a = genMatrix(LENGTH, LENGTH);
            int[][] b = genMatrix(LENGTH, LENGTH);
            MMult.mult(a, b);
        }
    }

    public static int[][] genMatrix(int numRows, int numCols) {

        int[][] m = new int[numRows][numCols];

        Random rand = new Random();

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {

                m[i][j] = rand.nextInt();
            }
        }
        return m;
    }

}