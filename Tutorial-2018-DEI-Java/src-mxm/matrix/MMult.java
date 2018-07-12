package matrix;

/**
 * Copyright 2014 Tiago D. R. Carvalho.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on
 * an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the
 * specific language governing permissions and limitations under the License. under the License.
 */

public class MMult {

    /**
     * Matrix multiplication
     * 
     * @param a
     *            input matrix A
     * @param b
     *            input matrix B
     * @return a vector of the operation: a*b
     */
    public static int[][] mult(int[][] a, int[][] b) {
        int rowsInA = a.length;
        int columnsInA = a[0].length; // same as rows in B
        int columnsInB = b[0].length;
        int[][] c = new int[rowsInA][columnsInB];

        for (int i = 0; i < rowsInA; i++) {
            for (int j = 0; j < columnsInB; j++) {

                for (int k = 0; k < columnsInA; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return c;

    }

}
