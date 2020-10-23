package linearalgebra;
import linearalgebra.Vector;

/**
 * The Matrix class stores a 2D array of doubles and provides
 * common linear algebra operations on that array, in addition
 * to producing special matrices and checking for specific conditions.
 */

public class Matrix {
   /**
    * matrix is a 2D array of doubles to hold the entries in the Matrix
    */
   private double[][] matrix;
   
   /**
    * accepts a 2D array of doubles and copies them into the matrix field.
    * @param matrix a 2D array of doubles
    */
   public Matrix(double[][] matrix) {
      this.matrix = new double[matrix.length][matrix[0].length];
      
      for (int row = 0; row < matrix.length; row++) {
         for (int col = 0; col < matrix[row].length; col++) {
            this.matrix[row][col] = matrix[row][col];
         }
      }
   }
   
   /**
    * copy constructor accepts a Matrix object and duplicates it.
    * @param m the Matrix object we want to copy
    */
   public Matrix(Matrix m) {
      this.matrix = new double[m.getNumRows()][m.getNumColumns()];
      
      for (int row = 0; row < m.matrix.length; row++) {
         for (int col = 0; col < m.matrix[row].length; col++) {
            this.matrix[row][col] = m.matrix[row][col];
         }
      }
   }

   /**
    * adds corresponding entries in the two matrices. Calling Matrix and Matrix b
    * must have the same dimensions to be compatible.
    * @param b an m-by-n Matrix object
    * @return an m-by-n Matrix object whose entries are sums of corresponding
    *         entries in the calling Matrix and b
    */   
   public Matrix add(Matrix b) {
      return Matrix.add(this, b);
   }
   
   /**
    * adds corresponding entries in the two matrices. Matrix a and Matrix b
    * must have the same dimensions to be compatible.
    * @param a an m-by-n Matrix object
    * @param b an m-by-n Matrix object
    * @return an m-by-n Matrix object whose entries are sums of corresponding
    *         entries in a and b
    */
   public static Matrix add(Matrix a, Matrix b) {
      if ( (a.getNumRows() != b.getNumRows()) ||
           (a.getNumColumns() != b.getNumColumns()) ) {
         throw new IllegalArgumentException("Matrix objects have different dimensions.");
      }
      
      double[][] entries = new double[a.getNumRows()][a.getNumColumns()];
      
      for (int row = 0; row < a.matrix.length; row++) {
         for (int col = 0; col < a.matrix[0].length; col++) {
            entries[row][col] = a.matrix[row][col] + b.matrix[row][col];
         }
      }
      
      return new Matrix(entries);
   }
   
   /**
    * returns the entry in the row-th row and col-th column of m.
    * @param row the row of the desired entry
    * @param col the column of the desired entry
    * @return the entry in the row-th row and col-th col of m
    */  
   public double getEntry(int row, int col) {
      return Matrix.getEntry(this, row, col);
   }
   
   /**
    * returns the entry in the row-th row and col-th column of m.
    * @param m a Matrix object
    * @param row the row of the desired entry
    * @param col the column of the desired entry
    * @return the entry in the row-th row and col-th col of m
    */
   public static double getEntry(Matrix m, int row, int col) {
      if (row < 0 || row > m.matrix.length) {
         throw new IllegalArgumentException("Invalid value for row.");
      }
      
      if (col < 0 || col > m.matrix[0].length) {
         throw new IllegalArgumentException("Invalid value for col.");
      }
      
      return m.matrix[row][col];  
   }
   
   /**
    * @return the number of columns in the matrix
    */
   public int getNumColumns() {
      return this.matrix[0].length;
   }

   /**
    * @return the number of rows in the matrix
    */
   public int getNumRows() {
      return this.matrix.length;
   }

   /**
    * accepts a column number (starting at 0) and returns the entries in
    * the column as a Vector object
    * @param col the number of the desired column 
    *        (starting at 0, up to this.matrix[0].length-1)
    * @return a Vector containing the entries in the desired column
    */   
   public Vector getColumn(int col) {
      return Matrix.getColumn(this, col);
   }

   /**
    * accepts a column number (starting at 0) and returns the entries in
    * the column as a Vector object
    * @param m the Matrix object we want the column from
    * @param col the number of the desired column 
    *        (starting at 0, up to this.matrix[0].length-1)
    * @return a Vector containing the entries in the desired column
    */   
   public static Vector getColumn(Matrix m, int col) {
      if (m == null) {
         throw new IllegalArgumentException("Matrix is null.");
      }
      if (col >= m.matrix[0].length || col < 0) {
         throw new IllegalArgumentException("Column does not exist in matrix.");
      } 
      
      double[] columnVector = new double[m.matrix.length];
      
      for (int row = 0; row < m.matrix.length; row++) {
         columnVector[row] = m.matrix[row][col];
      }  
      
      return new Vector(columnVector);  
   }
 
   /**
    * accepts a row number (starting at 0) and returns the entries in
    * the row as a Vector object
    * @param row the number of the desired row (starting at 0, up to this.matrix.length-1)
    * @return a Vector containing the entries in the desired row
    */
   public Vector getRow(int row) {
      return Matrix.getRow(this, row);
   }
   
   /**
    * accepts a row number (starting at 0) and returns the entries in
    * the row as a Vector object
    * @param m the Matrix object we want the row from
    * @param row the number of the desired row (starting at 0, up to this.matrix.length-1)
    * @return a Vector containing the entries in the desired row
    */
   public static Vector getRow(Matrix m, int row) {
      if (m == null) {
         throw new IllegalArgumentException("Matrix is null.");
      }
      if (row >= m.matrix.length || row < 0) {
         throw new IllegalArgumentException("Row does not exist in matrix.");
      }
      
      return new Vector(m.matrix[row]);
   }
   
   /**   
    * getIdentityMatrix returns an n by n matrix with all zeroes
    * and ones on the diagonal.
    * @param n an integer greater than or equal to 1
    * @return an n-by-n matrix with ones on the diagonal and zeroes
    *         everywhere else
    */
   public static Matrix identityMatrix(int n) {
      if (n < 1) {
         throw new IllegalArgumentException("n must be >= 1");
      }
      
      double[][] entries = new double[n][n];
      
      for (int i = 0; i < n; i++) {
         entries[i][i] = 1.0;
      }
      
      return new Matrix(entries);
   }  

   /**
    * multiplies each entry in the Matrix m by a real number.
    * @param x a real number (double)
    * @return a Matrix whose entries are the entries in the calling matrix, 
    *         multiplied by x
    */   
   public Matrix multiply(double x) {
      return Matrix.multiply(this, x);
   }
   
   /**
    * multiplies each entry in the Matrix m by a real number.
    * @param m a Matrix object
    * @param x a real number (double)
    * @return a Matrix whose entries are the entries in m, 
    *         multiplied by x
    */
   public static Matrix multiply(Matrix m, double x) {
      double[][] entries = new double[m.getNumRows()][m.getNumColumns()];
      
      for (int row = 0; row < m.getNumRows(); row++) {
         for (int col = 0; col < m.getNumColumns(); col++) {
            entries[row][col] = m.getEntry(row, col) * x;
         }
      }
      
      return new Matrix(entries);
   }
   
   /**
    * multplies the given vector by the calling matrix.
    * Number of entries in the vector must match the
    * number of columns in the matrix passed.
    *
    * M = | a b c |
    *     | d e f |
    * 
    * v = | x |
    *     | y |
    *     | z |
    *
    * Mv = | ax + by + cz |
    *      | dx + ey + fz |
    *
    * The n-th entry in the returned vector is the dot product of the nth-row of m
    * with v.
    *
    */
   public Vector multiply(Vector v) {
      return Matrix.multiply(this, v);
   }  
   
   /**
    * multiplies the given vector by the given matrix. this method treats the
    * vector as a column vector. Number of entries in the vector must match the
    * number of columns in the matrix passed.
    *
    * M = | a b c |
    *     | d e f |
    * 
    * v = | x |
    *     | y |
    *     | z |
    *
    * Mv = | ax + by + cz |
    *      | dx + ey + fz |
    *
    * The n-th entry in the returned vector is the dot product of the nth-row of m
    * with v.
    *
    * @param M the Matrix object
    * @param v the Vector object
    * @return a Vector object containing Mv
    */
   
   public static Vector multiply(Matrix m, Vector v) {
      if (v.length() != m.getNumRows()) {
         throw new IllegalArgumentException("Incompatible shapes.\n");
      }
   
      double[] result = new double[m.matrix.length];
      
      for (int i = 0; i < m.matrix.length; i++) {
         result[i] = m.getRow(i).dot(v);
      }
      
      return new Vector(result);
   }
   
   /**
    * If the calling Matrix is a, returns the product ab.
    * Number of columns of calling matrix must match number of
    * rows of passed Matrix.
    * @param b a Matrix object
    * @return the product ab, where a is the calling matrix
    */
   public Matrix multiply(Matrix b) {
      return Matrix.multiply(this, b);
   }
   
   /**
    * multiplies two matrices together via matrix multiplication.
    * the [i][j]-th entry is the dot product of row i from Matrix a
    * and column j from Matrix b.
    * if a is an m-by-n Matrix and b is an n-by-p Matrix, then the
    * returned matrix ab will be an m-by-p Matrix.
    * The number of columns in a must match the number of rows in b,
    * or else an IllegalArgumentException will be thrown.
    * @param a an m-by-n Matrix object
    * @param b an n-by-p Matrix object
    * @return ab: an m-by-p Matrix object
    */
   public static Matrix multiply(Matrix a, Matrix b) {
      if (a.getNumColumns() != b.getNumRows()) {
         throw new IllegalArgumentException("Matrix dimensions are incompatible.");
      }
      
      double[][] entries = new double[a.getNumRows()][b.getNumColumns()];
      
      for (int row = 0; row < a.getNumRows(); row++) {
         for (int col = 0; col < b.getNumColumns(); col++) {
            entries[row][col] = Vector.dotProduct(a.getRow(row), b.getColumn(col));
         }
      }
      
      return new Matrix(entries);  
   }

   /**
    * accepts a Vector object and a column index and overwrites the entries 
    * in the specified row with the entries in the Vector object
    * if the column length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param col an index for the column
    * @param v a Vector object
    * @return a Matrix object with the col-th row replaced with v
    */   
   public Matrix setColumn(int col, Vector v) {
      return Matrix.setColumn(this, col, v);
   }

   /**
    * accepts an array of doubles, and a column index and overwrites the entries 
    * in the specified column with the entries in the given array.
    * if the clumn length and the array length do not match, an
    * IllegalArgumentException will be thrown.
    * @param col an index for the column
    * @param v a an array of doubles
    * @return a Matrix object with the col-th column replaced with v
    */     
   public Matrix setColumn(int col, double[] v) {
      return Matrix.setColumn(this, col, v);
   }

   /**
    * accepts a Matrix object, a Vector object and a column index and 
    * overwrites the entries in the specified column with the entries in 
    * the given Vector.
    * if the column length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param col an index for the col
    * @param v a Vector object containing the entries we want
    * @return a Matrix object with the col-th column replaced with v
    */
   public static Matrix setColumn(Matrix m, int col, Vector v) {
      if (m.matrix.length != v.length()) {
         throw new IllegalArgumentException("Vector length and does not match column length.");
      }
      
      if (col < 0 || col >= m.getNumColumns()) {
         throw new IllegalArgumentException("Invalid column: " + col);
      }
      
      Matrix n = new Matrix(m);
      
      for (int i = 0; i < n.matrix.length; i++) {
         n.matrix[i][col] = v.get(i);
      }
      
      return n;
   }
   
   /**
    * accepts an array of doubles, a Vector object, and a column index and 
    * overwrites the entries in the specified column with the entries in 
    * the given array
    * if the column length and the array length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param col an index for the col
    * @param v an array of doubles containing the entries we want
    * @return a Matrix object with the col-th column replaced with v
    */
   public static Matrix setColumn(Matrix m, int col, double[] v) {
      if (m.matrix.length != v.length) {
         throw new IllegalArgumentException("Array length and does not match row length.");
      }
      
      return Matrix.setColumn(m, col, new Vector(v));
   }

   /**
    * accepts an Vector object and a row index and overwrites the entries 
    * in the specified row with the entries in the Vector object
    * if the row length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param row an index for the row
    * @param v a Vector object
    * @return a Matrix object with the row-th row replaced with v
    */   
   public Matrix setRow(int row, Vector v) {
      return Matrix.setRow(this, row, v);
   }

   /**
    * accepts an array of doubles, and a row index and overwrites the entries 
    * in the specified row with the entries in the given array.
    * if the row length and the array length do not match, an
    * IllegalArgumentException will be thrown.
    * @param row an index for the row
    * @param v a an array of doubles
    * @return a Matrix object with the row-th row replaced with v
    */     
   public Matrix setRow(int row, double[] v) {
      return Matrix.setRow(this, row, v);
   }
   
   /**
    * accepts a Matrix object, a Vector object and a row index and 
    * overwrites the entries in the specified row with the entries in 
    * the given Vector.
    * if the row length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param row an index for the row
    * @param v a Vector object containing the entries we want
    * @return a Matrix object with the row-th row replaced with v
    */
   public static Matrix setRow(Matrix m, int row, Vector v) {
      if (m.matrix[row].length != v.length()) {
         throw new IllegalArgumentException("Vector length and does not match row length.");
      }
      
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("Invalid row: " + row);
      }
      
      Matrix n = new Matrix(m);
      
      for (int i = 0; i < n.matrix[row].length; i++) {
         n.matrix[row][i] = v.get(i);
      }
      
      return n;
   }

   /**
    * accepts a Matrix object, an array of doubles, and a row index and 
    * overwrites the entries in the specified row with the entries in 
    * the given array.
    * if the row length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param m a Matrix object
    * @param row an index for the row
    * @param v a an array of doubles
    * @return a Matrix object with the row-th row replaced with v
    */
   public static Matrix setRow(Matrix m, int row, double[] v) {
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("Invalid row: " + row);
      }
      
      if (m.matrix[row].length != v.length) {
         throw new IllegalArgumentException("Array length and does not match row length.");
      }
      
      return Matrix.setRow(m, row, new Vector(v));
   }
   
   /**
    * returns an array of two ints containing the dimensions of the Matrix
    * @param m a Matrix object
    * @return an array containing {numRows, numCols}
    */
   public static int[] shape(Matrix m) {
      return new int[] {m.getNumRows(), m.getNumColumns()};
   }

   /**
    * swaps the entries in col1 with the entries in col2.
    * uses two intermediary Vector objects to hold the entries
    * @param col1 the index of the first column to be swapped
    * @param col2 the index of the second column to be swapped
    * @return a Matrix with col1 and col2 swapped
    */
   public Matrix swapColumns(int col1, int col2) {
      return Matrix.swapColumns(this, col1, col2);
   }

   /**
    * swaps the entries in col1 with the entries in col2.
    * uses two intermediary Vector objects to hold the entries
    * @param m a Matrix object
    * @param col1 the index of the first column to be swapped
    * @param col2 the index of the second column to be swapped
    * @return a Matrix with col1 and col2 swapped
    */
   public static Matrix swapColumns(Matrix m, int col1, int col2) {
      if (col1 < 0 || col1 >= m.getNumColumns()) {
         throw new IllegalArgumentException("Invalid column: " + col1);
      }
      
      if (col2 < 0 || col2 >= m.getNumColumns()) {
         throw new IllegalArgumentException("Invalid column: " + col2);
      }
      
      return new Matrix(m).setColumn(col1, m.getColumn(col2))
                          .setColumn(col2, m.getColumn(col1));
   }
   
   /**
    * swaps the entries in row1 with the entries in row2.
    * uses two intermediary Vector objects to hold the entries
    * @param row1 the index of the first row to be swapped
    * @param row2 the index of the second row to be swapped
    * @return a Matrix with row1 and row2 swapped
    */
   public Matrix swapRows(int row1, int row2) {
      return Matrix.swapRows(this, row1, row2);
   }
   
   /**
    * swaps the entries in row1 with the entries in row2.
    * uses two intermediary Vector objects to hold the entries
    * @param m a Matrix object
    * @param row1 the index of the first row to be swapped
    * @param row2 the index of the second row to be swapped
    * @return a Matrix with row1 and row2 swapped
    */
   public static Matrix swapRows(Matrix m, int row1, int row2) {
      if (row1 < 0 || row1 >= m.getNumRows()) {
         throw new IllegalArgumentException("Invalid row: " + row1);
      }
      
      if (row2 < 0 || row2 >= m.getNumRows()) {
         throw new IllegalArgumentException("Invalid row: " + row2);
      }
      
      return new Matrix(m).setRow(row1, m.getRow(row2))
                          .setRow(row2, m.getRow(row1));
   }
   
   /**
    * transposes the matrix. swaps the rows and columns.
    * @param m a Matrix object
    * @return the transpose of the Matrix object
    */
   public static Matrix transpose(Matrix m) {
      double[][] n = new double[m.getNumColumns()][m.getNumRows()];
      
      for (int row = 0; row < n.length; row++) {
         for (int col = 0; col < n[row].length; col++) {
            n[row][col] = m.matrix[col][row];
         }
      }
      
      return new Matrix(n);
   }
   
   /**
    * returns a String containing the matrix entries in the following format:
    * [[1.2, 3.4, ..., 8.9],
    *  [9.7, 7.5, ..., 3.1]]
    * @return a String representation of the Matrix's entries
    */
   public String toString() {
      String str = "[";
      
      for (int i = 0; i < this.matrix.length; i++) {
         str += Matrix.getRow(this, i);
         
         if (i < this.matrix.length - 1) {
            str += ",\n ";
         } else {
            str += "]";
         }
      }
      
      return str;
   }
}