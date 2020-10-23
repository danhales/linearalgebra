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
    * @return the number of rows in the matrix
    */
   public int getNumRows() {
      return this.matrix.length;
   }
   
   /**
    * @return the number of columns in the matrix
    */
   public int getNumColumns() {
      return this.matrix[0].length;
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
    * returns an array of two ints containing the dimensions of the Matrix
    * @param m a Matrix object
    * @return an array containing {numRows, numCols}
    */
   public static int[] shape(Matrix m) {
      return new int[] {m.getNumRows(), m.getNumColumns()};
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