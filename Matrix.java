package linearalgebra;
import linearalgebra.Vector;

/**
 * The Matrix class stores a 2D array of doubles and provides
 * common linear algebra operations on that array, in addition
 * to producing special matrices and checking for specific conditions.
 */

public class Matrix {
   /**
    * entries is a 2D array of doubles to hold the entries in the Matrix
    */
   private double[][] entries;
   
   /*
    * threshold for double comparisons
    */
   public static final double THRESHOLD = Double.MIN_VALUE * 1000;
   
   /**
    * accepts a 2D array of doubles and copies them into the entries field.
    * @param entries a 2D array of doubles
    */
   public Matrix(double[][] entries) {
      this.entries = new double[entries.length][entries[0].length];
      
      for (int row = 0; row < entries.length; row++) {
         for (int col = 0; col < entries[row].length; col++) {
            this.entries[row][col] = entries[row][col];
         }
      }
   }
   
   /**
    * copy constructor accepts a Matrix object and duplicates it.
    * @param m the Matrix object we want to copy
    */
   public Matrix(Matrix m) {
      this.entries = new double[m.getNumRows()][m.getNumColumns()];
      
      for (int row = 0; row < m.entries.length; row++) {
         for (int col = 0; col < m.entries[row].length; col++) {
            this.entries[row][col] = m.entries[row][col];
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
      
      for (int row = 0; row < a.entries.length; row++) {
         for (int col = 0; col < a.entries[0].length; col++) {
            entries[row][col] = a.entries[row][col] + b.entries[row][col];
         }
      }
      
      return new Matrix(entries);
   }

   /**
    * adds the passed Vector to the specified column of the calling Matrix.
    * throws an IllegalArgumentException if the dimensions don't line up.
    * @param u a Vector object
    * @param col the column in m to which we want to add u
    * @return a Matrix with the values in u added to the specified column
    */
   public Matrix addVectorToColumn(Vector u, int col) {
      return Matrix.addVectorToColumn(this, u, col);
   }  
   
   /**
    * adds the passed Vector to the specified column of the passed Matrix.
    * throws an IllegalArgumentException if the dimensions don't line up.
    * @param m a Matrix object
    * @param u a Vector object
    * @param col the column in m to which we want to add u
    * @return a Matrix with the values in u added to the specified column
    */
   public static Matrix addVectorToColumn(Matrix m, Vector u, int col) {
      if (col < 0 || col >= m.getNumColumns()) {
         throw new IllegalArgumentException("col is not in the correct range");
      }
   
      if (u.length() != m.getColumn(col).length()) {
         throw new IllegalArgumentException("Vector length does not match row length");
      }
      
      return m.setColumn(col, m.getColumn(col).add(u));
   }
   
   /**
    * adds the passed Vector to the specified row of the calling Matrix.
    * throws an IllegalArgumentException if the dimensions don't line up.
    * @param u a Vector object
    * @param row the row in m to which we want to add u
    * @return a Matrix with the values in u added to the specified row
    */
   public Matrix addVectorToRow(Vector u, int row) {
      return Matrix.addVectorToRow(this, u, row);
   }   
   
   /**
    * adds the passed Vector to the specified row of the passed Matrix.
    * throws an IllegalArgumentException if the dimensions don't line up.
    * @param m a Matrix object
    * @param u a Vector object
    * @param row the row in m to which we want to add u
    * @return a Matrix with the values in u added to the specified row
    */
   public static Matrix addVectorToRow(Matrix m, Vector u, int row) {
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("row is not in the correct range");
      }
   
      if (u.length() != m.getRow(row).length()) {
         throw new IllegalArgumentException("Vector length does not match row length");
      }
      
      return m.setRow(row, m.getRow(row).add(u));
   }

   /**
    * minorMatrix accepts a row and a column, and returns the Matrix object
    * with that row and column dropped.
    * @param row an int
    * @param col an int
    * @return a copy of m, with the specified row and column dropped
    */
   public Matrix minorMatrix(int row, int col) {
      return Matrix.minorMatrix(this, row, col);
   }
   
   /**
    * minorMatrix accepts a row and a column, and returns the Matrix object
    * with that row and column dropped.
    * @param m a Matrix object
    * @param row an int
    * @param col an int
    * @return a copy of m, with the specified row and column dropped
    */
   public static Matrix minorMatrix(Matrix m, int row, int col) {
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("row is not in the correct range");
      }
      if (col < 0 || col >= m.getNumColumns()) {
         throw new IllegalArgumentException("col is not in the correct range");
      }
      
      return m.dropRow(row).dropColumn(col);
   }
   
   /**
    * determinant is an instance wrapper for static determinant method.
    * @return the determinant of the calling matrix.
    */
   public double determinant() {
      return Matrix.determinant(this);
   }
   
   /**
    * determinant computes the determinant of a square Matrix object.
    * The determinant is not defined for non-square matrices.
    * Determinant is computed recursively via cofactor expansion across
    * the first row:
    * https://en.wikipedia.org/wiki/Determinant#Laplace's_expansion_and_the_adjugate_entries
    * For a 1x1 Matrix, the determinant is the only entry.
    * @param m a Matrix object
    * @return the determinant of the Matrix (a double) 
    */
   public static double determinant(Matrix m) {
      if (!Matrix.isSquare(m)) {
         throw new IllegalArgumentException("Determinant undefined for non-square matrices");
      }
      
      double determinant = 0;
            
      if (m.isDiagonal() || m.isUpperTriangular() || m.isLowerTriangular()) {
         determinant = 1.0;
         
         // determinant is the product of the entries on the diagonal
         for (int i = 0; i < m.getNumRows(); i++) {
            determinant *= m.getEntry(i,i);
         }
      } else {
         if (Matrix.getNumRows(m) == 1) {
            determinant = Matrix.getEntry(m, 0, 0);
         } else {
            for (int col = 0; col < m.getNumColumns(); col++) {
               determinant += Math.pow(-1, col) *  
                              m.getEntry(0, col) *
                              m.minorMatrix(0, col).determinant();
            }
         }
      }

      return determinant;
   }
   
   /**
    * dropColumn accepts a non-negative int corresponding to a column index, and 
    * returns the Matrix, with that column dropped.
    * It does this by passing all columns except the column to be dropped
    * to fromColumnVectors.
    * @param col a column index {@literal (0 < col < entries[0].length)}
    * @return a Matrix with the specified column dropped
    */  
   public Matrix dropColumn(int col) {
      return Matrix.dropColumn(this, col);
   }
   
   /**
    * dropColumn accepts a Matrix object and a non-negative int corresponding
    * to a column index, and returns the Matrix, with that column dropped.
    * It does this by passing all columns except the column to be dropped
    * to fromColumnVectors.
    * @param m a Matrix object
    * @param col a column index {@literal (0 < col < entries[0].length)}
    * @return a Matrix with the specified column dropped
    */
   public static Matrix dropColumn(Matrix m, int col) {
      if (col < 0 || col >= m.getNumColumns()) {
         System.out.println("col = " + col);
         System.out.println("num columns = " + m.getNumColumns());
         throw new IllegalArgumentException("col is out of range");
      }
      
      Vector[] columns = new Vector[m.getNumColumns() - 1];
      
      for (int i = 0; i < m.entries[0].length; i++) {
         // if we're to the left of the column to be dropped, columns[i] = m.getColumn(i)
         if (i < col) {
            columns[i] = m.getColumn(i);
         }
         
         // if we're to the right of the column to be dropped, columns[i-1] = m.getColumn(i)
         if (i > col) {
            columns[i-1] = m.getColumn(i);
         }
      }
      
      return Matrix.fromColumnVectors(columns);
   }
   
   /**
    * dropRow accepts a non-negative int corresponding to a row index, and 
    * returns the Matrix, with that row dropped.
    * It does this by passing all rows except the row to be dropped
    * to fromRowVectors.
    * @param row a row index {@literal (0 < row < entries.length)}
    * @return a Matrix with the specified row dropped
    */  
   public Matrix dropRow(int row) {
      return Matrix.dropRow(this, row);
   }
   
   /**
    * dropRow accepts a Matrix object and a non-negative int corresponding
    * to a row index, and returns the Matrix, with that row dropped. It does 
    * this by passing all row except the row to be dropped to fromRowVectors.
    * @param m a Matrix object
    * @param row a row {@literal index (0 < col < entries.length)}
    * @return a Matrix with the specified row dropped
    */
   public static Matrix dropRow(Matrix m, int row) {
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("row is out of range");
      }
      
      Vector[] u = new Vector[m.getNumRows()-1];
      
      for (int i = 0; i < m.getNumRows(); i++) {
         if (i < row) {
            u[i] = m.getRow(i);
         }
         
         if (i > row) {
            u[i-1] = m.getRow(i);
         }
      }
      
      return Matrix.fromRowVectors(u);
   }
    
   /**
    * Constructs a Matrix object from column vectors.
    * Vectors must all have the same length (the length of the first Vector in
    * the array), or an IllegalArgumentException will be thrown.
    * @param vectors an array of Vector objects
    * @return a Matrix whose rows are the passed vectors
    */ 
   public static Matrix fromColumnVectors(Vector ... vectors) {
      return Matrix.transpose(Matrix.fromRowVectors(vectors));  
   }
   
   /**
    * Constructs a Matrix object from row vectors.
    * Vectors must all have the same length (the length of the first Vector in
    * the array), or an IllegalArgumentException will be thrown.
    * @param vectors an array of Vector objects
    * @return a Matrix whose rows are the passed vectors
    */ 
   public static Matrix fromRowVectors(Vector ... vectors) {
      for (Vector u : vectors) {
         if (u.length() != vectors[0].length()) {
            throw new IllegalArgumentException("Vectors do not have the same length.");
         }
      }      
      
      double[][] entries = new double[vectors.length][vectors[0].length()];
      
      for (int row = 0; row < vectors.length; row++) {
         for (int col = 0; col < vectors[row].length(); col++) {
            entries[row][col] = vectors[row].get(col);
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
      if (row < 0 || row > m.entries.length) {
         throw new IllegalArgumentException("Invalid value for row.");
      }
      
      if (col < 0 || col > m.entries[0].length) {
         throw new IllegalArgumentException("Invalid value for col.");
      }
      
      return m.entries[row][col];  
   }
   
   /**
    * @return the number of columns in the Matrix
    */
   public int getNumColumns() {
      return Matrix.getNumColumns(this);
   }

   /**
    * @param m a Matrix object
    * @return the number of columns in the Matrix
    */   
   public static int getNumColumns(Matrix m) {
      return m.entries[0].length;
   }  

   /**
    * @return the number of rows in the Matrix
    */
   public int getNumRows() {
      return Matrix.getNumRows(this);
   }
   
   /**
    * @param m a Matrix object
    * @return the number of columns in the Matrix
    */
   public static int getNumRows(Matrix m) {
      return m.entries.length;
   } 

   /**
    * accepts a column number (starting at 0) and returns the entries in
    * the column as a Vector object
    * @param col the number of the desired column 
    *        (starting at 0, up to this.entries[0].length-1)
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
    *        (starting at 0, up to this.entries[0].length-1)
    * @return a Vector containing the entries in the desired column
    */   
   public static Vector getColumn(Matrix m, int col) {
      if (col >= m.entries[0].length || col < 0) {
         throw new IllegalArgumentException("Column is out of range");
      } 
      
      double[] columnVector = new double[m.entries.length];
      
      for (int row = 0; row < m.entries.length; row++) {
         columnVector[row] = m.entries[row][col];
      }  
      
      return new Vector(columnVector);  
   }
 
   /**
    * accepts a row number (starting at 0) and returns the entries in
    * the row as a Vector object
    * @param row the number of the desired row (starting at 0, up to this.entries.length-1)
    * @return a Vector containing the entries in the desired row
    */
   public Vector getRow(int row) {
      return Matrix.getRow(this, row);
   }
   
   /**
    * accepts a row number (starting at 0) and returns the entries in
    * the row as a Vector object
    * @param m the Matrix object we want the row from
    * @param row the number of the desired row (starting at 0, up to this.entries.length-1)
    * @return a Vector containing the entries in the desired row
    */
   public static Vector getRow(Matrix m, int row) {
      if (row >= m.entries.length || row < 0) {
         throw new IllegalArgumentException("Row does not exist in Matrix");
      }
      
      return new Vector(m.entries[row]);
   }
   
   /**   
    * getIdentityMatrix returns an n by n Matrix with all zeroes
    * and ones on the diagonal.
    * @param n an integer greater than or equal to 1
    * @return an n-by-n Matrix with ones on the diagonal and zeroes
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
    * isDiagonal checks to see if all values except for the diagonal
    * are zero.
    * @return true if only the diagonal contains non-zero values, 
    *         false otherwise
    */
   public boolean isDiagonal() {
      return Matrix.isDiagonal(this);
   }
   
   /**
    * isDiagonal checks to see if all values except for the diagonal
    * are zero.
    * @param m a square Matrix object
    * @return true if only the diagonal contains non-zero values, 
    *         false otherwise
    */
   public static boolean isDiagonal(Matrix m) {
      if (!Matrix.isSquare(m)) {
         throw new IllegalArgumentException("Matrix is not square.");
      }  
   
      for (int row = 0; row < m.entries.length; row++) {
         for (int col = 0; col < m.entries[row].length; col++) {
            if (row != col) { // if we're not on the diagonal
            
               // if the value is non-zero
               if (Math.abs(m.entries[row][col]) > Matrix.THRESHOLD) {
                  return false;
               }
            }
         }
      }
      
      return true;  
   }
   
   /**
    * isLowerTriangular checks to see if all values below the diagonal
    * are zero.
    * @return true if all entries below the diagonal are zero, false otherwise
    */
   public boolean isLowerTriangular() {
      return Matrix.isLowerTriangular(this);
   }   
   
   /**
    * isLowerTriangular checks to see if all values below the diagonal
    * are zero.
    * @param m a Matrix object
    * @return true if all entries below the diagonal are zero, false otherwise
    */
   public static boolean isLowerTriangular(Matrix m) {
      
      for (int row = 1; row < m.entries.length; row++) {
         for (int col = 0; col < row; col++) { // only below the diagonal
            if (Math.abs(m.entries[row][col]) > Matrix.THRESHOLD) {
               return false;
            }
         }
      }
      
      return true;
   }

   /**
    * checks to see if there is exactly one 1 in each row and column,
    * and only zeroes everywhere else. Matrix must be square.
    * an n-by-n permutation Matrix rearranges the entries of an n-by-1
    * Vector
    * @return true if Matrix is a permutation Matrix, false otherwise
    */
   public boolean isPermutationMatrix() {
      return Matrix.isPermutationMatrix(this);
   }
   
   /**
    * checks to see if there is exactly one 1 in each row and column,
    * and only zeroes everywhere else. Matrix must be square.
    * an n-by-n permutation Matrix rearranges the entries of an n-by-1
    * Vector
    * @param m a Matrix object
    * @return true if Matrix is a permutation Matrix, false otherwise
    */
   public static boolean isPermutationMatrix(Matrix m) {
      if (!Matrix.isSquare(m)) {
         throw new IllegalArgumentException("Matrix is not square");
      }
      
      for (int i = 0; i < m.getNumRows(); i++) {
         if (!(m.getRow(i).isCanonicalBasisVector() 
            && m.getColumn(i).isCanonicalBasisVector())) {
            return false;
         }
      }
      
      return true;
   }
   
   /**
    * isSparse checks to see if most of the Matrix entries are zero.
    * "most" is defined as no more than max(num columns, num rows)
    * @return true if the Matrix contains very few nonzero entries,
    *         false otherwise
    */
   public boolean isSparse() {
      return Matrix.isSparse(this);
   }
   
   /**
    * isSparse checks to see if most of the Matrix entries are zero.
    * "most" is defined as no more than max(num columns, num rows)
    * @param m a Matrix object
    * @return true if the Matrix contains very few nonzero entries,
    *         false otherwise
    */
   public static boolean isSparse(Matrix m) {
      int numNonzero = 0;
      
      for (int row = 0; row < m.entries.length; row++) {
         for (int col = 0; col < m.entries[row].length; col++) {
            if (Math.abs(m.entries[row][col]) > Matrix.THRESHOLD) {
               numNonzero++;
            }
         }
      }
      
      return (numNonzero <= Math.max(m.getNumRows(), m.getNumColumns()));
   }
   
   /**
    * isSparse checks to see if most of the Matrix entries are zero.
    * "most" is defined as no more the specified proportion
    * @param p the desired proportion {@literal (0 < p <= 1)}
    * @return true if the Matrix contains very few nonzero entries,
    *         false otherwise
    */
   public boolean isSparse(double p) {
      return Matrix.isSparse(this, p);
   }
   
   /**
    * isSparse checks to see if most of the Matrix entries are zero.
    * "most" is defined as no more the specified proportion
    * @param m a Matrix object
    * @param p the desired proportion {@literal (0 < p <= 1)}
    * @return true if the Matrix contains very few nonzero entries,
    *         false otherwise
    */
   public static boolean isSparse(Matrix m, double p) {
      if (p < Matrix.THRESHOLD || p > 1) {
         throw new IllegalArgumentException("p is not in the correct range (0,1]");
      }
   
      int numNonzero = 0;
      
      for (int row = 0; row < m.entries.length; row++) {
         for (int col = 0; col < m.entries[row].length; col++) {
            if (Math.abs(m.entries[row][col]) > Matrix.THRESHOLD) {
               numNonzero++;
            }
         }
      }
      
      return ((double)numNonzero / (m.getNumRows() * m.getNumColumns())) <= p;
   }
   
   /**
    * isSquare checks the dimensions of the Matrix to see if they are
    * equal. if so, the Matrix is square, and the method returns true.
    * @return true if the number of rows equals the number of columns,
    *         false otherwise
    */
   public boolean isSquare() {
      return Matrix.isSquare(this);
   }
   
   /**
    * isSquare checks the dimensions of the Matrix to see if they are
    * equal. if so, the Matrix is square, and the method returns true.
    * @param m a Matrix object
    * @return true if the number of rows equals the number of columns,
    *         false otherwise
    */  
   public static boolean isSquare(Matrix m) {
      return (m.getNumRows() == m.getNumColumns());
   }

   /**
    * isUpperTriangular checks to see if all values above the diagonal
    * are zero.
    * @return true if all entries above the diagonal are zero, false otherwise
    */   
   public boolean isUpperTriangular() {
      return Matrix.isUpperTriangular(this);
   }
   
   /**
    * isUpperTriangular checks to see if all values above the diagonal
    * are zero.
    * @param m a Matrix object
    * @return true if all entries above the diagonal are zero, false otherwise
    */
   public static boolean isUpperTriangular(Matrix m) {
      
      for (int row = 0; row < m.entries.length; row++) {
         for (int col = row+1; col < m.entries[row].length; col++) { // only below the diagonal
            if (Math.abs(m.entries[row][col]) > Matrix.THRESHOLD) {
               return false;
            }
         }
      }
      
      return true;
   }

   /**
    * multiplies each entry in the Matrix m by a real number.
    * @param x a real number (double)
    * @return a Matrix whose entries are the entries in the calling entries, 
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
    * multplies the given vector by the calling entries. Number of entries in the 
    * vector must match the number of columns in calling Matrix.
    *
    * M = | a b c |
    *     | d e f |
    * 
    * u = | x |
    *     | y |
    *     | z |
    *
    * Mu = | ax + by + cz |
    *      | dx + ey + fz |
    *
    * The n-th entry in the returned vector is the dot product of the nth-row of m
    * with u.
    * @param u a Vector object
    * @return a Vector object containing the linear combination mu
    */
   public Vector multiply(Vector u) {
      return Matrix.multiply(this, u);
   }  
   
   /**
    * multiplies the given vector by the given entries. this method treats the
    * vector as a column vector. Number of entries in the vector must match the
    * number of columns in the entries passed.
    *
    * M = | a b c |
    *     | d e f |
    * 
    * u = | x |
    *     | y |
    *     | z |
    *
    * Mu = | ax + by + cz |
    *      | dx + ey + fz |
    *
    * The n-th entry in the returned vector is the dot product of the nth-row of m
    * with u.
    *
    * @param m the Matrix object
    * @param u the Vector object
    * @return a Vector object containing the linear combination mu
    */
   
   public static Vector multiply(Matrix m, Vector u) {
      if (u.length() != m.getNumColumns()) {
         throw new IllegalArgumentException("Incompatible shapes.\n");
      }
   
      double[] result = new double[m.entries.length];
      
      for (int i = 0; i < m.entries.length; i++) {
         result[i] = m.getRow(i).dot(u);
      }
      
      return new Vector(result);
   }
   
   /**
    * If the calling Matrix is a, returns the product ab.
    * Number of columns of calling entries must match number of
    * rows of passed Matrix.
    * @param b a Matrix object
    * @return the product ab, where a is the calling entries
    */
   public Matrix multiply(Matrix b) {
      return Matrix.multiply(this, b);
   }
   
   /**
    * multiplies two matrices together via entries multiplication.
    * the [i][j]-th entry is the dot product of row i from Matrix a
    * and column j from Matrix b.
    * if a is an m-by-n Matrix and b is an n-by-p Matrix, then the
    * returned entries ab will be an m-by-p Matrix.
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
            entries[row][col] = Vector.dot(a.getRow(row), b.getColumn(col));
         }
      }
      
      return new Matrix(entries);  
   }

   /**
    * accepts a Vector object and a column index and overwrites the entries 
    * in the specified row with the entries in the Vector object
    * if the column length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param col an index for the column
    * @param u a Vector object
    * @return a Matrix object with the col-th row replaced with u
    */   
   public Matrix setColumn(int col, Vector u) {
      return Matrix.setColumn(this, col, u);
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
    * @param u a Vector object containing the entries we want
    * @return a Matrix object with the col-th column replaced with u
    */
   public static Matrix setColumn(Matrix m, int col, Vector u) {
      if (m.entries.length != u.length()) {
         throw new IllegalArgumentException("Vector length and does not match column length.");
      }
      
      if (col < 0 || col >= m.getNumColumns()) {
         throw new IllegalArgumentException("Invalid column: " + col);
      }
      
      Matrix n = new Matrix(m);
      
      for (int i = 0; i < n.entries.length; i++) {
         n.entries[i][col] = u.get(i);
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
      if (m.entries.length != v.length) {
         throw new IllegalArgumentException("Array length and does not match row length.");
      }
      
      return Matrix.setColumn(m, col, new Vector(v));
   }
   
   /**
    * Changes the entry in the row-th row and the col-th column to value.
    * @param m a Matrix object
    * @param value a double––the value we want to set
    * @param row the index of the row
    * @param col the index of the column
    * @return the Matrix m, modified so the entry at [row][col] is value
    */
   public static Matrix setEntry(Matrix m, double value, int row, int col) {
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("row is out of range");
      }
      
      if (col < 0 || col >= m.getNumColumns()) {
         throw new IllegalArgumentException("col is out of range");
      }
   
      Matrix n = new Matrix(m);
      n.entries[row][col] = value;
      return n;
   }
    

   /**
    * accepts an Vector object and a row index and overwrites the entries 
    * in the specified row with the entries in the Vector object
    * if the row length and the Vector length do not match, an
    * IllegalArgumentException will be thrown.
    * @param row an index for the row
    * @param u a Vector object
    * @return a Matrix object with the row-th row replaced with u
    */   
   public Matrix setRow(int row, Vector u) {
      return Matrix.setRow(this, row, u);
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
    * @param u a Vector object containing the entries we want
    * @return a Matrix object with the row-th row replaced with u
    */
   public static Matrix setRow(Matrix m, int row, Vector u) {
      if (m.entries[row].length != u.length()) {
         throw new IllegalArgumentException("Vector length and does not match row length.");
      }
      
      if (row < 0 || row >= m.getNumRows()) {
         throw new IllegalArgumentException("Invalid row: " + row);
      }
      
      Matrix n = new Matrix(m);
      
      for (int i = 0; i < n.entries[row].length; i++) {
         n.entries[row][i] = u.get(i);
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
      
      if (m.entries[row].length != v.length) {
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
    * Subtracts the entries of n from the corresponding entries in m, and
    * returns the matrix this-n. this and n must have the same shape, or an
    * IllegalArgumentException will be thrown.
    * @param n the matrix on the right of m-n
    * @return a Matrix where each entry is this[row][col] - n[row][col]
    */
   public Matrix subtract(Matrix n) {
      return Matrix.subtract(this, n);
   }  
   
   /**
    * Subtracts the entries of n from the corresponding entries in m, and
    * returns the matrix m-n. m and n must have the same shape, or an
    * IllegalArgumentException will be thrown.
    * @param m the matrix on the left of m-n
    * @param n the matrix on the right of m-n
    * @return a Matrix where each entry is m[row][col] - n[row][col]
    */
   public static Matrix subtract(Matrix m, Matrix n) {
      return Matrix.add(m, Matrix.multiply(n, -1));
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
    * decomposes the Matrix object into column vectors
    * @return an array of Vector objects
    */
   public Vector[] toColumnVectors() {
      return Matrix.toColumnVectors(this);
   }
   
   /**
    * decomposes the Matrix object into column vectors
    * @param m a Matrix object
    * @return an array of Vector objects
    */
   public static Vector[] toColumnVectors(Matrix m) {
      return Matrix.toRowVectors(m.transpose());
   }
   
   /**
    * decomposes the Matrix object into row vectors
    * @return an array of Vector objects
    */
   public Vector[] toRowVectors() {
      return Matrix.toRowVectors(this);
   }   
   
   /**
    * decomposes the Matrix object into row vectors
    * @param m a Matrix object
    * @return an array of Vector objects
    */
   public static Vector[] toRowVectors(Matrix m) {
      Vector[] v = new Vector[m.getNumRows()];
      
      for (int i = 0; i < m.getNumRows(); i++) {
         v[i] = new Vector(m.getRow(i));
      }
      
      return v;
   }
   
   
   /**
    * returns the sum of the entries on the diagonal of a square entries, which is
    * also the sum of the entries's eigenvalues
    * @return the sum of the entries on the diagonal
    */
   public double trace() {
      return Matrix.trace(this);
   }
   
   /**
    * returns the sum of the entries on the diagonal of a square entries, which is
    * also the sum of the entries's eigenvalues
    * @param m a square Matrix object
    * @return the sum of the entries on the diagonal
    */
   public static double trace(Matrix m) {
      if (!Matrix.isSquare(m)) {
         throw new IllegalArgumentException("Trace is not defined for non-square entries");
      }
      
      double trace = 0;
      
      for (int i = 0; i < Matrix.getNumRows(m); i++) {
         trace += m.getEntry(i, i);
      }
      
      return trace;
   }
   
   /**
    * transposes the entries. swaps the rows and columns.
    * @return the transpose of the Matrix object
    */
   public Matrix transpose() {
      return Matrix.transpose(this);
   }   
   
   /**
    * transposes the entries. swaps the rows and columns.
    * @param m a Matrix object
    * @return the transpose of the Matrix object
    */
   public static Matrix transpose(Matrix m) {
      double[][] n = new double[m.getNumColumns()][m.getNumRows()];
      
      for (int row = 0; row < n.length; row++) {
         for (int col = 0; col < n[row].length; col++) {
            n[row][col] = m.entries[col][row];
         }
      }
      
      return new Matrix(n);
   }
   
   /**
    * returns a String containing the entries entries in the following format:
    * [[1.2, 3.4, ..., 8.9],
    *  [9.7, 7.5, ..., 3.1]]
    * @return a String representation of the Matrix's entries
    */
   public String toString() {
      String str = "[";
      
      for (int i = 0; i < this.entries.length; i++) {
         str += Matrix.getRow(this, i);
         
         if (i < this.entries.length - 1) {
            str += ",\n ";
         } else {
            str += "]";
         }
      }
      
      return str;
   }
}
