package linearalgebra;
import linearalgebra.Vector;

/**
 * The Matrix class stores a 2D array of doubles and provides
 * common linear algebra operations on that array, in addition
 * to producing special matrices and checking for specific conditions.
 */

public class Matrix {
   private double[][] matrix;
   
   public Matrix(double[][] matrix) {
      this.matrix = new double[matrix.length][matrix[0].length];
      
      for (int row = 0; row < matrix.length; row++) {
         for (int col = 0; col < matrix[row].length; col++) {
            this.matrix[row][col] = matrix[row][col];
         }
      }
   } 
   
   public int getNumRows() {
      return this.matrix.length;
   }
   
   public int getNumColumns() {
      return this.matrix[0].length;
   }
   
   public static Vector getRow(Matrix m, int row) {
      if (m == null) {
         throw new IllegalArgumentException("Matrix is null.");
      }
      if (row >= m.matrix.length) {
         throw new IllegalArgumentException("Row does not exist in matrix.");
      }
      
      return new Vector(m.matrix[row]);
   } 
   
   public static Vector getColumn(Matrix m, int col) {
      if (m == null) {
         throw new IllegalArgumentException("Matrix is null.");
      }
      if (col >= m.matrix[0].length) {
         throw new IllegalArgumentException("Column does not exist in matrix.");
      } 
      
      double[] columnVector = new double[m.matrix.length];
      
      for (int row = 0; row < m.matrix.length; row++) {
         columnVector[row] = m.matrix[row][col];
      }  
      
      return new Vector(columnVector);  
   }  
   
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