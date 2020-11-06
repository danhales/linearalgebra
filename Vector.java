package linearalgebra;
import linearalgebra.Matrix;

/**
 * The Vector class provides basic vector operations for Euclidean vectors
 * represented as arrays of real numbers. 
 *
 * All operations between two vectors are designed for vectors of the same
 * length, and no checking is done. For loops are controlled by the length
 * of the first vector, so if the second vector is longer, an Exception may
 * not be thrown as expected.
 */
 
public class Vector {
   /**
    * entries contains the entries in the vector
    */
   private double[] entries;
   
   /*
    * threshold for double comparisons
    */
   public static final double THRESHOLD = Double.MIN_VALUE * 1000;
   
   /**
    * Constructor makes a copy of the array passed.
    * @param entries an array containing the entries in the vector
    */
   public Vector(double ... entries) {
      this.entries = new double[entries.length];
      for (int i = 0; i < entries.length; i++) {
         this.entries[i] = entries[i];      
      }
   }
   
   /**
    * copy constructor copies the entires in vect into v
    * @param u a Vector object
    */
   public Vector(Vector u) {
      this.entries = new double[u.entries.length];
      
      for (int i = 0; i < u.entries.length; i++) {
         this.entries[i] = u.entries[i];
      }
   }
   
   /**
    * add method accepts a vector and adds it to the current vector
    * @param u the vector to add onto the calling vector.
    * @return a Vector object whose entries are the element-wise sums
    *    of the calling vector and the argument
    */
   public Vector add(Vector u) {
      return Vector.add(this, u);
   }


   /**
    * add method accepts two vectors and returns their element-wise
    * sum in a new Vector object. Assumes v1 and v2 have the same 
    * length.
    * @param u1 a Vector object
    * @param u2 a Vector object
    * @return a Vector objects whose entries are the sums of corresponding
    *         entries in u1 and u2
    */
   public static Vector add(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
      
      double[] sums = new double[u1.length()];
      
      for (int i = 0; i < sums.length; i++) {
         sums[i] = u1.get(i) + u2.get(i);
      }
      
      return new Vector(sums);
   }
   
   /**
    * angleDegrees method computes the angle between the two vectors, 
    * computed as arccos(u1.dot(u2) / (u1.magnitude() * u2.magnitude())
    * @param u1 a Vector object
    * @param u2 a Vector object
    * @return the angle between u1 and u2 (in radians)
    */
   public static double angleDegrees(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
      return Vector.angleRadians(u1, u2) * 180 / Math.PI;
   }
   
   /**
    * angleRadians method computes the angle between the two vectors, 
    * computed as arccos(u1.dot(u2) / (u1.magnitude() * u2.magnitude())
    * @param u1 a Vector object
    * @param u2 a Vector object
    * @return the angle between u1 and u2 (in radians)
    */
   public static double angleRadians(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
      return Math.acos(Vector.dot(u1, u2) / (u1.magnitude() * u2.magnitude()));
   }
   
   /**
    * checkLengths method accepts two vectors and throws and 
    * IllegalArgumentException if they are not the same lengths.
    * @param u1 a Vector object
    * @param u2 a Vector object
    */
   public static void checkLengths(Vector u1, Vector u2) {
      if (u1.length() != u2.length()) {
         throw new IllegalArgumentException("Vectors are different lengths");
      }
   }
   
   /**
    * cross computes the cross product (this x u)
    * @param u the vector to cross the calling vector with
    * @return the cross product this x u
    */
   public Vector cross(Vector u) {
      return Vector.cross(this, u);
   }
   
   /**
    * cross method takes two vectors of length 3 and returns their
    * cross product. Note that this operation is anticommutative, so
    * cross(a, b) = -cross(b, a)
    * @param a the left vector Vector
    * @param b the right vector Vector
    * @return the cross product a X b
    */
   public static Vector cross(Vector a, Vector b) {
      // check to make sure both vectors are the right length
      if (a.length() != 3) {
         throw new IllegalArgumentException("Invalid vector length (first vector)");
      }
      if (a.length() != 3) {
         throw new IllegalArgumentException("Invalid vector length (second vector)");
      }
      Vector.checkLengths(a, b); // just in case      
      
      double[] entries = new double[] {
         a.entries[1] * b.entries[2] - a.entries[2] * b.entries[1],
         a.entries[2] * b.entries[0] - a.entries[0] * b.entries[2],
         a.entries[0] * b.entries[1] - a.entries[1] * b.entries[0]};
         
      return new Vector(entries);
   }
   
   /**
    * dot method computes the dot product of the calling vector and
    * the passed vectored.
    * assumes vectors have the same length.
    * @param u a Vector object
    * @return the sum of the products of corresponding elements
    */
   public double dot(Vector u) {
      return Vector.dot(this, u);
   }
   
   /**
    * dot method computes the dot product of two vectors.
    * assumes vectors have the same length.
    * @param u1 a Vector object
    * @param u2 a Vector object
    * @return the sum of the products of corresponding elements
    */
   public static double dot(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
   
      double sum = 0;
      
      for (int i = 0; i < u1.length(); i++) {
         sum += (u1.get(i) * u2.get(i));
      }
      
      return sum;
   }

   /**
    * Returns the entry in the specified position.
    * @param position the position to return
    * @return the value in entries[position]
    */
   public double get(int position) {
      return Vector.get(this, position);
   }
   
   /**
    * Returns the entry in the specified position.
    * @param u a Vector object
    * @param position the position to return
    * @return the value in u[position]
    */
   public static double get(Vector u, int position) {
      return u.entries[position];
   }

   /**
    * returns a copy of entries, not a reference to entries.
    * @return a copy of the array entries
    */
   public double[] getEntries() {
      double[] entries = new double[this.entries.length];
      
      for (int i = 0; i < this.entries.length; i++) {
         entries[i] = this.entries[i];
      }
      
      return entries;
   }

   
   /**
    * identityVector returns an additive identity vector (whose entries are
    * all zeros).
    * @param length the length of the vector
    * @return a vector with all zeros
    */
   public static Vector identityVector(int length) {
      return new Vector(new double[length]);
   }
   
   /**
    * inverseVector returns the additive inverse of the calling vector.
    * @return a Vector with the signs flipped on all entries
    */
   public Vector inverseVector() {
      return this.multiply(-1);
   }
   
   /**
    * inverseVector returns the additive inverse of the vector passed.
    * @param u a Vector
    * @return a Vector whose entries have the signs flipped
    */
   public static Vector inverseVector(Vector u) {
      return Vector.multiply(u, -1);
   }
   
   /**
    * checks to see if the Vector object is a canonical basis Vector,
    * i.e. it has a one in exactly one entry and zeroes everywhere
    * else.
    * @return true if Vector contains all zeroes and a single one
    */
   public boolean isCanonicalBasisVector() {
      return Vector.isCanonicalBasisVector(this);
   }
   
   /**
    * checks to see if the Vector object is a canonical basis Vector,
    * i.e. it has a one in exactly one entry and zeroes everywhere
    * else.
    * @param u a Vector object
    * @return true if Vector contains all zeroes and a single one
    */
   public static boolean isCanonicalBasisVector(Vector u) {
      int numOnes = 0;
      int numZeros = 0;
      
      for (int i = 0; i < u.length(); i++) {
         if (Math.abs(1 - u.get(i)) < Vector.THRESHOLD) {
            numOnes++;
         }
      }
      
      return numOnes == 1;
   }

   /**
    * isZero checks to see if all entries are zero.
    * @return true if all entries are zero, false otherwise
    */
   public boolean isZero() {
      return Vector.isZero(this);
   }   
   
   /**
    * isZero checks to see if all entries are zero.
    * @param u a Vector object
    * @return true if all entries in u are zero, false otherwise
    */
   public static boolean isZero(Vector u) {
      for (double entry : u.entries) {
         if (Math.abs(entry) > Vector.THRESHOLD) { // if a non-zero entry is found
            return false;
         }
      }
      
      return true;
   }

   /**
    * length method returns the number of entries in the 
    * vector.
    * @return the length of v
    */
   public int length() {
      return Vector.length(this);
   }

   /**
    * length method returns the number of entries in the 
    * vector.
    * @param u a Vector object
    * @return the length of u
    */
   public static int length(Vector u) {
      return u.entries.length;
   }

   
   /**
    * Creates a linear combination (weighted sum) of the Vector objects
    * and the weights. Throws IllegalArgumentException if the length of
    * the weights array does not match the length of the vectors array
    * @param vectors an array of Vector objects
    * @param weights an array of doubles to weight the sum
    * @return the linear combination of the vectors with the weights
    */
   public static Vector linearCombination(Vector[] vectors, double[] weights) {
      if (vectors.length != weights.length) {
         throw new IllegalArgumentException("Number of vectors does not match number of weights.");
      }
      
      // start the sum by weighting the first vector
      Vector sum = new Vector(vectors[0].multiply(weights[0]));
      
      // weight and add each Vector onto sum
      for (int i = 1; i < vectors.length; i++) {
         sum = sum.add(vectors[i].multiply(weights[i]));
      }
      
      return sum;
   }
   
   /**
    * magnitude method is a wrapper for pnorm, with p=2
    * @return the magnitude of the vector
    */
   public double magnitude() {
      return Vector.magnitude(this);
   }
   
   /**
    * magnitude method is a wrapper for pnorm, with p=2
    * @param u a Vector object
    * @return the magnitude of the vector
    */
   public static double magnitude(Vector u) {
      return Vector.pnorm(u, 2);
   }
   
   /**
    * multiply method accepts a scalar to and multiplies each element of
    * entries by that value.
    * @param scalar the real number to multiply the entries by
    * @return a Vector object whose entries are the element-wise sums
    *    of the calling vector and the argument
    */
   public Vector multiply(double scalar) {
      return Vector.multiply(this, scalar);
   }
   
   /**
    * multiply accepts a Vector object and a scalar and returns
    * a Vector whose entries are the entries of the Vector, multiplied
    * by the scalar.
    * @param u a Vector object
    * @param scalar a real number
    * @return the scalar product of the vector and the scalar
    */
   public static Vector multiply(Vector u, double scalar) {   
      double[] products = new double[u.length()];
      
      for (int i = 0; i < products.length; i++) {
         products[i] = scalar * u.get(i);
      }
      
      return new Vector(products);
   }

   /**
    * normalize scales the calling vector by dividing it by its 
    * magnitude. if the zero vector is passed, an IllegalArgumentException 
    * is thrown.
    * @return a Vector object
    */
   public Vector normalize() {
        return Vector.normalize(this);
   }  
   
   /**
    * normalize scales the passed vector by dividing it by its 
    * magnitude. if the zero vector is passed, an IllegalArgumentException 
    * is thrown.
    * @param u a Vector object
    * @return a Vector object
    */
   public static Vector normalize(Vector u) {
      if (u.isZero()) {
         throw new IllegalArgumentException();
      } else {
         return u.multiply(1.0/u.magnitude());
      }
   } 
   
   /**
    * scalarTripleProduct computes a.dot(b.cross(c))
    * @param a a Vector object
    * @param b a Vector object
    * @param c a Vector object
    * @return the scalar triple product a.dot(b.cross(c))
    */
   public static double scalarTripleProduct(Vector a, Vector b, Vector c) {
      return Vector.dot(a, Vector.cross(b, c));
   }
      
   /**
    * The outer product is matrix multiplication on this and
    * the transpose of u.
    * @param u a Vector
    * @return the outer product
    */
   public Matrix outerProduct(Vector u) {
      return Vector.outerProduct(this, u);
   }
   
   /**
    * The outer product is matrix multiplication on u1 and
    * the transpose of u2.
    * @param u1 a Vector
    * @param u2 a Vector
    * @return the outer product
    */
   public static Matrix outerProduct(Vector u1, Vector u2) {
      Matrix m1 = Matrix.fromColumnVectors(u1);
      Matrix m2 = Matrix.fromRowVectors(u2);
      return Matrix.multiply(m1, m2);
   }

   /**
    * an instance method that calls normL1 on the current object.
    * @param p a real number greater than or equal to 1
    * @return the L2 norm of the calling vector.
    */
   public double pnorm(double p) {
      return Vector.pnorm(this, p);
   }

   /**
    * pnorm accepts a Vector and a value for p and returns the Lp norm 
    * (the p-th root of the sum of the p-th power of the absolute value of 
    * the enttries
    * @param u a Vector object
    * @param p a real number greater than or equal to 1
    * @return the Lp norm of the vector
    */
   public static double pnorm(Vector u, double p) {
      if (p < 1) {
         throw new IllegalArgumentException("p must be >= 1");
      }
   
      double sum = 0;
      
      for (int i = 0; i < u.length(); i++) {
         sum += Math.pow(Math.abs(u.get(i)), p);
      }
      
      return Math.pow(sum, 1/p);
   }
   
   /**
    * projects the calling Vector onto the passed Vector
    * @param u the Vector we want to project this one onto
    * @return the orthogonal projection of this vector onto u
    */
   public Vector orthogonalProjection(Vector u) {
      return Vector.orthogonalProjection(this, u);
   }
   
   /**
    * projects u1 onto u2
    * @param u1 the Vector whose orthogonal projection we want
    * @param u2 the Vector we are projecting u1 onto
    * @return the orthogonal projection
    */
   public static Vector orthogonalProjection(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
      
      if (u2.isZero()) {
         throw new IllegalArgumentException("Cannot project onto zero vector");
      }
      
      return Vector.multiply(u2, u1.dot(u2) / u2.dot(u2));
   }
   
   /**
    * set method modifies the element at index to equal value.
    * @param index the index we want to modify
    * @param value the new value
    * @return a Vector with the value at index updated
    */
   public Vector set(int index, double value) {
      return Vector.set(this, index, value);
   }

   /**
    * set method modifies the element at index to equal value.
    * @param u a Vector object
    * @param index the index we want to modify
    * @param value the new value
    * @return a Vector with the value at index updated
    */   
   public static Vector set(Vector u, int index, double value) {
      if (index < 0 || index >= u.length()) {
         throw new IllegalArgumentException("Index is out of range");
      }
   
      double[] entries = new double[u.entries.length];
      entries[index] = value;
      
      return new Vector(entries);
   } 

   /**
    * Sets the values in the entries array.
    * @param entries an array of doubles
    */
   public void setEntries(double[] entries) {
      this.entries = new double[entries.length];
      
      for (int i = 0; i < entries.length; i++) {
         this.entries[i] = entries[i];
      }
   }

   /**
    * subtract method subtracts the passed Vector from the calling Vector.
    * @param u a Vector object
    * @return a Vector object whose entries are the difference of the
    *         entries in the calling Vector and the respective entries 
    *         in v
    */
   public Vector subtract(Vector u) {
      return Vector.subtract(this, u);
   }

   /**
    * subtract method returns the difference of two vectors. note
    * that difference is a special case of sum (v1 + (-1)*v2)
    * @param v1 a Vector object
    * @param v2 a Vector object
    * @return a new Vector object whose entries are the differences of
    *         the entries in v1 and v2 (v1 - v2)
    */
   public static Vector subtract(Vector v1, Vector v2) {
      return Vector.add(v1, v2.multiply(-1));
   }  
   
   /**
    * Return a String containing the vector represented as a row in brackets, e.g.
    * [1.0, 2.2, 3.1, 4.9, 5.7]
    * @return a String representation of the vector
    */
   @Override
   public String toString() {
      String str = "[";
      String sep = ", ";
      
      for (int i = 0; i < this.entries.length; i++) {
         str += this.entries[i];
         
         if (i < (this.entries.length - 1)) { // if we're not at the last entry
            str += sep;
         }
      }
      
      return str + "]";
   }
}