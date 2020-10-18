package linearalgebra;

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
    * v contains the entries in the vector
    */
   private double[] v;
   
   /**
    * Constructor makes a copy of the array passed.
    * @param v an array containing the entries in the vector
    */
   public Vector(double ... v) {
      this.v = new double[v.length];
      for (int i = 0; i < v.length; i++) {
         this.v[i] = v[i];      
      }
   }
   
   /**
    * add method accepts a vector and adds it to the current vector
    * @param u the vector to add onto the calling vector.
    * @return a Vector object whose entries are the element-wise sums
    *    of the calling vector and the argument
    */
   public Vector add(Vector u) {
      return Vector.sum(this, u);
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
      return Math.acos(Vector.dotProduct(u1, u2) / (u1.magnitude() * u2.magnitude()));
   }
   
   /**
    * checkLengths method accepts two vectors and throws and 
    * IllegalArgumentException if they are not the same lengths.
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
      return Vector.crossProduct(this, u);
   }
   
   /**
    * crossProduct method takes two vectors of length 3 and returns
    * their cross product. Note that this operation is anticommutative,
    * so crossProduct(a, b) = -crossProduct(b, a)
    * @param a the left vector Vector
    * @param b the right vector Vector
    * @return the cross product a X b
    */
   public static Vector crossProduct(Vector a, Vector b) {
      // check to make sure both vectors are the right length
      if (a.length() != 3) {
         throw new IllegalArgumentException("Invalid vector length (first vector)");
      }
      if (a.length() != 3) {
         throw new IllegalArgumentException("Invalid vector length (second vector)");
      }
      checkLengths(a, b); // just in case      
      
      double[] entries = new double[] {
         a.v[1] * b.v[2] - a.v[2] * b.v[1],
         a.v[2] * b.v[0] - a.v[0] * b.v[2],
         a.v[0] * b.v[1] - a.v[1] * b.v[0]};
         
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
      return Vector.dotProduct(this, u);
   }
   
   /**
    * dotProduct method computes the dot product of two vectors.
    * assumes vectors have the same length.
    * @param u1 a Vector object
    * @param u2 a Vector object
    * @return the sum of the products of corresponding elements
    */
   public static double dotProduct(Vector u1, Vector u2) {
      checkLengths(u1, u2);
   
      double sum = 0;
      
      for (int i = 0; i < u1.length(); i++) {
         sum += (u1.get(i) * u2.get(i));
      }
      
      return sum;
   }
   
   /**
    * identityVector returns an identity vector (whose entries are
    * all zeros).
    * @param length the length of the vector
    * @return a vector with all zeros
    */
   public static Vector identityVector(int length) {
      return new Vector(new double[length]);
   }
   
   /**
    * inverse returns the inverse of the calling vector.
    * @return a Vector with the signs flipped on all entries
    */
   public Vector inverse() {
      return this.multiply(-1);
   }
   
   /**
    * inverseVector returns the additive inverse of the vector passed.
    * @param u a Vector
    * @return a Vector whose entries have the signs flipped
    */
   public static Vector inverseVector(Vector u) {
      return Vector.product(u, -1);
   }
   
   /**
    * isZero checks to see if all entries are zero.
    * @return true if all entries in v are zero, false otherwise
    */
   public boolean isZero() {
      for (double entry : v) {
         if (entry != 0) { // if a non-zero entry is found
            return false;
         }
      }
      
      return true;
   }
   
   /**
    * magnitude method is a wrapper for normL2
    * @return the magnitude of the vector
    */
   public double magnitude() {
      return Vector.normL2(this);
   }
   
   /**
    * multiply method accepts a scalar to and multiplies each entry
    * in v by it.
    * @param u the vector to add onto the calling vector.
    * @return a Vector object whose entries are the element-wise sums
    *    of the calling vector and the argument
    */
   public Vector multiply(double scalar) {
      return Vector.product(this, scalar);
   }
   
   /**
    * normalize scales the passed vector by dividing it by its 
    * magnitude. if the zero vector is passed, an IllegalArgumentException 
    * is thrown.
    * @param u a Vector object
    * @return a Vector object
    */
   public static Vector normalize(Vector v) {
      if (v.isZero()) {
         throw new IllegalArgumentException();
      } else {
         return v.multiply(1.0/v.magnitude());
      }
   }
   
   /**
    * normalize scales the calling vector by dividing it by its 
    * magnitude. if the zero vector is passed, an IllegalArgumentException 
    * is thrown.
    * @return a Vector object
    */
   public Vector normalize() {
      if (this.isZero()) {
         throw new IllegalArgumentException();
      } else {
         return this.multiply(1.0/this.magnitude());
      }
   }   
   
   
   /**
    * normL1 accepts a Vector and returns the L1 norm (the sum
    * of the absolute values of the vector's entries)
    * @param u a Vector object
    * @return the L1 norm of the vector
    */
   public static double normL1(Vector u) {
      double sum = 0;
      
      for (int i = 0; i < u.length(); i++) {
         sum += Math.abs(u.get(i));
      }
      
      return sum;
   }
   
   /**
    * an instance method that calls normL1 on the current object.
    * @return the L1 norm of the calling vector.
    */
   public double normL1() {
      return Vector.normL1(this);
   }

   /**
    * normL2 accepts a Vector and returns the L1 norm (the sum
    * of the squares of the vector's entries)
    * @param u a Vector object
    * @return the L2 norm of the vector
    */
   public static double normL2(Vector u) {
      double sum = 0;
      
      for (int i = 0; i < u.length(); i++) {
         sum += Math.pow(u.get(i), 2);
      }
      
      return Math.sqrt(sum);
   }
   
   /**
    * an instance method that calls normL1 on the current object.
    * @return the L2 norm of the calling vector.
    */
   public double normL2() {
      return Vector.normL2(this);
   }

   /**
    * product accepts a Vector object and a scalar and returns
    * a Vector whose entries are the entries of the Vector, multiplied
    * by the scalar.
    * @param u a Vector object
    * @param scalar a real number
    * @return the scalar product of the vector and the scalar
    */
   public static Vector product(Vector u, double scalar) {   
      double[] products = new double[u.length()];
      
      for (int i = 0; i < products.length; i++) {
         products[i] = scalar * u.get(i);
      }
      
      return new Vector(products);
   }
   
   /**
    * scalarTripleProduct computes a.dot(b.cross(c))
    * @param a a Vector object
    * @param b a Vector object
    * @param c a Vector object
    * @return the scalar triple product a.dot(b.cross(c))
    */
   public static double scalarTripleProduct(Vector a, Vector b, Vector c) {
      return Vector.dotProduct(a, Vector.crossProduct(b, c));
   }
   
   /**
    * sum method accepts two vectors and returns their element-wise
    * sum in a new Vector object. Assumes v1 and v2 have the same 
    * length.
    * @param u1 a Vector object
    * @param u2 a Vector object
    */
   public static Vector sum(Vector u1, Vector u2) {
      checkLengths(u1, u2);
      
      double[] sums = new double[u1.length()];
      
      for (int i = 0; i < sums.length; i++) {
         sums[i] = u1.get(i) + u2.get(i);
      }
      
      return new Vector(sums);
   }
   
   // Setters, getters, and overridden methods.
   
   /**
    * Returns the entry in the passed position.
    * @param position the position to return
    * @return the value in v[position]
    */
   public double get(int position) {
      return v[position];
   }
   
   /**
    * getLength method returns the number of entries in the 
    * vector.
    * @return the length of v
    */
   public int length() {
      return this.v.length;
   }
      
   /**
    * Returns a copy of v, not a reference to v.
    * @return a copy of the array v
    */
   public double[] getV() {
      double[] v = new double[this.v.length];
      
      for (int i = 0; i < this.v.length; i++) {
         v[i] = this.v[i];
      }
      
      return v;
   }
   
   /**
    * Sets the values in the v array.
    * @param v an array of doubles
    */
   public void setV(double[] v) {
      this.v = new double[v.length];
      
      for (int i = 0; i < v.length; i++) {
         this.v[i] = v[i];
      }
   }
   
   /**
    * set method modifies the element at index to equal value.
    * @param index the index we want to modify
    * @param value the new value
    */
   public void set(int index, double value) {
      this.v[index] = value;
   }
   
   /**
    * Return a String containing the vector represented as a row in brackets.
    * @return a String representation of the vector
    */
   @Override
   public String toString() {
      String str = "[";
      String sep = "  ";
      
      for (int i = 0; i < this.v.length; i++) {
         str += this.v[i];
         
         if (i < (this.v.length - 1)) { // if we're not at the last entry
            str += sep;
         }
      }
      
      return str + "]";
   }
}