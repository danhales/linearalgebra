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
    * copy constructor copies the entires in vect into v
    * @param vect a Vector object
    */
   public Vector(Vector vect) {
      this.v = new double[vect.v.length];
      
      for (int i = 0; i < vect.v.length; i++) {
         this.v[i] = vect.v[i];
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
      Vector.checkLengths(a, b); // just in case      
      
      double[] entries = new double[] {
         a.v[1] * b.v[2] - a.v[2] * b.v[1],
         a.v[2] * b.v[0] - a.v[0] * b.v[2],
         a.v[0] * b.v[1] - a.v[1] * b.v[0]};
         
      return new Vector(entries);
   }
   
   /**
    * difference method returns the difference of two vectors. note
    * that difference is a special case of sum (v1 + (-1)*v2)
    * @param v1 a Vector object
    * @param v2 a Vector object
    * @return a new Vector object whose entries are the differences of
    *         the entries in v1 and v2 (v1 - v2)
    */
   public static Vector difference(Vector v1, Vector v2) {
      return Vector.sum(v1, v2.multiply(-1));
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
      Vector.checkLengths(u1, u2);
   
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
         if (Math.abs(entry) > 0.0000000001) { // if a non-zero entry is found
            return false;
         }
      }
      
      return true;
   }
   
   /**
    * Creates a linear combination (weighted sum) of the Vector objects
    * and the weights. Throws IllegalArgumentException if the length of
    * the weights array does not match the length of the vectors array
    * @param vectors an array of Vector objects
    * @param weights an array of doubles to weight the sum
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
    * magnitude method is a wrapper for normL2
    * @return the magnitude of the vector
    */
   public double magnitude() {
      return Vector.pnorm(this, 2);
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
        return Vector.normalize(this);
   }   

   /**
    * pnorm accepts a Vector and a value for p and returns the Lp norm 
    * (the p-th root of the sum of the p-th power of the absolute value of 
    * the enttries
    * @param u a Vector object
    * @param p a real number >= 1
    * @return the L2 norm of the vector
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
    * an instance method that calls normL1 on the current object.
    * @return the L2 norm of the calling vector.
    */
   public double pnorm(double p) {
      return Vector.pnorm(this, p);
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
    * subtract method subtracts the passed Vector from the calling Vector.
    * @param v a Vector object
    * @return a Vector object whose entries are the difference of the
    *         entries in the calling Vector and the respective entries 
    *         in v
    */
   public Vector subtract(Vector v) {
      return Vector.difference(this, v);
   }
   
   /**
    * sum method accepts two vectors and returns their element-wise
    * sum in a new Vector object. Assumes v1 and v2 have the same 
    * length.
    * @param u1 a Vector object
    * @param u2 a Vector object
    */
   public static Vector sum(Vector u1, Vector u2) {
      Vector.checkLengths(u1, u2);
      
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
    * Return a String containing the vector represented as a row in brackets, e.g.
    * [1.0, 2.2, 3.1, 4.9, 5.7]
    * @return a String representation of the vector
    */
   @Override
   public String toString() {
      String str = "[";
      String sep = ", ";
      
      for (int i = 0; i < this.v.length; i++) {
         str += this.v[i];
         
         if (i < (this.v.length - 1)) { // if we're not at the last entry
            str += sep;
         }
      }
      
      return str + "]";
   }
}