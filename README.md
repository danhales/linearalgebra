# linearalgebra
A linear algebra package in Java.

For flexibility, most methods are overloaded to include both `static` versions (which require the object as the first parameter) and instance methods. The `static` method does the heavy lifting, and the instance method is a wrapper for the `static` method using `this` as an argument (e.g. `Matrix.transpose(this)`).

No operations are carried out in-place––all operations that modify the object return a copy of the object with the modifications.

Current classes and capabilities (as of November 1, 2020)

## `Vector` class
-------------------------
Because linear algebra (as it applies to machine learning) largely revolves around lists of numbers, the Vector class is a fancy wrapper for an array of doubles, and provides common functionality for working with both that array, in addition to operations involving other Vector objects.

### Fields

* `entries` - array of doubles

### Constructors

* `Vector(double...)` - accepts either a list of doubles (e.g. `Vector(1,2,3)`), or an array of `double`s (e.g. `Vector(new double[] {1,2,3})`)

* `Vector(Vector)` - copies the entries in the passed `Vector` object's private field into a new Vector object

### Setters / Getters / Special calls to constructor

* `getEntries()` - returns a copy of (not a reference to) `entries`

* `get(int)` - returns the value in the specified position of the `Vector`, e.g. `u.get(3)` returns the entry in position 3

* `identityVector(int)` - returns the additive identity `Vector` with the specified length, which is just a list of zeros

* `set(int, double)` - sets the value at the specified index to be equal to the passed value

* `setV(double[])` - sets the values in the `v` field

* `toString()` - returns a Python-style representation of the `Vector`, e.g. `"[1, 2, 3]`"

### Unary Operations (operations on a single `Vector`)

* `inverseVector(Vector)` / `inverse()` - returns the additive inverse of the `Vector`, which is the same vector, with the signs flipped (e.g. {1, -2, 3} becomes {-1, 2, -3})

* `isCanonicalBasisVector(Vector)` - checks to see if `Vector` is all zeros, except for a single 1

* `isZero(Vector)` - checks to see if the `Vector` is essentially zero, which I've defined as smaller than Double.MIN_VALUE * 10.

* `length(Vector)` - returns the number of entries in the `Vector`, which is the length of `v`

* `magnitude(Vector)` - returns the <a href="https://en.wikipedia.org/wiki/Euclidean_space#Euclidean_norm">Euclidean norm</a> of the `Vector`

* `normalize()` - scales a `Vector` `v` with the factor `(1/v.magnitude())`. Throws an `IllegalArgumentException` if `v.isZero()` returns `true`

* `pnorm(double)` - computes the Lp norm with the given value of `p`

### Binary Operations (operations on two `Vector`s)

* `add(Vector)` / `add(Vector, Vector)` - adds corresponding entries in two `Vector` objects, provided they have the same length.

* `angleDegrees(Vector, Vector)` / `angleRadians(Vector, Vector)` - returns the angle between the two `Vector`s (arccos(u1.dot(u2) / (u1.magnitude() * u2.magnitude()))) in either degrees or radians

* `checkLengths(Vector, Vector)` checks to make sure the two `Vector` - objects have compatible lengths and throws an `IllegalArgumentException` if they have different lengths

* `cross(Vector)` / `cross(Vector, Vector)`- computes the cross product for two three-dimensional `Vector`s. The cross product is only defined for three-dimensional vectors, so an `IllegalArgumentException` will be thrown if a `Vector` any other dimension is passed.

* `dot(Vector)` / `dot(Vector, Vector)` - computes the dot product for two `Vector`s of the same length, which is the sum of products of corresponding entries, e.g. {1, 2, 3} dot {4, 5, 6} = (1)(4) + (2)(5) + (3)(6).

* `subtract(Vector)` / `subtract(Vector,Vector)` - subtracts the passed `Vector` from the calling `Vector`, or the second `Vector` from the first

### Other Operations

* `linearCombination(Vector[], double[])` - returns a weighted sum by multiplying the `Vector` objects in the `Vector` array by the weights in the `double` array.

* `multiply(Vector, double)` / `multiply(double)` - computes the scalar product by multiplying each entry in the `Vector` by the specified `double`

* `scalarTripleProduct(Vector, Vector, Vector)` - computes the scalar triple product of the three vectors, e.g. a.dot(b.cros(c)). This operation is NOT commutative.

## `Matrix` class
-------------------------
The `Matrix` class takes advantage of many of the properties of the `Vector` class in order to work with a 2D array of `double`s. Although there are many use cases possible for a `Matrix` class, I am focusing functionality in two primary areas before diving into more specific uses: machine learning and linear algebra operations. For machine learning purposes, we'll mostly be needing to multiply a `Vector` object by a `Matrix` object, pass the output through an activation function of some kind, and update the entries in the `Matrix` according to some kind of loss function. To accomplish this, we need a variety of linear algebra tricks, such as grabbing a row or column as a `Vector`, updating a single row or column, adding, subtracting, or multiplying `Matrix` objects, multiplying `Matrix` objects, and performing elementary operations such as adding a multiple of one row onto another row.

### Fields

* `entries` - a 2D `double` array. It is assumed that the array is rectangular.

### Constructors

* `Matrix(double[][])` - accepts the 2D array directly and copies its entries into the private field `entries`

* `Matrix(Matrix)` - copy constructor

### Setters / Getters / Special calls to constructor

I'm using the terms "setters" and "getters" very loosely here in order to organize this `README`, because methods like `getColumn` do not return the value in a private field––they return a subset of entries in `entries` as a `Vector` object.

* `minorMatrix(Matrix, int, int)` - returns a copy of the `Matrix` object with the indicated row and column dropped.

* `fromColumnVectors(Vector...)` - accepts either an array of `Vector` objects, or a list of `Vector` objects as parameters, and returns a `Matrix` where those `Vector` objects make up the columns. `Vector` objects must all have the same length, or an `IllegalArgumentException` will be thrown. Implemented by returning the `transpose` of the value returned by `fromRowVectors`

* `fromRowVectors(Vector...)` - accepts an array of `Vector` objects, or a list of `Vector` objects as parameters, and returns a `Matrix` where the `Vector` objects make up the rows. `Vector` objects must all have the same length, or an `IllegalArgumentException` will be thrown.

* `dropColumn(Matrix, int)` / `dropRow(Matrix, int)` - returns a copy of the `Matrix` object with the specified row/column removed. Note that this operation does change the index of the entries after the dropped column/row.

* `getColumn(Matrix,m int)` / `getRow(Matrix, int)` - accepts an `int` for the index of the desired column/row, and returns the entries in that column as a `Vector` object. If the column/row index is out of range, an `IllegalArgumentException` is thrown

* `getEntry(int, int)` - accepts an `int` for the row index and an `int` for the column index, and returns the entry at `entries[row][col]`. If either of the values are out of range, an `IllegalArgumentException` is thrown.

* `getNumColumns(Matrix)` / `getNumRows(Matrix)` - returns the length of the first row in `entries` (columns), i.e. `entries[0].length`, or the length of `entries`, i.e. `entries.length`. We assume the array is rectangular, and not ragged.

* `identityMatrix(int)` returns an n-by-n Matrix (where `n` is specified as a parameter) which has ones on the diagonal, and zeros elsewhere.

* `minorMatrix(Matrix, int, int)` - returns a copy of the matrix with the specified row (first int) and specified column (second int) dropped

* `setColumn(Matrix, int, double[])` / `setColumn(int, Vector)` / `setRow(Matrix, int, double[])` / `setRow(Matrix, int, Vector)` - accepts an int for the index of the column/row we want to set and either an array of new entries or a `Vector` of new entries, and replaces the entries in that column/row with the new entries. If the length of the new entries is not equal to the length of the column/row being set, or if the column/row index specified is out of range, an `IllegalArgumentException` will be thrown.

* `toColumnVectors(Matrix)` - returns the columns of the `Matrix` as an array of `Vector` objects

* `toRowVectors(Matrix)` - returns the rows of the `Matrix` as an array of `Vector` objects

* `toString()` - represents entries Python-style as a stack of row vectors, e.g.

`[[1, 2, 3],
  [4, -5, 6],
  [-7, 8, 9]]
`

### Unary Operations (operations on a single `Matrix`)

* `determinant(Matrix)` - recursively computes the determinant of a square matrix using the Laplace expansion, O(n!). Not recommended for use.

* `isDiagonal(Matrix)` - returns `true` if the `Matrix` is square and all nonzero entries (threshold-checked) are located on the diagonal

* `isLowerTriangular(Matrix)` - returns `true` if all values below the diagonal are zero

* `isPermutationMatrix(Matrix)` - returns `true` if the Matrix is all zeros, except for a single 1 in each row and each column

* `isSparse(Matrix)` / `isSparse(Matrix, double)` - counts the number of (threshold-checked) zero entries. No-arg version returns `true` if there are no more nonzero entries than max{number of rows, number of columns}. Version that accepts a `double` allows user to specify the proportion explicitly. As of now, this is simply a check, and no optimizations for sparse matrices are carried out.

* `isSquare(Matrix)` - checks to see if the number of rows equals the number of columns. All methods that operate only on square matrices call this method immediately and throw an `IllegalArgumentException` if it returns `false`.

* `isUpperTriangular(Matrix)` - checks to see if all entries above the diagonal are (threshold-checked) zero

* `shape(Matrix)` - returns the shape of the `Matrix` as an array of length two: {numRows, numCols}

* `swapColumns(Matrix, int, int)` - accepts two column indices, and swaps the columns at those indices. Throws an `IllegalArgumentException` if either index is out of range.

* `swapRows(Matrix, int, int)` - accepts two row indices, and swaps the rows at those indices. Throws an `IllegalArgumentException` if either index is out of range.

* `trace(Matrix)` - returns the sum of the values on the diagonal of a square Matrix

* `transpose(Matrix)` - returns the transpose of the `Matrix`, which swaps rows and columns. That is, row `i` becomes column `i` and column `j` becomes row `j`.

### Binary Operations (operations on two `Matrix` objects)

* `add(Matrix, Matrix)` - returns a new `Matrix` whose entries are sums of the corresponding entries of two `Matrix` objects.

* `multiply(Matrix, Matrix)` - returns a new `Matrix` which is the product of matrix multiplication. The number of columns in the left `Matrix` must match the number of rows in the right `Matrix`, or an `IllegalArgumentException` will be thrown. Matrix multiplication is NOT commutative, so `multiply(a, b)` will generally not equal `multiply(b, a)`, and `multiply(a, b)` being defined does not imply `multiply(b, a)` is even defined. In the instance method, the calling `Matrix` is on the left in the multiplication, i.e. `Matrix.multiply(this, m)`.

* `subtract(Matrix, Matrix)` - subtracts the second matrix from the first (or the passed matrix from the calling matrix).

### Other Operations - operations between Matrix objects and other types of objects

* `addVectorToColumn(Matrix, Vector, int)` - adds the entries in the passed `Vector` object to the specified column of the `Matrix` object. Length of `Vector` must match the length of the column or an `IllegalArgumentException` will be thrown.

* `addVectorToRow(Matrix, Vector, int)` - adds the entries in the passed `Vector` object to the specified row of the `Matrix` object. Length of `Vector` must match the length of the row or an `IllegalArgumentException` will be thrown.

* `multiply(Matrix, double)` - multiplies every entry in the `Matrix` by the passed `double`.

* `multiply(Matrix, Vector)` - multiplies the `Vector` object by the `Matrix` object and returns the resultant `Vector`.