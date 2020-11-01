# linearalgebra
A linear algebra package in Java.

For flexibility, most methods are overloaded to include both static versions (which require the object as the first parameter) and instance methods (which call the static method with the `this` keyword).

Current classes and capabilities (as of November 1, 2020)

Vector
------
Because linear algebra (as it applies to machine learning) largely revolves around lists of numbers, the Vector class is a fancy wrapper for an array of doubles, and provides common functionality for working with both that array, in addition to operations involving other Vector objects.

## Fields

### `v` - array of doubles

## Constructors

### `Vector(double...)` accepts either a list of doubles (e.g. Vector(1,2,3)), or an array of doubles (e.g. Vector(new double[] {1,2,3}))

### `Vector(Vector)` copies the entries in the passed Vector object's private field into a new Vector object

## Setters / Getters / Special calls to constructor

### `getV()` returns a copy of `v` (not a reference to v)

### `get(int)` returns the value in the specified position of the `Vector`, e.g. `u.get(3)` returns the entry in position 3

### `identityVector(int)` returns the additive identity `Vector` with the specified length, which is just a list of zeros

### `set(int, double)` sets the value at the specified index to be equal to the passed value

### `setV(double[])` sets the values in the `v` field

### `toString()` returns a Python-style representation of the `Vector`, e.g. `"[1, 2, 3]`"

## Unary Operations (operations on a single `Vector`)

### `inverseVector(Vector)` / `inverse()` returns the additive inverse of the vector, which is the same vector, with the signs flipped (e.g. {1, -2, 3} becomes {-1, 2, -3})

### `isZero()` checks to see if the `Vector` is essentially zero, which I've defined as smaller than Double.MIN_VALUE * 10.

### `length()` returns the number of entries in the `Vector`, which is the length of `v`

### `magnitude()` returns the <a href="https://en.wikipedia.org/wiki/Euclidean_space#Euclidean_norm">Euclidean norm</a> of the `Vector`

### `normalize()` scales a `Vector` `v` with the factor `(1/v.magnitude())`. Throws an `IllegalArgumentException` if `v.isZero()` returns `true`

### `pnorm(double)` computes the Lp norm with the given value of `p`

## Binary Operations (operations on two `Vector`s)

### `add(Vector)` / `sum(Vector,Vector)` adds corresponding entries in two `Vector` objects, provided they have the same length. For readability, I've verbified the static-method `sum` into the instnace method `add`, so to compute the sum of vectors, you call `Vector.sum(v1, v2)`, but to add `v2` onto `v1`, you simply call `v1.add(v2)`.

### `angleDegrees(Vector, Vector)` / `angleRadians(Vector, Vector)` returns the angle between the two `Vector`s (arccos(u1.dot(u2) / (u1.magnitude() * u2.magnitude()))) in either degrees or radians

### `checkLengths(Vector, Vector)` checks to make sure the two `Vector` objects have compatible lengths and throws an `IllegalArgumentException` if they have different lengths

### `cross` / `crossProduct(Vector, Vector)` computes the cross product for two three-dimensional `Vector`s. The cross product is only defined for three-dimensional vectors, so an `IllegalArgumentException` will be thrown if a `Vector` any other dimension is passed. Same naming logic as `add` / `sum`. This operation is NOT commutative.

### `dot` / `dotProduct` computes the dot product for two `Vector`s of the same length, which is the sum of products of corresponding entries, e.g. {1, 2, 3} dot {4, 5, 6} = (1)(4) + (2)(5) + (3)(6). Same naming logic as `add` / `sum`.

### `subtract(Vector)` / `difference(Vector,Vector)` subtracts the second `Vector` from the first. Same naming logic as `add` / `sum`.

## Other Operations

### `linearCombination(Vector[], double[])` returns a weighted sum by multiplying the `Vector` objects in the `Vector` array by the weights in the `double` array.

### `product(Vector, double)` / `multiply(double)` computes the scalar product by multiplying each entry in the `Vector` by the specified `double`

### `scalarTripleProduct(Vector, Vector, Vector)` computes the scalar triple product of the three vectors, e.g. a.dot(b.cros(c)). This operation is NOT commutative.

-- Notes on `Matrix.java` coming soon!
