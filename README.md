#   SoftBLAS

A BLAS/LAPACK implementation using [Berkeley SoftFloat](http://www.jhauser.us/arithmetic/SoftFloat.html) rather than hardware acceleration.

![](./img/logo.jpg)

Following SoftFloat 3e and requiring a 64-bit OS, all quantities are passed by value. For relevant extra-reading on floating point arithmetic and information about how SoftBLAS is useful to Urbit, check out [~lagrev-nocfep's USTJ Volume 1 article.](https://urbitsystems.tech/article/v01-i01/the-desert-of-the-reals-floating-point-arithmetic-on-deterministic-systems). If you're interested in finding out more / contributing, join the Numerics Tlon group: `~mopfel-winrux/numeric-computation-and-machine-learning`.

**Status ~2024.5.24:  `REAL`-valued operations are “complete” (BLAS is a pseudo-standard).**

##  Releases

- 1.0.0 (commit `7d05697aea5363dcf5f877a9c8b464e9c352d3d4`).  Basic suite of `REAL`-valued operations suitable for use with half, single, double, and quadruple precision floating-point numbers.
- 1.1.0 (commit `f94acbcfd26cebd8e135ad9e8c7caa156fcc4ac9`).  Errors changed to `exit(-1)` instead of `return`.

### Planned work

- [ ] Complete complex-valued functions (in progress).
- [ ] Use a linter.
- [ ] Add (kelvin) versioning, at least on interface.

##  BLAS

[BLAS naming conventions](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/naming-conventions-for-blas-routines.html) are followed for 32/64 bits, but extensions to the prefix scheme are necessary for 16/128 bit widths; we use:

| Bit Width | Real | Real Form | Bit Width | Complex | Complex Form |
| --- | --- | --- | --- | --- | --- |
|  16 | `h` |  `float16_t` |  16 | `i` |  `complex16_t` |
|  32 | `s` |  `float32_t` |  32 | `c` |  `complex32_t` |
|  64 | `d` |  `float64_t` |  64 | `z` |  `complex64_t` |
| 128 | `q` | `float128_t` | 128 | `v` | `complex128_t` |

The rounding mode should be set via the global variable `softblas_roundingMode` (a `char` `typedef`), defined in `softblas.h`.  Valid values are:

- `'n'` - round to nearest (even)
- `'z'` - round towards zero
- `'u'` - round up
- `'d'` - round down
- `'a'` - round away from zero

Per Wikipedia:

> BLAS functionality is categorized into three sets of routines called "levels", which correspond to both the chronological order of definition and publication, as well as the degree of the polynomial in the complexities of algorithms; Level 1 BLAS operations typically take linear time, $O(n)$, Level 2 operations quadratic time and Level 3 operations cubic time.  Modern BLAS implementations typically provide all three levels.
>
> ### Level 1
> This level consists of all the routines described in the original presentation of BLAS (1979), which defined only ''vector operations'' on strided arrays: dot products, vector norms, a generalized vector addition of the form
>
> $$
> \boldsymbol{y} \leftarrow \alpha \boldsymbol{x} + \boldsymbol{y}
> $$
>
> (called `axpy`, "`a` `x` plus `y`") and several other operations.
>
> ### Level 2
> This level contains matrix-vector operations including, among other things, a *ge*neralized *m*atrix-*v*ector multiplication (`gemv`):
>
> $$
> \boldsymbol{y} \leftarrow \alpha \boldsymbol{A} \boldsymbol{x} + \beta \boldsymbol{y}
> $$
>
> as well as a solver for $\boldsymbol{x}$ in the linear equation
>
> $$
> \boldsymbol{T} \boldsymbol{x} = \boldsymbol{y}
> $$
>
> with $\boldsymbol{T}$ being triangular. … The Level 2 subroutines are especially intended to improve performance of programs using BLAS on vector processors, where Level 1 BLAS are suboptimal "because they hide the matrix-vector nature of the operations from the compiler."
>
> ### Level 3
> This level … contains ''matrix-matrix operations'', including a "general matrix multiplication" (`gemm`), of the form
>
> $$
> \boldsymbol{C} \leftarrow \alpha \boldsymbol{A} \boldsymbol{B} + \beta \boldsymbol{C}\textrm{,}
> $$
>
> where $\boldsymbol{A}$ and $\boldsymbol{B}$ can optionally be transposed or hermitian-conjugated inside the routine, and all three matrices may be strided. The ordinary matrix multiplication $\boldsymbol{A B}$ can be performed by setting $α$ to one and $\boldsymbol{C}$ to an all-zeros matrix of the appropriate size.
>
> Also included in Level 3 are routines for computing
>
> $$
> \boldsymbol{B} \leftarrow \alpha \boldsymbol{T}^{-1} \boldsymbol{B}\textrm{,}
> $$
>
> where $\boldsymbol{T}$ is a triangular matrix, among other functionality.

![](./img/logo2.jpg)

### Level 1 Functions

- `asum` - ||**x**||₁, where **x** is a vector (sum of absolute values)
    - Implemented: `h`, `s`, `d`, `q`, `c`
    - Yet to implement: `i`, `z`, `v`
- `axpy` - **y** ← a**x** + **y**, where a is a scalar and **x**, **y** are vectors
    - Implemented: `h`, `s`, `d`, `q`, `c`
    - Yet to implement: `i`, `z`, `v`
- `copy` - **y** ← **x**, where **x**, **y** are vectors
    - Implemented: `h`, `s`, `d`, `q`, `c`
    - Yet to implement: `i`, `z`, `v`
- `dot` - **x**•**y**, where **x**, **y** are vectors
    - Implemented: `h`, `s`, `d`, `q`
    - Yet to implement: `i`, `c`, `z`, `v`
- `dotc` - <**x**, **y**> = **x***•**y**, where **x**, **y** are complex vectors (complex inner product)
    - Implemented: `c`
    - Yet to implement: `i`, `z`, `v`
- `nrm2` - ||**x**||₂, where **x** is a vector (Pythagorean distance)
    - Implemented: `h`, `s`, `d`, `q`, `c`
    - Yet to implement: `i`, `z`, `v`
- `scal` - **x** ← a**x**, where a is a scalar and **x** is a vector
    - Implemented: `h`, `s`, `d`, `q`
    - Yet to implement: `i`, `c`, `z`, `v`
- `swap` - (**x**, **y**) ← (**y**, **x**), where **x** and **y** are vectors
    - Implemented: `h`, `s`, `d`, `q`
    - Yet to implement: `i`, `c`, `z`, `v`
- `iamax` - argmaxᵢ(|**x**ᵢ|), where **x** is a real vector (index of max absolute value)
    - Implemented: `h`, `s`, `d`, `q`

#### Extras:
- `sdot` - a+**x**•**y**, where a is a scalar and **x**, **y** are vectors and intermediate calculations use a higher precision float
    - Implemented: `h`, `s`
- `rot` - (**x**ᵢ // **y**ᵢ) ← R(**x**ᵢ // **y**ᵢ), where R is a rotation matrix (c s // -s c), and **x**, **y** are vectors
    - Implemented: `h`, `s`, `d`, `c`
- `rotg` - Finds the Rotation matrix R (among other things), such that R(x // y) = (r // 0), where x, y are scalars
    - Implemented: `h`, `s`, `d`
- `rotm` - Performs a modified Givens rotation
    - Implemented: `h`, `s`, `d`
- `rotmg` - Set up a modified Givens rotation
    - Implemented: `h`, `s`, `d`

### Level 2 Functions

- `gemv` - **y** ← A**x**+b**y**, where b is a scalar, **x**, **y** are vectors and A is matrix (possibly transposed)
    - Implemented: `h`, `s`, `d`, `q`
    - Yet to implement: `i`, `c`, `z`, `v`

### Level 3 Functions

- `gemm` - C ← aAB+bC, where a, b are scalars, and A, B, C are matrices (possibly transposed).
    - Implemented: `h`, `s`, `d`, `q`
    - Yet to implement: `i`, `c`, `z`, `v`

#### Pending

- `?gemm3m` accelerated variants which calculate gemm for complex matrices with [25% fewer multiplications.](https://eprints.maths.manchester.ac.uk/348/1/0613043.pdf)

![](./img/logo3.jpg)

### Auxiliary Macros

These are provided as convenient extensions of SoftFloat definitions.

- `f16_ge(x,y)` → `( !(f16_lt( x,y ) ))
- `f16_gt(x,y)` → `( !(f16_le( x,y ) ))
- `f16_ne(x,y)` → `( !(f16_eq( x,y ) ))
- `f32_ge(x,y)` → `( !(f32_lt( x,y ) ))
- `f32_gt(x,y)` → `( !(f32_le( x,y ) ))
- `f32_ne(x,y)` → `( !(f32_eq( x,y ) ))
- `f64_ge(x,y)` → `( !(f64_lt( x,y ) ))
- `f64_gt(x,y)` → `( !(f64_le( x,y ) ))
- `f64_ne(x,y)` → `( !(f64_eq( x,y ) ))`
- `f128M_ge(x,y)` → `( !(f128M_lt( x,y ) ))`
- `f128M_gt(x,y)` → `( !(f128M_le( x,y ) ))`
- `f128M_ne(x,y)` → `( !(f128M_eq( x,y ) ))`

- `ABS(x)` → `( (x) >= 0 ? (x) : -(x) )`
- `f16_abs(x)` → `(float16_t){ ( (uint16_t)(x.v) & 0x7fff ) }`
- `f32_abs(x)` → `(float32_t){ ( (uint32_t)(x.v) & 0x7fffffff ) }`
- `f64_abs(x)` → `(float64_t){ ( (uint64_t)(x.v) & 0x7fffffffffffffff ) }`
- `f128_abs(x)` → `(float128_t){ ( (uint64_t)(x.v[0]) & 0x7fffffffffffffff ), ( (uint64_t)(x.v[1]) & 0x7fffffffffffffff ) }`
- `f16M_abs(x)` → `(float16_t*){ ( (uint16_t)(x->v) & 0x7fff ) }`
- `f32M_abs(x)` → `(float32_t*){ ( (uint32_t)(x->v) & 0x7fffffff ) }`
- `f64M_abs(x)` → `(float64_t*){ ( (uint64_t)(x->v) & 0x7fffffffffffffff ) }`
- `f128M_abs(x)` → `(float128_t*){ ( (uint64_t)(x->v[0]) & 0x7fffffffffffffff ), ( (uint64_t)(x->v[1]) & 0x7fffffffffffffff ) }`

- `MAX(x, y)` → `( (x) > (y) ? (x) : (y) )`
- `f16_max(x, y)` → `( f16_gt( (x) , (y) ) ? (x) : (y) )`
- `f32_max(x, y)` → `( f32_gt( (x) , (y) ) ? (x) : (y) )`
- `f64_max(x, y)` → `( f64_gt( (x) , (y) ) ? (x) : (y) )`
- `f128M_max(x, y)` → `( f128M_gt( (x) , (y) ) ? (x) : (y) )`

- `MIN(x, y)` → `( (x) > (y) ? (y) : (x) )`
- `f16_min(x, y)` → `( f16_gt( (x) , (y) ) ? (y) : (x) )`
- `f32_min(x, y)` → `( f32_gt( (x) , (y) ) ? (y) : (x) )`
- `f64_min(x, y)` → `( f64_gt( (x) , (y) ) ? (y) : (x) )`
- `f128_min(x, y)` → `( f128_gt( (x) , (y) ) ? (y) : (x) )`

- `f16_ceil(a)` → `f16_roundToInt( a, softfloat_round_max, false )`
- `f32_ceil(a)` → `f32_roundToInt( a, softfloat_round_max, false )`
- `f64_ceil(a)` → `f64_roundToInt( a, softfloat_round_max, false )`
- `f128M_ceil(a, b)` → `f128M_roundToInt( a, softfloat_round_max, false, b )`

- `f16_floor(a)` → `f16_roundToInt( a, softfloat_round_min, false )`
- `f32_floor(a)` → `f32_roundToInt( a, softfloat_round_min, false )`
- `f64_floor(a)` → `f64_roundToInt( a, softfloat_round_min, false )`
- `f128M_floor(a, b)` → `f128M_roundToInt( a, softfloat_round_min, false, b )`


##  References

The following were found to be particularly helpful in composing these functions.

- [BLAS](https://www.netlib.org/blas/)
- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS)
- [ATLAS](https://math-atlas.sourceforge.net/faq.html) LAPACK routimes
- [Intel MKL](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/)
- [IBM Engineering and Scientific Subroutines](https://www.ibm.com/docs/en/essl/6.3?topic=subroutines-matrix-operation)


##  Testing

![](./img/logo4.jpg)

To run the test suite, you need to make sure the SoftFloat and Munit Git submodules are properly fetched.
```sh
$ git submodule init
$ git submodule update
```

Then build the Linux x86_64 SoftFloat version:
```sh
$ cd SoftFloat/build/Linux-x86_64-GCC
$ make
```

(If your setup is different, build the relevant SoftFloat version and modify `LDFLAGS` in the Makefile.)

```sh
$ make tests
$ ./test_all
Running test suite with seed 0xa623450c...
/blas/level1/test_sasum_0            [ OK    ] [ 0.00000478 / 0.00000284 CPU ]
/blas/level1/test_sasum_12345        [ OK    ] [ 0.00000800 / 0.00000671 CPU ]
/blas/level1/test_saxpy_0            [ OK    ] [ 0.00002413 / 0.00002288 CPU ]
/blas/level1/test_saxpy_sum          [ OK    ] [ 0.00002411 / 0.00002285 CPU ]
/blas/level1/test_saxpy_stride       [ OK    ] [ 0.00002446 / 0.00002332 CPU ]
/blas/level1/test_saxpy_neg_stride   [ OK    ] [ 0.00002927 / 0.00002805 CPU ]
* * *
132 of 132 (100%) tests successful, 0 (0%) test skipped.
```


##  Versioning

When the initial complete release is complete, SoftBLAS will be [kelvin versioned](https://docs.urbit.org/glossary/kelvin).  We will start relatively hot, at room temperature or so.
