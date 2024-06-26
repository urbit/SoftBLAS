#   SoftBLAS

A BLAS/LAPACK implementation using [Berkeley SoftFloat](http://www.jhauser.us/arithmetic/SoftFloat.html) rather than hardware acceleration.

![](./img/logo.jpg)

Following SoftFloat 3e and requiring a 64-bit OS, all quantities are passed by value.

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

#### `s` 32-Bit Single Precision

- `sasum` - sum of absolute values
- `saxpy` - y = a*x + y
- `scopy` - copy x into y
- `sdot` - dot product
- ~~`sdsdot` - dot product with extended precision accumulation (returns `float64_t`)~~
- `snrm2` - Euclidean norm
- ~~`srot` - apply Givens rotation~~
- ~~`srotg` - set up Givens rotation~~
- ~~`srotm` - apply modified Givens rotation~~
- ~~`srotmg` - set up modified Givens rotation~~
- `sscal` - x = a*x
- `sswap` - swap x and y
- `isamax` - index of max abs value

#### `d` 64-Bit Double Precision

- `dasum` - sum of absolute values
- `daxpy` - y = a*x + y
- `dcopy` - copy x into y
- `ddot` - dot product
- ~~`dsdot` - dot product with extended precision accumulation (returns `float64_t`)~~
- `dnrm2` - Euclidean norm
- ~~`drot` - apply Givens rotation~~
- ~~`drotg` - set up Givens rotation~~
- ~~`drotm` - apply modified Givens rotation~~
- ~~`drotmg` - set up modified Givens rotation~~
- `dscal` - x = a*x
- `dswap` - swap x and y
- `idamax` - index of max abs value

#### `h` 16-Bit Half Precision

- `hasum` - sum of absolute values
- `haxpy` - y = a*x + y
- `hcopy` - copy x into y
- `hdot` - dot product
- `hnrm2` - Euclidean norm
- ~~`hrot` - apply Givens rotation~~
- ~~`hrotg` - set up Givens rotation~~
- ~~`hrotm` - apply modified Givens rotation~~
- ~~`hrotmg` - set up modified Givens rotation~~
- `hscal` - x = a*x
- `hswap` - swap x and y
- `ihamax` - index of max abs value

#### `q` 128-Bit Quadruple Precision

- `qasum` - sum of absolute values
- `qaxpy` - y = a*x + y
- `qcopy` - copy x into y
- `qdot` - dot product
- `qnrm2` - Euclidean norm
- ~~`qrot` - apply Givens rotation~~
- ~~`qrotg` - set up Givens rotation~~
- ~~`qrotm` - apply modified Givens rotation~~
- ~~`qrotmg` - set up modified Givens rotation~~
- `qscal` - x = a*x
- `qswap` - swap x and y
- `iqamax` - index of max abs value

#### `c` 32-Bit Single Precision (Complex)

- `casum` - sum of absolute values
- `caxpy` - y = a*x + y
- `ccopy` - copy x into y
- `cdot` - dot product
- ~~`cdsdot` - dot product with extended precision accumulation (returns `float64_t`)~~
- `cnrm2` - Euclidean norm
- ~~`crot` - apply Givens rotation~~
- ~~`crotg` - set up Givens rotation~~
- ~~`crotm` - apply modified Givens rotation~~
- ~~`crotmg` - set up modified Givens rotation~~
- `cscal` - x = a*x
- `cswap` - swap x and y

#### `z` 64-Bit Double Precision (Complex)

- `zasum` - sum of absolute values
- `zaxpy` - y = a*x + y
- `zcopy` - copy x into y
- `zdot` - dot product
- ~~`zsdot` - dot product with extended precision accumulation (returns `float64_t`)~~
- `znrm2` - Euclidean norm
- ~~`zrot` - apply Givens rotation~~
- ~~`zrotg` - set up Givens rotation~~
- ~~`zrotm` - apply modified Givens rotation~~
- ~~`zrotmg` - set up modified Givens rotation~~
- `zscal` - x = a*x
- `zswap` - swap x and y

#### `i` 16-Bit Half Precision (Complex)

- `iasum` - sum of absolute values
- `iaxpy` - y = a*x + y
- `icopy` - copy x into y
- `idot` - dot product
- `inrm2` - Euclidean norm
- ~~`irot` - apply Givens rotation~~
- ~~`irotg` - set up Givens rotation~~
- ~~`irotm` - apply modified Givens rotation~~
- ~~`irotmg` - set up modified Givens rotation~~
- `iscal` - x = a*x
- `iswap` - swap x and y

#### `v` 128-Bit Quadruple Precision (Complex)

- `vasum` - sum of absolute values
- `vaxpy` - y = a*x + y
- `vcopy` - copy x into y
- `vdot` - dot product
- `vnrm2` - Euclidean norm
- ~~`vrot` - apply Givens rotation~~
- ~~`vrotg` - set up Givens rotation~~
- ~~`vrotm` - apply modified Givens rotation~~
- ~~`vrotmg` - set up modified Givens rotation~~
- `vscal` - x = a*x
- `vswap` - swap x and y

#### Pending

- `?rot??` Givens rotation functions in Level 1
- Various OpenBLAS extended functions like `i?amin`

### Level 2 Functions

#### `s` 32-Bit Single Precision

- `sgemv` - computes a matrix-vector product using a general matrix

#### `d` 64-Bit Single Precision

- `dgemv` - computes a matrix-vector product using a general matrix

#### `h` 16-Bit Single Precision

- `hgemv` - computes a matrix-vector product using a general matrix

#### `q` 128-Bit Single Precision

- `qgemv` - computes a matrix-vector product using a general matrix

#### `s` 32-Bit Single Precision (Complex)

- `cgemv` - computes a matrix-vector product using a general matrix

#### `d` 64-Bit Single Precision (Complex)

- `zgemv` - computes a matrix-vector product using a general matrix

#### `h` 16-Bit Single Precision (Complex)

- `igemv` - computes a matrix-vector product using a general matrix

#### `q` 128-Bit Single Precision (Complex)

- `vgemv` - computes a matrix-vector product using a general matrix

### Level 3 Functions

#### `s` 32-Bit Single Precision

- `sgemm` - computes a matrix-matrix product using a general matrix

#### `d` 64-Bit Single Precision

- `dgemm` - computes a matrix-matrix product using a general matrix

#### `h` 16-Bit Single Precision

- `hgemm` - computes a matrix-matrix product using a general matrix

#### `h` 128-Bit Single Precision

- `qgemm` - computes a matrix-matrix product using a general matrix

#### `s` 32-Bit Single Precision (Complex)

- `cgemm` - computes a matrix-matrix product using a general matrix

#### `d` 64-Bit Single Precision (Complex)

- `zgemm` - computes a matrix-matrix product using a general matrix

#### `h` 16-Bit Single Precision (Complex)

- `igemm` - computes a matrix-matrix product using a general matrix

#### `h` 128-Bit Single Precision (Complex)

- `vgemm` - computes a matrix-matrix product using a general matrix

#### Pending

- `?gemm3m` accelerated variants

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


##  References

The following were found to be particularly helpful in composing these functions.

- [BLAS](https://www.netlib.org/blas/)
- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS)
- [ATLAS](https://math-atlas.sourceforge.net/faq.html) LAPACK routimes
- [Intel MKL](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/)
- [IBM Engineering and Scientific Subroutines](https://www.ibm.com/docs/en/essl/6.3?topic=subroutines-matrix-operation)


##  Testing

![](./img/logo4.jpg)

To run the test suite, you need a build of [SoftFloat 3e](http://www.jhauser.us/arithmetic/SoftFloat-3/doc/SoftFloat.html) and [µnit](https://nemequ.github.io/munit/) in the `SoftBLAS/` directory.

- SoftFloat 3e at `SoftBLAS/SoftFloat/build/Linux-x86_64-GCC/`
- µnit at `SoftBLAS/munit/`

If your setup is different, modify `LDFLAGS` in the Makefile.

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
