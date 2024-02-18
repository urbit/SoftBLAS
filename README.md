# SoftBLAS

A BLAS/LAPACK implmentation using [Berkeley SoftFloat](http://www.jhauser.us/arithmetic/SoftFloat.html) rather than hardware acceleration.

Following SoftFloat 3e and requiring a 64-bit OS, all quantities are passed by value.

TODO:

- compare function signatures against reference CBLAS (`const` &c.)
- finish L2, L3

[BLAS naming conventions](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/naming-conventions-for-blas-routines.html) are followed for 32/64 bits, but extensions to the prefix scheme are necessary for other bit widths:

| Bits | Real | Complex |
| --- | --- | --- |
|  16 | `h` | `i` |
|  32 | `s` | `c` |
|  64 | `d` | `z` |
| 128 | `q` | `v` |

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

### Level 1 Functions

#### `s` 32-Bit Single Precision

- `sasum` - sum of absolute values
- `saxpy` - y = a*x + y
- `scopy` - copy x into y
- `sdot` - dot product
- `sdsdot` - dot product with extended precision accumulation (returns `float64_t`)
- `snrm2` - Euclidean norm
- `srot` - apply Givens rotation
- `srotg` - set up Givens rotation
- `srotm` - apply modified Givens rotation
- `srotmg` - set up modified Givens rotation
- `sscal` - x = a*x
- `sswap` - swap x and y
- `isamax` - index of max abs value

#### `d` 64-Bit Double Precision

- `dasum` - sum of absolute values
- `daxpy` - y = a*x + y
- `dcopy` - copy x into y
- `ddot` - dot product
- `dnrm2` - Euclidean norm
- `dsdot` - dot product with extended precision accumulation (returns `float64_t`)
- `drot` - apply Givens rotation
- `drotg` - set up Givens rotation
- `drotm` - apply modified Givens rotation
- `drotmg` - set up modified Givens rotation
- `dscal` - x = a*x
- `dswap` - swap x and y
- `idamax` - index of max abs value

#### `h` 16-Bit Half Precision

- `hasum` - sum of absolute values
- `haxpy` - y = a*x + y
- `hcopy` - copy x into y
- `hdot` - dot product
- `hnrm2` - Euclidean norm
- `hsdot` - dot product with extended precision accumulation (returns `float32_t`)
- `hrot` - apply Givens rotation
- `hrotg` - set up Givens rotation
- `hrotm` - apply modified Givens rotation
- `hrotmg` - set up modified Givens rotation
- `hscal` - x = a*x
- `hswap` - swap x and y
- `ihamax` - index of max abs value

#### Pending

- `{hi|sc|dz|qv}nrm` - Euclidean norm, complex
- `q*` - quadruple-precision functions

### Level 2 Functions

#### `s` 32-Bit Single Precision

- `sgemv` - computes a matrix-vector product using a general matrix

#### `d` 64-Bit Single Precision

- `dgemv` - computes a matrix-vector product using a general matrix

#### `h` 16-Bit Single Precision

- `hgemv` - computes a matrix-vector product using a general matrix

### Level 3 Functions

#### `s` 32-Bit Single Precision

- `sgemm` - computes a matrix-matrix product using a general matrix

#### `d` 64-Bit Single Precision

- `dgemm` - computes a matrix-matrix product using a general matrix

#### `h` 16-Bit Single Precision

- `hgemm` - computes a matrix-matrix product using a general matrix

TODO:  `gemm3m` variants

### Auxiliary Macros

- `f16_ge` → `!f16_lt`
- `f16_gt` → `!f16_le`
- `f16_ne` → `!f16_eq`
- `f32_ge` → `!f32_lt`
- `f32_gt` → `!f32_le`
- `f32_ne` → `!f32_eq`
- `f64_ge` → `!f64_lt`
- `f64_gt` → `!f64_le`
- `f64_ne` → `!f64_eq`
- `f128_ge` → `!f128_lt`
- `f128_gt` → `!f128_le`
- `f128_ne` → `!f128_eq`

- `ABS(x)` → `( (x) >= 0 ? (x) : -(x) )`
- `f16_abs(x)` → `( (x) & 0x7fff )`
- `f32_abs(x)` → `( (x) & 0x7fffffff )`
- `f64_abs(x)` → `( (x) & 0x7fffffffffffffff )`

- `MAX(x, y)` → `( (x) > (y) ? (x) : (y) )`
- `f16_max(x, y)` → `( f16_gt( (x) , (y) ) ? (x) : (y) )`
- `f32_max(x, y)` → `( f32_gt( (x) , (y) ) ? (x) : (y) )`
- `f64_max(x, y)` → `( f64_gt( (x) , (y) ) ? (x) : (y) )`

- `MIN(x, y)` → `( (x) > (y) ? (y) : (x) )`
- `f16_min(x, y)` → `( f16_gt( (x) , (y) ) ? (y) : (x) )`
- `f32_min(x, y)` → `( f32_gt( (x) , (y) ) ? (y) : (x) )`
- `f64_min(x, y)` → `( f64_gt( (x) , (y) ) ? (y) : (x) )`

---

- sasum simplified logic
- saxpy simplified logic, kept neg behavior
- scopy simplified logic, kept neg behavior
- sdot simplified logic, kept neg behavior

##  Testing

To run the test suite, you need a copy of the source of [SoftFloat 3e](http://www.jhauser.us/arithmetic/SoftFloat-3/doc/SoftFloat.html) and [µnit](https://nemequ.github.io/munit/) in the `SoftBLAS/` directory.

```sh
gcc -o test_sasum -I.. -I../../../include  -I../../../SoftFloat/source/include -I../../../munit sasum.c ../../../src/blas/level1/sasum.c ../../../munit/munit.c
```
