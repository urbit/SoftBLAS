# SoftBLAS

A BLAS/LAPACK implmentation using [Berkeley SoftFloat](http://www.jhauser.us/arithmetic/SoftFloat.html) rather than hardware acceleration.

[BLAS naming conventions](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/naming-conventions-for-blas-routines.html) are followed for 32/64 bits, but extensions to the prefix scheme are necessary for other bit widths:

| Bits | Real | Complex |
| --- | --- | --- |
| 16 | `h` | `i` |
| 32 | `s` | `c` |
| 64 | `d` | `z` |
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
