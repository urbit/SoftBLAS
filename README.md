# SoftBLAS

A BLAS/LAPACK implmentation using [Berkeley SoftFloat](http://www.jhauser.us/arithmetic/SoftFloat.html) rather than hardware acceleration.

[BLAS naming conventions](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2023-1/naming-conventions-for-blas-routines.html) are followed for 32/64 bits, but extensions to the prefix scheme are necessary for other bit widths:

| Bits | Real | Complex |
| --- | --- | --- |
| 16 | `h` | `i` |
| 32 | `s` | `c` |
| 64 | `d` | `z` |
| 128 | `q` | `v` |
