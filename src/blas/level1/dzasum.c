#include "softblas.h"

//  Sum of |Re| + |Im| over a complex vector (the BLAS SCASUM convention).
float64_t dzasum(uint64_t N, const complex64_t *CX, int64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t stemp = { SB_REAL64_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) {
        return nan_unify_d(stemp);
    }
    for (uint64_t i = 0; i < N; i++) {
        complex64_t z = CX[i * incX];
        stemp = f64_add(stemp, f64_add(f64_abs(z.real), f64_abs(z.imag)));
    }

    return nan_unify_d(stemp);
}
