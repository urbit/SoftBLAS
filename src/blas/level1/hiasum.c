#include "softblas.h"

//  Sum of |Re| + |Im| over a complex vector (the BLAS SCASUM convention).
float16_t hiasum(uint64_t N, const complex16_t *CX, int64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float16_t stemp = { SB_REAL16_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) {
        return nan_unify_h(stemp);
    }
    for (uint64_t i = 0; i < N; i++) {
        complex16_t z = CX[i * incX];
        stemp = f16_add(stemp, f16_add(f16_abs(z.real), f16_abs(z.imag)));
    }

    return nan_unify_h(stemp);
}
