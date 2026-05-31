#include "softblas.h"

//  Sum of |Re| + |Im| over a complex vector (the BLAS SCASUM convention).
float32_t scasum(uint64_t N, const complex32_t *CX, int64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float32_t stemp = { SB_REAL32_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) {
        return nan_unify_s(stemp);
    }
    for (uint64_t i = 0; i < N; i++) {
        complex32_t z = CX[i * incX];
        stemp = f32_add(stemp, f32_add(f32_abs(z.real), f32_abs(z.imag)));
    }

    return nan_unify_s(stemp);
}
