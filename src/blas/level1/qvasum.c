#include "softblas.h"

//  Sum of |Re| + |Im| over a quad-complex vector (the BLAS ?asum convention).
float128_t qvasum(uint64_t N, const complex128_t *CX, int64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float128_t qtemp = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) { nan_unify_q(&qtemp); return qtemp; }
    for (uint64_t i = 0; i < N; i++) {
        complex128_t z = CX[i * incX];
        qtemp = f128_add(qtemp, f128_add(f128_abs(z.real), f128_abs(z.imag)));
    }

    nan_unify_q(&qtemp);
    return qtemp;
}
