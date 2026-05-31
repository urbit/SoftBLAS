#include "softblas.h"

float64_t dasum(uint64_t N, const float64_t *DX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t dtemp = { SB_REAL64_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) return nan_unify_d(dtemp);   // reject 0/garbage stride or an overflowing one (B5) (cf. dnrm2)
    for (uint64_t i = 0; i < N; i++) {
        dtemp = f64_add(dtemp, f64_abs(DX[i*incX]));
    }

    return nan_unify_d(dtemp);
}