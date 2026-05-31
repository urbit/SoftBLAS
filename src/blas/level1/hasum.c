#include "softblas.h"

float16_t hasum(uint64_t N, const float16_t *HX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float16_t htemp = { SB_REAL16_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) return nan_unify_h(htemp);   // reject 0/garbage stride or an overflowing one (B5) (cf. hnrm2)
    for (uint64_t i = 0; i < N; i++) {
        htemp = f16_add(htemp, f16_abs(HX[i*incX]));
    }

    return nan_unify_h(htemp);
}
