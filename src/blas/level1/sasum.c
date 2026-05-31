#include "softblas.h"

float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float32_t stemp = { SB_REAL32_ZERO };

    if (incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) return nan_unify_s(stemp);   // reject 0/garbage stride or an overflowing one (B5) (cf. snrm2)
    for (uint64_t i = 0; i < N; i++) {
        stemp = f32_add(stemp, f32_abs(SX[i*incX]));
    }

    return nan_unify_s(stemp);
}
