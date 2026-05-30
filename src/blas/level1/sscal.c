#include "softblas.h"

void sscal(uint64_t N, float32_t SA, float32_t *SX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        SX[iX] = nan_unify_s(f32_mul(SA, SX[iX]));
        iX += incX;
    }
}
