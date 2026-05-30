#include "softblas.h"

//  Apply a plane rotation: [x_i; y_i] <- [c s; -s c] [x_i; y_i].
void hrot(const uint64_t N, float16_t *X, const uint64_t incX, float16_t *Y, const uint64_t incY, const float16_t c, const float16_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        float16_t tmp  = f16_add(f16_mul(c, X[iX]), f16_mul(s, Y[iY]));
        Y[iY] = nan_unify_h(f16_sub(f16_mul(c, Y[iY]), f16_mul(s, X[iX])));
        X[iX] = nan_unify_h(tmp);
        iX += incX;
        iY += incY;
    }
}
