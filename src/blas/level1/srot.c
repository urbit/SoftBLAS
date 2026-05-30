#include "softblas.h"

//  Apply a plane rotation: [x_i; y_i] <- [c s; -s c] [x_i; y_i].
void srot(const uint64_t N, float32_t *X, const uint64_t incX, float32_t *Y, const uint64_t incY, const float32_t c, const float32_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        float32_t tmp  = f32_add(f32_mul(c, X[iX]), f32_mul(s, Y[iY]));
        Y[iY] = nan_unify_s(f32_sub(f32_mul(c, Y[iY]), f32_mul(s, X[iX])));
        X[iX] = nan_unify_s(tmp);
        iX += incX;
        iY += incY;
    }
}
