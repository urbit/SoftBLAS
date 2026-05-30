#include "softblas.h"

//  Apply a plane rotation: [x_i; y_i] <- [c s; -s c] [x_i; y_i].
void drot(const uint64_t N, float64_t *X, const uint64_t incX, float64_t *Y, const uint64_t incY, const float64_t c, const float64_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        float64_t tmp  = f64_add(f64_mul(c, X[iX]), f64_mul(s, Y[iY]));
        Y[iY] = nan_unify_d(f64_sub(f64_mul(c, Y[iY]), f64_mul(s, X[iX])));
        X[iX] = nan_unify_d(tmp);
        iX += incX;
        iY += incY;
    }
}
