#include "softblas.h"

float64_t ddot(const uint64_t N, const float64_t *X, const int64_t incX, const float64_t *Y, const int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t ddot = { SB_REAL64_ZERO };

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        ddot = f64_add(ddot, f64_mul(X[iX], Y[iY]));
        iX += incX;
        iY += incY;
    }

    return ddot;
}
