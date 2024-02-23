#include "softblas.h"

float16_t hdot(const uint64_t N, const float16_t *X, const int64_t incX, const float16_t *Y, const int64_t incY) {
    float16_t hdot = { SB_REAL16_ZERO };

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        hdot = f16_add(hdot, f16_mul(X[iX], Y[iY]));
        iX += incX;
        iY += incY;
    }

    return hdot;
}