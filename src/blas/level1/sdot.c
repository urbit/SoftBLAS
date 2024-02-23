#include "softblas.h"

float32_t sdot(const uint64_t N, const float32_t *X, const int64_t incX, const float32_t *Y, const int64_t incY) {
    float32_t sdot = { SB_REAL32_ZERO };

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        sdot = f32_add(sdot, f32_mul(X[iX], Y[iY]));
        iX += incX;
        iY += incY;
    }

    return sdot;
}
