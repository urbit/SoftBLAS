#include "softblas.h"

float32_t sdot(const uint64_t N, const float32_t *X, const int64_t incX, const float32_t *Y, const int64_t incY) {
    float32_t sdot = SB_REAL32_ZERO;

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-n + 1) * incX;
    if (incY < 0) iy = (-n + 1) * incY;
    for (uint64_t i = 0; i < n; i++) {
        sdot = f32_add(sdot, f32_mul(sx[ix], sy[iy]));
        ix += incX;
        iy += incY;
    }

    return sdot;
}
