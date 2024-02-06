#include "softblas.h"

float32_t sdsdot(const uint64_t N, const float32_t alpha, const float32_t *SX, const int64_t incX, const float32_t *SY, const int64_t incY) {
    float64_t dsdot = f32_to_f64(alpha);

    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-n + 1) * incX;
    if (incY < 0) iy = (-n + 1) * incY;
    for (uint64_t i = 0; i < n; i++) {
        dsdot = f64_add(dsdot, f64_mul(f32_to_f64(SX[ix]), f32_to_f64(SY[iy])));
        ix += incX;
        iy += incY;
    }

    return(f64_to_f32(dot));
}
