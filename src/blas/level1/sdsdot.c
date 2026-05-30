#include "softblas.h"

//  sdsdot: dot product of single-precision vectors accumulated in double
//  precision, plus the single-precision bias alpha. Returns single precision.
float32_t sdsdot(const uint64_t N, const float32_t alpha, const float32_t *SX, const int64_t incX, const float32_t *SY, const int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t dsdot = f32_to_f64(alpha);
    int64_t ix = 0, iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        dsdot = f64_add(dsdot, f64_mul(f32_to_f64(SX[ix]), f32_to_f64(SY[iy])));
        ix += incX;
        iy += incY;
    }
    return nan_unify_s(f64_to_f32(dsdot));
}
