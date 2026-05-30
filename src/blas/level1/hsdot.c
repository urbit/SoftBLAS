#include "softblas.h"

//  hsdot: dot product of half-precision vectors accumulated in single
//  precision, plus the single-precision bias (alpha widened from half).
float32_t hsdot(const uint64_t N, const float16_t alpha, const float16_t *X, const int64_t incX, const float16_t *Y, const int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float32_t dot = f16_to_f32(alpha);
    int64_t ix = 0, iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        dot = f32_add(dot, f32_mul(f16_to_f32(X[ix]), f16_to_f32(Y[iy])));
        ix += incX;
        iy += incY;
    }
    return nan_unify_s(dot);
}
