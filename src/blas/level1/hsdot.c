#include "softblas.h"

float32_t hsdot(const uint16_t N, const float16_t alpha, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY) {
    float32_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f32_add(dot, f32_mul(f16_to_f32(*X), f16_to_f32(*Y)));
    }
    return(dot);
}
