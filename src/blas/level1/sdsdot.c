#include "softblas.h"

float64_t sdsdot(const uint32_t N, const float32_t alpha, const float32_t *X, const uint32_t incX, const float32_t *Y, const uint32_t incY) {
    float64_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f64_add(dot, f64_mul(f32_to_f64(*X), f32_to_f64(*Y)));
    }
    return(dot);
}
