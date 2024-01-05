#include "softblas.h"

float64_t dsdot(const uint64_t N, const float64_t alpha, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY) {
    float64_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f64_add(dot, f64_mul(*X, *Y));
    }
    return dot;
}
