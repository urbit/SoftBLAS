#include "softfloat.h"

float64_t sdsdot( const uint64_t   N,       // length
                  const float32_t  alpha,   // accumulator
                  const float32_t *X,       // vector
                  const uint64_t   incX,    // stride
                  const float32_t *Y,       // vector
                  const uint64_t   incY     // stride
                ) {
    float64_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f64_add(dot, f64_mul(float32_to_f64(*X), float32_to_f64(*Y)));
    }
    return(dot);
}
