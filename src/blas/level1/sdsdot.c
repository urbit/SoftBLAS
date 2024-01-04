#include "softfloat.h"

float32_t sdsdot( const uint64_t   N,       // length
                  const float32_t  alpha,   // accumulator
                  const float32_t *X,       // vector
                  const uint64_t   incX,    // stride
                  const float32_t *Y,       // vector
                  const uint64_t   incY     // stride
                ) {
    float32_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f32_add(dot, f32_mul(*X, *Y));
    }
    return(dot);
}
