#include "softfloat.h"

float64_t dsdsdot( const uint64_t   N,      // length
                   const float64_t  alpha,  // accumulator
                   const float64_t *X,      // vector
                   const uint64_t   incX,   // stride
                   const float64_t *Y,      // vector
                   const uint64_t   incY    // stride
                 ) {
    float64_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f64_add(dot, f64_mul(*X, *Y));
    }
    return dot;
}
