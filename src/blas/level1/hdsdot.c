#include "softfloat.h"

float16_t hdsdot( const uint64_t   N,       // length
                  const float16_t  alpha,   // accumulator
                  const float16_t *X,       // vector
                  const uint64_t   incX,    // stride
                  const float16_t *Y,       // vector
                  const uint64_t   incY     // stride
                ) {
    float16_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f16_add(dot, f16_mul(*X, *Y));
    }
    return(dot);
}
