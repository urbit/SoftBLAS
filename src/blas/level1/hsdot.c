#include "softfloat.h"

float32_t hsdot( const uint64_t   N,       // length
                  const float16_t  alpha,   // accumulator
                  const float16_t *X,       // vector
                  const uint64_t   incX,    // stride
                  const float16_t *Y,       // vector
                  const uint64_t   incY     // stride
                ) {
    float32_t dot = alpha;
    for (uint64_t i = N; i; i--, X += incX, Y += incY) {
        dot = f32_add(dot, f32_mul(float16_to_f32(*X), float16_to_f32(*Y)));
    }
    return(dot);
}
