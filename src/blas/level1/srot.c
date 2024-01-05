#include "softblas.h"

void srot(const uint32_t N, float32_t *X, const uint32_t  incX, float32_t *Y, const uint32_t incY, const float32_t c, const float32_t s) {
    float32_t tmp;

    if (c != SB_REAL32_ONE || s != SB_REAL32_ZERO) {
        if (incX == 1 && incY == 1) {
            for (uint64_t i=0; i != N; i++) {
                tmp = f32_add(f32_mul(c, X[i]), f32_mul(s, Y[i]));
                Y[i] = f32_sub(f32_mul(c, Y[i]), f32_mul(s, X[i]));
                X[i] = tmp;
            }
        }
        else
        {
            for (uint64_t i=N; i; i--, Y += incY, X += incX)
            {
                tmp = f32_add(f32_mul(c, *X), f32_mul(s, *Y));
                *Y = f32_sub(f32_mul(c, *Y), f32_mul(s, *X));
                *X = tmp;
            }
        }
    }
}
