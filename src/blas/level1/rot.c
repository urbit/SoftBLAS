#include "softblas.h"
#include "softfloat.h"

void srot(const uint64_t  N,       // length
             float32_t      *X,       // vector
             const uint64_t  incX,    // stride
             float32_t      *Y,       // vector
             const uint64_t  incY,    // stride
             const float32_t c,       // cos(theta)
             const float32_t s        // sin(theta)
          ) {
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

void drot(const uint64_t  N,       // length
             float64_t      *X,       // vector
             const uint64_t  incX,    // stride
             float64_t      *Y,       // vector
             const uint64_t  incY,    // stride
             const float64_t c,       // cos(theta)
             const float64_t s        // sin(theta)
          ) {
    float64_t tmp;

    if (c != SB_REAL64_ONE || s != SB_REAL64_ZERO) {
        if (incX == 1 && incY == 1) {
            for (uint64_t i=0; i != N; i++) {
                tmp = f64_add(f64_mul(c, X[i]), f64_mul(s, Y[i]));
                Y[i] = f64_sub(f64_mul(c, Y[i]), f64_mul(s, X[i]));
                X[i] = tmp;
            }
        }
        else
        {
            for (uint64_t i=N; i; i--, Y += incY, X += incX)
            {
                tmp = f64_add(f64_mul(c, *X), f64_mul(s, *Y));
                *Y = f64_sub(f64_mul(c, *Y), f64_mul(s, *X));
                *X = tmp;
            }
        }
    }
}

void hrot(const uint64_t  N,       // length
             float16_t      *X,       // vector
             const uint64_t  incX,    // stride
             float16_t      *Y,       // vector
             const uint64_t  incY,    // stride
             const float16_t c,       // cos(theta)
             const float16_t s        // sin(theta)
          ) {
    float16_t tmp;

    if (c != SB_REAL32_ONE || s != SB_REAL32_ZERO) {
        if (incX == 1 && incY == 1) {
            for (uint64_t i=0; i != N; i++) {
                tmp = f16_add(f16_mul(c, X[i]), f16_mul(s, Y[i]));
                Y[i] = f16_sub(f16_mul(c, Y[i]), f16_mul(s, X[i]));
                X[i] = tmp;
            }
        }
        else
        {
            for (uint64_t i=N; i; i--, Y += incY, X += incX)
            {
                tmp = f16_add(f16_mul(c, *X), f16_mul(s, *Y));
                *Y = f16_sub(f16_mul(c, *Y), f16_mul(s, *X));
                *X = tmp;
            }
        }
    }
}
