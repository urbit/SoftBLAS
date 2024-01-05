#include "softfloat.h"

float32_t sdot( const uint64_t   N,       // length
                const float32_t *X,       // vector
                const uint64_t   incX,    // stride
                const float32_t *Y,       // vector
                const uint64_t   incY     // stride
              ) {
    float32_t sdot = SB_REAL32_ZERO;

    if (incX == 1 && incY == 1) {
        uint64_t m = n % 5;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                sdot = f32_add(sdot, f32_mul(sx[i], sy[i]));
            }
            if (n < 5) {
                return sdot;
            }
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < n; i += 5) {
            sdot = f32_add(sdot, f32_add(f32_add(f32_add(f32_add(f32_mul(sx[i], sy[i]), f32_mul(sx[i+1], sy[i+1])), f32_mul(sx[i+2], sy[i+2])), f32_mul(sx[i+3], sy[i+3])), f32_mul(sx[i+4], sy[i+4])));
        }
    } else {
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-n + 1) * incX;
        if (incY < 0) iy = (-n + 1) * incY;
        for (uint64_t i = 0; i < n; i++) {
            sdot = f32_add(sdot, f32_mul(sx[ix], sy[iy]));
            ix += incX;
            iy += incY;
        }
    }
    return sdot;
}
