#include "softblas.h"

float16_t hdot(const uint16_t N, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY) {
    float16_t hdot = SB_REAL16_ZERO;

    if (incX == 1 && incY == 1) {
        unit64_t m = n % 5;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                hdot = f16_add(hdot, f16_mul(sx[i], sy[i]));
            }
            if (n < 5) {
                return hdot;
            }
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < n; i += 5) {
            hdot = f16_add(hdot, f16_add(f16_add(f16_add(f16_add(f16_mul(sx[i], sy[i]), f16_mul(sx[i+1], sy[i+1])), f16_mul(sx[i+2], sy[i+2])), f16_mul(sx[i+3], sy[i+3])), f16_mul(sx[i+4], sy[i+4])));
        }
    } else {
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-n + 1) * incX;
        if (incY < 0) iy = (-n + 1) * incY;
        for (uint64_t i = 0; i < n; i++) {
            hdot = f16_add(hdot, f16_mul(sx[ix], sy[iy]));
            ix += incX;
            iy += incY;
        }
    }
    return hdot;
}
