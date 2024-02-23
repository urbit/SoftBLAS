#include "softblas.h"

float16_t hnrm2(uint64_t N, const float16_t *X, uint64_t incX) {
    float16_t norm = { SB_REAL16_ZERO };
    
    if (N < 1 || incX < 1) {
        return norm;
    } else if (N == 1) {
        norm = f16_abs(X[0]);
        return norm;
    } else {
        float16_t scale = { SB_REAL16_ZERO };
        float16_t ssq = { SB_REAL16_ONE };
        float16_t absXI;

        for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
            if (f16_ne(X[ix], (float16_t){ SB_REAL16_ZERO })) {
                absXI = f16_abs(X[ix]);
                if (f16_lt(scale, absXI)) {
                    ssq = f16_add((float16_t){ SB_REAL16_ONE }, f16_mul(ssq, f16_div(f16_mul(scale, scale), f16_mul(absXI, absXI))));
                    scale = absXI;
                } else {
                    ssq = f16_add(ssq, f16_div(f16_mul(absXI, absXI), f16_mul(scale, scale)));
                }
            }
        }

        norm = f16_mul(scale, f16_sqrt(ssq));

        return norm;
    }
}
