#include "softblas.h"

//  Euclidean norm of a complex vector: sqrt( sum |z_i|^2 ), computed as a
//  scaled sum of squares over the 2N real components (Re and Im) to avoid
//  overflow/underflow.
float16_t hinrm2(uint64_t N, const complex16_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float16_t ZERO = { SB_REAL16_ZERO };

    if (N < 1 || incX < 1 || SB_STRIDE_OVERFLOWS(N, incX)) {
        return ZERO;
    }

    float16_t scale = { SB_REAL16_ZERO };
    float16_t ssq   = { SB_REAL16_ONE };

    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float16_t parts[2] = { CX[ix].real, CX[ix].imag };
        for (int p = 0; p < 2; p++) {
            if (f16_ne(parts[p], ZERO)) {
                float16_t absXI = f16_abs(parts[p]);
                if (f16_lt(scale, absXI)) {
                    float16_t r = f16_div(scale, absXI);
                    ssq = f16_add((float16_t){ SB_REAL16_ONE }, f16_mul(ssq, f16_mul(r, r)));
                    scale = absXI;
                } else {
                    float16_t r = f16_div(absXI, scale);
                    ssq = f16_add(ssq, f16_mul(r, r));
                }
            }
        }
    }

    return nan_unify_h(f16_mul(scale, f16_sqrt(ssq)));
}
