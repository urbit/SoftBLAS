#include "softblas.h"

//  Euclidean norm of a complex vector: sqrt( sum |z_i|^2 ), computed as a
//  scaled sum of squares over the 2N real components (Re and Im) to avoid
//  overflow/underflow.
float32_t scnrm2(uint64_t N, const complex32_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float32_t ZERO = { SB_REAL32_ZERO };

    if (N < 1 || incX < 1) {
        return ZERO;
    }

    float32_t scale = { SB_REAL32_ZERO };
    float32_t ssq   = { SB_REAL32_ONE };

    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float32_t parts[2] = { CX[ix].real, CX[ix].imag };
        for (int p = 0; p < 2; p++) {
            if (f32_ne(parts[p], ZERO)) {
                float32_t absXI = f32_abs(parts[p]);
                if (f32_lt(scale, absXI)) {
                    float32_t r = f32_div(scale, absXI);
                    ssq = f32_add((float32_t){ SB_REAL32_ONE }, f32_mul(ssq, f32_mul(r, r)));
                    scale = absXI;
                } else {
                    float32_t r = f32_div(absXI, scale);
                    ssq = f32_add(ssq, f32_mul(r, r));
                }
            }
        }
    }

    return nan_unify_s(f32_mul(scale, f32_sqrt(ssq)));
}
