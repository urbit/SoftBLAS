#include "softblas.h"

//  Euclidean norm of a complex vector: sqrt( sum |z_i|^2 ), computed as a
//  scaled sum of squares over the 2N real components (Re and Im) to avoid
//  overflow/underflow.
float64_t dznrm2(uint64_t N, const complex64_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float64_t ZERO = { SB_REAL64_ZERO };

    if (N < 1 || incX < 1) {
        return ZERO;
    }

    float64_t scale = { SB_REAL64_ZERO };
    float64_t ssq   = { SB_REAL64_ONE };

    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float64_t parts[2] = { CX[ix].real, CX[ix].imag };
        for (int p = 0; p < 2; p++) {
            if (f64_ne(parts[p], ZERO)) {
                float64_t absXI = f64_abs(parts[p]);
                if (f64_lt(scale, absXI)) {
                    float64_t r = f64_div(scale, absXI);
                    ssq = f64_add((float64_t){ SB_REAL64_ONE }, f64_mul(ssq, f64_mul(r, r)));
                    scale = absXI;
                } else {
                    float64_t r = f64_div(absXI, scale);
                    ssq = f64_add(ssq, f64_mul(r, r));
                }
            }
        }
    }

    return nan_unify_d(f64_mul(scale, f64_sqrt(ssq)));
}
