#include "softblas.h"

//  Euclidean norm of a quad-complex vector: sqrt( sum |z_i|^2 ), via a scaled
//  sum of squares over the 2N real components to avoid overflow/underflow.
float128_t qvnrm2(uint64_t N, const complex128_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float128_t ZERO = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    if (N < 1 || incX < 1) {
        return ZERO;
    }

    float128_t scale = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };
    float128_t ssq   = { SB_REAL128L_ONE,  SB_REAL128U_ONE };

    for (uint64_t k = 0, ix = 0; k < N; k++, ix += incX) {
        float128_t parts[2] = { CX[ix].real, CX[ix].imag };
        for (int p = 0; p < 2; p++) {
            if (f128_ne(parts[p], ZERO)) {
                float128_t absXI = f128_abs(parts[p]);
                if (f128_lt(scale, absXI)) {
                    float128_t r = f128_div(scale, absXI);
                    ssq = f128_add((float128_t){ SB_REAL128L_ONE, SB_REAL128U_ONE }, f128_mul(ssq, f128_mul(r, r)));
                    scale = absXI;
                } else {
                    float128_t r = f128_div(absXI, scale);
                    ssq = f128_add(ssq, f128_mul(r, r));
                }
            }
        }
    }

    float128_t result = f128_mul(scale, f128_sqrt(ssq));
    nan_unify_q(&result);
    return result;
}
