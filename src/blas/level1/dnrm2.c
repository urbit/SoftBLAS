#include "softblas.h"

float64_t dnrm2(uint64_t N, const float64_t *X, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float64_t norm = { SB_REAL64_ZERO };
    
    if (N < 1 || incX < 1) {
        return norm;
    } else if (N == 1) {
        norm = f64_abs(X[0]);
        return norm;
    } else {
        float64_t scale = { SB_REAL64_ZERO };
        float64_t ssq = { SB_REAL64_ONE };
        float64_t absXI;

        for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
            if (f64_ne(X[ix], (float64_t){ SB_REAL64_ZERO })) {
                absXI = f64_abs(X[ix]);
                if (f64_lt(scale, absXI)) {
                    ssq = f64_add((float64_t){ SB_REAL64_ONE }, f64_mul(ssq, f64_div(f64_mul(scale, scale), f64_mul(absXI, absXI))));
                    scale = absXI;
                } else {
                    ssq = f64_add(ssq, f64_div(f64_mul(absXI, absXI), f64_mul(scale, scale)));
                }
            }
        }

        norm = f64_mul(scale, f64_sqrt(ssq));

        return norm;
    }
}
