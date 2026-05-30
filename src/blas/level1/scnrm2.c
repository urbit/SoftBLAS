#include "softblas.h"

float32_t scnrm2(uint64_t N, const complex32_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
      float32_t norm, scale, ssq, temp;
      
      if (N < 1 || incX < 1) {
            norm = SB_REAL32_ZERO;
      } else {
            scale = SB_REAL32_ZERO;
            ssq = SB_REAL32_ONE;

        for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
            if (f32_ne(X[ix], (complex32_t){ SB_REAL32_ZERO })) {
                absXI = f32_abs(X[ix]);
                if (f32_lt(scale, absXI)) {
                    ssq = f32_add((complex32_t){ SB_REAL32_ONE }, f32_mul(ssq, f32_div(f32_mul(scale, scale), f32_mul(absXI, absXI))));
                    scale = absXI;
                } else {
                    ssq = f32_add(ssq, f32_div(f32_mul(absXI, absXI), f32_mul(scale, scale)));
                }
            }
        }

        norm = f32_mul(scale, f32_sqrt(ssq));

        return norm;
    }
}
