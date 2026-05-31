#include "softblas.h"

//  Index of the first element of maximal magnitude |Re| + |Im| (the BLAS
//  ICAMAX / scabs1 convention), returned 0-based to match isamax. Returns 0
//  for N < 1, a single element, or a non-positive stride.
uint64_t icamax(uint64_t N, const complex32_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 0;

    uint64_t imax = 0;
    uint64_t ix = 0;
    complex32_t z = CX[0];
    float32_t smax = f32_add(f32_abs(z.real), f32_abs(z.imag));
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        z = CX[ix];
        float32_t mag = f32_add(f32_abs(z.real), f32_abs(z.imag));
        if (f32_gt(mag, smax)) {
            imax = i;
            smax = mag;
        }
        ix += incX;
    }

    return imax;
}
