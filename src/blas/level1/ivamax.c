#include "softblas.h"

//  Index of the first element of maximal magnitude |Re| + |Im| over a
//  quad-complex vector, 0-based. Returns 0 for N < 1, a single element, or a
//  non-positive stride.
uint64_t ivamax(uint64_t N, const complex128_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 0;

    uint64_t imax = 0;
    uint64_t ix = 0;
    complex128_t z = CX[0];
    float128_t smax = f128_add(f128_abs(z.real), f128_abs(z.imag));
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        z = CX[ix];
        float128_t mag = f128_add(f128_abs(z.real), f128_abs(z.imag));
        if (f128_gt(mag, smax)) {
            imax = i;
            smax = mag;
        }
        ix += incX;
    }

    return imax;
}
