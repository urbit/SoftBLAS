#include "softblas.h"

float32_t scasum(uint64_t N, const complex32_t *CX, uint64_t incX) {
    float32_t stemp = { SB_REAL32_ZERO };
    
    for (uint64_t i = 0; i < N; i++) {
        stemp = f32_add(stemp, f32_add(f32_abs(CX[i*incX].real), f32_abs(CX[i*incX].imag)));
    }

    return nan_unify_s(stemp);
}
