#include "softblas.h"

float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX) {
    float32_t stemp = { SB_REAL32_ZERO };
    
    for (uint64_t i = 0; i < N; i++) {
        stemp = f32_add(stemp, f32_abs(SX[i*incX]));
    }

    return nan_unify_s(stemp);
}
