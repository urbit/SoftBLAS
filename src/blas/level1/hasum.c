#include "softblas.h"

float16_t hasum(uint64_t N, const float16_t *HX, uint64_t incX) {
    float16_t htemp = { SB_REAL16_ZERO };
    
    for (uint64_t i = 0; i < N; i++) {
        htemp = f16_add(htemp, f16_abs(HX[i*incX]));
    }

    return nan_unify_h(htemp);
}
