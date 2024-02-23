#include "softblas.h"

float64_t dasum(uint64_t N, const float64_t *DX, uint64_t incX) {
    float64_t dtemp = { SB_REAL64_ZERO };
    
    for (uint64_t i = 0; i < N; i++) {
        dtemp = f64_add(dtemp, f64_abs(DX[i*incX]));
    }

    return nan_unify_d(dtemp);
}