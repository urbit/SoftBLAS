#include "softblas.h"

uint64_t idamax(uint64_t N, const float64_t *DX, uint64_t incX) {
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 1;

    uint64_t idamax = 0;
    uint64_t ix = 0;
    float64_t dmax = f64_abs(DX[0]);
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        if (f64_gt(f64_abs(DX[ix]), dmax)) {
            idamax = i;
            dmax = f64_abs(DX[ix]);
        }
        ix += incX;
    }
    
    return idamax;
}
