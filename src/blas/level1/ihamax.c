#include "softblas.h"

uint64_t ihamax(uint64_t N, const float16_t *HX, uint64_t incX) {
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 1;

    uint64_t ihamax = 0;
    uint64_t ix = 0;
    float16_t hmax = f16_abs(HX[0]);
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        if (f16_gt(f16_abs(HX[ix]), hmax)) {
            ihamax = i;
            hmax = f16_abs(HX[ix]);
        }
        ix += incX;
    }
    
    return ihamax;
}
