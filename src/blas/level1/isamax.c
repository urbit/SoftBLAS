#include "softblas.h"

uint64_t isamax(uint64_t N, const float32_t *SX, uint64_t incX) {
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 1;

    uint64_t isamax = 0;
    uint64_t ix = 0;
    float32_t smax = f32_abs(SX[0]);
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        if (f32_gt(f32_abs(SX[ix]), smax)) {
            isamax = i;
            smax = f32_abs(SX[ix]);
        }
        ix += incX;
    }
    
    return isamax;
}
