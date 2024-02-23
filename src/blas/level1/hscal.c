#include "softblas.h"

void hscal(uint64_t N, float16_t HA, float16_t *HX, uint64_t incX) {
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        HX[iX] = f16_mul(HA, HX[iX]);
        iX += incX;
    }
}
