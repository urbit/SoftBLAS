#include "softblas.h"

void sscal(uint64_t N, float32_t SA, float32_t *SX, uint64_t incX) {
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        SX[iX] = f32_mul(SA, SX[iX]);
        iX += incX;
    }
}
