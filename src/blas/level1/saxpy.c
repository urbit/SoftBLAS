#include "softblas.h"

void saxpy(uint64_t N, float32_t SA, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY) {
    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        SY[iY] = f32_add(SY[iY], f32_mul(SA, SX[iX]));
        iX += incX;
        iY += incY;
    }
}
