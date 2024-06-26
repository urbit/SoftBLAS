#include "softblas.h"

void sswap(uint64_t N, float32_t *SX, uint64_t incX, float32_t *SY, uint64_t incY) {
    float32_t stemp;

    if (N == 0) exit(-1);

    uint64_t iX = 0;
    uint64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        stemp = SX[iX];
        SX[iX] = SY[iY];
        SY[iY] = stemp;
        iX += incX;
        iY += incY;
    }
}
