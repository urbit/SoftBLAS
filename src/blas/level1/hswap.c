#include "softblas.h"

void hswap(uint64_t N, float16_t *HX, uint64_t incX, float16_t *HY, uint64_t incY) {
    float16_t htemp;

    if (N == 0) exit(-1);

    uint64_t iX = 0;
    uint64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        htemp = HX[iX];
        HX[iX] = HY[iY];
        HY[iY] = htemp;
        iX += incX;
        iY += incY;
    }
}
