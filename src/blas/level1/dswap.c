#include "softblas.h"

void dswap(uint64_t N, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY) {
    float64_t dtemp;

    if (N == 0) exit(-1);

    uint64_t iX = 0;
    uint64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        dtemp = DX[iX];
        DX[iX] = DY[iY];
        DY[iY] = dtemp;
        iX += incX;
        iY += incY;
    }
}
