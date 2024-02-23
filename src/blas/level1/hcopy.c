#include "softblas.h"

void hcopy(uint64_t N, const float16_t *HX, int64_t incX, float16_t *HY, int64_t incY) {
    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        HY[iY] = HX[iX];
        iX += incX;
        iY += incY;
    }
}