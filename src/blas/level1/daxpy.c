#include "softblas.h"

void daxpy(uint64_t N, float64_t DA, float64_t *DX, int64_t incX, float64_t *DY, int64_t incY) {
    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        DY[iY] = f64_add(DY[iY], f64_mul(DA, DX[iX]));
        iX += incX;
        iY += incY;
    }
}
