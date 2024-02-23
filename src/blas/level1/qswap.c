#include "softblas.h"

void qswap(uint64_t N, float128_t *QX, uint64_t incX, float128_t *QY, uint64_t incY) {
    float128_t qtemp;

    if (N == 0) return;

    uint64_t iX = 0;
    uint64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        qtemp = QX[iX];
        QX[iX] = QY[iY];
        QY[iY] = qtemp;
        iX += incX;
        iY += incY;
    }
}
