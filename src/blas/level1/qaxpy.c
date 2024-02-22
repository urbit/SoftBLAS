#include "softblas.h"

void qaxpy(uint64_t N, float128_t QA, float128_t *QX, int64_t incX, float128_t *QY, int64_t incY) {
    int64_t iX = 0;
    int64_t iY = 0;
    float128_t *QT;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        f128M_mul(&QA, &(QX[iX]), QT);
        f128M_add(&(QY[iY]), QT, QT);
        QY[iY] = *QT;
        iX += incX;
        iY += incY;
    }
}
