#include "softblas.h"

void qswap(uint64_t N, float128_t *QX, int64_t incX, float128_t *QY, int64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float128_t qtemp;

    if (N == 0) return;

    int64_t iX = 0;
    int64_t iY = 0;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        qtemp = QX[iX];
        QX[iX] = QY[iY];
        nan_unify_q(&(QX[iX]));
        QY[iY] = qtemp;
        nan_unify_q(&(QY[iY]));
        iX += incX;
        iY += incY;
    }
}
