#include "softblas.h"

float128_t qdot(const uint64_t N, const float128_t *X, const int64_t incX, const float128_t *Y, const int64_t incY) {
    // XXX note that 0,1 is lo,hi
    float128_t qdot = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    int64_t iX = 0;
    int64_t iY = 0;
    float128_t QT;
    if (incX < 0) iX = (-N + 1) * incX;
    if (incY < 0) iY = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        f128M_mul(&(Y[iY]), &(X[iX]), &QT);
        f128M_add(&qdot, &QT, &qdot);
        iX += incX;
        iY += incY;
    }

    return qdot;
}
