#include "softblas.h"

//  Apply a plane rotation: [x_i; y_i] <- [c s; -s c] [x_i; y_i].
void qrot(const uint64_t N, float128_t *X, const uint64_t incX, float128_t *Y, const uint64_t incY, const float128_t c, const float128_t s, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    float128_t cx, sy, cy, sx, newx, newy;
    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        f128M_mul(&c, &(X[iX]), &cx);
        f128M_mul(&s, &(Y[iY]), &sy);
        f128M_add(&cx, &sy, &newx);   // newx = c*x + s*y
        f128M_mul(&c, &(Y[iY]), &cy);
        f128M_mul(&s, &(X[iX]), &sx);
        f128M_sub(&cy, &sx, &newy);   // newy = c*y - s*x
        Y[iY] = newy;
        nan_unify_q(&(Y[iY]));
        X[iX] = newx;
        nan_unify_q(&(X[iX]));
        iX += incX;
        iY += incY;
    }
}
