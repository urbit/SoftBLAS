#include "softblas.h"

//  Apply a modified Givens rotation H (selected by the flag P[0]) to the
//  pair of vectors X, Y. H is built from P[1..4] as in the reference SROTM.
//  Quad-precision port of srotm.
void qrotm(const uint64_t N, float128_t *X, const uint64_t incX, float128_t *Y, const uint64_t incY, const float128_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float128_t flag   = P[0];
    const float128_t NEGTWO = { SB_REAL128L_NEGTWO, SB_REAL128U_NEGTWO };
    const float128_t NEGONE = { SB_REAL128L_NEGONE, SB_REAL128U_NEGONE };
    const float128_t ZERO   = { SB_REAL128L_ZERO,   SB_REAL128U_ZERO };
    const float128_t ONE    = { SB_REAL128L_ONE,    SB_REAL128U_ONE };

    if (f128M_eq(&flag, &NEGTWO)) return;

    float128_t h11, h21, h12, h22;
    if (f128M_eq(&flag, &NEGONE)) {
        h11 = P[1]; h21 = P[2]; h12 = P[3]; h22 = P[4];
    } else if (f128M_eq(&flag, &ZERO)) {
        h11 = ONE;  h21 = P[2]; h12 = P[3]; h22 = ONE;
    } else { /* flag == 1 */
        h11 = P[1]; h21 = NEGONE; h12 = ONE; h22 = P[4];
    }

    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        float128_t w = X[iX], z = Y[iY];
        float128_t a, b, nx, ny;
        f128M_mul(&w, &h11, &a);
        f128M_mul(&z, &h12, &b);
        f128M_add(&a, &b, &nx);       // nx = w*h11 + z*h12
        f128M_mul(&w, &h21, &a);
        f128M_mul(&z, &h22, &b);
        f128M_add(&a, &b, &ny);       // ny = w*h21 + z*h22
        X[iX] = nx;
        nan_unify_q(&(X[iX]));
        Y[iY] = ny;
        nan_unify_q(&(Y[iY]));
        iX += incX;
        iY += incY;
    }
}
