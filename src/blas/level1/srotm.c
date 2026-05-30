#include "softblas.h"

//  Apply a modified Givens rotation H (selected by the flag P[0]) to the
//  pair of vectors X, Y.  H is built from P[1..4] as in the reference SROTM.
void srotm(const uint64_t N, float32_t *X, const uint64_t incX, float32_t *Y, const uint64_t incY, const float32_t *P, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const float32_t flag  = P[0];
    const float32_t NEGTWO = { SB_REAL32_NEGTWO };
    const float32_t NEGONE = { SB_REAL32_NEGONE };
    const float32_t ZERO   = { SB_REAL32_ZERO };
    const float32_t ONE    = { SB_REAL32_ONE };

    if (f32_eq(flag, NEGTWO)) return;

    float32_t h11, h21, h12, h22;
    if (f32_eq(flag, NEGONE)) {
        h11 = P[1]; h21 = P[2]; h12 = P[3]; h22 = P[4];
    } else if (f32_eq(flag, ZERO)) {
        h11 = ONE;  h21 = P[2]; h12 = P[3]; h22 = ONE;
    } else { /* flag == 1 */
        h11 = P[1]; h21 = NEGONE; h12 = ONE; h22 = P[4];
    }

    uint64_t iX = 0, iY = 0;
    for (uint64_t i = 0; i < N; i++) {
        float32_t w = X[iX], z = Y[iY];
        X[iX] = nan_unify_s(f32_add(f32_mul(w, h11), f32_mul(z, h12)));
        Y[iY] = nan_unify_s(f32_add(f32_mul(w, h21), f32_mul(z, h22)));
        iX += incX;
        iY += incY;
    }
}
