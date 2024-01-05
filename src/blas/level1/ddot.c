#include "softblas.h"

ddot(const uint64_t N, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY) {
    float64_t ddot = SB_REAL64_ZERO;

    if (incX == 1 && incY == 1) {
        uint64_t m = N % 5;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                ddot = f64_add(ddot, f64_mul(X[i], Y[i]));
            }
            if (N < 5) {
                return ddot;
            }
        }
        uint64_t mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 5) {
            ddot = f64_add(ddot, f64_add(f64_add(f64_add(f64_add(f64_mul(X[i], Y[i]), f64_mul(X[i+1], Y[i+1])), f64_mul(X[i+2], Y[i+2])), f64_mul(X[i+3], Y[i+3])), f64_mul(X[i+4], Y[i+4])));
        }
    } else {
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-N + 1) * incX;
        if (incY < 0) iy = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            ddot = f64_add(ddot, f64_mul(X[ix], Y[iy]));
            ix += incX;
            iy += incY;
        }
    }
    return ddot;
}
