#include "softblas.h"

uint64_t iqamax(uint64_t N, const float128_t *QX, uint64_t incX) {
    if (N < 1 || incX <= 0) return 0;
    if (N == 1) return 1;

    uint64_t iqamax = 0;
    uint64_t ix = 0;
    float128_t qmax = f128_abs(QX[0]);
    ix += incX;
    for (uint64_t i = 1; i < N; i++) {
        if (f128M_gt(&f128_abs(QX[ix]), &qmax)) {
            iqamax = i;
            qmax = f128_abs(QX[ix]);
        }
        ix += incX;
    }
    
    return iqamax;
}
