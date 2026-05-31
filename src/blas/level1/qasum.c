#include "softblas.h"

float128_t qasum(uint64_t N, const float128_t *QX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    //  Stack accumulator: a BLAS kernel should not touch the heap.
    float128_t qtemp = {{ SB_REAL128L_ZERO, SB_REAL128U_ZERO }};
    if (incX < 1) { nan_unify_q(&qtemp); return qtemp; }   // reject 0/garbage stride (cf. qnrm2)
    for (uint64_t i = 0; i < N; i++) {
        f128M_add(&qtemp, &(f128_abs(QX[i*incX])), &qtemp);
    }

    nan_unify_q(&qtemp);
    return qtemp;
}
