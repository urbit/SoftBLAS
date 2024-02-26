#include "softblas.h"

float128_t qasum(uint64_t N, const float128_t *QX, uint64_t incX) {
    float128_t *qtemp = qvec((float128_pair_t[]){
            {.hi = SB_REAL128U_ZERO, .lo = SB_REAL128L_ZERO}},
        1);
    for (uint64_t i = 0; i < N; i++) {
        f128M_add(qtemp, &(f128_abs(QX[i*incX])), qtemp);
    }

    nan_unify_q(qtemp);
    return *qtemp;
}
