#include "softblas.h"

void qscal(uint64_t N, float128_t QA, float128_t *QX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        f128M_mul(&QA, &(QX[iX]), &(QX[iX]));
        nan_unify_q(&(QX[iX]));
        iX += incX;
    }
}
