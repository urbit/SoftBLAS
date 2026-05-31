#include "softblas.h"

//  Scale a complex vector in place: CX <- CA * CX (BLAS CSCAL).
void zscal(uint64_t N, complex64_t CA, complex64_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        CX[iX] = nan_unify_z(c64_mul(CA, CX[iX]));
        iX += incX;
    }
}
