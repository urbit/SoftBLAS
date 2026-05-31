#include "softblas.h"

//  Scale a complex vector in place: CX <- CA * CX (quad complex).
void vscal(uint64_t N, complex128_t CA, complex128_t *CX, uint64_t incX, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        CX[iX] = c128_mul(CA, CX[iX]);
        nan_unify_v(&(CX[iX]));
        iX += incX;
    }
}
