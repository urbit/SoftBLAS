#include "softblas.h"

void dscal(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX) {
    uint64_t iX = 0;
    for (uint64_t i = 0; i < N; i++) {
        DX[iX] = f64_mul(DA, DX[iX]);
        iX += incX;
    }
}
