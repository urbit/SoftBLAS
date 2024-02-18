#include "softblas.h"

void csrot(const uint64_t N, complex32_t *CX, const uint64_t incX, complex32_t *CY, const uint64_t incY, const complex32_t c, const complex32_t s) {
    complex32_t temp;
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        temp = c32_add(c32_mul(c, CX[ix]), c32_mul(s, CY[iy]));
        CY[iy] = c32_sub(c32_mul(c, CY[iy]), c32_mul(s, CX[ix]));
        CX[ix] = temp;
        ix += incX;
        iy += incY;
    }
}
