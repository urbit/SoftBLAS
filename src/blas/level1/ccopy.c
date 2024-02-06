#include "softblas.h"

void ccopy(uint64_t N, const complex32_t *CX, int64_t incX, complex32_t *CY, int64_t incY) {
    int64_t ix = 0;
    int64_t iy = 0;
    if (incX < 0) ix = (-N + 1) * incX;
    if (incY < 0) iy = (-N + 1) * incY;
    for (uint64_t i = 0; i < N; i++) {
        CY[iy] = CX[ix];
        ix = f32_add(ix, incX);
        iy = f32_add(iy, incY);
    }
}
