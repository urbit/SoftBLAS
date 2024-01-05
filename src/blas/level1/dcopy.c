#include "softblas.h"

void dcopy(uint64_t N,
         const float64_t *DX,
         uint64_t incX,
         float64_t *DY,
         uint64_t incY
        ) {
    if (incX == 1 && incY == 1) {
        // code for both increments equal to 1

        // clean-up loop
        uint64_t m = N % 7;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                DY[i] = DX[i];
            }
            if (N < 7) return;
        }
        uint64_t mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 7) {
            DY[i] = DX[i];
            DY[i + 1] = DX[i + 1];
            DY[i + 2] = DX[i + 2];
            DY[i + 3] = DX[i + 3];
            DY[i + 4] = DX[i + 4];
            DY[i + 5] = DX[i + 5];
            DY[i + 6] = DX[i + 6];
        }
    } else {
        // code for unequal increments or equal increments not equal to 1
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-N + 1) * incX;
        if (incY < 0) iy = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            DY[iy] = DX[ix];
            ix = f64_add(ix, incX);
            iy = f64_add(iy, incY);
        }
    }
}
