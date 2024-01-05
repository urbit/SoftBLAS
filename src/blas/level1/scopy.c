#include "softblas.h"

void scopy(uint32_t N, const float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY) {
    if (incX == 1 && incY == 1) {
        // code for both increments equal to 1

        // clean-up loop
        uint64_t m = N % 7;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                SY[i] = SX[i];
            }
            if (N < 7) return;
        }
        uint64_t mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 7) {
            SY[i] = SX[i];
            SY[i + 1] = SX[i + 1];
            SY[i + 2] = SX[i + 2];
            SY[i + 3] = SX[i + 3];
            SY[i + 4] = SX[i + 4];
            SY[i + 5] = SX[i + 5];
            SY[i + 6] = SX[i + 6];
        }
    } else {
        // code for unequal increments or equal increments not equal to 1
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-N + 1) * incX;
        if (incY < 0) iy = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            SY[iy] = SX[ix];
            ix = f32_add(ix, incX);
            iy = f32_add(iy, incY);
        }
    }
}
