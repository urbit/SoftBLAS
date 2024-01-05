#include "softblas.h"

void hcopy(uint16_t N, const float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY) {
    if (incX == 1 && incY == 1) {
        // code for both increments equal to 1

        // clean-up loop
        uint64_t m = N % 7;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                HY[i] = HX[i];
            }
            if (N < 7) return;
        }
        uint64_t mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 7) {
            HY[i] = HX[i];
            HY[i + 1] = HX[i + 1];
            HY[i + 2] = HX[i + 2];
            HY[i + 3] = HX[i + 3];
            HY[i + 4] = HX[i + 4];
            HY[i + 5] = HX[i + 5];
            HY[i + 6] = HX[i + 6];
        }
    } else {
        // code for unequal increments or equal increments not equal to 1
        uint64_t ix = 0;
        uint64_t iy = 0;
        if (incX < 0) ix = (-N + 1) * incX;
        if (incY < 0) iy = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            HY[iy] = HX[ix];
            ix = f16_add(ix, incX);
            iy = f16_add(iy, incY);
        }
    }
}
