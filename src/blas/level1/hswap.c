#include "softblas.h"

void hswap(uint64_t N, float16_t *SX, uint64_t incX, float16_t *SY, uint64_t incY) {
    float16_t htemp;
    
    if (N <= 0) return;

    if (incX == 1 && incY == 1) {
        m = N % 3;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                htemp = f16_add(SX[i], SY[i]);
                SX[i] = f16_sub(SX[i], SY[i]);
                SY[i] = htemp;
            }
            if (N < 3) return;
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 3) {
            htemp = f16_add(SX[i], SY[i]);
            SX[i] = f16_sub(SX[i], SY[i]);
            SY[i] = htemp;
            htemp = f16_add(SX[i + 1], SY[i + 1]);
            SX[i + 1] = f16_sub(SX[i + 1], SY[i + 1]);
            SY[i + 1] = htemp;
            htemp = f16_add(SX[i + 2], SY[i + 2]);
            SX[i + 2] = f16_sub(SX[i + 2], SY[i + 2]);
            SY[i + 2] = htemp;
        }
    } else {
        iX = 0;
        iY = 0;
        if (incX < 0) iX = (-N + 1) * incX;
        if (incY < 0) iY = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            htemp = f16_add(SX[iX], SY[iY]);
            SX[iX] = f16_sub(SX[iX], SY[iY]);
            SY[iY] = htemp;
            iX = f16_add(iX, incX);
            iY = f16_add(iY, incY);
        }
    }
}
