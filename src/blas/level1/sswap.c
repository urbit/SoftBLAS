#include "softblas.h"

void sswap(uint32_t N, float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY) {
    float32_t stemp;
    
    if (N <= 0) return;

    if (incX == 1 && incY == 1) {
        m = N % 3;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                stemp = f32_add(SX[i], SY[i]);
                SX[i] = f32_sub(SX[i], SY[i]);
                SY[i] = stemp;
            }
            if (N < 3) return;
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 3) {
            stemp = f32_add(SX[i], SY[i]);
            SX[i] = f32_sub(SX[i], SY[i]);
            SY[i] = stemp;
            stemp = f32_add(SX[i + 1], SY[i + 1]);
            SX[i + 1] = f32_sub(SX[i + 1], SY[i + 1]);
            SY[i + 1] = stemp;
            stemp = f32_add(SX[i + 2], SY[i + 2]);
            SX[i + 2] = f32_sub(SX[i + 2], SY[i + 2]);
            SY[i + 2] = stemp;
        }
    } else {
        iX = 0;
        iY = 0;
        if (incX < 0) iX = (-N + 1) * incX;
        if (incY < 0) iY = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            stemp = f32_add(SX[iX], SY[iY]);
            SX[iX] = f32_sub(SX[iX], SY[iY]);
            SY[iY] = stemp;
            iX = f32_add(iX, incX);
            iY = f32_add(iY, incY);
        }
    }
}
