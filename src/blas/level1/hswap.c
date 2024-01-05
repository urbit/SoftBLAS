#include "softblas.h"

void hswap(uint16_t N, float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY) {
    float16_t htemp;
    
    if (N <= 0) return;

    if (incX == 1 && incY == 1) {
        m = N % 3;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                htemp = f16_add(HX[i], HY[i]);
                HX[i] = f16_sub(HX[i], HY[i]);
                HY[i] = htemp;
            }
            if (N < 3) return;
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 3) {
            htemp = f16_add(HX[i], HY[i]);
            HX[i] = f16_sub(HX[i], HY[i]);
            HY[i] = htemp;
            htemp = f16_add(HX[i + 1], HY[i + 1]);
            HX[i + 1] = f16_sub(HX[i + 1], HY[i + 1]);
            HY[i + 1] = htemp;
            htemp = f16_add(HX[i + 2], HY[i + 2]);
            HX[i + 2] = f16_sub(HX[i + 2], HY[i + 2]);
            HY[i + 2] = htemp;
        }
    } else {
        iX = 0;
        iY = 0;
        if (incX < 0) iX = (-N + 1) * incX;
        if (incY < 0) iY = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            htemp = f16_add(HX[iX], HY[iY]);
            HX[iX] = f16_sub(HX[iX], HY[iY]);
            HY[iY] = htemp;
            iX = f16_add(iX, incX);
            iY = f16_add(iY, incY);
        }
    }
}
