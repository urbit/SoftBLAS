#include "softblas.h"

void dswap(uint64_t N, float64_t *SX, uint64_t incX, float64_t *SY, uint64_t incY) {
    float64_t dtemp;
    
    if (N <= 0) return;

    if (incX == 1 && incY == 1) {
        m = N % 3;
        if (m != 0) {
            for (uint64_t i = 0; i < m; i++) {
                dtemp = f64_add(SX[i], SY[i]);
                SX[i] = f64_sub(SX[i], SY[i]);
                SY[i] = dtemp;
            }
            if (N < 3) return;
        }
        mp1 = m + 1;
        for (uint64_t i = mp1; i < N; i += 3) {
            dtemp = f64_add(SX[i], SY[i]);
            SX[i] = f64_sub(SX[i], SY[i]);
            SY[i] = dtemp;
            dtemp = f64_add(SX[i + 1], SY[i + 1]);
            SX[i + 1] = f64_sub(SX[i + 1], SY[i + 1]);
            SY[i + 1] = dtemp;
            dtemp = f64_add(SX[i + 2], SY[i + 2]);
            SX[i + 2] = f64_sub(SX[i + 2], SY[i + 2]);
            SY[i + 2] = dtemp;
        }
    } else {
        iX = 0;
        iY = 0;
        if (incX < 0) iX = (-N + 1) * incX;
        if (incY < 0) iY = (-N + 1) * incY;
        for (uint64_t i = 0; i < N; i++) {
            dtemp = f64_add(SX[iX], SY[iY]);
            SX[iX] = f64_sub(SX[iX], SY[iY]);
            SY[iY] = dtemp;
            iX = f64_add(iX, incX);
            iY = f64_add(iY, incY);
        }
    }
}
