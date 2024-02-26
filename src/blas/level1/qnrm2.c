#include "softblas.h"

float128_t qnrm2(uint64_t N, const float128_t *X, uint64_t incX) {
    float128_t norm = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };
    
    if (N < 1 || incX < 1) {
        return norm;
    } else if (N == 1) {
        norm = f128_abs(X[0]);
        return norm;
    } else {
        float128_t scale = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };
        float128_t ssq = { SB_REAL128L_ONE, SB_REAL128U_ONE };
        float128_t absXI;
        float128_t *QT = malloc(sizeof(float128_t));
        float128_t *QU = malloc(sizeof(float128_t));
        float128_t *QV = malloc(sizeof(float128_t));
        float128_t *QW = malloc(sizeof(float128_t));
        float128_t *QX = malloc(sizeof(float128_t));
        float128_t *QY = malloc(sizeof(float128_t));
        const float128_t QZERO = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };
        const float128_t QONE  = { SB_REAL128L_ONE,  SB_REAL128U_ONE };

        for (uint64_t ix = 0; ix < 1 + (N - 1) * incX; ix += incX) {
            if (f128M_ne(&(X[ix]), &QZERO)) {
                absXI = f128_abs(X[ix]);
                if (f128M_lt(&scale, &absXI)) {
                    f128M_mul(&absXI, &absXI, QT);
                    f128M_mul(&scale, &scale, QU);
                    f128M_div(QU, QT, QW);
                    f128M_mul(&ssq, QW, QX);
                    f128M_add(&QONE, QX, &ssq);
                    scale = absXI;
                } else {
                    f128M_mul(&absXI, &absXI, QT);
                    f128M_mul(&scale, &scale, QU);
                    f128M_div(QT, QU, QW);
                    f128M_add(&ssq, QW, &ssq);
                }
            }
        }

        f128M_sqrt(&ssq, QY);
        f128M_mul(&scale, QY, &norm);

        free(QT);
        free(QU);
        free(QV);
        free(QW);
        free(QX);
        free(QY);

        return norm;
    }
}
