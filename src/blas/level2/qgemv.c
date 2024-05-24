#include "softblas.h"

void qgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float128_t alpha, const float128_t *A, const uint64_t lda, const float128_t *X, const int64_t incX, const float128_t beta, float128_t *Y, const uint64_t incY) {
    const float128_t ZERO = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    if (Layout != 'C' && Layout != 'c' && Layout != 'R' && Layout != 'r') {
        // Invalid order parameter
        exit(-1);
    }

    if (Trans != 'N' && Trans != 'n' && Trans != 'T' && Trans != 't') {
        // Invalid trans parameter
        exit(-1);
    }

    if (Layout == 'C' || Layout == 'c') {
        // Column-major order
        if (Trans == 'N' || Trans == 'n') {
            // No transpose
            float128_t qtemp;
            float128_t qtemp2;
            for (uint64_t i = 0; i < M; ++i) {
                float128_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    f128M_mul(&(A[(i*incY) + (j*incX) * lda]), &(X[j * incX]), &qtemp);
                    f128M_add(&dotProduct, &qtemp, &dotProduct);
                }
                f128M_mul(&alpha, &dotProduct, &qtemp);
                f128M_mul(&beta, &(Y[i * incY]), &qtemp2);
                f128M_add(&qtemp, &qtemp2, &(Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            float128_t qtemp;
            float128_t qtemp2;
            for (uint64_t i = 0; i < N; ++i) {
                float128_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    f128M_mul(&(A[(j*incX) + (i*incY) * lda]), &(X[j * incX]), &qtemp);
                    f128M_add(&dotProduct, &qtemp, &dotProduct);
                }
                f128M_mul(&alpha, &dotProduct, &qtemp);
                f128M_mul(&beta, &(Y[i * incY]), &qtemp2);
                f128M_add(&qtemp, &qtemp2, &(Y[i * incY]));
            }
        }
    } else if (Layout == 'R' || Layout == 'r') {
        // Row-major order
        if (Trans == 'N' || Trans == 'n') {
            // No transpose
            float128_t qtemp;
            float128_t qtemp2;
            for (uint64_t i = 0; i < M; ++i) {
                float128_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    f128M_mul(&(A[(i*incY) * lda + (j*incX)]), &(X[j * incX]), &qtemp);
                    f128M_add(&dotProduct, &qtemp, &dotProduct);
                }
                f128M_mul(&alpha, &dotProduct, &qtemp);
                f128M_mul(&beta, &(Y[i * incY]), &qtemp2);
                f128M_add(&qtemp, &qtemp2, &(Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            float128_t qtemp;
            float128_t qtemp2;
            for (uint64_t i = 0; i < N; ++i) {
                float128_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    f128M_mul(&(A[(j*incX) * lda + (i*incY)]), &(X[j * incX]), &qtemp);
                    f128M_add(&dotProduct, &qtemp, &dotProduct);
                }
                f128M_mul(&alpha, &dotProduct, &qtemp);
                f128M_mul(&beta, &(Y[i * incY]), &qtemp2);
                f128M_add(&qtemp, &qtemp2, &(Y[i * incY]));
            }
        } 
    }
}
