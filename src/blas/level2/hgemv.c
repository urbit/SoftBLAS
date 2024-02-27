#include "softblas.h"
#include <stdio.h>

void hgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *X, const int64_t incX, const float16_t beta, float16_t *Y, const uint64_t incY) {
    const float16_t ZERO = { SB_REAL16_ZERO };

    if (Layout != 'C' && Layout != 'c' && Layout != 'R' && Layout != 'r') {
        printf("Invalid order parameter\n");
        return;
    }

    if (Trans != 'N' && Trans != 'n' && Trans != 'T' && Trans != 't') {
        printf("Invalid trans parameter\n");
        return;
    }

    if (Layout == 'C' || Layout == 'c') {
        // Column-major order
        if (Trans == 'N' || Trans == 'n') {
            // No transpose
            for (uint64_t i = 0; i < M; ++i) {
                float16_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    dotProduct = f16_add(dotProduct, f16_mul(A[(i*incY) + (j*incX) * lda], X[j * incX]));
                }
                Y[i * incY] = f16_add(f16_mul(alpha, dotProduct), f16_mul(beta, Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            for (uint64_t i = 0; i < N; ++i) {
                float16_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    dotProduct = f16_add(dotProduct, f16_mul(A[(j*incX) + (i*incY) * lda], X[j * incX]));
                }
                Y[i * incY] = f16_add(f16_mul(alpha, dotProduct), f16_mul(beta, Y[i * incY]));
            }
        }
    } else if (Layout == 'R' || Layout == 'r') {
        // Row-major order
        if (Trans == 'N' || Trans == 'n') {
            // No transpose
            for (uint64_t i = 0; i < M; ++i) {
                float16_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    dotProduct = f16_add(dotProduct, f16_mul(A[(i*incY) * lda + (j*incX)], X[j * incX]));
                }
                Y[i * incY] = f16_add(f16_mul(alpha, dotProduct), f16_mul(beta, Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            for (uint64_t i = 0; i < N; ++i) {
                float16_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    dotProduct = f16_add(dotProduct, f16_mul(A[(j*incX) * lda + (i*incY)], X[j * incX]));
                }
                Y[i * incY] = f16_add(f16_mul(alpha, dotProduct), f16_mul(beta, Y[i * incY]));
            }
        } 
    }
}
