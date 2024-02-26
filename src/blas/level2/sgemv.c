#include "softblas.h"

void sgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *X, const int64_t incX, const float32_t beta, float32_t *Y, const int64_t incY) {
    const float32_t ZERO = { SB_REAL32_ZERO };

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
                float32_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    dotProduct = f32_add(dotProduct, f32_mul(A[i + j * lda], X[j * incX]));
                }
                Y[i * incY] = f32_add(f32_mul(alpha, dotProduct), f32_mul(beta, Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            for (uint64_t i = 0; i < N; ++i) {
                float32_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    dotProduct = f32_add(dotProduct, f32_mul(A[j + i * lda], X[j * incX]));
                }
                Y[i * incY] = f32_add(f32_mul(alpha, dotProduct), f32_mul(beta, Y[i * incY]));
            }
        }
    } else if (Layout == 'R' || Layout == 'r') {
        // Row-major order
        if (Trans == 'N' || Trans == 'n') {
            // No transpose
            for (uint64_t i = 0; i < M; ++i) {
                float32_t dotProduct = ZERO;
                for (uint64_t j = 0; j < N; ++j) {
                    dotProduct = f32_add(dotProduct, f32_mul(A[i * lda + j], X[j * incX]));
                }
                Y[i * incY] = f32_add(f32_mul(alpha, dotProduct), f32_mul(beta, Y[i * incY]));
            }
        } else if (Trans == 'T' || Trans == 't') {
            // Transpose
            for (uint64_t i = 0; i < N; ++i) {
                float32_t dotProduct = ZERO;
                for (uint64_t j = 0; j < M; ++j) {
                    dotProduct = f32_add(dotProduct, f32_mul(A[j * lda + i], X[j * incX]));
                }
                Y[i * incY] = f32_add(f32_mul(alpha, dotProduct), f32_mul(beta, Y[i * incY]));
            }
        } 
    }
}
