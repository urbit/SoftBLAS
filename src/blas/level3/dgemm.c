#include "softblas.h"
#include <stdio.h>

void dgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *B, const uint64_t ldb, const float64_t beta, float64_t *C, const uint64_t ldc) {
    const float64_t ZERO = {SB_REAL64_ZERO};

    if (transA != 'N' && transA != 'n' && transA != 'T' && transA != 't') {
        printf("Invalid transA parameter\n");
        return;
    }

    if (transB != 'N' && transB != 'n' && transB != 'T' && transB != 't') {
        printf("Invalid transB parameter\n");
        return;
    }

    for (uint64_t i = 0; i < M; i++) {
        for (uint64_t j = 0; j < P; j++) {
            float64_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldc : k + j * ldc;
                float64_t a = A[indexA];
                float64_t b = B[indexB];
                dotProduct = f64_add(dotProduct, f64_mul(a, b));
            }
            C[j + i * ldc] = f64_add(f64_mul(alpha, dotProduct), f64_mul(beta, C[j + i * ldc]));
        }
    }
}
