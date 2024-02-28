#include "softblas.h"

void sgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *B, const uint64_t ldb, const float32_t beta, float32_t *C, const uint64_t ldc) {
    const float32_t ZERO = {SB_REAL32_ZERO};

    if (transA != 'N' && transA != 'n' && transA != 'T' && transA != 't') {
        printf("Invalid transA parameter\n");
        return;
    }

    if (transB != 'N' && transB != 'n' && transB != 'T' && transB != 't') {
        printf("Invalid transB parameter\n");
        return;
    }

    printf("\n");
    for (uint64_t i = 0; i < M; i++) {
        for (uint64_t j = 0; j < P; j++) {
            float32_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldc : k + j * ldc;
                float32_t a = A[indexA];
                float32_t b = B[indexB];
                dotProduct = f32_add(dotProduct, f32_mul(a, b));
            }
            C[j + i * ldc] = f32_add(f32_mul(alpha, dotProduct), f32_mul(beta, C[j + i * ldc]));
        }
    }
}
