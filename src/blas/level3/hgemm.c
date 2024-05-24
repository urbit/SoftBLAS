#include "softblas.h"

void hgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *B, const uint64_t ldb, const float16_t beta, float16_t *C, const uint64_t ldc) {
    const float16_t ZERO = {SB_REAL16_ZERO};

    if (transA != 'N' && transA != 'n' && transA != 'T' && transA != 't') {
        // Invalid transA parameter
        exit(-1);
    }

    if (transB != 'N' && transB != 'n' && transB != 'T' && transB != 't') {
        // Invalid transB parameter
        exit(-1);
    }

    for (uint64_t i = 0; i < M; i++) {
        for (uint64_t j = 0; j < P; j++) {
            float16_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldc : k + j * ldc;
                float16_t a = A[indexA];
                float16_t b = B[indexB];
                dotProduct = f16_add(dotProduct, f16_mul(a, b));
            }
            C[j + i * ldc] = f16_add(f16_mul(alpha, dotProduct), f16_mul(beta, C[j + i * ldc]));
        }
    }
}
