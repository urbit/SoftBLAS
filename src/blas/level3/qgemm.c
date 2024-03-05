#include "softblas.h"
#include <stdio.h>

void qgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float128_t alpha, const float128_t *A, const uint64_t lda, const float128_t *B, const uint64_t ldb, const float128_t beta, float128_t *C, const uint64_t ldc) {
    const float128_t ZERO = { SB_REAL128L_ZERO, SB_REAL128U_ZERO };

    if (transA != 'N' && transA != 'n' && transA != 'T' && transA != 't') {
        printf("Invalid transA parameter\n");
        return;
    }

    if (transB != 'N' && transB != 'n' && transB != 'T' && transB != 't') {
        printf("Invalid transB parameter\n");
        return;
    }

    float128_t qtemp;
    float128_t qtemp2;
    for (uint64_t i = 0; i < M; i++) {
        for (uint64_t j = 0; j < P; j++) {
            float128_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldc : k + j * ldc;
                float128_t a = A[indexA];
                float128_t b = B[indexB];
                f128M_mul(&a, &b, &qtemp);
                f128M_add(&dotProduct, &qtemp, &dotProduct);
            }
            f128M_mul(&alpha, &dotProduct, &qtemp);
            f128M_mul(&beta, &C[j + i * ldc], &qtemp2);
            f128M_add(&qtemp, &qtemp2, &(C[j + i * ldc]));
        }
    }
}
