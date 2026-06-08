#include "softblas.h"

void cgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const complex32_t alpha, const complex32_t *A, const uint64_t lda, const complex32_t *B, const uint64_t ldb, const complex32_t beta, complex32_t *C, const uint64_t ldc, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const complex32_t ZERO = SB_COMPLEX32_ZERO;

    if (transA != 'N' && transA != 'n' && transA != 'T' && transA != 't' && transA != 'C' && transA != 'c') {
        // Invalid transA parameter
        return;
    }

    if (transB != 'N' && transB != 'n' && transB != 'T' && transB != 't' && transB != 'C' && transB != 'c') {
        // Invalid transB parameter
        return;
    }

    for (uint64_t i = 0; i < M; i++) {
        for (uint64_t j = 0; j < P; j++) {
            complex32_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldb : k + j * ldb;
                complex32_t a = A[indexA];
                complex32_t b = B[indexB];
                //  'C' is the conjugate transpose: same index as 'T', conjugated.
                if (transA == 'C' || transA == 'c') a = c32_conj(a);
                if (transB == 'C' || transB == 'c') b = c32_conj(b);
                dotProduct = c32_add(dotProduct, c32_mul(a, b));
            }
            C[j + i * ldc] = nan_unify_c(c32_add(c32_mul(alpha, dotProduct), c32_mul(beta, C[j + i * ldc])));
        }
    }
}
