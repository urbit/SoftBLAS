#include "softblas.h"

void igemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const complex16_t alpha, const complex16_t *A, const uint64_t lda, const complex16_t *B, const uint64_t ldb, const complex16_t beta, complex16_t *C, const uint64_t ldc, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const complex16_t ZERO = SB_COMPLEX16_ZERO;

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
            complex16_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldb : k + j * ldb;
                complex16_t a = A[indexA];
                complex16_t b = B[indexB];
                //  'C' is the conjugate transpose: same index as 'T', conjugated.
                if (transA == 'C' || transA == 'c') a = c16_conj(a);
                if (transB == 'C' || transB == 'c') b = c16_conj(b);
                dotProduct = c16_add(dotProduct, c16_mul(a, b));
            }
            C[j + i * ldc] = nan_unify_i(c16_add(c16_mul(alpha, dotProduct), c16_mul(beta, C[j + i * ldc])));
        }
    }
}
