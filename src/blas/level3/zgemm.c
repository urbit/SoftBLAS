#include "softblas.h"

void zgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const complex64_t alpha, const complex64_t *A, const uint64_t lda, const complex64_t *B, const uint64_t ldb, const complex64_t beta, complex64_t *C, const uint64_t ldc, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const complex64_t ZERO = SB_COMPLEX64_ZERO;

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
            complex64_t dotProduct = ZERO;
            for (uint64_t k = 0; k < N; k++) {
                uint64_t indexA = (transA == 'N' || transA == 'n') ? k + i * lda : i + k * lda;
                uint64_t indexB = (transB == 'N' || transB == 'n') ? j + k * ldb : k + j * ldb;
                complex64_t a = A[indexA];
                complex64_t b = B[indexB];
                //  'C' is the conjugate transpose: same index as 'T', conjugated.
                if (transA == 'C' || transA == 'c') a = c64_conj(a);
                if (transB == 'C' || transB == 'c') b = c64_conj(b);
                dotProduct = c64_add(dotProduct, c64_mul(a, b));
            }
            C[j + i * ldc] = nan_unify_z(c64_add(c64_mul(alpha, dotProduct), c64_mul(beta, C[j + i * ldc])));
        }
    }
}
