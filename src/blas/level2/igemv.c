#include "softblas.h"

void igemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const complex16_t alpha, const complex16_t *A, const uint64_t lda, const complex16_t *X, const int64_t incX, const complex16_t beta, complex16_t *Y, const uint64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const complex16_t ZERO = SB_COMPLEX16_ZERO;

    if (Layout != 'C' && Layout != 'c' && Layout != 'R' && Layout != 'r') {
        // Invalid order parameter
        return;
    }

    if (Trans != 'N' && Trans != 'n' && Trans != 'T' && Trans != 't' && Trans != 'C' && Trans != 'c') {
        // Invalid trans parameter
        return;
    }

    //  Reject non-positive strides: incX is signed, and a negative value in the
    //  unsigned index expression X[j * incX] would wrap to a huge OOB offset.
    if (incX < 1 || incY < 1) {
        return;
    }

    const bool colmajor = (Layout == 'C' || Layout == 'c');
    const bool trans = (Trans != 'N' && Trans != 'n');
    //  'C' is the conjugate transpose: index like 'T', but conjugate A.
    const bool conj = (Trans == 'C' || Trans == 'c');

    if (!trans) {
        //  y <- alpha * A * x + beta * y   (length M, inner length N)
        for (uint64_t i = 0; i < M; ++i) {
            complex16_t dotProduct = ZERO;
            for (uint64_t j = 0; j < N; ++j) {
                uint64_t idx = colmajor ? (i + j * lda) : (i * lda + j);
                dotProduct = c16_add(dotProduct, c16_mul(A[idx], X[j * incX]));
            }
            Y[i * incY] = nan_unify_i(c16_add(c16_mul(alpha, dotProduct), c16_mul(beta, Y[i * incY])));
        }
    } else {
        //  y <- alpha * op(A) * x + beta * y   (length N, inner length M)
        for (uint64_t i = 0; i < N; ++i) {
            complex16_t dotProduct = ZERO;
            for (uint64_t j = 0; j < M; ++j) {
                uint64_t idx = colmajor ? (j + i * lda) : (j * lda + i);
                complex16_t a = A[idx];
                if (conj) a = c16_conj(a);
                dotProduct = c16_add(dotProduct, c16_mul(a, X[j * incX]));
            }
            Y[i * incY] = nan_unify_i(c16_add(c16_mul(alpha, dotProduct), c16_mul(beta, Y[i * incY])));
        }
    }
}
