#include "softblas.h"

void zgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const complex64_t alpha, const complex64_t *A, const uint64_t lda, const complex64_t *X, const int64_t incX, const complex64_t beta, complex64_t *Y, const uint64_t incY, const uint_fast8_t rndMode) {
    _set_rounding(rndMode);
    const complex64_t ZERO = SB_COMPLEX64_ZERO;

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
            complex64_t dotProduct = ZERO;
            for (uint64_t j = 0; j < N; ++j) {
                uint64_t idx = colmajor ? (i + j * lda) : (i * lda + j);
                dotProduct = c64_add(dotProduct, c64_mul(A[idx], X[j * incX]));
            }
            Y[i * incY] = nan_unify_z(c64_add(c64_mul(alpha, dotProduct), c64_mul(beta, Y[i * incY])));
        }
    } else {
        //  y <- alpha * op(A) * x + beta * y   (length N, inner length M)
        for (uint64_t i = 0; i < N; ++i) {
            complex64_t dotProduct = ZERO;
            for (uint64_t j = 0; j < M; ++j) {
                uint64_t idx = colmajor ? (j + i * lda) : (j * lda + i);
                complex64_t a = A[idx];
                if (conj) a = c64_conj(a);
                dotProduct = c64_add(dotProduct, c64_mul(a, X[j * incX]));
            }
            Y[i * incY] = nan_unify_z(c64_add(c64_mul(alpha, dotProduct), c64_mul(beta, Y[i * incY])));
        }
    }
}
