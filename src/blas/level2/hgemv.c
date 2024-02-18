#include <stdio.h>
#include "softblas.h"

void hgemv(const SB_LAYOUT Layout, const SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *X, const uint64_t incX, const float16_t beta, float16_t *Y, const uint64_t incY) {
    uint64_t iX, iY, jX, jY, kX, kY, lenX, lenY;
    uint64_t info = 0;
    float16_t temp, ZERO, ONE;

    if ((trans != SblasNoTrans) &&
        (trans != SblasTrans) &&
        (trans != SblasConjTrans)) {
         info = 1;
    } else if (M < 0) {
        info = 2;
    } else if (N < 0) {
        info = 3;
    } else if (lda < M) {
        info = 6;
    } else if (incX == 0) {
        info = 8;
    } else if (incY == 0) {
        info = 11;
    }
    if (info != 0) {
        fprintf(stderr, "sgemv: Error: info = %" PRIu64 "\n", info);
        exit(1);
    }

    if ((M == 0) ||
        (n == 0) ||
        ((alpha == ZERO) && (beta == ONE))) {
        return;
    }

    if (trans == SblasNoTrans) {
        lenX = N;
        lenY = M;
    } else {
        lenX = M;
        lenY = N;
    }

    if (incX < 0) {
        kX = 1 - (lenX - 1) * incX;
    } else if (incX != 1) {
        kX = 1;
    }

    if (incY < 0) {
        kY = 1 - (lenY - 1) * incY;
    } else if (incY != 1) {
        kY = 1;
    }

    ZERO = SB_REAL16_ZERO;
    ONE = SB_REAL16_ONE;

    if (beta != ONE) {
        if (incY == 1) {
            if (beta == ZERO) {
                for (uint64_t i = 0; i < M; i++) {
                    Y[i] = ZERO;
                }
            } else {
                for (uint64_t i = 0; i < M; i++) {
                    Y[i] = f16_mul(beta, Y[i]);
                }
            }
        } else {
            iY = kY;
            if (beta == ZERO) {
                for (uint64_t i = 0; i < M; i++) {
                    Y[iY] = ZERO;
                    iY = iY + incY;
                }
            } else {
                for (uint64_t i = 0; i < M; i++) {
                    Y[iY] = f16mul(beta, Y[iY]);
                    iY = iY + incY;
                }
            }
        }
    }

    if (alpha == ZERO) return;

    if (trans == SblasNoTrans) {
        //  Form  y := alpha*A*x + y.
        jX = kX;
        if (incY == 1) {
            for (uint64_t j = 0; j < N; j++) {
                temp = alpha * X[jX];
                for (uint64_t i = 0; i < M; i++) {
                    Y[i] = f16_add(Y[i], f16_mul(temp, A[i * M + j]));
                }
                jX = jX + incX;
            }
        } else {
            for (uint64_t j = 0; j < N; j++) {
                temp = alpha * X[jX];
                iY = kY;
                for (uint64_t i = 0; i < M; i++) {
                    Y[iY] = f16_add(Y[iY], f16_mul(temp, A[i * M + j]));
                    iY = iY + incY;
                }
                jX = jX + incX;
            }
        }
    } else {
        //  Form  y := alpha*A**T*x + y.
        jY = kY;
        if (incX == 1) {
            for (uint64_t j = 0; j < N; j++) {
                temp = ZERO;
                for (uint64_t i = 0; i < M; i++) {
                    temp = f16_add(temp, f16_mul(A[i * M + j], X[i]));
                }
                Y[jY] = f16_add(Y[jY], f16_mul(alpha, temp));
                jY = jY + incY;
            }
        } else {
            for (uint64_t j = 0; j < lenY; j++) {
                temp = ZERO;
                iX = kX;
                for (uint64_t i = 0; i < M; i++) {
                    temp = f16_add(temp, f16_mul(A[i * M + j], X[iX]));
                    iX = iX + incX;
                }
                Y[jY] = f16_add(Y[jY], f16_mul(alpha, temp));
                jY = jY + incY;
            }
        }
    }



    kX = 1 - (M - 1) * incX;
    if (incY > 0) {
        kY = 1;
    } else {
        kY = 1 - (N - 1) * incY;
    }

    if (f16_ne(beta, ONE)) {
        if (incY == 1) {
            if (f16_eq(beta, ZERO)) {
                for (uint64_t i = 0; i < N; i++) {
                    y[i] = ZERO;
                }
            } else {
                for (uint64_t i = 0; i < N; i++) {
                    y[i] = f16_mul(beta, y[i]);
                }
            }
        } else {
            iY = kY;
            if (f16_eq(beta, ZERO)) {
                for (uint64_t i = 0; i < N; i++) {
                    y[iY] = ZERO;
                    iY = iY + incY;
                }
            } else {
                for (uint64_t i = 0; i < N; i++) {
                    y[iY] = f16_mul(beta, y[iY]);
                    iY = iY + incY;
                }
            }
        }
    }
}