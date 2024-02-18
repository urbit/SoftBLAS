#include <stdio.h>
#include "softblas.h"

void hgemm(const SB_LAYOUT Layout, const SB_TRANSPOSE transA, const SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *B, const uint64_t ldb, const float16_t beta, float16_t *C, const uint64_t ldc) {
    uint64_t nColA, nRowA, nRowB;
    bool64_t notA, notB;
    uint64_t info = 0;
    float16_t temp, ZERO, ONE;

    ZERO = SB_REAL16_ZERO;
    ONE = SB_REAL16_ONE;

    notA = (transA == SblasNoTrans);
    notB = (transB == SblasNoTrans);
    if notA {
        nRowA = M;
        nColA = K;
    } else {
        nRowA = K;
        nColA = M;
    }
    if notB {
        nRowB = K;
    } else {
        nRowB = N;
    }

    if (!notA &&
        (transA != SblasTrans) &&
        (transA != SblasConjTrans)) {
        info = 1;
    } else if (!notB &&
               (transB != SblasTrans) &&
               (transB != SblasConjTrans)) {
        info = 2;
    } else if M < 0 {
        info = 3;
    } else if N < 0 {
        info = 4;
    } else if K < 0 {
        info = 5;
    } else if lda < nRowA {
        info = 8;
    } else if ldb < nRowB {
        info = 10;
    } else if ldc < M {
        info = 13;
    }
    if (info != 0) {
        fprintf(stderr, "sgemm: Error: info = %" PRIu64 "\n", info);
        exit(1);
    }

    if (M == 0 ||
        N == 0 ||
        ((f16_eq(alpha, ZERO) || f16_eq(K, ZERO)) && f16_eq(beta, ONE))) {
        return;
    }

    if f16_eq(alpha, ZERO) {
        if f16_eq(beta, ZERO) {
            for (uint64_t j = 0; j < N; j++) {
                for (uint64_t i = 0; i < M; i++) {
                    C[i * ldc + j] = ZERO;
                }
            }
        } else {
            for (uint64_t j = 0; j < N; j++) {
                for (uint64_t i = 0; i < M; i++) {
                    C[i * ldc + j] = f16_mul(beta, C[i * ldc + j]);
                }
            }
        }
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

    if f16_eq(alpha, ZERO) {
        if f16_eq(beta, ZERO) {
            for (uint64_t i = 0; i < M; i++) {
                C[i*M + j] = ZERO;
            }
        } else {
            for (uint64_t i = 0; i < M; i++) {
                C[i*M + j] = f16_mul(beta, Y[i*M + j]);
            }
        }
    }

    if notB {
        if notA {
            // Form  C := alpha*A*B + beta*C.
            for (uint64_t j = 1; j < N; j++) {
                if f16_eq(beta, ZERO) {
                    for (uint64_t i = 0; i < M; i++) {
                        C[i*M + j] = ZERO;
                    }
                } else if f16_ne(beta, ONE) {
                    for (uint64_t i = 0; i < M; i++) {
                        C[i*M + j] = f16_mul(beta, C[i*M + j]);
                    }
                }
                for (uint64_t l = 1; l < K; l++) {
                    temp = f16_mul(alpha, B[l * ldb + j]);
                    for (uint64_t i = 0; i < M; i++) {
                        C[i*M + j] = f16_add(C[i*M + j], f16_mul(temp, A[i * lda + l]));
                    }
                }
            }
        } else {
            //  Form  C := alpha*A**T*B + beta*C
            for (uint64_t j = 0; j < N; j++) {
                for (uint64_t i = 0; i < M; i++) {
                    temp = ZERO;
                    for (uint64_t l = 0; l < K; l++) {
                        temp = f16_add(temp, f16_mul(A[l * lda + i], B[l * ldb + j]));
                    }
                    if f16_eq(beta, ZERO) {
                        C[i * ldc + j] = f16_mul(alpha, temp);
                    } else {
                        C[i * ldc + j] = f16_add(f16_mul(alpha, temp), f16_mul(beta, C[i * ldc + j]));
                    }
                }
            }
        }
    } else {
        if notA {
            //  Form  C := alpha*A*B**T + beta*C
            for (uint64_t j = 0; j < N; j++) {
                if f16_eq(beta, ZERO) {
                    for (uint64_t i = 0; i < M; i++) {
                        C[i * ldc + j] = ZERO;
                    }
                } else if f16_ne(beta, ONE) {
                    for (uint64_t i = 0; i < M; i++) {
                        C[i * ldc + j] = f16_mul(beta, C[i * ldc + j]);
                    }
                }
                for (uint64_t l = 0; l < K; l++) {
                    temp = f16_mul(alpha, B[j * ldb + l]);
                    for (uint64_t i = 0; i < M; i++) {
                        C[i * ldc + j] = f16_add(C[i * ldc + j], f16_mul(temp, A[i * lda + l]));
                    }
                }
            }
        } else {
            //  Form  C := alpha*A**T*B**T + beta*C
            for (uint64_t j = 0; j < N; j++) {
                for (uint64_t i = 0; i < M; i++) {
                    temp = ZERO;
                    for (uint64_t l = 0; l < K; l++) {
                        temp = f16_add(temp, f16_mul(A[l * lda + i], B[j * ldb + l]));
                    }
                    if f16_eq(beta, ZERO) {
                        C[i * ldc + j] = f16_mul(alpha, temp);
                    } else {
                        C[i * ldc + j] = f16_add(f16_mul(alpha, temp), f16_mul(beta, C[i * ldc + j]));
                    }
                }
            }
        }
    }
}