#include "test.h"

//  C4 (transpose) + C5 (padded lda) coverage for gemv in double/half/quad,
//  mirroring test_sgemv_trans / test_sgemv_padlda. Expected results are built
//  with the same vec helper, so only the input encodings are precision-specific.

//  Transpose: A is 2x3 row-major [[1,2,3],[4,5,6]], x=[1,1] (length M=2).
//  y = A^T x has length N=3: [1+4, 2+5, 3+6] = [5,7,9].
MunitResult test_dgemv_trans(const MunitParameter params[], void* u) {
    const float64_t alpha = { SB_REAL64_ONE }, beta = { SB_REAL64_ZERO };
    float64_t* A = dvec((double[]){1,2,3, 4,5,6}, 6);
    float64_t* X = dvec((double[]){1, 1}, 2);
    float64_t* Y = dvec((double[]){0, 0, 0}, 3);
    dgemv('R', 'T', 2, 3, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float64_t* RY = dvec((double[]){5, 7, 9}, 3);
    for (uint64_t i = 0; i < 3; i++) assert_ullong(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
MunitResult test_hgemv_trans(const MunitParameter params[], void* u) {
    const float16_t alpha = { SB_REAL16_ONE }, beta = { SB_REAL16_ZERO };
    float16_t* A = hvec((uint16_t[]){0x3c00,0x4000,0x4200, 0x4400,0x4500,0x4600}, 6);  // 1..6
    float16_t* X = hvec((uint16_t[]){0x3c00, 0x3c00}, 2);
    float16_t* Y = hvec((uint16_t[]){0, 0, 0}, 3);
    hgemv('R', 'T', 2, 3, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float16_t* RY = hvec((uint16_t[]){0x4500, 0x4700, 0x4880}, 3);  // 5,7,9
    for (uint64_t i = 0; i < 3; i++) assert_ushort(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
MunitResult test_qgemv_trans(const MunitParameter params[], void* u) {
    const float128_t alpha = {{0x0ull,0x3fff000000000000ull}}, beta = {{0x0ull,0x0ull}};
    float128_t* A = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0},{.hi=0x4000000000000000,.lo=0},{.hi=0x4000800000000000,.lo=0},
        {.hi=0x4001000000000000,.lo=0},{.hi=0x4001400000000000,.lo=0},{.hi=0x4001800000000000,.lo=0}}, 6);
    float128_t* X = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0},{.hi=0x3fff000000000000,.lo=0}}, 2);
    float128_t* Y = qvec((float128_pair_t[]){{.hi=0,.lo=0},{.hi=0,.lo=0},{.hi=0,.lo=0}}, 3);
    qgemv('R', 'T', 2, 3, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float128_t* RY = qvec((float128_pair_t[]){
        {.hi=0x4001400000000000,.lo=0},{.hi=0x4001c00000000000,.lo=0},{.hi=0x4002200000000000,.lo=0}}, 3); // 5,7,9
    for (uint64_t i = 0; i < 3; i++) assert_ullong(Y[i].v[1], ==, RY[i].v[1]);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}

//  Padded lda: 2x2 logical A=[[1,2],[3,4]] stored with a pad column (lda=3),
//  x=[1,1]. y = A x = [1+2, 3+4] = [3,7]; the pad column (99) is skipped.
MunitResult test_dgemv_padlda(const MunitParameter params[], void* u) {
    const float64_t alpha = { SB_REAL64_ONE }, beta = { SB_REAL64_ZERO };
    float64_t* A = dvec((double[]){1,2,99, 3,4,99}, 6);
    float64_t* X = dvec((double[]){1, 1}, 2);
    float64_t* Y = dvec((double[]){0, 0}, 2);
    dgemv('R', 'N', 2, 2, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float64_t* RY = dvec((double[]){3, 7}, 2);
    for (uint64_t i = 0; i < 2; i++) assert_ullong(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
MunitResult test_hgemv_padlda(const MunitParameter params[], void* u) {
    const float16_t alpha = { SB_REAL16_ONE }, beta = { SB_REAL16_ZERO };
    float16_t* A = hvec((uint16_t[]){0x3c00,0x4000,0x5630, 0x4200,0x4400,0x5630}, 6);  // 1,2,99,3,4,99
    float16_t* X = hvec((uint16_t[]){0x3c00, 0x3c00}, 2);
    float16_t* Y = hvec((uint16_t[]){0, 0}, 2);
    hgemv('R', 'N', 2, 2, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float16_t* RY = hvec((uint16_t[]){0x4200, 0x4700}, 2);  // 3,7
    for (uint64_t i = 0; i < 2; i++) assert_ushort(Y[i].v, ==, RY[i].v);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
MunitResult test_qgemv_padlda(const MunitParameter params[], void* u) {
    const float128_t alpha = {{0x0ull,0x3fff000000000000ull}}, beta = {{0x0ull,0x0ull}};
    float128_t* A = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0},{.hi=0x4000000000000000,.lo=0},{.hi=0x4058c00000000000,.lo=0},
        {.hi=0x4000800000000000,.lo=0},{.hi=0x4001000000000000,.lo=0},{.hi=0x4058c00000000000,.lo=0}}, 6); // 1,2,99,3,4,99
    float128_t* X = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0},{.hi=0x3fff000000000000,.lo=0}}, 2);
    float128_t* Y = qvec((float128_pair_t[]){{.hi=0,.lo=0},{.hi=0,.lo=0}}, 2);
    qgemv('R', 'N', 2, 2, alpha, A, 3, X, 1, beta, Y, 1, 'n');
    float128_t* RY = qvec((float128_pair_t[]){
        {.hi=0x4000800000000000,.lo=0},{.hi=0x4001c00000000000,.lo=0}}, 2);  // 3,7
    for (uint64_t i = 0; i < 2; i++) assert_ullong(Y[i].v[1], ==, RY[i].v[1]);
    free(A); free(X); free(Y); free(RY);
    return MUNIT_OK;
}
