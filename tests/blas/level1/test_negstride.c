#include "test.h"

//  C2: negative-stride coverage for the swap and copy families (dot/axpy
//  already have it). A negative incX uses the offset start (-N+1)*incX, so a
//  stride of -1 reverses the indexed vector. Inputs 1,2,3 / 4,5,6 in each
//  precision; expected results are built with the same vec helper.

//  swap(3, X, -1, Y, 1): X is indexed in reverse, Y forward -> X gets reverse
//  of Y and Y gets reverse of X. X=[1,2,3],Y=[4,5,6] -> X=[6,5,4],Y=[3,2,1].
MunitResult test_dswap_negstride(const MunitParameter params[], void* u) {
    float64_t* X = dvec((double[]){1.0, 2.0, 3.0}, 3);
    float64_t* Y = dvec((double[]){4.0, 5.0, 6.0}, 3);
    dswap(3, X, -1, Y, 1, 'n');
    float64_t* RX = dvec((double[]){6.0, 5.0, 4.0}, 3);
    float64_t* RY = dvec((double[]){3.0, 2.0, 1.0}, 3);
    for (uint64_t i = 0; i < 3; i++) {
        assert_ullong(X[i].v, ==, RX[i].v);
        assert_ullong(Y[i].v, ==, RY[i].v);
    }
    free(X); free(Y); free(RX); free(RY);
    return MUNIT_OK;
}
MunitResult test_hswap_negstride(const MunitParameter params[], void* u) {
    float16_t* X = hvec((uint16_t[]){0x3c00, 0x4000, 0x4200}, 3);  // 1,2,3
    float16_t* Y = hvec((uint16_t[]){0x4400, 0x4500, 0x4600}, 3);  // 4,5,6
    hswap(3, X, -1, Y, 1, 'n');
    float16_t* RX = hvec((uint16_t[]){0x4600, 0x4500, 0x4400}, 3); // 6,5,4
    float16_t* RY = hvec((uint16_t[]){0x4200, 0x4000, 0x3c00}, 3); // 3,2,1
    for (uint64_t i = 0; i < 3; i++) {
        assert_ushort(X[i].v, ==, RX[i].v);
        assert_ushort(Y[i].v, ==, RY[i].v);
    }
    free(X); free(Y); free(RX); free(RY);
    return MUNIT_OK;
}
MunitResult test_qswap_negstride(const MunitParameter params[], void* u) {
    float128_t* X = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0}, {.hi=0x4000000000000000,.lo=0},
        {.hi=0x4000800000000000,.lo=0}}, 3);  // 1,2,3
    float128_t* Y = qvec((float128_pair_t[]){
        {.hi=0x4001000000000000,.lo=0}, {.hi=0x4001400000000000,.lo=0},
        {.hi=0x4001800000000000,.lo=0}}, 3);  // 4,5,6
    qswap(3, X, -1, Y, 1, 'n');
    float128_t* RX = qvec((float128_pair_t[]){
        {.hi=0x4001800000000000,.lo=0}, {.hi=0x4001400000000000,.lo=0},
        {.hi=0x4001000000000000,.lo=0}}, 3);  // 6,5,4
    float128_t* RY = qvec((float128_pair_t[]){
        {.hi=0x4000800000000000,.lo=0}, {.hi=0x4000000000000000,.lo=0},
        {.hi=0x3fff000000000000,.lo=0}}, 3);  // 3,2,1
    for (uint64_t i = 0; i < 3; i++) {
        assert_ullong(X[i].v[1], ==, RX[i].v[1]);
        assert_ullong(Y[i].v[1], ==, RY[i].v[1]);
    }
    free(X); free(Y); free(RX); free(RY);
    return MUNIT_OK;
}

//  copy(3, X, -1, Y, 1): Y = reverse(X). X=[1,2,3] -> Y=[3,2,1].
MunitResult test_scopy_negstride(const MunitParameter params[], void* u) {
    float32_t* X = svec((float[]){1.0f, 2.0f, 3.0f}, 3);
    float32_t* Y = svec((float[]){0.0f, 0.0f, 0.0f}, 3);
    scopy(3, X, -1, Y, 1, 'n');
    float32_t* R = svec((float[]){3.0f, 2.0f, 1.0f}, 3);
    for (uint64_t i = 0; i < 3; i++) assert_ulong(Y[i].v, ==, R[i].v);
    free(X); free(Y); free(R);
    return MUNIT_OK;
}
MunitResult test_dcopy_negstride(const MunitParameter params[], void* u) {
    float64_t* X = dvec((double[]){1.0, 2.0, 3.0}, 3);
    float64_t* Y = dvec((double[]){0.0, 0.0, 0.0}, 3);
    dcopy(3, X, -1, Y, 1, 'n');
    float64_t* R = dvec((double[]){3.0, 2.0, 1.0}, 3);
    for (uint64_t i = 0; i < 3; i++) assert_ullong(Y[i].v, ==, R[i].v);
    free(X); free(Y); free(R);
    return MUNIT_OK;
}
MunitResult test_hcopy_negstride(const MunitParameter params[], void* u) {
    float16_t* X = hvec((uint16_t[]){0x3c00, 0x4000, 0x4200}, 3);  // 1,2,3
    float16_t* Y = hvec((uint16_t[]){0, 0, 0}, 3);
    hcopy(3, X, -1, Y, 1, 'n');
    float16_t* R = hvec((uint16_t[]){0x4200, 0x4000, 0x3c00}, 3);  // 3,2,1
    for (uint64_t i = 0; i < 3; i++) assert_ushort(Y[i].v, ==, R[i].v);
    free(X); free(Y); free(R);
    return MUNIT_OK;
}
MunitResult test_qcopy_negstride(const MunitParameter params[], void* u) {
    float128_t* X = qvec((float128_pair_t[]){
        {.hi=0x3fff000000000000,.lo=0}, {.hi=0x4000000000000000,.lo=0},
        {.hi=0x4000800000000000,.lo=0}}, 3);  // 1,2,3
    float128_t* Y = qvec((float128_pair_t[]){
        {.hi=0,.lo=0}, {.hi=0,.lo=0}, {.hi=0,.lo=0}}, 3);
    qcopy(3, X, -1, Y, 1, 'n');
    float128_t* R = qvec((float128_pair_t[]){
        {.hi=0x4000800000000000,.lo=0}, {.hi=0x4000000000000000,.lo=0},
        {.hi=0x3fff000000000000,.lo=0}}, 3);  // 3,2,1
    for (uint64_t i = 0; i < 3; i++) assert_ullong(Y[i].v[1], ==, R[i].v[1]);
    free(X); free(Y); free(R);
    return MUNIT_OK;
}
