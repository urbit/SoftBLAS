#include "test.h"

MunitResult test_dnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){0.0, 0.0}, 2);

    float64_t D = (float64_t) dnrm2(2, (float64_t*)DX, 1);
    float64_t R = (float64_t){0.0};

    assert_ullong(D.v, ==, R.v);
    
    free(DX);

    return MUNIT_OK;
}

MunitResult test_dnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);

    float64_t D = (float64_t) dnrm2(5, (float64_t*)DX, 1);
    float64_t R = {*(uint64_t*)&(double){7.416198487095663}};  // sqrt(55.0)

    assert_ullong(D.v, ==, R.v);
    
    free(DX);

    return MUNIT_OK;
}

MunitResult test_dnrm2_stride(const MunitParameter params[],
                              void *user_data) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0}, 9);

    float64_t D = (float64_t) dnrm2(5, (float64_t*)DX, 2);
    float64_t R = {*(uint64_t*)&(double){2.23606797749979}};  // sqrt(5.0)

    assert_ullong(D.v, ==, R.v);

    free(DX);

    return MUNIT_OK;
}
