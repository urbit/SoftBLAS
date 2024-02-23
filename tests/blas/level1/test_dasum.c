#include "test.h"

MunitResult test_dasum_0(const MunitParameter params[],
                         void* user_data_or_fixture) {
    const uint64_t N = 2;
    float64_t* DX = dvec((double[]){0.0, 0.0}, N);

    float64_t D = (float64_t) dasum(N, (float64_t*)DX, 1);
    float64_t R = (float64_t){0.0};

    assert_int(D.v, ==, R.v);
    
    return MUNIT_OK;
}

MunitResult test_dasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture) {
    const uint64_t N = 5;
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, N);

    float64_t D = (float64_t) dasum(N, (float64_t*)DX, 1);
    float64_t R = {*(uint64_t*)&(double){15.0}};
    
    assert_int(D.v, ==, R.v);
    
    return MUNIT_OK;
}

// TODO test stride
