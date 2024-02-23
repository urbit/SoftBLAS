#include "test.h"

MunitResult test_idamax_0(const MunitParameter params[],
                          void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){0.0, 0.0, 0.0}, 3);

    uint64_t I = idamax(3, DX, 1);
    uint64_t R = 0;

    assert_int(I, ==, R);
    
    free(DX);

    return MUNIT_OK;
}

MunitResult test_idamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);

    uint64_t I = idamax(5, DX, 1);
    uint64_t R = 4;

    assert_int(I, ==, R);
    
    free(DX);

    return MUNIT_OK;
}

MunitResult test_idamax_stride(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);

    uint64_t I = idamax(2, DX, 3);
    // note that idamax counts against the stride index, not the overall index
    uint64_t R = 1;

    assert_int(I, ==, R);
    
    free(DX);

    return MUNIT_OK;
}

MunitResult test_idamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 3.0, 5.0, 4.0, 2.0}, 5);

    uint64_t I = idamax(5, DX, 1);
    uint64_t R = 2;

    assert_int(I, ==, R);
    
    free(DX);

    return MUNIT_OK;
}
