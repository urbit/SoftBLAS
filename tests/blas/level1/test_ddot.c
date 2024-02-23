#include "test.h"

MunitResult test_ddot_0(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 1);
    float64_t RZ = {*(uint64_t*)&(double){0.0}};

    assert_ullong(DZ.v, ==, RZ.v);

    free(DX);
    free(DY);

    return MUNIT_OK;
}

MunitResult test_ddot_12345(const MunitParameter params[],
                            void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){5.0, 4.0, 3.0, 2.0, 1.0}, 5);

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 1);
    float64_t RZ = {*(uint64_t*)&(double){35.0}};

    assert_ullong(DZ.v, ==, RZ.v);

    free(DX);
    free(DY);

    return MUNIT_OK;
}

MunitResult test_ddot_stride(const MunitParameter params[],
                             void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0}, 9);

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 2);
    float64_t RZ = {*(uint64_t*)&(double){15.0}};

    assert_ullong(DZ.v, ==, RZ.v);

    free(DX);
    free(DY);

    return MUNIT_OK;
}

MunitResult test_ddot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){5.0, 4.0, 3.0, 2.0, 1.0}, 5);

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, -1);
    float64_t RZ = {*(uint64_t*)&(double){55.0}};

    assert_ullong(DZ.v, ==, RZ.v);

    free(DX);
    free(DY);

    return MUNIT_OK;
}