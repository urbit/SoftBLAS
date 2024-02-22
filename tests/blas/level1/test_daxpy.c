#include "test.h"

MunitResult test_daxpy_0(const MunitParameter params[], void *user_data) {
    const uint64_t N = 5;
    const float64_t DA = {*(uint64_t*)&(double){0.0}};
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, 5);
    float64_t* DY = dvec((double[]){10.0, 20.0, 30.0, 40.0, 50.0}, 5);
    float64_t* RY = dvec((double[]){10.0, 20.0, 30.0, 40.0, 50.0}, 5);

    daxpy(N, DA, DX, 1, DY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_daxpy_sum(const MunitParameter params[],
                           void *user_data) {
    const uint64_t N = 5;
    const float64_t DA = {*(uint64_t*)&(double){2.0}};
    float64_t* DX = dvec((double[]){1.0, -2.0, 3.0, -4.0, 5.0}, N);
    float64_t* DY = dvec((double[]){10.0, 20.0, 30.0, 40.0, 50.0}, N);
    float64_t* RY = dvec((double[]){12.0, 16.0, 36.0, 32.0, 60.0}, N);

    daxpy(N, DA, DX, 1, DY, 1);

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_daxpy_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float64_t DA = {*(uint64_t*)&(double){2.0}};
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, N);
    float64_t* DY = dvec((double[]){1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0}, 9);
    float64_t* RY = dvec((double[]){3.0, 2.0, 5.0, 2.0, 7.0, 2.0, 9.0, 2.0, 11.0}, 9);

    daxpy(N, DA, DX, 1, DY, 2);

    for (uint64_t i = 0; i < 9; i++) {
        assert_ulong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}

MunitResult test_daxpy_neg_stride(const MunitParameter params[],
                              void *user_data) {
    const uint64_t N = 5;
    const float64_t DA = {*(uint64_t*)&(double){2.0}};
    float64_t* DX = dvec((double[]){1.0, 2.0, 3.0, 4.0, 5.0}, N);
    float64_t* DY = dvec((double[]){5.0, 4.0, 3.0, 2.0, 1.0}, N);
    float64_t* RY = dvec((double[]){15.0, 12.0, 9.0, 6.0, 3.0}, N);

    daxpy(N, DA, DX, 1, DY, -1);

    for (uint64_t i = 0; i < N; i++) {
        assert_ulong(DY[i].v, ==, RY[i].v);
    }

    free(DX);
    free(DY);
    free(RY);

    return MUNIT_OK;
}
