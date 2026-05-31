#include "test.h"

MunitResult test_ddot_0(const MunitParameter params[],
                           void* user_data_or_fixture) {
    float64_t* DX = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);
    float64_t* DY = dvec((double[]){0.0, 0.0, 0.0, 0.0, 0.0}, 5);

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 1, 'n');
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

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 1, 'n');
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

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, 2, 'n');
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

    float64_t DZ = ddot(5, (float64_t*)DX, 1, DY, -1, 'n');
    float64_t RZ = {*(uint64_t*)&(double){55.0}};

    assert_ullong(DZ.v, ==, RZ.v);

    free(DX);
    free(DY);

    return MUNIT_OK;
}

//  C6: all five rounding modes, double precision. dot([1, 2^-53],[1, 1]) =
//  1 + 2^-53, an exact tie; the negative case distinguishes z/d/u/a.
MunitResult test_ddot_rounding_modes(const MunitParameter params[], void* u) {
    const uint64_t ONE = 0x3ff0000000000000ull, UP = 0x3ff0000000000001ull,
                   NEG_ONE = 0xbff0000000000000ull, NEG_UP = 0xbff0000000000001ull;
    struct { char m; uint64_t want; } pos[] = {
        {'n', ONE}, {'z', ONE}, {'d', ONE}, {'u', UP}, {'a', UP} };
    for (uint64_t k = 0; k < 5; k++) {
        float64_t* X = dvec((double[]){ 1.0, 0x1p-53 }, 2);
        float64_t* Y = dvec((double[]){ 1.0, 1.0 }, 2);
        assert_ullong(ddot(2, X, 1, Y, 1, pos[k].m).v, ==, pos[k].want);
        free(X); free(Y);
    }
    struct { char m; uint64_t want; } neg[] = {
        {'n', NEG_ONE}, {'z', NEG_ONE}, {'u', NEG_ONE}, {'d', NEG_UP}, {'a', NEG_UP} };
    for (uint64_t k = 0; k < 5; k++) {
        float64_t* X = dvec((double[]){ -1.0, -0x1p-53 }, 2);
        float64_t* Y = dvec((double[]){ 1.0, 1.0 }, 2);
        assert_ullong(ddot(2, X, 1, Y, 1, neg[k].m).v, ==, neg[k].want);
        free(X); free(Y);
    }
    return MUNIT_OK;
}