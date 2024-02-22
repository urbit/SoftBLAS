#include "test.h"

int main(int argc, char* argv[MUNIT_ARRAY_PARAM(argc + 1)]) {
    MunitSuite suite = {
        "/blas/level1",
        (MunitTest[]){
            {"/test_sasum_0", test_sasum_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sasum_12345", test_sasum_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_0", test_saxpy_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_sum", test_saxpy_sum, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_stride", test_saxpy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_neg_stride", test_saxpy_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}
        },
        NULL,
        1,
        MUNIT_SUITE_OPTION_NONE
    };

    return munit_suite_main(&suite, NULL, argc, argv);
}
