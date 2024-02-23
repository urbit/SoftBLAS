#include "test.h"

int main(int argc, char* argv[MUNIT_ARRAY_PARAM(argc + 1)]) {
    MunitSuite suite = {
        "/blas/level1",
        (MunitTest[]){
            {"/test_sasum_0", test_sasum_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sasum_12345", test_sasum_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sasum_stride", test_sasum_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_0", test_saxpy_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_sum", test_saxpy_sum, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_stride", test_saxpy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_saxpy_neg_stride", test_saxpy_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_scopy_all", test_scopy_all, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_scopy_stride", test_scopy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sdot_0", test_sdot_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sdot_12345", test_sdot_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sdot_stride", test_sdot_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_sdot_neg_stride", test_sdot_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            
            {"/test_dasum_0", test_dasum_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_dasum_12345", test_dasum_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_dasum_stride", test_dasum_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_daxpy_0", test_daxpy_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_daxpy_sum", test_daxpy_sum, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_daxpy_stride", test_daxpy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_daxpy_neg_stride", test_daxpy_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_dcopy_all", test_dcopy_all, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_dcopy_stride", test_dcopy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_ddot_0", test_ddot_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_ddot_12345", test_ddot_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_ddot_stride", test_ddot_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_ddot_neg_stride", test_ddot_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            
            {"/test_hasum_0", test_hasum_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hasum_12345", test_hasum_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hasum_stride", test_hasum_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_haxpy_0", test_haxpy_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_haxpy_sum", test_haxpy_sum, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_haxpy_stride", test_haxpy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_haxpy_neg_stride", test_haxpy_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hcopy_all", test_hcopy_all, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hcopy_stride", test_hcopy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hdot_0", test_hdot_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hdot_12345", test_hdot_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hdot_stride", test_hdot_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_hdot_neg_stride", test_hdot_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},

            {"/test_qasum_0", test_qasum_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qasum_12345", test_qasum_12345, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qasum_stride", test_qasum_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qaxpy_0", test_qaxpy_0, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qaxpy_sum", test_qaxpy_sum, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qaxpy_stride", test_qaxpy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qaxpy_neg_stride", test_qaxpy_neg_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qcopy_all", test_qcopy_all, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},
            {"/test_qcopy_stride", test_qcopy_stride, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL},

            {NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL}
        },
        NULL,
        1,
        MUNIT_SUITE_OPTION_NONE
    };

    return munit_suite_main(&suite, NULL, argc, argv);
}
