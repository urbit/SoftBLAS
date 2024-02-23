#ifndef TEST_H
#define TEST_H

#include <softfloat.h>
#include <softblas.h>

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

//  Test function prototypes
MunitResult test_haxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_haxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_haxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_haxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);

MunitResult test_sasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_sasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_saxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_saxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_saxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_saxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);

MunitResult test_dasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_dasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_daxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_daxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_daxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_daxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);

MunitResult test_qasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_qasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_qaxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_qaxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_qaxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qaxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);

#endif // TEST_H
