#ifndef TEST_H
#define TEST_H

#include <softfloat.h>
#include <softblas.h>

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

//  Test function prototypes for BLAS Level 1
MunitResult test_sasum_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_sasum_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sasum_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_saxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_saxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_saxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_saxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);
MunitResult test_scopy_all(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_scopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_sdot_0(const MunitParameter params[],
                        void* user_data_or_fixture);
MunitResult test_sdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_sdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture);
MunitResult test_snrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_snrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_snrm2_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_sscal_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_sscal_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sscal_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_sswap_two(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_sswap_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_isamax_0(const MunitParameter params[],
                          void* user_data_or_fixture);
MunitResult test_isamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_isamax_stride(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_isamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture);

MunitResult test_dasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_dasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_dasum_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_daxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_daxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_daxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_daxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);
MunitResult test_dcopy_all(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_dcopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_ddot_0(const MunitParameter params[],
                        void* user_data_or_fixture);
MunitResult test_ddot_12345(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_ddot_stride(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_ddot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture);
MunitResult test_dnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_dnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dnrm2_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_dscal_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_dscal_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dscal_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_dswap_two(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_dswap_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_idamax_0(const MunitParameter params[],
                          void* user_data_or_fixture);
MunitResult test_idamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_idamax_stride(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_idamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture);

MunitResult test_hasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_hasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_hasum_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_haxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_haxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_haxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_haxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);
MunitResult test_hcopy_all(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_hcopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_hdot_0(const MunitParameter params[],
                        void* user_data_or_fixture);
MunitResult test_hdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_hdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture);
MunitResult test_hnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_hnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hnrm2_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_hscal_0(const MunitParameter params[],
                        void* user_data_or_fixture);
MunitResult test_hscal_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hscal_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_hswap_two(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_hswap_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_ihamax_0(const MunitParameter params[],
                          void* user_data_or_fixture);
MunitResult test_ihamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_ihamax_stride(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_ihamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture);

MunitResult test_qasum_0(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_qasum_12345(const MunitParameter params[],
                       void* user_data_or_fixture);
MunitResult test_qasum_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qaxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_qaxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_qaxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qaxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);
MunitResult test_qcopy_all(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qcopy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qdot_0(const MunitParameter params[],
                        void* user_data_or_fixture);
MunitResult test_qdot_12345(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_qdot_stride(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qdot_neg_stride(const MunitParameter params[],
                                 void* user_data_or_fixture);
MunitResult test_qnrm2_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_qnrm2_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qnrm2_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qscal_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_qscal_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qscal_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qswap_two(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_qswap_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_iqamax_0(const MunitParameter params[],
                          void* user_data_or_fixture);
MunitResult test_iqamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_iqamax_stride(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_iqamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture);

//  Test function prototypes for BLAS Level 2
MunitResult test_sgemv_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemv_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemv_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemv_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_sgemv_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);


#endif // TEST_H
