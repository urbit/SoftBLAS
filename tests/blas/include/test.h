#ifndef TEST_H
#define TEST_H

#include <softfloat.h>
#include <softblas.h>

#define MUNIT_ENABLE_ASSERT_ALIASES
#include "munit.h"

//  Test-only helpers: build float{16,32,64,128}_t arrays from literals.
//  (Heap-allocating; never used by the SoftBLAS library itself.)
/*  Convert an array of half-precision floats to an array of float16_t.
    The array of floats is assumed to be in little-endian format.

    Example usage:
        float16_t* HX = hvec((uint16_t[]){
            0x3C00,  // 1.0
            0xBC00,  // -1.0
            0x4248,  // 3.14
        }, 3);
        free(HX);
 */
static inline float16_t* hvec(uint16_t values[], uint64_t size) {
    float16_t* result = malloc(size * sizeof(float16_t));
    if (result == NULL) abort();
    
    for (uint64_t i = 0; i < size; i++) {
        result[i].v = values[i];
    }
    
    return result;
}

/*  Convert an array of single-precision floats to an array of float32_t.
    The array of floats is assumed to be in little-endian format. 

    Example usage:
        float32_t* SX = svec((float[]){
            1.0f,
            -1.0f,
            3.14159274f
        }, 3);
        free(SX);
 */
static inline float32_t* svec(float values[], uint64_t size) {
    float32_t* result = malloc(size * sizeof(float32_t));
    if (result == NULL) abort();
    
    for (uint64_t i = 0; i < size; i++) {
        result[i].v = *(uint32_t*)&values[i];
    }
    
    return result;
}

/*  Convert an array of double-precision floats to an array of float64_t.
    The array of floats is assumed to be in little-endian format. 

    Example usage:
        float64_t* DX = dvec((double[]){
            1.0,
            -1.0,
            3.1415926535897931
        }, 3);
        free(DX);
 */
static inline float64_t* dvec(double values[], uint64_t size) {
    float64_t* result = malloc(size * sizeof(float64_t));
    if (result == NULL) abort();
    
    for (uint64_t i = 0; i < size; i++) {
        result[i].v = *(uint64_t*)&values[i];
    }
    
    return result;
}

typedef struct {
    uint64_t hi;
    uint64_t lo;
} float128_pair_t;

/*  Convert an array of quad-precision float pairs to an array of float128_t.
    The array of floats is assumed to be in little-endian format. 

    Example usage:
        float128_t* QX = qvec((float128_pair_t[]){
            {.hi = 0x3FF0000000000000, .lo = 0x0000000000000000},  // 1.0
            {.hi = 0xBFF0000000000000, .lo = 0x0000000000000000},  // -1.0
            {.hi = 0x400921FB54442D18, .lo = 0x469898cc51701b80},  // 3.141592653589793238462643383279502797479068098137295573004504331874296718662975536062731407582759857177734375
        }, 3);
        free(QX);
 */
static inline float128_t* qvec(float128_pair_t pairs[], uint64_t size) {
    float128_t* result = malloc(size * sizeof(float128_t));
    if (result == NULL) abort();
    
    for (uint64_t i = 0; i < size; i++) {
        result[i].v[0] = pairs[i].lo;
        result[i].v[1] = pairs[i].hi;
    }
    
    return result;
}

//  Build an array of complex32_t from interleaved (re, im) float literals;
//  size is the number of complex elements.
static inline complex32_t* cvec(float vals[], uint64_t size) {
    complex32_t* result = malloc(size * sizeof(complex32_t));
    if (result == NULL) abort();
    for (uint64_t i = 0; i < size; i++) {
        result[i].real.v = *(uint32_t*)&vals[2*i];
        result[i].imag.v = *(uint32_t*)&vals[2*i+1];
    }
    return result;
}

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
MunitResult test_saxpy_rounding(const MunitParameter params[],
                                void* user_data_or_fixture);
MunitResult test_sasum_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_sdot_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_snrm2_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_saxpy_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_scopy_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_sscal_zero(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_saxpy_rounding_modes(const MunitParameter params[],
                 void* user_data_or_fixture);
MunitResult test_qasum_layout(const MunitParameter params[],
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
MunitResult test_sswap_zero(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_sswap_negstride(const MunitParameter params[],
                                 void* user_data_or_fixture);
MunitResult test_isamax_0(const MunitParameter params[],
                          void* user_data_or_fixture);
MunitResult test_isamax_12345(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_isamax_stride(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_isamax_13542(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_isamax_one(const MunitParameter params[],
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
MunitResult test_qasum_negpi(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qaxpy_0(const MunitParameter params[],
                         void* user_data_or_fixture);
MunitResult test_qaxpy_sum(const MunitParameter params[],
                           void* user_data_or_fixture);
MunitResult test_qaxpy_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_qaxpy_neg_stride(const MunitParameter params[],
                                  void* user_data_or_fixture);
MunitResult test_qaxpy_nan_unify(const MunitParameter params[],
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
MunitResult test_sgemv_incx(const MunitParameter params[],
                            void* user_data_or_fixture);
MunitResult test_sgemv_badlayout(const MunitParameter params[],
                                 void* user_data_or_fixture);

MunitResult test_dgemv_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemv_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemv_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemv_stride(const MunitParameter params[],
                              void* user_data_or_fixture);
MunitResult test_dgemv_incx(const MunitParameter params[],
                            void* user_data_or_fixture);

MunitResult test_hgemv_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hgemv_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hgemv_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hgemv_stride(const MunitParameter params[],
                              void* user_data_or_fixture);

MunitResult test_qgemv_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qgemv_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qgemv_12345(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qgemv_stride(const MunitParameter params[],
                              void* user_data_or_fixture);

//  Test function prototypes for BLAS Level 3
MunitResult test_sgemm_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemm_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemm_3x2x1_0(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_sgemm_5x4x3(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_sgemm_ldb(const MunitParameter params[],
                           void* user_data_or_fixture);

MunitResult test_dgemm_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemm_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemm_3x2x1_0(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_dgemm_5x4x3(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_dgemm_ldb(const MunitParameter params[],
                           void* user_data_or_fixture);

MunitResult test_hgemm_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hgemm_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_hgemm_3x2x1_0(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_hgemm_5x4x3(const MunitParameter params[],
                             void* user_data_or_fixture);

MunitResult test_qgemm_0_row(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qgemm_0_col(const MunitParameter params[],
                             void* user_data_or_fixture);
MunitResult test_qgemm_3x2x1_0(const MunitParameter params[],
                               void* user_data_or_fixture);
MunitResult test_qgemm_5x4x3(const MunitParameter params[],
                             void* user_data_or_fixture);


//  NaN canonicalization (test_nan.c)
MunitResult test_saxpy_nan(const MunitParameter params[], void* u);
MunitResult test_sdot_nan(const MunitParameter params[], void* u);
MunitResult test_snrm2_nan(const MunitParameter params[], void* u);
MunitResult test_sscal_nan(const MunitParameter params[], void* u);
MunitResult test_scopy_nan(const MunitParameter params[], void* u);
MunitResult test_sswap_nan(const MunitParameter params[], void* u);
MunitResult test_sgemv_nan(const MunitParameter params[], void* u);
MunitResult test_sgemm_nan(const MunitParameter params[], void* u);
MunitResult test_daxpy_nan(const MunitParameter params[], void* u);
MunitResult test_haxpy_nan(const MunitParameter params[], void* u);
//  C1 expansion: Inf arithmetic + asum/nrm2_B + d/h/q reductions
MunitResult test_saxpy_inf_minus_inf(const MunitParameter params[], void* u);
MunitResult test_sscal_zero_times_inf(const MunitParameter params[], void* u);
MunitResult test_sasum_nan(const MunitParameter params[], void* u);
MunitResult test_dasum_nan(const MunitParameter params[], void* u);
MunitResult test_hasum_nan(const MunitParameter params[], void* u);
MunitResult test_qasum_nan(const MunitParameter params[], void* u);
MunitResult test_snrm2_B_nan(const MunitParameter params[], void* u);
MunitResult test_ddot_nan(const MunitParameter params[], void* u);
MunitResult test_dnrm2_nan(const MunitParameter params[], void* u);
MunitResult test_dscal_nan(const MunitParameter params[], void* u);
MunitResult test_hdot_nan(const MunitParameter params[], void* u);
MunitResult test_hnrm2_nan(const MunitParameter params[], void* u);
MunitResult test_hscal_nan(const MunitParameter params[], void* u);
MunitResult test_qnrm2_nan(const MunitParameter params[], void* u);
MunitResult test_qdot_nan(const MunitParameter params[], void* u);
MunitResult test_qscal_nan(const MunitParameter params[], void* u);

//  Transpose & padding coverage (C4/C5)
MunitResult test_sgemv_trans(const MunitParameter params[], void* u);
MunitResult test_sgemv_padlda(const MunitParameter params[], void* u);
MunitResult test_sgemm_transA(const MunitParameter params[], void* u);
MunitResult test_sgemm_transB(const MunitParameter params[], void* u);
MunitResult test_hgemm_ldb(const MunitParameter params[], void* u);
MunitResult test_qgemm_ldb(const MunitParameter params[], void* u);
//  Rounding-mode coverage (C6): all five modes through a 2nd/3rd routine.
MunitResult test_sdot_rounding_modes(const MunitParameter params[], void* u);
MunitResult test_ddot_rounding_modes(const MunitParameter params[], void* u);
MunitResult test_sgemv_rounding_modes(const MunitParameter params[], void* u);
//  Negative-stride coverage (C2): swap + copy families.
MunitResult test_dswap_negstride(const MunitParameter params[], void* u);
MunitResult test_hswap_negstride(const MunitParameter params[], void* u);
MunitResult test_qswap_negstride(const MunitParameter params[], void* u);
MunitResult test_scopy_negstride(const MunitParameter params[], void* u);
MunitResult test_dcopy_negstride(const MunitParameter params[], void* u);
MunitResult test_hcopy_negstride(const MunitParameter params[], void* u);
MunitResult test_qcopy_negstride(const MunitParameter params[], void* u);
//  N==0 no-op sweeps (C3): double/half/quad L1.
MunitResult test_dn0(const MunitParameter params[], void* u);
MunitResult test_hn0(const MunitParameter params[], void* u);
MunitResult test_qn0(const MunitParameter params[], void* u);
//  gemv transpose (C4) + padded lda (C5), double/half/quad.
MunitResult test_dgemv_trans(const MunitParameter params[], void* u);
MunitResult test_hgemv_trans(const MunitParameter params[], void* u);
MunitResult test_qgemv_trans(const MunitParameter params[], void* u);
MunitResult test_dgemv_padlda(const MunitParameter params[], void* u);
MunitResult test_hgemv_padlda(const MunitParameter params[], void* u);
MunitResult test_qgemv_padlda(const MunitParameter params[], void* u);
//  gemm TT combo (C4): completes single-precision {N,T}^2.
MunitResult test_sgemm_transAB(const MunitParameter params[], void* u);

//  Complex GEMM (test_cgemm.c): i/c/z/v gemm, incl. 'C' conjugate transpose.
MunitResult test_cgemm_basic(const MunitParameter params[], void* u);
MunitResult test_cgemm_conjtrans(const MunitParameter params[], void* u);
MunitResult test_cgemm_alphabeta(const MunitParameter params[], void* u);
MunitResult test_cgemm_ldb(const MunitParameter params[], void* u);
MunitResult test_zgemm_basic(const MunitParameter params[], void* u);
MunitResult test_zgemm_conjtrans(const MunitParameter params[], void* u);
MunitResult test_igemm_basic(const MunitParameter params[], void* u);
MunitResult test_vgemm_basic(const MunitParameter params[], void* u);

//  Complex GEMV (test_cgemv.c): i/c/z/v gemv, layouts + 'C' conjugate transpose.
MunitResult test_cgemv_basic(const MunitParameter params[], void* u);
MunitResult test_cgemv_colmajor(const MunitParameter params[], void* u);
MunitResult test_cgemv_conjtrans(const MunitParameter params[], void* u);
MunitResult test_zgemv_basic(const MunitParameter params[], void* u);
MunitResult test_igemv_basic(const MunitParameter params[], void* u);
MunitResult test_vgemv_basic(const MunitParameter params[], void* u);

//  Complex Level-1 (test_complex.c)
MunitResult test_caxpy_basic(const MunitParameter params[], void* u);
MunitResult test_ccopy_basic(const MunitParameter params[], void* u);
MunitResult test_cdotc_conj(const MunitParameter params[], void* u);
//  Unconjugated complex dot product (cdotu / zdotu / idotu / vdotu).
MunitResult test_cdotu_basic(const MunitParameter params[], void* u);
MunitResult test_cdotu_vs_cdotc(const MunitParameter params[], void* u);
MunitResult test_cdotu_negstride(const MunitParameter params[], void* u);
MunitResult test_zdotu_basic(const MunitParameter params[], void* u);
MunitResult test_idotu_basic(const MunitParameter params[], void* u);
MunitResult test_vdotu_basic(const MunitParameter params[], void* u);
MunitResult test_scasum_basic(const MunitParameter params[], void* u);
MunitResult test_scnrm2_basic(const MunitParameter params[], void* u);
MunitResult test_csrot_basic(const MunitParameter params[], void* u);
//  Rotation (test_rot.c)
MunitResult test_srot_basic(const MunitParameter params[], void* u);
MunitResult test_srotg_basic(const MunitParameter params[], void* u);
MunitResult test_srotg_345(const MunitParameter params[], void* u);
MunitResult test_srotm_basic(const MunitParameter params[], void* u);


//  Modified Givens setup (test_rotmg.c)
MunitResult test_srotmg_basic(const MunitParameter params[], void* u);
MunitResult test_srotmg_flag_neg2(const MunitParameter params[], void* u);
MunitResult test_drotmg_basic(const MunitParameter params[], void* u);
//  Extended-precision dot (test_sdsdot.c)
MunitResult test_sdsdot_basic(const MunitParameter params[], void* u);
MunitResult test_hsdot_basic(const MunitParameter params[], void* u);

MunitResult test_srotmg_singular(const MunitParameter params[], void* u);
MunitResult test_srotmg_nan_d1(const MunitParameter params[], void* u);
MunitResult test_srotmg_inf_d1(const MunitParameter params[], void* u);
MunitResult test_hrotmg_inf_d1(const MunitParameter params[], void* u);
//  Quad rotation + complex Level-1 routines (qrot/qrotg/qrotm/qrotmg, c*).
MunitResult test_qrotmg_basic(const MunitParameter params[], void* u);
MunitResult test_qrotmg_inf_d1(const MunitParameter params[], void* u);
MunitResult test_qrot_basic(const MunitParameter params[], void* u);
MunitResult test_qrotg_basic(const MunitParameter params[], void* u);
MunitResult test_qrotm_basic(const MunitParameter params[], void* u);
MunitResult test_cscal_basic(const MunitParameter params[], void* u);
MunitResult test_cswap_basic(const MunitParameter params[], void* u);
MunitResult test_icamax_basic(const MunitParameter params[], void* u);
//  Complex coverage (C8): NaN canonicalization + negative stride.
MunitResult test_caxpy_nan(const MunitParameter params[], void* u);
MunitResult test_cscal_nan(const MunitParameter params[], void* u);
MunitResult test_ccopy_nan(const MunitParameter params[], void* u);
MunitResult test_caxpy_negstride(const MunitParameter params[], void* u);
MunitResult test_ccopy_negstride(const MunitParameter params[], void* u);
//  Huge-stride overflow guard (B5).
MunitResult test_nrm2_huge_stride(const MunitParameter params[], void* u);
//  Half/double/quad complex Level-1 (i / z / v families).
MunitResult test_iaxpy_basic(const MunitParameter params[], void* u);
MunitResult test_icopy_basic(const MunitParameter params[], void* u);
MunitResult test_idotc_basic(const MunitParameter params[], void* u);
MunitResult test_iscal_basic(const MunitParameter params[], void* u);
MunitResult test_iswap_basic(const MunitParameter params[], void* u);
MunitResult test_ihrot_basic(const MunitParameter params[], void* u);
MunitResult test_hiasum_basic(const MunitParameter params[], void* u);
MunitResult test_hinrm2_basic(const MunitParameter params[], void* u);
MunitResult test_iiamax_basic(const MunitParameter params[], void* u);
MunitResult test_zaxpy_basic(const MunitParameter params[], void* u);
MunitResult test_zcopy_basic(const MunitParameter params[], void* u);
MunitResult test_zdotc_basic(const MunitParameter params[], void* u);
MunitResult test_zscal_basic(const MunitParameter params[], void* u);
MunitResult test_zswap_basic(const MunitParameter params[], void* u);
MunitResult test_zdrot_basic(const MunitParameter params[], void* u);
MunitResult test_dzasum_basic(const MunitParameter params[], void* u);
MunitResult test_dznrm2_basic(const MunitParameter params[], void* u);
MunitResult test_izamax_basic(const MunitParameter params[], void* u);
MunitResult test_vaxpy_basic(const MunitParameter params[], void* u);
MunitResult test_vcopy_basic(const MunitParameter params[], void* u);
MunitResult test_vdotc_basic(const MunitParameter params[], void* u);
MunitResult test_vscal_basic(const MunitParameter params[], void* u);
MunitResult test_vswap_basic(const MunitParameter params[], void* u);
MunitResult test_vqrot_basic(const MunitParameter params[], void* u);
MunitResult test_qvasum_basic(const MunitParameter params[], void* u);
MunitResult test_qvnrm2_basic(const MunitParameter params[], void* u);
MunitResult test_ivamax_basic(const MunitParameter params[], void* u);
//  Complex i/z/v NaN canonicalization + negative stride; complex B5 guard.
MunitResult test_iaxpy_nan(const MunitParameter params[], void* u);
MunitResult test_iaxpy_negstride(const MunitParameter params[], void* u);
MunitResult test_icopy_negstride(const MunitParameter params[], void* u);
MunitResult test_zaxpy_nan(const MunitParameter params[], void* u);
MunitResult test_zaxpy_negstride(const MunitParameter params[], void* u);
MunitResult test_zcopy_negstride(const MunitParameter params[], void* u);
MunitResult test_vaxpy_nan(const MunitParameter params[], void* u);
MunitResult test_vaxpy_negstride(const MunitParameter params[], void* u);
MunitResult test_vcopy_negstride(const MunitParameter params[], void* u);
MunitResult test_complex_huge_stride(const MunitParameter params[], void* u);


//  Fast (Anderson) nrm2 (test_nrm2_B.c)
MunitResult test_snrm2_B_med(const MunitParameter params[], void* u);
MunitResult test_snrm2_B_big(const MunitParameter params[], void* u);
MunitResult test_snrm2_B_small(const MunitParameter params[], void* u);
MunitResult test_dnrm2_B_med(const MunitParameter params[], void* u);
MunitResult test_dnrm2_B_big(const MunitParameter params[], void* u);
MunitResult test_hnrm2_B_med(const MunitParameter params[], void* u);
MunitResult test_hnrm2_B_big(const MunitParameter params[], void* u);
MunitResult test_qnrm2_B_med(const MunitParameter params[], void* u);

#endif // TEST_H
