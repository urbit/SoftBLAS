#include "test.h"

//  NaN determinism: SoftFloat propagates a (quieted) NaN payload, so every
//  value-producing routine must canonicalize its output to the per-precision
//  canonical NaN (HALFNAN/SINGNAN/DOUBNAN/QUADNAN). Inputs below carry a
//  non-canonical payload (…0001); outputs must come back canonical.

#define NCNAN_S 0x7fc00001u            // non-canonical single NaN
#define NCNAN_D 0x7ff0000000000001ull  // non-canonical double NaN
#define NCNAN_H 0x7e01                 // non-canonical half NaN

MunitResult test_saxpy_nan(const MunitParameter params[], void* u) {
    const float32_t SA = { SB_REAL32_ONE };
    float32_t SX[1] = {{ NCNAN_S }}, SY[1] = {{ 0 }};
    saxpy(1, SA, SX, 1, SY, 1, 'n');
    assert_ulong(SY[0].v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_sdot_nan(const MunitParameter params[], void* u) {
    float32_t SX[1] = {{ NCNAN_S }}, SY[1] = {{ 0x3f800000u }};
    float32_t r = sdot(1, SX, 1, SY, 1, 'n');
    assert_ulong(r.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_snrm2_nan(const MunitParameter params[], void* u) {
    float32_t SX[1] = {{ NCNAN_S }};
    float32_t r = snrm2(1, SX, 1, 'n');
    assert_ulong(r.v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_sscal_nan(const MunitParameter params[], void* u) {
    const float32_t SA = { SB_REAL32_ONE };
    float32_t SX[1] = {{ NCNAN_S }};
    sscal(1, SA, SX, 1, 'n');
    assert_ulong(SX[0].v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_scopy_nan(const MunitParameter params[], void* u) {
    float32_t SX[1] = {{ NCNAN_S }}, SY[1] = {{ 0 }};
    scopy(1, SX, 1, SY, 1, 'n');
    assert_ulong(SY[0].v, ==, (uint32_t)SINGNAN);   // a copied NaN is canonicalized too
    return MUNIT_OK;
}
MunitResult test_sswap_nan(const MunitParameter params[], void* u) {
    float32_t SX[1] = {{ NCNAN_S }}, SY[1] = {{ 0x3f800000u }};
    sswap(1, SX, 1, SY, 1, 'n');
    assert_ulong(SY[0].v, ==, (uint32_t)SINGNAN);   // the NaN moved into SY
    assert_ulong(SX[0].v, ==, 0x3f800000u);         // the finite value moved into SX
    return MUNIT_OK;
}
MunitResult test_sgemv_nan(const MunitParameter params[], void* u) {
    const float32_t alpha = { SB_REAL32_ONE }, beta = { SB_REAL32_ZERO };
    float32_t A[1] = {{ NCNAN_S }}, X[1] = {{ 0x3f800000u }}, Y[1] = {{ 0 }};
    sgemv('R', 'N', 1, 1, alpha, A, 1, X, 1, beta, Y, 1, 'n');
    assert_ulong(Y[0].v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_sgemm_nan(const MunitParameter params[], void* u) {
    const float32_t alpha = { SB_REAL32_ONE }, beta = { SB_REAL32_ZERO };
    float32_t A[1] = {{ NCNAN_S }}, B[1] = {{ 0x3f800000u }}, C[1] = {{ 0 }};
    sgemm('N', 'N', 1, 1, 1, alpha, A, 1, B, 1, beta, C, 1, 'n');
    assert_ulong(C[0].v, ==, (uint32_t)SINGNAN);
    return MUNIT_OK;
}
MunitResult test_daxpy_nan(const MunitParameter params[], void* u) {
    const float64_t DA = { SB_REAL64_ONE };
    float64_t DX[1] = {{ NCNAN_D }}, DY[1] = {{ 0 }};
    daxpy(1, DA, DX, 1, DY, 1, 'n');
    assert_ullong(DY[0].v, ==, (uint64_t)DOUBNAN);
    return MUNIT_OK;
}
MunitResult test_haxpy_nan(const MunitParameter params[], void* u) {
    const float16_t HA = { SB_REAL16_ONE };
    float16_t HX[1] = {{ NCNAN_H }}, HY[1] = {{ 0 }};
    haxpy(1, HA, HX, 1, HY, 1, 'n');
    assert_ushort(HY[0].v, ==, (uint16_t)HALFNAN);
    return MUNIT_OK;
}
