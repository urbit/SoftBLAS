#ifndef SOFTBLAS_H
#define SOFTBLAS_H

#include "softfloat.h"

//  TYPES

enum SB_LAYOUT {
   SblasRowMajor=101,     // row-major arrays
   SblasColMajor=102};    // column-major arrays
enum SB_TRANSPOSE {
   SblasNoTrans=78,       // trans='N'
   SblasTrans=84,         // trans='T'
   SblasConjTrans=67};    // trans='C'
enum SB_UPLO {
   SblasUpper=85,         // uplo ='U'
   SblasLower=76};        // uplo ='L'
enum SB_DIAG {
   SblasNonUnit=78,       // diag ='N'
   SblasUnit=85};         // diag ='U'
enum SB_SIDE {
   SblasLeft=76,          // side ='L'
   SblasRight=82};        // side ='R'

typedef struct {
   float32_t real;
   float32_t imag;
} complex32_t;
typedef struct {
   float64_t real;
   float64_t imag;
} complex64_t;
typedef struct {
   float16_t real;
   float16_t imag;
} complex16_t;
typedef struct {
   float128_t real;
   float128_t imag;
} complex128_t;

//  CONSTANTS

#define SB_REAL16_ONE 0x3c00
#define SB_REAL32_ONE 1.0f
#define SB_REAL64_ONE 1.0
// #define SB_REAL128U_ONE 1.0
// #define SB_REAL128L_ONE 0.0
#define SB_REAL16_ZERO 0x0000
#define SB_REAL32_ZERO 0.0f
#define SB_REAL64_ZERO 0.0
// #define SB_REAL128U_ZERO 0.0
// #define SB_REAL128L_ZERO 0.0
#define SB_REAL16_NEGONE 0xbc00
#define SB_REAL32_NEGONE -1.0f
#define SB_REAL64_NEGONE -1.0
// #define SB_REAL128U_NEGONE -1.0
// #define SB_REAL128L_NEGONE 0.0
#define SB_REAL16_NEGTWO 0xc000
#define SB_REAL32_NEGTWO -2.0f
#define SB_REAL64_NEGTWO -2.0

#define SB_COMPLEX16_ZERO {SB_REAL16_ZERO, SB_REAL16_ZERO}
#define SB_COMPLEX32_ZERO {SB_REAL32_ZERO, SB_REAL32_ZERO}
#define SB_COMPLEX64_ZERO {SB_REAL64_ZERO, SB_REAL64_ZERO}
// #define SB_COMPLEX128_ZERO {SB_REAL128U_ZERO, SB_REAL128L_ZERO}

#define SB_COMPLEX16_ONE {SB_REAL16_ONE, SB_REAL16_ZERO}
#define SB_COMPLEX32_ONE {SB_REAL32_ONE, SB_REAL32_ZERO}
#define SB_COMPLEX64_ONE {SB_REAL64_ONE, SB_REAL64_ZERO}
// #define SB_COMPLEX128_ONE {SB_REAL128U_ONE, SB_REAL128L_ONE}

#define SB_COMPLEX16_NEGONE {SB_REAL16_NEGONE, SB_REAL16_ZERO}
#define SB_COMPLEX32_NEGONE {SB_REAL32_NEGONE, SB_REAL32_ZERO}
#define SB_COMPLEX64_NEGONE {SB_REAL64_NEGONE, SB_REAL64_ZERO}
// #define SB_COMPLEX128_NEGONE {SB_REAL128U_NEGONE, SB_REAL128L_NEGONE}

#define SB_COMPLEX16_NEGTWO {SB_REAL16_NEGTWO, SB_REAL16_ZERO}
#define SB_COMPLEX32_NEGTWO {SB_REAL32_NEGTWO, SB_REAL32_ZERO}
#define SB_COMPLEX64_NEGTWO {SB_REAL64_NEGTWO, SB_REAL64_ZERO}
// #define SB_COMPLEX128_NEGTWO {SB_REAL128U_NEGTWO, SB_REAL128L_ZERO}

#define SB_COMPLEX16_I {SB_REAL16_ZERO, SB_REAL16_ONE}
#define SB_COMPLEX32_I {SB_REAL32_ZERO, SB_REAL32_ONE}
#define SB_COMPLEX64_I {SB_REAL64_ZERO, SB_REAL64_ONE}
// #define SB_COMPLEX128_I {SB_REAL128L_ZERO, SB_REAL128U_ONE}

//  SOFTFLOAT FUNCTIONS

#define f16_ge !f16_lt
#define f16_gt !f16_le
#define f16_ne !f16_eq
#define f32_ge !f32_lt
#define f32_gt !f32_le
#define f32_ne !f32_eq
#define f64_ge !f64_lt
#define f64_gt !f64_le
#define f64_ne !f64_eq
#define f128_ge !f128_lt
#define f128_gt !f128_le
#define f128_ne !f128_eq

#define ABS(x) ( (x) >= 0 ? (x) : -(x) )
#define f16_abs(x) ( (x) & (float16_t)0x7fff )
#define f32_abs(x) ( (x) & (float32_t)0x7fffffff )
#define f64_abs(x) ( (x) & (float64_t)0x7fffffffffffffff )

#define MAX(x, y) ( (x) > (y) ? (x) : (y) )
#define f16_max(x, y) ( f16_gt( (x) , (y) ) ? (x) : (y) )
#define f32_max(x, y) ( f32_gt( (x) , (y) ) ? (x) : (y) )
#define f64_max(x, y) ( f64_gt( (x) , (y) ) ? (x) : (y) )

#define MIN(x, y) ( (x) > (y) ? (y) : (x) )
#define f16_min(x, y) ( f16_gt( (x) , (y) ) ? (y) : (x) )
#define f32_min(x, y) ( f32_gt( (x) , (y) ) ? (y) : (x) )
#define f64_min(x, y) ( f64_gt( (x) , (y) ) ? (y) : (x) )

#define c16_eq(a, b)  ( (a.real == b.real) && (a.imag == b.imag) )
#define c32_eq(a, b)  ( (a.real == b.real) && (a.imag == b.imag) )
#define c64_eq(a, b)  ( (a.real == b.real) && (a.imag == b.imag) )
#define c128_eq(a, b) ( (a.real == b.real) && (a.imag == b.imag) )
#define c16_ne  !c16_eq
#define c32_ne  !c32_eq
#define c64_ne  !c64_eq
#define c128_ne !c128_eq

#define c16_mag(a)  f16_sqrt( f16_add( f16_mul(a.real, a.real), f16_mul(a.imag, a.imag) ) )
#define c32_mag(a)  f32_sqrt( f32_add( f32_mul(a.real, a.real), f32_mul(a.imag, a.imag) ) )
#define c64_mag(a)  f64_sqrt( f64_add( f64_mul(a.real, a.real), f64_mul(a.imag, a.imag) ) )
#define c128_mag(a) f128_sqrt( f128_add( f128_mul(a.real, a.real), f128_mul(a.imag, a.imag) ) )

#define REAL(a)  (a.real)
#define IMAG(a)  (a.imag)

#define c16_add(a, b) { f16_add(a.real, b.real), f16_add(a.imag, b.imag) }
#define c32_add(a, b) { f32_add(a.real, b.real), f32_add(a.imag, b.imag) }
#define c64_add(a, b) { f64_add(a.real, b.real), f64_add(a.imag, b.imag) }
#define c128_add(a, b) { f128_add(a.real, b.real), f128_add(a.imag, b.imag) }

#define c16_sub(a, b) { f16_sub(a.real, b.real), f16_sub(a.imag, b.imag) }
#define c32_sub(a, b) { f32_sub(a.real, b.real), f32_sub(a.imag, b.imag) }
#define c64_sub(a, b) { f64_sub(a.real, b.real), f64_sub(a.imag, b.imag) }
#define c128_sub(a, b) { f128_sub(a.real, b.real), f128_sub(a.imag, b.imag) }

#define c16_mul(a, b) { f16_sub( f16_mul(a.real, b.real), f16_mul(a.imag, b.imag) ), f16_add( f16_mul(a.real, b.imag), f16_mul(a.imag, b.real) ) }
#define c32_mul(a, b) { f32_sub( f32_mul(a.real, b.real), f32_mul(a.imag, b.imag) ), f32_add( f32_mul(a.real, b.imag), f32_mul(a.imag, b.real) ) }
#define c64_mul(a, b) { f64_sub( f64_mul(a.real, b.real), f64_mul(a.imag, b.imag) ), f64_add( f64_mul(a.real, b.imag), f64_mul(a.imag, b.real) ) }
#define c128_mul(a, b) { f128_sub( f128_mul(a.real, b.real), f128_mul(a.imag, b.imag) ), f128_add( f128_mul(a.real, b.imag), f128_mul(a.imag, b.real) ) }

#define c16_div(a, b) { f16_div( f16_add( f16_mul(a.real, b.real), f16_mul(a.imag, b.imag) ), f16_add( f16_mul(b.real, b.real), f16_mul(b.imag, b.imag) ) ), f16_div( f16_sub( f16_mul(a.imag, b.real), f16_mul(a.real, b.imag) ), f16_add( f16_mul(b.real, b.real), f16_mul(b.imag, b.imag) ) ) }
#define c32_div(a, b) { f32_div( f32_add( f32_mul(a.real, b.real), f32_mul(a.imag, b.imag) ), f32_add( f32_mul(b.real, b.real), f32_mul(b.imag, b.imag) ) ), f32_div( f32_sub( f32_mul(a.imag, b.real), f32_mul(a.real, b.imag) ), f32_add( f32_mul(b.real, b.real), f32_mul(b.imag, b.imag) ) ) }
#define c64_div(a, b) { f64_div( f64_add( f64_mul(a.real, b.real), f64_mul(a.imag, b.imag) ), f64_add( f64_mul(b.real, b.real), f64_mul(b.imag, b.imag) ) ), f64_div( f64_sub( f64_mul(a.imag, b.real), f64_mul(a.real, b.imag) ), f64_add( f64_mul(b.real, b.real), f64_mul(b.imag, b.imag) ) ) }
#define c128_div(a, b) { f128_div( f128_add( f128_mul(a.real, b.real), f128_mul(a.imag, b.imag) ), f128_add( f128_mul(b.real, b.real), f128_mul(b.imag, b.imag) ) ), f128_div( f128_sub( f128_mul(a.imag, b.real), f128_mul(a.real, b.imag) ), f128_add( f128_mul(b.real, b.real), f128_mul(b.imag, b.imag) ) ) }

#define CONJ(a) { a.real, FNEG(a.imag) }
#define FNEG(x) ( x ^ (1 << (sizeof(x) * 8 - 1)) )

//  SOFTBLAS PROTOTYPES

//  Level 1
//    Single-precision
float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX);
void saxpy(uint64_t N, float32_t SA, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
void scopy(uint64_t N, const float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
float32_t sdot(const uint64_t N, const float32_t *X, const int64_t incX, const float32_t *Y, const int64_t incY);

float32_t sdsdot(const uint64_t N, const float32_t alpha, const float32_t *X, const uint64_t incX, const float32_t *Y, const uint64_t incY);
float32_t snrm2(uint64_t N, const float32_t *X, uint64_t incX);
void srot(const uint64_t N, float32_t *X, const uint64_t  incX, float32_t *Y, const uint64_t incY, const float32_t c, const float32_t s);
void srotg(float32_t *a, float32_t *b, float32_t *c, float32_t *s);
void srotm(const uint64_t N, float32_t *X, const uint64_t incX, float32_t *Y, const uint64_t incY, const float32_t *P);
void srotmg(float32_t *D1, float32_t *D2, float32_t *X1, const float32_t y1, float32_t *P);
void sscal(uint64_t N, float32_t SA, float32_t *SX, uint64_t incX);
void sswap(uint64_t N, float32_t *SX, uint64_t incX, float32_t *SY, uint64_t incY);
uint64_t isamax(uint64_t N, const float32_t *SX, uint64_t incX);

//    Double-precision
float64_t dasum(uint64_t N, const float64_t *DX, uint64_t incX);
void daxpy(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
void dcopy(uint64_t N, const float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
float64_t ddot(const uint64_t N, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY);
float64_t dnrm2(uint64_t N, const float64_t *X, uint64_t incX);
void drot(const uint64_t N, float64_t *X, const uint64_t  incX, float64_t *Y, const uint64_t incY, const float64_t c, const float64_t s);
void drotg(float64_t *a, float64_t *b, float64_t *c, float64_t *s);
void drotm(const uint64_t N, float64_t *X, const uint64_t incX, float64_t *Y, const uint64_t incY, const float64_t *P);
void drotmg(float64_t *D1, float64_t *D2, float64_t *X1, const float64_t y1, float64_t *P);
void dscal(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX);
float64_t dsdot(const uint64_t N, const float64_t alpha, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY);
void dswap(uint64_t N, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
uint64_t idamax(uint64_t N, const float64_t *DX, uint64_t incX);

//    Half-precision
float16_t hasum(uint64_t N, const float16_t *HX, uint64_t incX);
void haxpy(uint64_t N, float16_t HA, float16_t *HX, uint64_t incX, float16_t *HY, uint64_t incY);
void hcopy(uint64_t N, const float16_t *HX, uint64_t incX, float16_t *HY, uint64_t incY);
float16_t hdot(const uint64_t N, const float16_t *X, const uint64_t incX, const float16_t *Y, const uint64_t incY);
float16_t hnrm2(uint64_t N, const float16_t *X, uint64_t incX);
void hrot(const uint64_t N, float16_t *X, const uint64_t  incX, float16_t *Y, const uint64_t incY, const float16_t c, const float16_t s);
void hrotg(float16_t *a, float16_t *b, float16_t *c, float16_t *s);
void hrotm(const uint64_t N, float16_t *X, const uint64_t incX, float16_t *Y, const uint64_t incY, const float16_t *P);
void hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P);
void hscal(uint64_t N, float16_t HA, float16_t *HX, uint64_t incX);
float16_t hsdot(const uint64_t N, const float16_t alpha, const float16_t *X, const uint64_t incX, const float16_t *Y, const uint64_t incY);
void hswap(uint64_t N, float16_t *HX, uint64_t incX, float16_t *HY, uint64_t incY);
uint64_t ihamax(uint64_t N, const float16_t *HX, uint64_t incX);

//    Quad-precision

//    Complex single-precision
float32_t scasum(uint64_t N, const complex32_t *SX, int64_t incX);
void caxpy(uint64_t N, complex32_t SA, complex32_t *HX, int64_t incX, complex32_t *HY, int64_t incY);
void ccopy(uint64_t N, const complex32_t *CX, int64_t incX, complex32_t *CY, int64_t incY);
complex32_t cdotc(uint64_t N, const complex32_t *CX, int64_t incX, const complex32_t *CY, int64_t incY);
float32_t scnrm2(uint64_t N, const complex32_t *CX, uint64_t incX);
void csrot(const uint64_t N, complex32_t *CX, const uint64_t incX, complex32_t *CY, const uint64_t incY, const complex32_t c, const complex32_t s);

void cscal(uint64_t N, complex32_t SA, complex32_t *SX, uint64_t incX);
void cswap(uint64_t N, complex32_t *SX, uint64_t incX, complex32_t *SY, uint64_t incY);
uint64_t icamax(uint64_t N, const complex32_t *SX, uint64_t incX);

//    Complex double-precision

//    Complex half-precision

//    Complex quad-precision

// Level 2

void sgemv(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *X, const uint64_t incX, const float32_t beta, float32_t *Y, const uint64_t incY);
void dgemv(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *X, const uint64_t incX, const float64_t beta, float64_t *Y, const uint64_t incY);
void hgemv(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *X, const uint64_t incX, const float16_t beta, float16_t *Y, const uint64_t incY);

// Level 3

void sgemm(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE transA, const enum SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *B, const uint64_t ldb, const float32_t beta, float32_t *C, const uint64_t ldc);
void dgemm(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE transA, const enum SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *B, const uint64_t ldb, const float64_t beta, float64_t *C, const uint64_t ldc);
void hgemm(const enum SB_LAYOUT Layout, const enum SB_TRANSPOSE transA, const enum SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *B, const uint64_t ldb, const float16_t beta, float16_t *C, const uint64_t ldc);

// NAN unification

#define HALFNAN 0x7e00
#define SINGNAN 0x7fc00000
#define DOUBNAN 0x7ff8000000000000
#define QUADNAN 0x7fff800000000000

static inline bool nan_test_h(float16_t a) {
    return f16_ne(a, a);
}

static inline float16_t nan_unify_h(float16_t a) {
    if ( nan_test_h(a) ) {
        *(uint16_t*)(&a) = HALFNAN;
    }
    return a;
}

static inline bool nan_test_s(float32_t a) {
    return f32_ne(a, a);
}

static inline float32_t nan_unify_s(float32_t a) {
    if ( nan_test_s(a) ) {
        *(uint32_t*)(&a) = SINGNAN;
    }
    return a;
}

static inline bool nan_test_d(float64_t a) {
    return f64_ne(a, a);
}

static inline float64_t nan_unify_d(float64_t a) {
    if ( nan_test_d(a) ) {
        *(uint64_t*)(&a) = DOUBNAN;
    }
    return a;
}

static inline bool nan_test_q(float128_t* a) {
    return !f128M_eq(a, a);
}

static inline void nan_unify_q(float128_t* a) {
    if ( nan_test_q(a) ) {
        *( (uint64_t*)a)    = 0;
        *(((uint64_t*)a)+1) = QUADNAN;
    }
}

#endif // SOFTBLAS_H
