#ifndef SOFTBLAS_H
#define SOFTBLAS_H

#include "softfloat.h"
#include <stdlib.h>

//  TYPES

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
#define SB_REAL32_ONE 0x3f800000
#define SB_REAL64_ONE 0x3ff0000000000000
#define SB_REAL128U_ONE 0x3fff000000000000
#define SB_REAL128L_ONE 0x0000000000000000
#define SB_REAL16_ZERO 0x0000
#define SB_REAL32_ZERO 0x00000000
#define SB_REAL64_ZERO 0x0000000000000000
#define SB_REAL128U_ZERO 0x0000000000000000
#define SB_REAL128L_ZERO 0x0000000000000000
#define SB_REAL16_NEGONE 0xbc00
#define SB_REAL32_NEGONE 0xbf800000
#define SB_REAL64_NEGONE 0xbff0000000000000
#define SB_REAL128U_NEGONE 0xbfff000000000000
#define SB_REAL128L_NEGONE 0x0000000000000000
#define SB_REAL16_NEGTWO 0xc000
#define SB_REAL32_NEGTWO 0xc0000000
#define SB_REAL64_NEGTWO 0xc000000000000000
#define SB_REAL128U_NEGTWO 0xc000000000000000
#define SB_REAL128L_NEGTWO 0x0000000000000000

#define SB_COMPLEX16_ZERO {SB_REAL16_ZERO, SB_REAL16_ZERO}
#define SB_COMPLEX32_ZERO {SB_REAL32_ZERO, SB_REAL32_ZERO}
#define SB_COMPLEX64_ZERO {SB_REAL64_ZERO, SB_REAL64_ZERO}
#define SB_COMPLEX128_ZERO {SB_REAL128L_ZERO, SB_REAL128U_ZERO}

#define SB_COMPLEX16_ONE {SB_REAL16_ONE, SB_REAL16_ZERO}
#define SB_COMPLEX32_ONE {SB_REAL32_ONE, SB_REAL32_ZERO}
#define SB_COMPLEX64_ONE {SB_REAL64_ONE, SB_REAL64_ZERO}
#define SB_COMPLEX128_ONE {SB_REAL128L_ONE, SB_REAL128U_ONE}

#define SB_COMPLEX16_NEGONE {SB_REAL16_NEGONE, SB_REAL16_ZERO}
#define SB_COMPLEX32_NEGONE {SB_REAL32_NEGONE, SB_REAL32_ZERO}
#define SB_COMPLEX64_NEGONE {SB_REAL64_NEGONE, SB_REAL64_ZERO}
#define SB_COMPLEX128_NEGONE {SB_REAL128L_NEGONE, SB_REAL128U_NEGONE}

#define SB_COMPLEX16_NEGTWO {SB_REAL16_NEGTWO, SB_REAL16_ZERO}
#define SB_COMPLEX32_NEGTWO {SB_REAL32_NEGTWO, SB_REAL32_ZERO}
#define SB_COMPLEX64_NEGTWO {SB_REAL64_NEGTWO, SB_REAL64_ZERO}
#define SB_COMPLEX128_NEGTWO {SB_REAL128L_NEGTWO, SB_REAL128U_NEGTWO}

#define SB_COMPLEX16_I {SB_REAL16_ZERO, SB_REAL16_ONE}
#define SB_COMPLEX32_I {SB_REAL32_ZERO, SB_REAL32_ONE}
#define SB_COMPLEX64_I {SB_REAL64_ZERO, SB_REAL64_ONE}
// #define SB_COMPLEX128_I {SB_REAL128L_ZERO, SB_REAL128U_ONE}

//  SOFTFLOAT FUNCTIONS

#define f16_ge(x,y) ( !(f16_lt( x,y ) ))
#define f16_gt(x,y) ( !(f16_le( x,y ) ))
#define f16_ne(x,y) ( !(f16_eq( x,y ) ))
#define f32_ge(x,y) ( !(f32_lt( x,y ) ))
#define f32_gt(x,y) ( !(f32_le( x,y ) ))
#define f32_ne(x,y) ( !(f32_eq( x,y ) ))
#define f64_ge(x,y) ( !(f64_lt( x,y ) ))
#define f64_gt(x,y) ( !(f64_le( x,y ) ))
#define f64_ne(x,y) ( !(f64_eq( x,y ) ))
#define f128M_ge(x,y) ( !(f128M_lt( x,y ) ))
#define f128M_gt(x,y) ( !(f128M_le( x,y ) ))
#define f128M_ne(x,y) ( !(f128M_eq( x,y ) ))

#define ABS(x) ( (x) >= 0 ? (x) : -(x) )
#define f16_abs(x) (float16_t){ ( (uint16_t)(x.v) & 0x7fff ) }
#define f32_abs(x) (float32_t){ ( (uint32_t)(x.v) & 0x7fffffff ) }
#define f64_abs(x) (float64_t){ ( (uint64_t)(x.v) & 0x7fffffffffffffff ) }
#define f128_abs(x) (float128_t){ ( (uint64_t)(x.v[0]) & 0x7fffffffffffffff ), ( (uint64_t)(x.v[1]) & 0x7fffffffffffffff ) }
#define f16M_abs(x) (float16_t*){ ( (uint16_t)(x->v) & 0x7fff ) }
#define f32M_abs(x) (float32_t*){ ( (uint32_t)(x->v) & 0x7fffffff ) }
#define f64M_abs(x) (float64_t*){ ( (uint64_t)(x->v) & 0x7fffffffffffffff ) }
#define f128M_abs(x) (float128_t*){ ( (uint64_t)(x->v[0]) & 0x7fffffffffffffff ), ( (uint64_t)(x->v[1]) & 0x7fffffffffffffff ) }

#define MAX(x, y) ( (x) > (y) ? (x) : (y) )
#define f16_max(x, y) ( f16_gt( (x) , (y) ) ? (x) : (y) )
#define f32_max(x, y) ( f32_gt( (x) , (y) ) ? (x) : (y) )
#define f64_max(x, y) ( f64_gt( (x) , (y) ) ? (x) : (y) )
#define f128M_max(x, y) ( f128M_gt( (x) , (y) ) ? (x) : (y) )

#define MIN(x, y) ( (x) > (y) ? (y) : (x) )
#define f16_min(x, y) ( f16_gt( (x) , (y) ) ? (y) : (x) )
#define f32_min(x, y) ( f32_gt( (x) , (y) ) ? (y) : (x) )
#define f64_min(x, y) ( f64_gt( (x) , (y) ) ? (y) : (x) )
#define f128_min(x, y) ( f128_gt( (x) , (y) ) ? (y) : (x) )

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
#define FNEG(x) ( x.v ^ (1 << (sizeof(x.v) * 8 - 1)) )

//  SOFTBLAS PROTOTYPES

//  Level 1
//    Single-precision
float32_t sasum(uint64_t N, const float32_t *SX, uint64_t incX);
void saxpy(uint64_t N, float32_t SA, float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
void scopy(uint64_t N, const float32_t *SX, int64_t incX, float32_t *SY, int64_t incY);
float32_t sdot(const uint64_t N, const float32_t *X, const int64_t incX, const float32_t *Y, const int64_t incY);
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
void daxpy(uint64_t N, float64_t DA, float64_t *DX, int64_t incX, float64_t *DY, int64_t incY);
void dcopy(uint64_t N, const float64_t *DX, int64_t incX, float64_t *DY, int64_t incY);
float64_t ddot(const uint64_t N, const float64_t *X, const int64_t incX, const float64_t *Y, const int64_t incY);
float64_t dnrm2(uint64_t N, const float64_t *X, uint64_t incX);
void drot(const uint64_t N, float64_t *X, const uint64_t  incX, float64_t *Y, const uint64_t incY, const float64_t c, const float64_t s);
void drotg(float64_t *a, float64_t *b, float64_t *c, float64_t *s);
void drotm(const uint64_t N, float64_t *X, const uint64_t incX, float64_t *Y, const uint64_t incY, const float64_t *P);
void drotmg(float64_t *D1, float64_t *D2, float64_t *X1, const float64_t y1, float64_t *P);
void dscal(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX);
void dswap(uint64_t N, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
uint64_t idamax(uint64_t N, const float64_t *DX, uint64_t incX);

//    Half-precision
float16_t hasum(uint64_t N, const float16_t *HX, uint64_t incX);
void haxpy(uint64_t N, float16_t HA, float16_t *HX, int64_t incX, float16_t *HY, int64_t incY);
void hcopy(uint64_t N, const float16_t *HX, int64_t incX, float16_t *HY, int64_t incY);
float16_t hdot(const uint64_t N, const float16_t *X, const int64_t incX, const float16_t *Y, const int64_t incY);
float16_t hnrm2(uint64_t N, const float16_t *X, uint64_t incX);
void hrot(const uint64_t N, float16_t *X, const uint64_t  incX, float16_t *Y, const uint64_t incY, const float16_t c, const float16_t s);
void hrotg(float16_t *a, float16_t *b, float16_t *c, float16_t *s);
void hrotm(const uint64_t N, float16_t *X, const uint64_t incX, float16_t *Y, const uint64_t incY, const float16_t *P);
void hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P);
void hscal(uint64_t N, float16_t HA, float16_t *HX, uint64_t incX);
void hswap(uint64_t N, float16_t *HX, uint64_t incX, float16_t *HY, uint64_t incY);
uint64_t ihamax(uint64_t N, const float16_t *HX, uint64_t incX);

//    Quad-precision
float128_t qasum(uint64_t N, const float128_t *QX, uint64_t incX);
void qaxpy(uint64_t N, float128_t QA, float128_t *QX, int64_t incX, float128_t *QY, int64_t incY);
void qcopy(uint64_t N, const float128_t *QX, int64_t incX, float128_t *QY, int64_t incY);
float128_t qdot(const uint64_t N, const float128_t *X, const int64_t incX, const float128_t *Y, const int64_t incY);
float128_t qnrm2(uint64_t N, const float128_t *X, uint64_t incX);
void qrot(const uint64_t N, float16_t *X, const uint64_t  incX, float16_t *Y, const uint64_t incY, const float16_t c, const float16_t s);
void qrotg(float128_t *a, float128_t *b, float128_t *c, float128_t *s);
void qrotm(const uint64_t N, float128_t *X, const uint64_t incX, float128_t *Y, const uint64_t incY, const float128_t *P);
void qrotmg(float128_t *D1, float128_t *D2, float128_t *X1, const float128_t y1, float128_t *P);
void qscal(uint64_t N, float128_t QA, float128_t *QX, uint64_t incX);
void qswap(uint64_t N, float128_t *QX, uint64_t incX, float128_t *QY, uint64_t incY);
uint64_t iqamax(uint64_t N, const float128_t *QX, uint64_t incX);

//    Complex single-precision
float32_t scasum(uint64_t N, const complex32_t *CX, int64_t incX);
void caxpy(uint64_t N, complex32_t CA, complex32_t *CX, int64_t incX, complex32_t *HY, int64_t incY);
void ccopy(uint64_t N, const complex32_t *CX, int64_t incX, complex32_t *CY, int64_t incY);
complex32_t cdotc(uint64_t N, const complex32_t *CX, int64_t incX, const complex32_t *CY, int64_t incY);
float32_t scnrm2(uint64_t N, const complex32_t *CX, uint64_t incX);
void csrot(const uint64_t N, complex32_t *CX, const uint64_t incX, complex32_t *CY, const uint64_t incY, const complex32_t c, const complex32_t s);

void cscal(uint64_t N, complex32_t CA, complex32_t *CX, uint64_t incX);
void cswap(uint64_t N, complex32_t *CX, uint64_t incX, complex32_t *CY, uint64_t incY);
uint64_t icamax(uint64_t N, const complex32_t *CX, uint64_t incX);

//    Complex double-precision

//    Complex half-precision

//    Complex quad-precision

// Level 2

void sgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *X, const int64_t incX, const float32_t beta, float32_t *Y, const uint64_t incY);
void dgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *X, const int64_t incX, const float64_t beta, float64_t *Y, const uint64_t incY);
void hgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *X, const int64_t incX, const float16_t beta, float16_t *Y, const uint64_t incY);
void qgemv(const char Layout, const char Trans, const uint64_t M, const uint64_t N, const float128_t alpha, const float128_t *A, const uint64_t lda, const float128_t *X, const int64_t incX, const float128_t beta, float128_t *Y, const uint64_t incY);

// Level 3

void sgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *B, const uint64_t ldb, const float32_t beta, float32_t *C, const uint64_t ldc);
void dgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *B, const uint64_t ldb, const float64_t beta, float64_t *C, const uint64_t ldc);
void hgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *B, const uint64_t ldb, const float16_t beta, float16_t *C, const uint64_t ldc);
void qgemm(const char transA, const char transB, const uint64_t M, const uint64_t N, const uint64_t P, const float128_t alpha, const float128_t *A, const uint64_t lda, const float128_t *B, const uint64_t ldb, const float128_t beta, float128_t *C, const uint64_t ldc);

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
        *( (uint64_t*)a)    = 0x0000000000000000;
        *(((uint64_t*)a)+1) = QUADNAN;
    }
}

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
    
    for (uint64_t i = 0; i < size; i++) {
        result[i].v[0] = pairs[i].lo;
        result[i].v[1] = pairs[i].hi;
    }
    
    return result;
}

#endif // SOFTBLAS_H
