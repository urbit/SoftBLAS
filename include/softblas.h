#include "softfloat.h"

//  TYPES

enum SB_LAYOUT {
   SblasRowMajor=101,     // row-major arrays
   SblasColMajor=102};    // column-major arrays
enum SBLAS_TRANSPOSE {
   SblasNoTrans=78,       // trans='N'
   SblasTrans=84,         // trans='T'
   SblasConjTrans=67};    // trans='C'
enum SBLAS_UPLO {
   SblasUpper=85,         // uplo ='U'
   SblasLower=76};        // uplo ='L'
enum SBLAS_DIAG {
   SblasNonUnit=78,       // diag ='N'
   SblasUnit=85};         // diag ='U'
enum SBLAS_SIDE {
   SblasLeft=76,          // side ='L'
   SblasRight=82};        // side ='R'

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
#define f16_abs(x) ( (x) & 0x7fff )
#define f32_abs(x) ( (x) & 0x7fffffff )
#define f64_abs(x) ( (x) & 0x7fffffffffffffff )

#define MAX(x, y) ( (x) > (y) ? (x) : (y) )
#define f16_max(x, y) ( f16_gt( (x) , (y) ) ? (x) : (y) )
#define f32_max(x, y) ( f32_gt( (x) , (y) ) ? (x) : (y) )
#define f64_max(x, y) ( f64_gt( (x) , (y) ) ? (x) : (y) )

#define MIN(x, y) ( (x) > (y) ? (y) : (x) )
#define f16_min(x, y) ( f16_gt( (x) , (y) ) ? (y) : (x) )
#define f32_min(x, y) ( f32_gt( (x) , (y) ) ? (y) : (x) )
#define f64_min(x, y) ( f64_gt( (x) , (y) ) ? (y) : (x) )

//  SOFTBLAS PROTOTYPES

//  Level 1
float32_t sasum(uint32_t N, const float32_t *SX, uint32_t incX);
void saxpy(uint32_t N, float32_t SA, float32_t *HX, uint32_t incX, float32_t *HY, uint32_t incY);
void scopy(uint32_t N, const float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY);
float32_t sdot(const uint32_t N, const float32_t *X, const uint32_t incX, const float32_t *Y, const uint32_t incY);
float32_t sdsdot(const uint32_t N, const float32_t alpha, const float32_t *X, const uint32_t incX, const float32_t *Y, const uint32_t incY);
float32_t snrm2(uint64_t N, const float32_t *X, uint64_t incX);
void srot(const uint32_t N, float32_t *X, const uint32_t  incX, float32_t *Y, const uint32_t incY, const float32_t c, const float32_t s);
void srotg(float32_t *a, float32_t *b, float32_t *c, float32_t *s);
void srotm(const uint32_t N, float32_t *X, const uint32_t incX, float32_t *Y, const uint32_t incY, const float32_t *P);
void srotmg(float32_t *D1, float32_t *D2, float32_t *X1, const float32_t y1, float32_t *P);
void sscal(uint32_t N, float32_t SA, float32_t *SX, uint32_t incX);
void sswap(uint32_t N, float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY);
uint64_t isamax(uint32_t N, const float32_t *SX, uint32_t incX);

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

float16_t hasum(uint16_t N, const float16_t *HX, uint16_t incX);
void haxpy(uint16_t N, float16_t HA, float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY);
void hcopy(uint16_t N, const float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY);
float16_t hdot(const uint16_t N, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY);
float16_t hnrm2(uint64_t N, const float16_t *X, uint64_t incX);
void hrot(const uint16_t N, float16_t *X, const uint16_t  incX, float16_t *Y, const uint16_t incY, const float16_t c, const float16_t s);
void hrotg(float16_t *a, float16_t *b, float16_t *c, float16_t *s);
void hrotm(const uint16_t N, float16_t *X, const uint16_t incX, float16_t *Y, const uint16_t incY, const float16_t *P);
void hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P);
void hscal(uint16_t N, float16_t HA, float16_t *HX, uint16_t incX);
float16_t hsdot(const uint16_t N, const float16_t alpha, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY);
void hswap(uint16_t N, float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY);
uint64_t ihamax(uint16_t N, const float16_t *HX, uint16_t incX);

// Level 2

void sgemv(const SB_LAYOUT Layout, const SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *X, const uint64_t incX, const float32_t beta, float32_t *Y, const uint64_t incY);
void dgemv(const SB_LAYOUT Layout, const SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *X, const uint64_t incX, const float64_t beta, float64_t *Y, const uint64_t incY);
void hgemv(const SB_LAYOUT Layout, const SB_TRANSPOSE trans, const uint64_t M, const uint64_t N, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *X, const uint64_t incX, const float16_t beta, float16_t *Y, const uint64_t incY);

// Level 3

void sgemm(const SB_LAYOUT Layout, const SB_TRANSPOSE transA, const SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float32_t alpha, const float32_t *A, const uint64_t lda, const float32_t *B, const uint64_t ldb, const float32_t beta, float32_t *C, const uint64_t ldc);
void dgemm(const SB_LAYOUT Layout, const SB_TRANSPOSE transA, const SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float64_t alpha, const float64_t *A, const uint64_t lda, const float64_t *B, const uint64_t ldb, const float64_t beta, float64_t *C, const uint64_t ldc);
void hgemm(const SB_LAYOUT Layout, const SB_TRANSPOSE transA, const SB_TRANSPOSE transB, const uint64_t M, const uint64_t N, const uint64_t K, const float16_t alpha, const float16_t *A, const uint64_t lda, const float16_t *B, const uint64_t ldb, const float16_t beta, float16_t *C, const uint64_t ldc);
