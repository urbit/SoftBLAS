#include "softfloat.h"

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
sasum(uint32_t N, const float32_t *SX, uint32_t incX);
saxpy(uint32_t N, float32_t SA, float32_t *HX, uint32_t incX, float32_t *HY, uint32_t incY);
scopy(uint32_t N, const float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY);
sdot(const uint32_t N, const float32_t *X, const uint32_t incX, const float32_t *Y, const uint32_t incY);
srot(const uint32_t N, float32_t *X, const uint32_t  incX, float32_t *Y, const uint32_t incY, const float32_t c, const float32_t s);
srotg(float32_t *a, float32_t *b, float32_t *c, float32_t *s);
srotm(const uint32_t N, float32_t *X, const uint32_t incX, float32_t *Y, const uint32_t incY, const float32_t *P);
srotmg(float32_t *D1, float32_t *D2, float32_t *X1, const float32_t y1, float32_t *P)
sscal(uint32_t N, float32_t SA, float32_t *SX, uint32_t incX);
ssdot( const uint32_t N, const float32_t alpSA, const float32_t *X, const uint32_t incX, const float32_t *Y, const uint32_t incY);
sswap(uint32_t N, float32_t *SX, uint32_t incX, float32_t *SY, uint32_t incY)
isamax(uint32_t N, const float32_t *SX, uint32_t incX);

dasum(uint64_t N, const float64_t *DX, uint64_t incX);
daxpy(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
dcopy(uint64_t N, const float64_t *DX, uint64_t incX, float64_t *DY, uint64_t incY);
ddot(const uint64_t N, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY);
drot(const uint64_t N, float64_t *X, const uint64_t  incX, float64_t *Y, const uint64_t incY, const float64_t c, const float64_t s);
drotg(float64_t *a, float64_t *b, float64_t *c, float64_t *s);
drotm(const uint64_t N, float64_t *X, const uint64_t incX, float64_t *Y, const uint64_t incY, const float64_t *P);
drotmg(float64_t *D1, float64_t *D2, float64_t *X1, const float64_t y1, float64_t *P)
dscal(uint64_t N, float64_t DA, float64_t *DX, uint64_t incX);
dsdot( const uint64_t N, const float64_t alpha, const float64_t *X, const uint64_t incX, const float64_t *Y, const uint64_t incY);
dswap(uint64_t N, float64_t *DX, uint64_t incX, float64_t *SY, uint64_t incY)
idamax(uint64_t N, const float64_t *DX, uint64_t incX);

hasum(uint16_t N, const float16_t *HX, uint16_t incX);
haxpy(uint16_t N, float16_t HA, float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY);
hcopy(uint16_t N, const float16_t *HX, uint16_t incX, float16_t *HY, uint16_t incY);
hdot(const uint16_t N, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY);
hrot(const uint16_t N, float16_t *X, const uint16_t  incX, float16_t *Y, const uint16_t incY, const float16_t c, const float16_t s);
hrotg(float16_t *a, float16_t *b, float16_t *c, float16_t *s);
hrotm(const uint16_t N, float16_t *X, const uint16_t incX, float16_t *Y, const uint16_t incY, const float16_t *P);
hrotmg(float16_t *D1, float16_t *D2, float16_t *X1, const float16_t y1, float16_t *P)
hscal(uint16_t N, float16_t HA, float16_t *HX, uint16_t incX);
hsdot( const uint16_t N, const float16_t alpha, const float16_t *X, const uint16_t incX, const float16_t *Y, const uint16_t incY);
hswap(uint16_t N, float16_t *SX, uint16_t incX, float16_t *SY, uint16_t incY)
ihamax(uint16_t N, const float16_t *HX, uint16_t incX);
