#include "softfloat.h"

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
#define MAX(x, y) ( (x) > (y) ? (x) : (y) )
#define MIN(x, y) ( (x) > (y) ? (y) : (x) )
