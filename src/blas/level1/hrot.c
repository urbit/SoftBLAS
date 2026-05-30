void hrot(const uint16_t N, float16_t *X, const uint16_t  incX, float16_t *Y, const uint16_t incY, const float16_t c, const float16_t s) {
    _set_rounding(rndMode);
    float16_t tmp;

    if (c != SB_REAL32_ONE || s != SB_REAL32_ZERO) {
        if (incX == 1 && incY == 1) {
            for (uint64_t i=0; i != N; i++) {
                tmp = f16_add(f16_mul(c, X[i]), f16_mul(s, Y[i]));
                Y[i] = f16_sub(f16_mul(c, Y[i]), f16_mul(s, X[i]));
                X[i] = tmp;
            }
        }
        else
        {
            for (uint64_t i=N; i; i--, Y += incY, X += incX)
            {
                tmp = f16_add(f16_mul(c, *X), f16_mul(s, *Y));
                *Y = f16_sub(f16_mul(c, *Y), f16_mul(s, *X));
                *X = tmp;
            }
        }
    }
}
