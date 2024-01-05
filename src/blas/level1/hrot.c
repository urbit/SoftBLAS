void hrot(const uint64_t  N,       // length
             float16_t      *X,       // vector
             const uint64_t  incX,    // stride
             float16_t      *Y,       // vector
             const uint64_t  incY,    // stride
             const float16_t c,       // cos(theta)
             const float16_t s        // sin(theta)
          ) {
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
