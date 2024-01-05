void drot(const uint64_t N, float64_t *X, const uint64_t  incX, float64_t *Y, const uint64_t incY, const float64_t c, const float64_t s) {
    float64_t tmp;

    if (c != SB_REAL64_ONE || s != SB_REAL64_ZERO) {
        if (incX == 1 && incY == 1) {
            for (uint64_t i=0; i != N; i++) {
                tmp = f64_add(f64_mul(c, X[i]), f64_mul(s, Y[i]));
                Y[i] = f64_sub(f64_mul(c, Y[i]), f64_mul(s, X[i]));
                X[i] = tmp;
            }
        }
        else
        {
            for (uint64_t i=N; i; i--, Y += incY, X += incX)
            {
                tmp = f64_add(f64_mul(c, *X), f64_mul(s, *Y));
                *Y = f64_sub(f64_mul(c, *Y), f64_mul(s, *X));
                *X = tmp;
            }
        }
    }
}