#define TYPE \
    fmpz

#define SET(x, y) \
    fmpz_set((x), (y))

#define ADD(x, y, z)                  \
    do {                              \
        fmpz_add((x), (y), (z));      \
        if (fmpz_cmpabs((x), p) >= 0) \
            fmpz_sub((x), (x), p);    \
    } while (0)

#define SUB(x, y, z)               \
    do {                           \
        fmpz_sub((x), (y), (z));   \
        if (fmpz_sgn(x) < 0)       \
            fmpz_add((x), (x), p); \
    } while (0)

#define NEG(x, y)                  \
    do {                           \
        if (*(y))                  \
            fmpz_sub((x), p, (y)); \
        else                       \
            fmpz_zero(y);          \
    } while (0)

#define VEC_INIT(len) \
    _fmpz_vec_init(len)

#define VEC_CLEAR(v, len) \
    _fmpz_vec_clear((v), (len))

#define VEC_SCALAR_MUL(v, w, len, c) \
    _fmpz_mod_poly_scalar_mul_fmpz((v), (w), (len), (c), p)

#define POLY_MUL(rop, op1, len1, op2, len2) \
    _fmpz_mod_poly_mul((rop), (op1), (len1), (op2), (len2), p)

#define POLY_EVALUATE(rop, op, len, c) \
    _fmpz_mod_poly_evaluate_fmpz((rop), (op), (len), (c), p)

