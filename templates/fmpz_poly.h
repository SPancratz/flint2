#define TYPE \
    fmpz

#define SET(x, y) \
    fmpz_set((x), (y))

#define ADD(x, y, z) \
    fmpz_add((x), (y), (z))

#define SUB(x, y, z) \
    fmpz_add((x), (y), (z))

#define NEG(x, y) \
    fmpz_neg((x), (y))

#define VEC_INIT(len) \
    _fmpz_vec_init(len)

#define VEC_CLEAR(v, len) \
    _fmpz_vec_clear((v), (len))

#define VEC_ADD(v, v1, v2, len) \
    _fmpz_vec_add((v), (v1), (v2), (len))

#define VEC_SUB(v, v1, v2, len) \
    _fmpz_vec_sub((v), (v1), (v2), (len))

#define VEC_SCALAR_MUL(v, w, len, c) \
    _fmpz_vec_scalar_mul_fmpz((v), (w), (len), (c))

#define POLY_MUL(rop, op1, len1, op2, len2) \
    _fmpz_poly_mul((rop), (op1), (len1), (op2), (len2))

#define POLY_EVALUATE(rop, op, len, c) \
    _fmpz_poly_evaluate_fmpz((rop), (op), (len), (c))

