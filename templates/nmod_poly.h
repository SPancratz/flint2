#define TYPE \
    mp_limb_t

#define SET(x, y) \
    *(x) = *(y)

#define ADD(x, y, z) \
    *(x) = n_addmod(*(y), *(z), mod.n)

#define SUB(x, y, z) \
    *(x) = n_submod(*(y), *(z), mod.n)

#define NEG(x, y) \
    *(x)= nmod_neg(*(y), mod)

#define VEC_INIT(len) \
    _nmod_vec_init(len)

#define VEC_CLEAR(v, len) \
    _nmod_vec_clear(v)

#define VEC_ADD(v, v1, v2, len) \
    _nmod_vec_add((v), (v1), (v2), (len), mod)

#define VEC_SUB(v, v1, v2, len) \
    _nmod_vec_sub((v), (v1), (v2), (len), mod)

#define VEC_SCALAR_MUL(v, w, len, c) \
    _nmod_vec_scalar_mul_nmod((v), (w), (len), *(c), mod)

#define POLY_MUL(rop, op1, len1, op2, len2) \
    _nmod_poly_mul((rop), (op1), (len1), (op2), (len2), mod)

#define POLY_EVALUATE(rop, op, len, c) \
    *(rop) = _nmod_poly_evaluate_nmod((op), (len), *(c), mod)

