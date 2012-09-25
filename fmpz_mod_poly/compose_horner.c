/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "fmpz_mod_poly.h"

#define TYPE \
    fmpz
#define SET(x, y) \
    fmpz_set((x), (y))
#define ADD(x, y, z) \
    do {                              \
        fmpz_add((x), (y), (z));      \
        if (fmpz_cmpabs((x), p) >= 0) \
            fmpz_sub((x), (x), p);    \
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

void _fmpz_mod_poly_compose_horner(fmpz *res, const fmpz *poly1, long len1, 
                                              const fmpz *poly2, long len2, 
                                              const fmpz_t p)
{
    #include "generics/poly_compose_horner.in"
}

void fmpz_mod_poly_compose_horner(fmpz_mod_poly_t res, 
                                  const fmpz_mod_poly_t poly1, 
                                  const fmpz_mod_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;

    if (len1 == 0)
    {
        fmpz_mod_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fmpz_mod_poly_set_fmpz(res, poly1->coeffs);
    }
    else
    {
        const long lenr = (len1 - 1) * (len2 - 1) + 1;

        if ((res != poly1) && (res != poly2))
        {
            fmpz_mod_poly_fit_length(res, lenr);
            _fmpz_mod_poly_compose_horner(res->coeffs, poly1->coeffs, len1, 
                                                       poly2->coeffs, len2,
                                                       &(res->p));
        }
        else
        {
            fmpz *t = _fmpz_vec_init(lenr);

            _fmpz_mod_poly_compose_horner(t, poly1->coeffs, len1,
                                             poly2->coeffs, len2, &(res->p));
            _fmpz_vec_clear(res->coeffs, res->alloc);
            res->coeffs = t;
            res->alloc  = lenr;
            res->length = lenr;
        }

        _fmpz_mod_poly_set_length(res, lenr);
        _fmpz_mod_poly_normalise(res);
    }
}
