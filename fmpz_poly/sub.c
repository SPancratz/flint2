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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Sebastian Pancratz

******************************************************************************/

#include "fmpz_poly.h"

#define SET(x, y) \
    fmpz_set((x), (y))
#define NEG(x, y) \
    fmpz_neg((x), (y))
#define VEC_SUB(v, v1, v2, len) \
    _fmpz_vec_sub((v), (v1), (v2), (len))

void _fmpz_poly_sub(fmpz * res, const fmpz * poly1, long len1, 
                                const fmpz * poly2, long len2)
{
    #include "generics/poly_sub.in"
}

void
fmpz_poly_sub(fmpz_poly_t res, const fmpz_poly_t poly1,
              const fmpz_poly_t poly2)
{
    long max = FLINT_MAX(poly1->length, poly2->length);

    fmpz_poly_fit_length(res, max);

    _fmpz_poly_sub(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs,
                   poly2->length);

    _fmpz_poly_set_length(res, max);
    _fmpz_poly_normalise(res);  /* there may have been cancellation */
}
