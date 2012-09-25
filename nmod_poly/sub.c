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

    Copyright (C) 2007, David Howden
    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
   
******************************************************************************/

#include "nmod_poly.h"

#define SET(x, y) \
    *(x) = *(y)
#define NEG(x, y) \
    *(x)= nmod_neg(*(y), mod)
#define VEC_SUB(v, v1, v2, len) \
    _nmod_vec_sub((v), (v1), (v2), (len), mod)

void _nmod_poly_sub(mp_ptr res, mp_srcptr poly1, long len1, 
                                mp_srcptr poly2, long len2, nmod_t mod)
{
    #include "generics/poly_sub.in"
}

void nmod_poly_sub(nmod_poly_t res, const nmod_poly_t poly1, 
                                    const nmod_poly_t poly2)
{
    const long max = FLINT_MAX(poly1->length, poly2->length);

    nmod_poly_fit_length(res, max);

    _nmod_poly_sub(res->coeffs, poly1->coeffs, poly1->length, 
                                poly2->coeffs, poly2->length, poly1->mod);

    res->length = max;
    _nmod_poly_normalise(res);  /* there may have been cancellation */
}
