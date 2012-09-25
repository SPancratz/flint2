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

    Copyright (C) 2010 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz

******************************************************************************/

#include "nmod_poly.h"

void
_nmod_poly_compose(mp_ptr res, mp_srcptr poly1, long len1, 
                               mp_srcptr poly2, long len2, nmod_t mod)
{
    if (len1 == 1)
        res[0] = poly1[0];
    else if (len2 == 1)
        res[0] = _nmod_poly_evaluate_nmod(poly1, len1, poly2[0], mod);
    else if (len1 <= 7)
        _nmod_poly_compose_horner(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_compose_divconquer(res, poly1, len1, poly2, len2, mod);
}

void nmod_poly_compose(nmod_poly_t res, 
                       const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    
    if (len1 == 0)
    {
        nmod_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        nmod_poly_fit_length(res, 1);
        res->coeffs[0] = poly1->coeffs[0];
        res->length = (res->coeffs[0] != 0);
    }
    else
    {
        const long lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            nmod_poly_fit_length(res, lenr);
            _nmod_poly_compose_horner(res->coeffs, poly1->coeffs, len1, 
                                                   poly2->coeffs, len2, poly1->mod);
        }
        else
        {
            nmod_poly_t t;
            nmod_poly_init2(t, poly1->mod.n, lenr);
            _nmod_poly_compose_horner(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, poly1->mod);
            nmod_poly_swap(res, t);
            nmod_poly_clear(t);
        }

        res->length = lenr;
        _nmod_poly_normalise(res);
    }
}
