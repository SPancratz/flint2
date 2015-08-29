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

    Copyright (C) 2015 Sebastian Pancratz

******************************************************************************/

#include "fmpz_poly.h"

/*
    Assumes that seq is an array of length n + 1 of initialized 
    polynomials.
 */
void fmpz_poly_sturm_sequence(fmpz_poly_t *seq, const fmpz_poly_t poly)
{
    slong i;
    slong n = poly->length;
    ulong d;
    fmpz_t c;

    fmpz_init(c);

    fmpz_poly_set(seq[0], poly);
    fmpz_poly_derivative(seq[1], poly);

    for (i = 2; i < n + 1 && !fmpz_poly_is_zero(seq[i-1]); ++i)
    {
        fmpz_poly_content(c, seq[i-1]);
        if (!fmpz_is_one(c))
            fmpz_poly_scalar_divexact_fmpz(seq[i-1], seq[i-1], c);

        /*
            Pseudo division over Z[X]:
            ell^d seq[i-2] = Q seq[i-1] + R 
         */
        fmpz_poly_pseudo_rem(seq[i], &d, seq[i-2], seq[i-1]);

        /*
            If the remainder were given over Q[X], we would 
            reverse the sign.  Here, using pseudo division 
            over Z[X], we carefully consider the sign of 
            ell^d, too
         */
        if (d % 2 == 0 || fmpz_sgn(fmpz_poly_lead(seq[i-1])) > 0)
            fmpz_poly_neg(seq[i], seq[i]);
    }

    for ( ; i < n + 1; ++i)
        fmpz_poly_zero(seq[i]);

    fmpz_clear(c);
}

