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

    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_mat.h"

int
fmpz_mat_get_d_mat(d_mat_t B, const fmpz_mat_t A)
{
    slong i, j;
    fmpz_t dmax;

    fmpz_init(dmax);
    fmpz_set_d(dmax, DBL_MAX);

    for (i = 0; i < A->r; i++)
    {
        for (j = 0; j < A->c; j++)
        {
            if (fmpz_cmpabs(fmpz_mat_entry(A, i, j), dmax) > 0)
            {
                fmpz_clear(dmax);
                return -1;
            }
            d_mat_entry(B, i, j) = fmpz_get_d(fmpz_mat_entry(A, i, j));
        }
    }
    fmpz_clear(dmax);
    return 0;
}
