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
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "mpf_mat.h"

void
mpf_mat_randtest(mpf_mat_t mat, flint_rand_t state, mp_bitcnt_t bits)
{
    slong r, c, i, j;

    r = mat->r;
    c = mat->c;

    _flint_rand_init_gmp(state);

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            mpf_urandomb(mpf_mat_entry(mat, i, j), state->gmp_state, bits);
}
