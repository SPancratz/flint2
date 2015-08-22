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

#include "mpf_vec.h"

int
_mpf_vec_dot2(mpf_t res, const mpf * vec1, const mpf * vec2, slong len2,
              mp_bitcnt_t prec)
{
    slong i;
    int r = 0;
    mpf_t tmp, tmp2;
    mpf_init2(tmp, prec);
    mpf_init2(tmp2, prec);

    flint_mpf_set_ui(res, 0);
    for (i = 0; i < len2; i++)
    {
        mpf_mul(tmp, vec1 + i, vec2 + i);
        mpf_add(res, res, tmp);
    }

    _mpf_vec_norm(tmp, vec1, len2);
    _mpf_vec_norm(tmp2, vec2, len2);
    mpf_mul(tmp, tmp, tmp2);
    mpf_div_2exp(tmp, tmp, prec);
    mpf_mul(tmp2, res, res);

    if (mpf_cmp(tmp2, tmp) > 0)
        r = 1;

    mpf_clear(tmp);
    mpf_clear(tmp2);

    return r;
}
