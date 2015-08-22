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

#include <math.h>
#include "fmpz_vec.h"

slong
_fmpz_vec_get_d_vec_2exp(double *appv, const fmpz * vec, slong len)
{
    slong *exp, i, maxexp = 0L;
    exp = (slong *) malloc(len * sizeof(slong));

    for (i = 0; i < len; i++)
    {
        appv[i] = fmpz_get_d_2exp(&exp[i], vec + i);
        if (exp[i] > maxexp)
            maxexp = exp[i];
    }

    for (i = 0; i < len; i++)
        appv[i] = ldexp(appv[i], exp[i] - maxexp);

    free(exp);
    return maxexp;
}
