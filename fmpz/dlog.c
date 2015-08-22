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

******************************************************************************/

#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

double
fmpz_dlog(const fmpz_t x)
{
    if (!COEFF_IS_MPZ(*x))
    {
        return log(*x);
    }
    else
    {
        double s;

#if defined(__MPIR_VERSION)
        slong e;
#else
        long e;
#endif        

        s = mpz_get_d_2exp(&e, COEFF_TO_PTR(*x));
       
        return log(s) + e * 0.69314718055994530942;
    }
}
