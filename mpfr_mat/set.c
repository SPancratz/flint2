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

    Copyright (C) 2008-2009 William Hart
    Copyright (C) 2010 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_mat.h"

void
mpfr_mat_set(mpfr_mat_t mat1, const mpfr_mat_t mat2)
{
    if (mat1 != mat2)
    {
        slong i;

        if (mat2->r && mat2->c)
            for (i = 0; i < mat2->r; i++)
                _mpfr_vec_set(mat1->rows[i], mat2->rows[i], mat2->c);
    }
}
