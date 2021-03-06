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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen
 
******************************************************************************/

#include <stdio.h>
#include <string.h>

#include "fmpz_vec.h"
#include "padic.h"
#include "qadic.h"

int flint_conway_polynomials [] = {
#if !defined( _MSC_VER)
#include "../build/CPimport.h"
#else
#include "CPimport.h"
#endif
  0
};

void qadic_ctx_init_conway(qadic_ctx_t ctx,
                           const fmpz_t p, slong d, slong min, slong max, 
                           const char *var, enum padic_print_mode mode)
{
    unsigned int position;

    if (fmpz_cmp_ui(p, 109987) > 0)
    {
        flint_printf("Exception (qadic_ctx_init_conway).  Conway polynomials \n");
        flint_printf("are only available for primes up to 109987.\n");
        abort();
    }

    for (position = 0; flint_conway_polynomials[position] != 0; position += 3+flint_conway_polynomials[position+1])
    {
        /* Different prime? */
        if (fmpz_cmp_ui(p, flint_conway_polynomials[position]))
            continue;

        /* Same degree? */
        if (d == flint_conway_polynomials[position+1])
        {
            slong i, j;

            /* Find number of non-zero coefficients */
            ctx->len = 1;

            for (i = 0; i < d; i++)
            {
                if (flint_conway_polynomials[position+2+i])
                    ctx->len ++;
            }

            ctx->a = _fmpz_vec_init(ctx->len);
            ctx->j = flint_malloc(ctx->len * sizeof(slong));

            /* Copy the polynomial */
            j = 0;

            for (i = 0; i < d; i++)
            {
                int coeff = flint_conway_polynomials[position+2+i];

                if (coeff)
                {
                    fmpz_set_ui(ctx->a + j, coeff);
                    ctx->j[j] = i;
                    j++;
                }
            }

            fmpz_set_ui(ctx->a + j, 1);
            ctx->j[j] = d;

            /* Complete the initialisation of the context */
            padic_ctx_init(&ctx->pctx, p, min, max, mode);

            ctx->var = flint_malloc(strlen(var) + 1);
            strcpy(ctx->var, var);

            return;
        }
    }

    flint_printf("Exception (qadic_ctx_init_conway).  The polynomial for \n");
    flint_printf("(p,d) = (%wd,%wd) is not present in the database.\n", *p, d);
    abort();
}

