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

/*
    Demo FLINT program for Sturm sequences
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "fmpz_poly.h"

void print_time(double s)
{
    if (s >= 1.)
    {
        printf("%fs", s);
    }
    else if (s >= 1.e-3)
    {
        printf("%fms", 1.e3 * s);
    }
    else if (s >= 1.e-6)
    {
        printf("%f\xC2\xB5s", 1.e6 * s);
    }
    else
    {
        printf("%fns", 1.e9 * s);
    }
}

long sign_changes(fmpz_poly_t *s, long n, const fmpz_t x)
{
    long c = 0, i;
    int prev = 0, curr;

    fmpz_t y;
    fmpz_init(y);

    for (i = 0; (i < n) && !prev; ++i)
    {
        fmpz_poly_evaluate_fmpz(y, s[i], x);
        prev = fmpz_sgn(y);
    }

    for ( ; i < n; ++i)
    {
        fmpz_poly_evaluate_fmpz(y, s[i], x);
        curr = fmpz_sgn(y);

        if (curr == 0)
            continue;

        if (curr != prev)
        {
            ++c;
            prev = curr;
        }
    }

    fmpz_clear(y);

    return c;
}

long sign_changes_pos_infty(fmpz_poly_t *s, long n)
{
    long c = 0, i = 0;
    int prev, curr;

    prev = fmpz_sgn(s[i]->coeffs + s[i]->length - 1);

    for (++i; i < n; ++i)
    {
        curr = fmpz_sgn(s[i]->coeffs + s[i]->length - 1);
        if (curr != prev)
        {
            ++c;
            prev = curr;
        }
    }

    return c;
}

long sign_changes_neg_infty(fmpz_poly_t *s, long n)
{
    long c = 0, d, i = 0;
    int prev, curr;

    d    = s[i]->length - 1;
    prev = fmpz_sgn(s[i]->coeffs + d) * (d % 2L == 0 ? 1 : -1);

    for (++i; i < n; ++i)
    {
        d    = s[i]->length - 1;
        curr = fmpz_sgn(s[i]->coeffs + d) * (d % 2L == 0 ? 1 : -1);

        if (curr != prev)
        {
            ++c;
            prev = curr;
        }
    }

    return c;
}

int main(int argc, char* argv[])
{
    /* Timing */
    clock_t c0, c1;
    double c;
    slong N = 10000;

    fmpz_poly_t f;
    fmpz_poly_t *seq;
    slong i, n;

    fmpz_poly_init(f);
    {
        /* fmpz_poly_set_str(f, "5  -1 -1 0 1 1"); */
        fmpz_poly_set_str(f, "11  0 7 27 -29 -82 35 75 -15 -26 2 3");
    }

    n = f->length;

    seq = flint_malloc((n+1) * sizeof(fmpz_poly_t));
    for (i = 0; i < (n+1); i++)
        fmpz_poly_init(seq[i]);

    c = 0.;
    for (i = 0; i < N; i++)
    {
        c0 = clock();
    
        fmpz_poly_sturm_sequence(seq, f);
    
        c1 = clock();
        c += (double) (c1 - c0) / CLOCKS_PER_SEC;
    }
    c /= (double) N;

    printf("Time: ");
    print_time(c);
    printf("\n");

    for (i = 0; i < (n+1); i++)
    {
        printf("%ld: ", i);
        fmpz_poly_print(seq[i]);
        printf("\n");
    }

    {
        fmpz_t lb, ub;
        long s, t;

        fmpz_init_set_si(lb, -2);
        fmpz_init_set_si(ub,  2);
        s = sign_changes(seq, f->length, lb);
        t = sign_changes(seq, f->length, ub);

        printf("sigma("); fmpz_print(lb); printf(") = %ld\n", s);
        printf("sigma("); fmpz_print(ub); printf(") = %ld\n", t);
        printf("Number of roots in ("); fmpz_print(lb); printf(", "); fmpz_print(ub); printf("] = %ld\n", s - t);

        fmpz_clear(lb);
        fmpz_clear(ub);
    }

    for (i = 0; i < (n+1); i++)
        fmpz_poly_clear(seq[i]);
    free(seq);

    return EXIT_SUCCESS;
}

