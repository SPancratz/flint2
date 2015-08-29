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
    Demo FLINT program
    
    Determines whether all roots of an integer polynomial lie 
    in the interval [a, b]
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

int main(int argc, char* argv[])
{

    fmpz_poly_t f;
    fmpz_t a, b;
    slong nMax = 11;
    int r;

    fmpz *w = _fmpz_vec_init(3 * nMax + 8);

    fmpz_poly_init(f);
    fmpz_init(a);
    fmpz_init(b);

    {
        fmpz_poly_set_str(f, "11  0 7 27 -29 -82 35 75 -15 -26 2 3");
        fmpz_set_si(a, -2);
        fmpz_set_si(b,  2);

        printf("f(x) = ");
        fmpz_poly_print_pretty(f, "x");
        printf("\n");

        r = _fmpz_poly_all_roots_in_interval(f->coeffs, f->length, a, b, w);
        
        printf("All roots in [");
        fmpz_print(a);
        printf(", ");
        fmpz_print(b);
        printf("]? %d\n", r);
    }

    {
        /* Timing */
        clock_t c0, c1;
        double c;
        slong i, N = 10000;

        c = 0.;
        for (i = 0; i < N; i++)
        {
            c0 = clock();

            r = _fmpz_poly_all_roots_in_interval(f->coeffs, f->length, a, b, w);
        
            c1 = clock();
            c += (double) (c1 - c0) / CLOCKS_PER_SEC;
        }
        c /= (double) N;

        printf("Time: ");
        print_time(c);
        printf("\n");
    }

    _fmpz_vec_clear(w, 3 * nMax + 8);
    fmpz_poly_clear(f);
    fmpz_clear(a);
    fmpz_clear(b);

    return EXIT_SUCCESS;
}

