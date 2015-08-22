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

    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
   int i, result;
   FLINT_TEST_INIT(state);
   
   flint_printf("powmod2_ui_preinv....");
   fflush(stdout);

   

   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t a, d, r1, r2, dinv;
      mpz_t a_m, d_m, r2_m;
      mp_limb_t exp;

      mpz_init(a_m);
      mpz_init(d_m);
      mpz_init(r2_m);
      
      d = n_randtest_not_zero(state);
      do
      {
         a = n_randtest(state) % d;
      } while (n_gcd(d, a) != UWORD(1));
      exp = n_randtest(state);
      
      dinv = n_preinvert_limb(d);
      r1 = n_powmod2_ui_preinv(a, exp, d, dinv);

      flint_mpz_set_ui(a_m, a);
      flint_mpz_set_ui(d_m, d);
      flint_mpz_powm_ui(r2_m, a_m, exp, d_m);      
      r2 = flint_mpz_get_ui(r2_m);
      
      result = (r1 == r2);
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("a = %wu, exp = %wu, d = %wu\n", a, exp, d); 
         flint_printf("r1 = %wu, r2 = %wu\n", r1, r2);
         abort();
      }

      mpz_clear(a_m);
      mpz_clear(d_m);
      mpz_clear(r2_m);
   }

   /* check 0^0 = 1 */
   for (i = 0; i < 10000 * flint_test_multiplier(); i++)
   {
      mp_limb_t d, r, dinv;
      
      d = n_randtest_not_zero(state);
      
      dinv = n_preinvert_limb(d);
      r = n_powmod2_ui_preinv(0, 0, d, dinv);

      result = (r == 1 || (d == 1 && r == 0));
      if (!result)
      {
         flint_printf("FAIL:\n");
         flint_printf("0^0 != 1 mod %wd\n", d);
         abort();
      }
   }

   FLINT_TEST_CLEANUP(state);
   
   flint_printf("PASS\n");
   return 0;
}
