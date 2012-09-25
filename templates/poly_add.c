/*
    Implements polynomial addition:

        POLY_ADD(res, poly1, len1, poly2, len2)

    Requires the following templates:

        * SET(x, y)
        * VEC_ADD(v, v1, v2, len)
 */

long i, min = FLINT_MIN(len1, len2);

VEC_ADD(res, poly1, poly2, min);

if (poly1 != res)
    for (i = min; i < len1; i++)
        SET(res + i, poly1 + i);

if (poly2 != res)
    for (i = min; i < len2; i++)
        SET(res + i, poly2 + i);
