/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*! @file izmax1.c
 * \brief Finds the index of the element whose real part has maximum absolute value
 *
 * <pre>
 *     -- LAPACK auxiliary routine (version 2.0) --
 *     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
 *     Courant Institute, Argonne National Lab, and Rice University
 *     October 31, 1992
 * </pre>
 */
#include <math.h>
#include "slu_dcomplex.h"
#include "slu_Cnames.h"

/*! \brief

<pre>
    Purpose
    =======

    IZMAX1 finds the index of the element whose real part has maximum
    absolute value.

    Based on IZAMAX from Level 1 BLAS.
    The change is to use the 'genuine' absolute value.

    Contributed by Nick Higham for use with ZLACON.

    Arguments
    =========

    N       (input) INT
            The number of elements in the vector CX.

    CX      (input) COMPLEX*16 array, dimension (N)
            The vector whose elements will be summed.

    INCX    (input) INT
            The spacing between successive values of CX.  INCX >= 1.

   =====================================================================
</pre>
*/

int izmax1_slu(int* n, doublecomplex* cx, int* incx) {
    /* System generated locals */
    int ret_val;
    double d__1;

    /* Local variables */
    double smax;
    int i, ix;

#define CX(I) cx[(I) - 1]

    ret_val = 0;
    if(*n < 1) {
        return ret_val;
    }
    ret_val = 1;
    if(*n == 1) {
        return ret_val;
    }
    if(*incx == 1) {
        goto L30;
    }

    /*     CODE FOR INCREMENT NOT EQUAL TO 1 */

    ix = 1;
    smax = (d__1 = CX(1).r, fabs(d__1));
    ix += *incx;
    for(i = 2; i <= *n; ++i) {
        if((d__1 = CX(ix).r, fabs(d__1)) <= smax) {
            goto L10;
        }
        ret_val = i;
        smax = (d__1 = CX(ix).r, fabs(d__1));
    L10:
        ix += *incx;
        /* L20: */
    }
    return ret_val;

    /*     CODE FOR INCREMENT EQUAL TO 1 */

L30:
    smax = (d__1 = CX(1).r, fabs(d__1));
    for(i = 2; i <= *n; ++i) {
        if((d__1 = CX(i).r, fabs(d__1)) <= smax) {
            goto L40;
        }
        ret_val = i;
        smax = (d__1 = CX(i).r, fabs(d__1));
    L40:;
    }
    return ret_val;

    /*     End of IZMAX1 */

} /* izmax1_slu */
