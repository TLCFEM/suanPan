/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file cpivotL.c
 * \brief Performs numerical pivoting
 *
 * <pre>
 * -- SuperLU routine (version 7.0.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
 * August 2024
 *
 * Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 *
 * THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
 * EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 *
 * Permission is hereby granted to use or copy this program for any
 * purpose, provided the above notices are retained on all copies.
 * Permission to modify the code and to distribute modified code is
 * granted, provided the above notices are retained, and a notice that
 * the code was modified is included with the above copyright notice.
 * </pre>
 */

#include <math.h>
#include <stdlib.h>
#include "slu_cdefs.h"

#undef DEBUG

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *   Performs the numerical pivoting on the current column of L,
 *   and the CDIV operation.
 *
 *   Pivot policy:
 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
 *           pivot row = k;
 *       ELSE IF abs(A_jj) >= thresh THEN
 *           pivot row = j;
 *       ELSE
 *           pivot row = m;
 *
 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
 *
 *   Return value: 0      success;
 *                 i > 0  U(i,i) is exactly zero.
 * </pre>
 */

int cpivotL(
    const int jcol,     /* in */
    const double u,     /* in - diagonal pivoting threshold */
    int* usepr,         /* re-use the pivot sequence given by perm_r/iperm_r */
    int* perm_r,        /* may be modified */
    int* iperm_r,       /* in - inverse of perm_r */
    int* iperm_c,       /* in - used to find diagonal of Pc*A*Pc' */
    int* pivrow,        /* out */
    GlobalLU_t* Glu,    /* modified - global LU data structures */
    SuperLUStat_t* stat /* output */
) {
    singlecomplex one = {1.0, 0.0};
    int fsupc;  /* first column in the supernode */
    int nsupc;  /* no of columns in the supernode */
    int nsupr;  /* no of rows in the supernode */
    int_t lptr; /* points to the starting subscript of the supernode */
    int pivptr, old_pivptr, diag, diagind;
    float pivmax, rtemp, thresh;
    singlecomplex temp;
    singlecomplex* lu_sup_ptr;
    singlecomplex* lu_col_ptr;
    int_t* lsub_ptr;
    int_t isub, icol, k, itemp;
    int_t *lsub, *xlsub;
    singlecomplex* lusup;
    int_t* xlusup;
    flops_t* ops = stat->ops;

    /* Initialize pointers */
    lsub = Glu->lsub;
    xlsub = Glu->xlsub;
    lusup = (singlecomplex*)Glu->lusup;
    xlusup = Glu->xlusup;
    fsupc = (Glu->xsup)[(Glu->supno)[jcol]];
    nsupc = jcol - fsupc; /* excluding jcol; nsupc >= 0 */
    lptr = xlsub[fsupc];
    nsupr = xlsub[fsupc + 1] - lptr;
    lu_sup_ptr = &lusup[xlusup[fsupc]]; /* start of the current supernode */
    lu_col_ptr = &lusup[xlusup[jcol]];  /* start of jcol in the supernode */
    lsub_ptr = &lsub[lptr];             /* start of row indices of the supernode */

#ifdef DEBUG
    if(jcol == MIN_COL) {
        printf("Before cdiv: col %d\n", jcol);
        for(k = nsupc; k < nsupr; k++)
            printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
    }
#endif

    /* Determine the largest abs numerical value for partial pivoting;
       Also search for user-specified pivot, and diagonal element. */
    if(*usepr) *pivrow = iperm_r[jcol];
    diagind = iperm_c[jcol];
    pivmax = 0.0;
    pivptr = nsupc;
    diag = SLU_EMPTY;
    old_pivptr = nsupc;
    for(isub = nsupc; isub < nsupr; ++isub) {
        rtemp = c_abs1(&lu_col_ptr[isub]);
        if(rtemp > pivmax) {
            pivmax = rtemp;
            pivptr = isub;
        }
        if(*usepr && lsub_ptr[isub] == *pivrow) old_pivptr = isub;
        if(lsub_ptr[isub] == diagind) diag = isub;
    }

    /* Test for singularity */
    if(pivmax == 0.0) {
#if 0
        // There is no valid pivot.
        // jcol represents the rank of U, 
        // report the rank, let dgstrf handle the pivot
	*pivrow = lsub_ptr[pivptr];
	perm_r[*pivrow] = jcol;
#endif
        *usepr = 0;
        return (jcol + 1);
    }

    thresh = u * pivmax;

    /* Choose appropriate pivotal element by our policy. */
    if(*usepr) {
        rtemp = c_abs1(&lu_col_ptr[old_pivptr]);
        if(rtemp != 0.0 && rtemp >= thresh)
            pivptr = old_pivptr;
        else
            *usepr = 0;
    }
    if(*usepr == 0) {
        /* Use diagonal pivot? */
        if(diag >= 0) { /* diagonal exists */
            rtemp = c_abs1(&lu_col_ptr[diag]);
            if(rtemp != 0.0 && rtemp >= thresh) pivptr = diag;
        }
        *pivrow = lsub_ptr[pivptr];
    }

    /* Record pivot row */
    perm_r[*pivrow] = jcol;

    /* Interchange row subscripts */
    if(pivptr != nsupc) {
        itemp = lsub_ptr[pivptr];
        lsub_ptr[pivptr] = lsub_ptr[nsupc];
        lsub_ptr[nsupc] = itemp;

        /* Interchange numerical values as well, for the whole snode, such
         * that L is indexed the same way as A.
         */
        for(icol = 0; icol <= nsupc; icol++) {
            itemp = pivptr + icol * nsupr;
            temp = lu_sup_ptr[itemp];
            lu_sup_ptr[itemp] = lu_sup_ptr[nsupc + icol * nsupr];
            lu_sup_ptr[nsupc + icol * nsupr] = temp;
        }
    } /* if */

    /* cdiv operation */
    ops[FACT] += 10 * (nsupr - nsupc);

    c_div(&temp, &one, &lu_col_ptr[nsupc]);
    for(k = nsupc + 1; k < nsupr; k++)
        cc_mult(&lu_col_ptr[k], &lu_col_ptr[k], &temp);

    return 0;
}
