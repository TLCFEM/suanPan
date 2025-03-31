/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*
 * -- SuperLU MT routine (version 3.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley,
 * and Xerox Palo Alto Research Center.
 * September 10, 2007
 *
 */
#include <stdlib.h>
#include "slu_mt_cdefs.h"

/*
 * Generate a banded square matrix A, with dimension n and semi-bandwidth b.
 */
void cband(int_t n, int_t b, int_t nonz, complex** nzval, int_t** rowind, int_t** colptr) {
    int iseed[] = {1992, 1993, 1994, 1995};
    register int_t i, j, ub, lb, ilow, ihigh, lasta = 0;
    complex* a;
    int_t *asub, *xa;
    complex* val;
    int_t* row;
    extern double dlaran_();

    printf("A banded matrix.");
    callocateA(n, nonz, nzval, rowind, colptr); /* Allocate storage */
    a = *nzval;
    asub = *rowind;
    xa = *colptr;
    ub = lb = b;

    for(i = 0; i < 4; ++i) iseed[i] = abs(iseed[i]) % 4096;
    if(iseed[3] % 2 != 1) ++iseed[3];

    for(j = 0; j < n; ++j) {
        xa[j] = lasta;
        val = &a[lasta];
        row = &asub[lasta];
        ilow = SUPERLU_MAX(0, j - ub);
        ihigh = SUPERLU_MIN(n - 1, j + lb);
        for(i = ilow; i <= ihigh; ++i) {
            val[i - ilow].r = dlaran_(iseed);
            row[i - ilow] = i;
        }
        lasta += ihigh - ilow + 1;
    } /* for j ... */
    xa[n] = lasta;
}

/*
 * Generate a block diagonal matrix A.
 */
void cblockdiag(int_t nb, /* number of blocks */
                int_t bs, /* block size */
                int_t nonz,
                complex** nzval,
                int_t** rowind,
                int_t** colptr) {
    int iseed[] = {1992, 1993, 1994, 1995};
    register int_t i, j, b, n, lasta = 0, cstart, rstart;
    complex* a;
    int_t *asub, *xa;
    complex* val;
    int_t* row;
    extern double dlaran_();

    n = bs * nb;
    printf("A block diagonal matrix: nb " IFMT ", bs " IFMT ", n " IFMT "\n", nb, bs, n);
    callocateA(n, nonz, nzval, rowind, colptr); /* Allocate storage */
    a = *nzval;
    asub = *rowind;
    xa = *colptr;

    for(i = 0; i < 4; ++i) iseed[i] = abs(iseed[i]) % 4096;
    if(iseed[3] % 2 != 1) ++iseed[3];

    for(b = 0; b < nb; ++b) {
        cstart = b * bs; /* start of the col # of the current block */
        rstart = b * bs; /* start of the row # of the current block */
        for(j = cstart; j < cstart + bs; ++j) {
            xa[j] = lasta;
            val = &a[lasta];
            row = &asub[lasta];
            for(i = 0; i < bs; ++i) {
                val[i].r = dlaran_(iseed);
                row[i] = i + rstart;
            }
            lasta += bs;
        } /* for j ... */
    } /* for b ... */

    xa[n] = lasta;
}
