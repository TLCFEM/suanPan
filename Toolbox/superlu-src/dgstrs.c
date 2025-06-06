/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file dgstrs.c
 * \brief Solves a system using LU factorization
 *
 * <pre>
 * -- SuperLU routine (version 3.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * October 15, 2003
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

#include "slu_ddefs.h"

/*! \brief
 *
 * <pre>
 * Purpose
 * =======
 *
 * DGSTRS solves a system of linear equations A*X=B or A'*X=B
 * with A sparse and B dense, using the LU factorization computed by
 * DGSTRF.
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * trans   (input) trans_t
 *          Specifies the form of the system of equations:
 *          = NOTRANS: A * X = B  (No transpose)
 *          = TRANS:   A'* X = B  (Transpose)
 *          = CONJ:    A**H * X = B  (Conjugate transpose)
 *
 * L       (input) SuperMatrix*
 *         The factor L from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use compressed row subscripts storage for supernodes,
 *         i.e., L has types: Stype = SLU_SC, Dtype = SLU_D, Mtype = SLU_TRLU.
 *
 * U       (input) SuperMatrix*
 *         The factor U from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use column-wise storage scheme, i.e., U has types:
 *         Stype = SLU_NC, Dtype = SLU_D, Mtype = SLU_TRU.
 *
 * perm_c  (input) int*, dimension (L->ncol)
 *	   Column permutation vector, which defines the
 *         permutation matrix Pc; perm_c[i] = j means column i of A is
 *         in position j in A*Pc.
 *
 * perm_r  (input) int*, dimension (L->nrow)
 *         Row permutation vector, which defines the permutation matrix Pr;
 *         perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = SLU_DN, Dtype = SLU_D, Mtype = SLU_GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * stat     (output) SuperLUStat_t*
 *          Record the statistics on runtime and floating-point operation count.
 *          See util.h for the definition of 'SuperLUStat_t'.
 *
 * info    (output) int*
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 * </pre>
 */

void dgstrs(trans_t trans, SuperMatrix* L, SuperMatrix* U, const int* perm_c, const int* perm_r, SuperMatrix* B, SuperLUStat_t* stat, int* info) {
#ifdef _CRAY
    _fcd ftcs1, ftcs2, ftcs3, ftcs4;
#endif
#ifdef USE_VENDOR_BLAS
    double alpha = 1.0, beta = 1.0;
    double* work_col;
#endif
    DNformat* Bstore;
    double* Bmat;
    SCformat* Lstore;
    NCformat* Ustore;
    double *Lval, *Uval;
    int fsupc, nrow, nsupr, nsupc, irow;
    int_t i, j, k, luptr, istart, iptr;
    int jcol, n, ldb, nrhs;
    double *work, *rhs_work, *soln;
    flops_t solve_ops;
    void dprint_soln(int n, int nrhs, const double* soln);

    /* Test input parameters ... */
    *info = 0;
    Bstore = B->Store;
    ldb = Bstore->lda;
    nrhs = B->ncol;
    if(trans != NOTRANS && trans != TRANS && trans != CONJ) *info = -1;
    else if(L->nrow != L->ncol || L->nrow < 0 ||
            L->Stype != SLU_SC || L->Dtype != SLU_D || L->Mtype != SLU_TRLU)
        *info = -2;
    else if(U->nrow != U->ncol || U->nrow < 0 ||
            U->Stype != SLU_NC || U->Dtype != SLU_D || U->Mtype != SLU_TRU)
        *info = -3;
    else if(ldb < SUPERLU_MAX(0, L->nrow) ||
            B->Stype != SLU_DN || B->Dtype != SLU_D || B->Mtype != SLU_GE)
        *info = -6;
    if(*info) {
        int ii = -(*info);
        input_error("dgstrs", &ii);
        return;
    }

    n = L->nrow;
    work = doubleCalloc((size_t)n * (size_t)nrhs);
    if(!work) ABORT("Malloc fails for local work[].");
    soln = doubleMalloc((size_t)n);
    if(!soln) ABORT("Malloc fails for local soln[].");

    Bmat = Bstore->nzval;
    Lstore = L->Store;
    Lval = Lstore->nzval;
    Ustore = U->Store;
    Uval = Ustore->nzval;
    solve_ops = 0;

    if(trans == NOTRANS) {
        /* Permute right hand sides to form Pr*B */
        for(i = 0; i < nrhs; i++) {
            rhs_work = &Bmat[(size_t)i * (size_t)ldb];
            for(k = 0; k < n; k++) soln[perm_r[k]] = rhs_work[k];
            for(k = 0; k < n; k++) rhs_work[k] = soln[k];
        }

        /* Forward solve PLy=Pb. */
        for(k = 0; k <= Lstore->nsuper; k++) {
            fsupc = L_FST_SUPC(k);
            istart = L_SUB_START(fsupc);
            nsupr = L_SUB_START(fsupc + 1) - istart;
            nsupc = L_FST_SUPC(k + 1) - fsupc;
            nrow = nsupr - nsupc;

            solve_ops += nsupc * (nsupc - 1) * nrhs;
            solve_ops += 2 * nrow * nsupc * nrhs;

            if(nsupc == 1) {
                for(j = 0; j < nrhs; j++) {
                    rhs_work = &Bmat[(size_t)j * (size_t)ldb];
                    luptr = L_NZ_START(fsupc);
                    for(iptr = istart + 1; iptr < L_SUB_START(fsupc + 1); iptr++) {
                        irow = L_SUB(iptr);
                        ++luptr;
                        rhs_work[irow] -= rhs_work[fsupc] * Lval[luptr];
                    }
                }
            }
            else {
                luptr = L_NZ_START(fsupc);
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
                ftcs1 = _cptofcd("L", strlen("L"));
                ftcs2 = _cptofcd("N", strlen("N"));
                ftcs3 = _cptofcd("U", strlen("U"));
                STRSM(ftcs1, ftcs1, ftcs2, ftcs3, &nsupc, &nrhs, &alpha, &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);

                SGEMM(ftcs2, ftcs2, &nrow, &nrhs, &nsupc, &alpha, &Lval[luptr + nsupc], &nsupr, &Bmat[fsupc], &ldb, &beta, &work[0], &n);
#else
                dtrsm_("L", "L", "N", "U", &nsupc, &nrhs, &alpha, &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);

                dgemm_("N", "N", &nrow, &nrhs, &nsupc, &alpha, &Lval[luptr + nsupc], &nsupr, &Bmat[fsupc], &ldb, &beta, &work[0], &n);
#endif
                for(j = 0; j < nrhs; j++) {
                    rhs_work = &Bmat[(size_t)j * (size_t)ldb];
                    work_col = &work[(size_t)j * (size_t)n];
                    iptr = istart + nsupc;
                    for(i = 0; i < nrow; i++) {
                        irow = L_SUB(iptr);
                        rhs_work[irow] -= work_col[i]; /* Scatter */
                        work_col[i] = 0.0;
                        iptr++;
                    }
                }
#else
                for(j = 0; j < nrhs; j++) {
                    rhs_work = &Bmat[(size_t)j * (size_t)ldb];
                    dlsolve(nsupr, nsupc, &Lval[luptr], &rhs_work[fsupc]);
                    dmatvec(nsupr, nrow, nsupc, &Lval[luptr + nsupc], &rhs_work[fsupc], &work[0]);

                    iptr = istart + nsupc;
                    for(i = 0; i < nrow; i++) {
                        irow = L_SUB(iptr);
                        rhs_work[irow] -= work[i];
                        work[i] = 0.0;
                        iptr++;
                    }
                }
#endif
            } /* else ... */
        } /* for L-solve */

#if (DEBUGlevel >= 2)
        printf("After L-solve: y=\n");
        dprint_soln(n, nrhs, Bmat);
#endif

        /*
         * Back solve Ux=y.
         */
        for(k = Lstore->nsuper; k >= 0; k--) {
            fsupc = L_FST_SUPC(k);
            istart = L_SUB_START(fsupc);
            nsupr = L_SUB_START(fsupc + 1) - istart;
            nsupc = L_FST_SUPC(k + 1) - fsupc;
            luptr = L_NZ_START(fsupc);

            solve_ops += nsupc * (nsupc + 1) * nrhs;

            if(nsupc == 1) {
                rhs_work = &Bmat[0];
                for(j = 0; j < nrhs; j++) {
                    rhs_work[fsupc] /= Lval[luptr];
                    rhs_work += ldb;
                }
            }
            else {
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
                ftcs1 = _cptofcd("L", strlen("L"));
                ftcs2 = _cptofcd("U", strlen("U"));
                ftcs3 = _cptofcd("N", strlen("N"));
                STRSM(ftcs1, ftcs2, ftcs3, ftcs3, &nsupc, &nrhs, &alpha, &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#else
                dtrsm_("L", "U", "N", "N", &nsupc, &nrhs, &alpha, &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#endif
#else
                for(j = 0; j < nrhs; j++)
                    dusolve(nsupr, nsupc, &Lval[luptr], &Bmat[(size_t)fsupc + (size_t)j * (size_t)ldb]);
#endif
            }

            for(j = 0; j < nrhs; ++j) {
                rhs_work = &Bmat[(size_t)j * (size_t)ldb];
                for(jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
                    solve_ops += 2 * (U_NZ_START(jcol + 1) - U_NZ_START(jcol));
                    for(i = U_NZ_START(jcol); i < U_NZ_START(jcol + 1); i++) {
                        irow = U_SUB(i);
                        rhs_work[irow] -= rhs_work[jcol] * Uval[i];
                    }
                }
            }

        } /* for U-solve */

#if (DEBUGlevel >= 2)
        printf("After U-solve: x=\n");
        dprint_soln(n, nrhs, Bmat);
#endif

        /* Compute the final solution X := Pc*X. */
        for(i = 0; i < nrhs; i++) {
            rhs_work = &Bmat[(size_t)i * (size_t)ldb];
            for(k = 0; k < n; k++) soln[k] = rhs_work[perm_c[k]];
            for(k = 0; k < n; k++) rhs_work[k] = soln[k];
        }

        stat->ops[SOLVE] = solve_ops;
    }
    else { /* Solve A'*X=B or CONJ(A)*X=B */
        /* Permute right hand sides to form Pc'*B. */
        for(i = 0; i < nrhs; i++) {
            rhs_work = &Bmat[(size_t)i * (size_t)ldb];
            for(k = 0; k < n; k++) soln[perm_c[k]] = rhs_work[k];
            for(k = 0; k < n; k++) rhs_work[k] = soln[k];
        }

        stat->ops[SOLVE] = 0;
        for(k = 0; k < nrhs; ++k) {
            /* Multiply by inv(U'). */
            sp_dtrsv("U", "T", "N", L, U, &Bmat[(size_t)k * (size_t)ldb], stat, info);

            /* Multiply by inv(L'). */
            sp_dtrsv("L", "T", "U", L, U, &Bmat[(size_t)k * (size_t)ldb], stat, info);
        }
        /* Compute the final solution X := Pr'*X (=inv(Pr)*X) */
        for(i = 0; i < nrhs; i++) {
            rhs_work = &Bmat[(size_t)i * (size_t)ldb];
            for(k = 0; k < n; k++) soln[k] = rhs_work[perm_r[k]];
            for(k = 0; k < n; k++) rhs_work[k] = soln[k];
        }
    }

    SUPERLU_FREE(work);
    SUPERLU_FREE(soln);
}

/*
 * Diagnostic print of the solution vector
 */
void dprint_soln(int n, int nrhs, const double* soln) {
    int i;

    for(i = 0; i < n; i++)
        printf("\t%d: %.4f\n", i, soln[i]);
}
