/*
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/
/*
 * -- SuperLU routine (version 4.1) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
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
 */
/*! \file
 * \brief Utility functions
 *
 * \ingroup Common
 */

#include <math.h>
#include "slu_ddefs.h"

/*! \brief Print message to error stream and exit program
 *
 * \param[in] mgs Message that is printed to error stream.
 */
void superlu_abort_and_exit(const char* msg) {
    fprintf(stderr, "%s", msg);
    exit(-1);
}

/*! \brief Set the default values for the options argument.
 *
 * \param[in,out] options Options struct that is filled with default values.
 */
void set_default_options(superlu_options_t* options) {
    options->Fact = DOFACT;
    options->Equil = YES;
    options->ColPerm = COLAMD;
    options->Trans = NOTRANS;
    options->IterRefine = NOREFINE;
    options->DiagPivotThresh = 1.0;
    options->SymmetricMode = NO;
    options->PivotGrowth = NO;
    options->ConditionNumber = NO;
    options->PrintStat = YES;
}

/*! \brief Set the default values for the options argument for ILU.
 *
 * \param[out] options Options struct that is filled with default values.
 */
void ilu_set_default_options(superlu_options_t* options) {
    set_default_options(options);

    /* further options for incomplete factorization */
    options->DiagPivotThresh = 0.1;
    options->RowPerm = LargeDiag_MC64;
    options->ILU_DropRule = DROP_BASIC | DROP_AREA;
    options->ILU_DropTol = 1e-4;
    options->ILU_FillFactor = 10.0;
    options->ILU_Norm = INF_NORM;
    options->ILU_MILU = SILU;
    options->ILU_MILU_Dim = 3.0; /* -log(n)/log(h) is perfect */
    options->ILU_FillTol = 1e-2;
}

/*! \brief Print the options setting.
 *
 * \param[in] options Options struct that is printed.
 */
void print_options(const superlu_options_t* options) {
    printf(".. options:\n");
    printf("\tFact\t %8d\n", options->Fact);
    printf("\tEquil\t %8d\n", options->Equil);
    printf("\tColPerm\t %8d\n", options->ColPerm);
    printf("\tDiagPivotThresh %8.4f\n", options->DiagPivotThresh);
    printf("\tTrans\t %8d\n", options->Trans);
    printf("\tIterRefine\t%4d\n", options->IterRefine);
    printf("\tSymmetricMode\t%4d\n", options->SymmetricMode);
    printf("\tPivotGrowth\t%4d\n", options->PivotGrowth);
    printf("\tConditionNumber\t%4d\n", options->ConditionNumber);
    printf("..\n");
}

/*! \brief Print the options setting.
 *
 * \param[in] options Options struct that is printed.
 */
void print_ilu_options(const superlu_options_t* options) {
    printf(".. ILU options:\n");
    printf("\tDiagPivotThresh\t%6.2e\n", options->DiagPivotThresh);
    printf("\ttau\t%6.2e\n", options->ILU_DropTol);
    printf("\tgamma\t%6.2f\n", options->ILU_FillFactor);
    printf("\tDropRule\t%0x\n", options->ILU_DropRule);
    printf("\tMILU\t%d\n", options->ILU_MILU);
    printf("\tMILU_ALPHA\t%6.2e\n", MILU_ALPHA);
    printf("\tDiagFillTol\t%6.2e\n", options->ILU_FillTol);
    printf("..\n");
}

/*! \brief Deallocate SuperMatrix
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix.
 */
void Destroy_SuperMatrix_Store(SuperMatrix* A) {
    SUPERLU_FREE(A->Store);
}

/*! \brief Deallocate SuperMatrix of type NC
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix of type NC.
 */
void Destroy_CompCol_Matrix(SuperMatrix* A) {
    SUPERLU_FREE(((NCformat*)A->Store)->rowind);
    SUPERLU_FREE(((NCformat*)A->Store)->colptr);
    SUPERLU_FREE(((NCformat*)A->Store)->nzval);
    SUPERLU_FREE(A->Store);
}

/*! \brief Deallocate SuperMatrix of type NR
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix of type NR.
 */
void Destroy_CompRow_Matrix(SuperMatrix* A) {
    SUPERLU_FREE(((NRformat*)A->Store)->colind);
    SUPERLU_FREE(((NRformat*)A->Store)->rowptr);
    SUPERLU_FREE(((NRformat*)A->Store)->nzval);
    SUPERLU_FREE(A->Store);
}

/*! \brief Deallocate SuperMatrix of type SC
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix of type SC.
 */
void Destroy_SuperNode_Matrix(SuperMatrix* A) {
    SUPERLU_FREE(((SCformat*)A->Store)->rowind);
    SUPERLU_FREE(((SCformat*)A->Store)->rowind_colptr);
    SUPERLU_FREE(((SCformat*)A->Store)->nzval);
    SUPERLU_FREE(((SCformat*)A->Store)->nzval_colptr);
    SUPERLU_FREE(((SCformat*)A->Store)->col_to_sup);
    SUPERLU_FREE(((SCformat*)A->Store)->sup_to_col);
    SUPERLU_FREE(A->Store);
}

/*! \brief Deallocate SuperMatrix of type NCP
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix of type NCP.
 */
void Destroy_CompCol_Permuted(SuperMatrix* A) {
    SUPERLU_FREE(((NCPformat*)A->Store)->colbeg);
    SUPERLU_FREE(((NCPformat*)A->Store)->colend);
    SUPERLU_FREE(A->Store);
}

/*! \brief Deallocate SuperMatrix of type DN
 *
 * Deallocate the structure pointing to the actual storage of the matrix.
 *
 * \param[in,out] A Deallocate all memory of this SuperMatrix of type DN.
 */
void Destroy_Dense_Matrix(SuperMatrix* A) {
    DNformat* Astore = A->Store;
    SUPERLU_FREE(Astore->nzval);
    SUPERLU_FREE(A->Store);
}

/*! \brief Reset repfnz[] for the current column
 */
void resetrep_col(const int nseg, const int* segrep, int* repfnz) {
    int_t i, irep;

    for(i = 0; i < nseg; i++) {
        irep = segrep[i];
        repfnz[irep] = SLU_EMPTY;
    }
}

/*! \brief Count the total number of nonzeros in factors L and U, and in the symmetrically reduced L.
 */
void countnz(const int n, int_t* xprune, int_t* nnzL, int_t* nnzU, GlobalLU_t* Glu) {
    int nsuper, fsupc, i, j;
    int_t jlen;
#if (DEBUGlevel >= 1)
    int_t irep = 0, nnzL0 = 0;
#endif
    int* xsup;
    int_t* xlsub;

    xsup = Glu->xsup;
    xlsub = Glu->xlsub;
    *nnzL = 0;
    *nnzU = (Glu->xusub)[n];
    nsuper = (Glu->supno)[n];

    if(n <= 0) return;

    /*
     * For each supernode
     */
    for(i = 0; i <= nsuper; i++) {
        fsupc = xsup[i];
        jlen = xlsub[fsupc + 1] - xlsub[fsupc];

        for(j = fsupc; j < xsup[i + 1]; j++) {
            *nnzL += jlen;
            *nnzU += j - fsupc + 1;
            jlen--;
        }
#if (DEBUGlevel >= 1)
        irep = xsup[i + 1] - 1;
        nnzL0 += xprune[irep] - xlsub[irep];
#endif
    }

#if (DEBUGlevel >= 1)
    printf("\tNo of nonzeros in symm-reduced L = %lld\n", (long long)nnzL0);
    fflush(stdout);
#endif
}

/*! \brief Count the total number of nonzeros in factors L and U.
 */
void ilu_countnz(const int n, int_t* nnzL, int_t* nnzU, GlobalLU_t* Glu) {
    int nsuper, fsupc, i, j;
    int jlen;
    int* xsup;
    int_t* xlsub;

    xsup = Glu->xsup;
    xlsub = Glu->xlsub;
    *nnzL = 0;
    *nnzU = (Glu->xusub)[n];
    nsuper = (Glu->supno)[n];

    if(n <= 0) return;

    /*
     * For each supernode
     */
    for(i = 0; i <= nsuper; i++) {
        fsupc = xsup[i];
        jlen = xlsub[fsupc + 1] - xlsub[fsupc];

        for(j = fsupc; j < xsup[i + 1]; j++) {
            *nnzL += jlen;
            *nnzU += j - fsupc + 1;
            jlen--;
        }
        // irep = xsup[i+1] - 1;
    }
}

/*! \brief Fix up the data storage lsub for L-subscripts.
 *
 * It removes the subscript sets for structural pruning,
 * and applies permutation to the remaining subscripts.
 */
void fixupL(const int n, const int* perm_r, GlobalLU_t* Glu) {
    int nsuper, fsupc, i, k;
    int_t nextl, j, jstrt;
    int* xsup;
    int_t *lsub, *xlsub;

    if(n <= 1) return;

    xsup = Glu->xsup;
    lsub = Glu->lsub;
    xlsub = Glu->xlsub;
    nextl = 0;
    nsuper = (Glu->supno)[n];

    /*
     * For each supernode ...
     */
    for(i = 0; i <= nsuper; i++) {
        fsupc = xsup[i];
        jstrt = xlsub[fsupc];
        xlsub[fsupc] = nextl;
        for(j = jstrt; j < xlsub[fsupc + 1]; j++) {
            lsub[nextl] = perm_r[lsub[j]]; /* Now indexed into P*A */
            nextl++;
        }
        for(k = fsupc + 1; k < xsup[i + 1]; k++)
            xlsub[k] = nextl; /* Other columns in supernode i */
    }

    xlsub[n] = nextl;
}

/*! \brief Diagnostic print of segment info after panel_dfs().
 */
void print_panel_seg(int n, int w, int jcol, int nseg, const int* segrep, const int* repfnz) {
    int j, k;

    for(j = jcol; j < jcol + w; j++) {
        printf("\tcol %d:\n", j);
        for(k = 0; k < nseg; k++)
            printf("\t\tseg %d, segrep %d, repfnz %d\n", k, segrep[k], repfnz[(j - jcol) * n + segrep[k]]);
    }
}

/*! \brief Initialize SuperLU stat
 *
 * \param[in,out] stat SuperLU stat that is initialized.
 */
void StatInit(SuperLUStat_t* stat) {
    register int i, w, panel_size, relax;

    panel_size = sp_ienv(1);
    relax = sp_ienv(2);
    w = SUPERLU_MAX(panel_size, relax);
    stat->panel_histo = int32Calloc(w + 1);
    stat->utime = (double*)SUPERLU_MALLOC(NPHASES * sizeof(double));
    if(!stat->utime) ABORT("SUPERLU_MALLOC fails for stat->utime");
    stat->ops = (flops_t*)SUPERLU_MALLOC(NPHASES * sizeof(flops_t));
    if(!stat->ops) ABORT("SUPERLU_MALLOC fails for stat->ops");
    for(i = 0; i < NPHASES; ++i) {
        stat->utime[i] = 0.;
        stat->ops[i] = 0.;
    }
    stat->TinyPivots = 0;
    stat->RefineSteps = 0;
    stat->expansions = 0;
#if (PRNTlevel >= 1)
    printf(".. parameters in sp_ienv():\n");
    printf("\t 1: panel size \t %4d \n"
           "\t 2: relax      \t %4d \n"
           "\t 3: max. super \t %4d \n"
           "\t 4: row-dim 2D \t %4d \n"
           "\t 5: col-dim 2D \t %4d \n"
           "\t 6: fill ratio \t %4d \n",
           sp_ienv(1), sp_ienv(2), sp_ienv(3), sp_ienv(4), sp_ienv(5), sp_ienv(6));
#endif
}

/*! \brief Display SuperLU stat
 *
 * Print content of SuperLU stat to output.
 *
 * \param[in] stat Display this SuperLU stat
 */
void StatPrint(SuperLUStat_t* stat) {
    double* utime;
    flops_t* ops;

    utime = stat->utime;
    ops = stat->ops;
    printf("Factor time  = %8.5f\n", utime[FACT]);
    if(utime[FACT] != 0.0)
        printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT], ops[FACT] * 1e-6 / utime[FACT]);

    printf("Solve time   = %8.4f\n", utime[SOLVE]);
    if(utime[SOLVE] != 0.0)
        printf("Solve flops = %e\tMflops = %8.2f\n", ops[SOLVE], ops[SOLVE] * 1e-6 / utime[SOLVE]);

    printf("Number of memory expansions: %d\n", stat->expansions);
}

/*! \brief Deallocate SuperLU stat
 *
 * Deallocate the structure pointing to the actual storage of SuperLU stat.
 *
 * \param[in,out] stat Deallocate all memory of this SuperLU stat
 */
void StatFree(SuperLUStat_t* stat) {
    SUPERLU_FREE(stat->panel_histo);
    SUPERLU_FREE(stat->utime);
    SUPERLU_FREE(stat->ops);
}

/*! \brief Get operations for LU factorization
 *
 * Read out number of operations (ops) needed for LU factorization.
 *
 * \param[in] stat SuperLU stat used to read out the opts.
 *
 * \return Number of operations needed for LU factorization.
 */
flops_t
LUFactFlops(SuperLUStat_t* stat) {
    return (stat->ops[FACT]);
}

/*! \brief Get operations for LU solve
 *
 * Read out number of operations (ops) needed for LU solve.
 *
 * \param[in] stat SuperLU stat used to read out the opts.
 *
 * \return Number of operations needed for LU solve.
 */
flops_t
LUSolveFlops(SuperLUStat_t* stat) {
    return (stat->ops[SOLVE]);
}

/*! \brief Fills an integer array with a given value.
 *
 * \param[in,out] a Integer array that is filled.
 * \param[in]     alen Length of integer array \a a.
 * \param[in]     ival Value to be filled in every element of \a a.
 */
void ifill(int* a, int alen, int ival) {
    register int i;
    for(i = 0; i < alen; i++) a[i] = ival;
}

/*! \brief Get the statistics of the supernodes
 */
#define NBUCKS 10

void super_stats(int nsuper, const int* xsup) {
    register int nsup1 = 0;
    int i, isize, whichb, bl, bh;
    int bucket[NBUCKS];
    int max_sup_size = 0;

    for(i = 0; i <= nsuper; i++) {
        isize = xsup[i + 1] - xsup[i];
        if(isize == 1) nsup1++;
        if(max_sup_size < isize) max_sup_size = isize;
    }

    printf("    Supernode statistics:\n\tno of super = %d\n", nsuper + 1);
    printf("\tmax supernode size = %d\n", max_sup_size);
    printf("\tno of size 1 supernodes = %d\n", nsup1);

    /* Histogram of the supernode sizes */
    ifill(bucket, NBUCKS, 0);

    for(i = 0; i <= nsuper; i++) {
        isize = xsup[i + 1] - xsup[i];
        whichb = (float)isize / max_sup_size * NBUCKS;
        if(whichb >= NBUCKS) whichb = NBUCKS - 1;
        bucket[whichb]++;
    }

    printf("\tHistogram of supernode sizes:\n");
    for(i = 0; i < NBUCKS; i++) {
        bl = (float)i * max_sup_size / NBUCKS;
        bh = (float)(i + 1) * max_sup_size / NBUCKS;
        printf("\tsnode: %d-%d\t\t%d\n", bl + 1, bh, bucket[i]);
    }
}

float SpaSize(int n, int np, float sum_npw) {
    return (sum_npw * 8 + np * 8 + n * 4) / 1024.;
}

float DenseSize(int n, float sum_nw) {
    return (sum_nw * 8 + n * 8) / 1024.;
}

/*! \brief Check whether repfnz[] == SLU_EMPTY after reset.
 */
void check_repfnz(int n, int w, int jcol, const int* repfnz) {
    int jj, k;

    for(jj = jcol; jj < jcol + w; jj++)
        for(k = 0; k < n; k++)
            if(repfnz[(jj - jcol) * n + k] != SLU_EMPTY) {
                fprintf(stderr, "col %d, repfnz_col[%d] = %d\n", jj, k, repfnz[(jj - jcol) * n + k]);
                ABORT("check_repfnz");
            }
}

/*! \brief Print a summary of the testing results.
 *
 * \param[in] type Array with three chars indicating the type for that the tests were run like "CGE" or "DGE".
 * \param[in] nfail Number of failed tests.
 * \param[in] nrun Number of tests run.
 * \param[in] nerrs Number of error messages recorded.
 */
void PrintSumm(const char* type, int nfail, int nrun, int nerrs) {
    if(nfail > 0)
        printf("%3s driver: %d out of %d tests failed to pass the threshold\n", type, nfail, nrun);
    else
        printf("All tests for %3s driver passed the threshold (%6d tests run)\n", type, nrun);

    if(nerrs > 0)
        printf("%6d error messages recorded\n", nerrs);
}

/*! \brief Print content of int array
 *
 * \param[in] what Vector name that is printed first.
 * \param[in] n Number of elements in array.
 * \param[in] vec Array of ints to be printed
 */
void print_int_vec(const char* what, int n, const int* vec) {
    int i;
    printf("%s\n", what);
    for(i = 0; i < n; ++i) printf("%d\t%d\n", i, vec[i]);
}

/*! \brief Print content of int array with index numbers after every tenth row
 *
 * Print all elements of an int array. After ten rows the index is printed to
 * make it more readable for humans.
 *
 * \param[in] name Vector name that is printed first.
 * \param[in] len Number of elements in array.
 * \param[in] x Array of ints to be printed.
 */
void slu_PrintInt10(const char* name, int len, const int* x) {
    register int i;

    printf("%10s:", name);
    for(i = 0; i < len; ++i) {
        if(i % 10 == 0) printf("\n\t[%2d-%2d]", i, i + 9);
        printf("%6d", x[i]);
    }
    printf("\n");
}

/*! \brief Validity check of a permutation
 *
 * \param[in] what String to be printed as part of displayed text for success or error.
 * \param[in] n Number of elements in permutation \a perm.
 * \param[in] perm Array describing the permutation.
 */
void check_perm(const char* what, int n, const int* perm) {
    register int i;
    int* marker;
    /*marker = (int *) calloc(n, sizeof(int));*/
    marker = int32Malloc(n);
    for(i = 0; i < n; ++i) marker[i] = 0;

    for(i = 0; i < n; ++i) {
        if(marker[perm[i]] == 1 || perm[i] >= n) {
            printf("%s: Not a valid PERM[%d] = %d\n", what, i, perm[i]);
            ABORT("check_perm");
        }
        else {
            marker[perm[i]] = 1;
        }
    }

    SUPERLU_FREE(marker);
    printf("check_perm: %s: n %d\n", what, n);
}
