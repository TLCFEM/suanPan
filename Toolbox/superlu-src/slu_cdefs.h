/*! \file
Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from U.S. Dept. of Energy)

All rights reserved.

The source code is distributed under BSD license, see the file License.txt
at the top-level directory.
*/

/*! @file slu_cdefs.h
 * \brief Header file for real operations
 *
 * <pre>
 * -- SuperLU routine (version 7.0.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November, 2010
 * August 2024
 *
 * Global data structures used in LU factorization -
 *
 *   nsuper: \#supernodes = nsuper + 1, numbered [0, nsuper].
 *   (xsup,supno): supno[i] is the supernode no to which i belongs;
 *	xsup(s) points to the beginning of the s-th supernode.
 *	e.g.   supno 0 1 2 2 3 3 3 4 4 4 4 4   (n=12)
 *	        xsup 0 1 2 4 7 12
 *	Note: dfs will be performed on supernode rep. relative to the new
 *	      row pivoting ordering
 *
 *   (xlsub,lsub): lsub[*] contains the compressed subscript of
 *	rectangular supernodes; xlsub[j] points to the starting
 *	location of the j-th column in lsub[*]. Note that xlsub
 *	is indexed by column.
 *	Storage: original row subscripts
 *
 *      During the course of sparse LU factorization, we also use
 *	(xlsub,lsub) for the purpose of symmetric pruning. For each
 *	supernode {s,s+1,...,t=s+r} with first column s and last
 *	column t, the subscript set
 *		lsub[j], j=xlsub[s], .., xlsub[s+1]-1
 *	is the structure of column s (i.e. structure of this supernode).
 *	It is used for the storage of numerical values.
 *	Furthermore,
 *		lsub[j], j=xlsub[t], .., xlsub[t+1]-1
 *	is the structure of the last column t of this supernode.
 *	It is for the purpose of symmetric pruning. Therefore, the
 *	structural subscripts can be rearranged without making physical
 *	interchanges among the numerical values.
 *
 *	However, if the supernode has only one column, then we
 *	only keep one set of subscripts. For any subscript interchange
 *	performed, similar interchange must be done on the numerical
 *	values.
 *
 *	The last column structures (for pruning) will be removed
 *	after the numerical LU factorization phase.
 *
 *   (xlusup,lusup): lusup[*] contains the numerical values of the
 *	rectangular supernodes; xlusup[j] points to the starting
 *	location of the j-th column in storage vector lusup[*]
 *	Note: xlusup is indexed by column.
 *	Each rectangular supernode is stored by column-major
 *	scheme, consistent with Fortran 2-dim array storage.
 *
 *   (xusub,ucol,usub): ucol[*] stores the numerical values of
 *	U-columns outside the rectangular supernodes. The row
 *	subscript of nonzero ucol[k] is stored in usub[k].
 *	xusub[i] points to the starting location of column i in ucol.
 *	Storage: new row subscripts; that is subscripts of PA.
 * </pre>
 */
#ifndef __SUPERLU_cSP_DEFS /* allow multiple inclusions */
#define __SUPERLU_cSP_DEFS

/*
 * File name:		csp_defs.h
 * Purpose:             Sparse matrix types and function prototypes
 * History:
 */

#ifdef _CRAY
#include <fortran.h>
#endif

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "slu_Cnames.h"
#include "superlu_config.h"
#include "supermatrix.h"
#include "slu_util.h"
#include "slu_scomplex.h"

/* -------- Prototypes -------- */

#ifdef __cplusplus
extern "C" {
#endif

/*! \brief Driver routines */
extern void
cgssv(superlu_options_t*, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int_t* info);
extern void
cgssvx(superlu_options_t*, SuperMatrix*, int*, int*, int*, char*, float*, float*, SuperMatrix*, SuperMatrix*, void*, int_t lwork, SuperMatrix*, SuperMatrix*, float*, float*, float*, float*, GlobalLU_t*, mem_usage_t*, SuperLUStat_t*, int_t* info);
/* ILU */
extern void
cgsisv(superlu_options_t*, SuperMatrix*, int*, int*, SuperMatrix*, SuperMatrix*, SuperMatrix*, SuperLUStat_t*, int*);
extern void
cgsisx(superlu_options_t* options, SuperMatrix* A, int* perm_c, int* perm_r, int* etree, char* equed, float* R, float* C, SuperMatrix* L, SuperMatrix* U, void* work, int_t lwork, SuperMatrix* B, SuperMatrix* X, float* recip_pivot_growth, float* rcond, GlobalLU_t* Glu, mem_usage_t* mem_usage, SuperLUStat_t* stat, int_t* info);

/*! \brief Supernodal LU factor related */
extern void
cCreate_CompCol_Matrix(SuperMatrix*, int, int, int_t, singlecomplex*, int_t*, int_t*, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_CompRow_Matrix(SuperMatrix*, int, int, int_t, singlecomplex*, int_t*, int_t*, Stype_t, Dtype_t, Mtype_t);
extern void cCompRow_to_CompCol(int, int, int_t, singlecomplex*, int_t*, int_t*, singlecomplex**, int_t**, int_t**);
extern void
cCopy_CompCol_Matrix(SuperMatrix*, SuperMatrix*);
extern void
cCreate_Dense_Matrix(SuperMatrix*, int, int, singlecomplex*, int, Stype_t, Dtype_t, Mtype_t);
extern void
cCreate_SuperNode_Matrix(SuperMatrix*, int, int, int_t, singlecomplex*, int_t*, int_t*, int_t*, int*, int*, Stype_t, Dtype_t, Mtype_t);
extern void
cCopy_Dense_Matrix(int, int, singlecomplex*, int, singlecomplex*, int);

extern void callocateA(int, int_t, singlecomplex**, int_t**, int_t**);
extern void cgstrf(superlu_options_t*, SuperMatrix*, int, int, int*, void*, int_t, int*, int*, SuperMatrix*, SuperMatrix*, GlobalLU_t*, SuperLUStat_t*, int_t* info);
extern int_t csnode_dfs(const int, const int, const int_t*, const int_t*, const int_t*, int_t*, int*, GlobalLU_t*);
extern int csnode_bmod(const int, const int, const int, singlecomplex*, singlecomplex*, GlobalLU_t*, SuperLUStat_t*);
extern void cpanel_dfs(const int, const int, const int, SuperMatrix*, int*, int*, singlecomplex*, int*, int*, int*, int_t*, int*, int*, int_t*, GlobalLU_t*);
extern void cpanel_bmod(const int, const int, const int, const int, singlecomplex*, singlecomplex*, int*, int*, GlobalLU_t*, SuperLUStat_t*);
extern int ccolumn_dfs(const int, const int, int*, int*, int*, int*, int*, int_t*, int*, int*, int_t*, GlobalLU_t*);
extern int ccolumn_bmod(const int, const int, singlecomplex*, singlecomplex*, int*, int*, int, GlobalLU_t*, SuperLUStat_t*);
extern int ccopy_to_ucol(int, int, int*, int*, int*, singlecomplex*, GlobalLU_t*);
extern int cpivotL(const int, const double, int*, int*, int*, int*, int*, GlobalLU_t*, SuperLUStat_t*);
extern void cpruneL(const int, const int*, const int, const int, const int*, const int*, int_t*, GlobalLU_t*);
extern void creadmt(int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void cGenXtrue(int, int, singlecomplex*, int);
extern void cFillRHS(trans_t, int, singlecomplex*, int, SuperMatrix*, SuperMatrix*);
extern void cgstrs(trans_t, SuperMatrix*, SuperMatrix*, const int*, const int*, SuperMatrix*, SuperLUStat_t*, int*);
/* ILU */
extern void cgsitrf(superlu_options_t*, SuperMatrix*, int, int, int*, void*, int_t, int*, int*, SuperMatrix*, SuperMatrix*, GlobalLU_t*, SuperLUStat_t*, int_t* info);
extern int cldperm(int, int, int_t, int_t[], int_t[], singlecomplex[], int[], float[], float[]);
extern int ilu_csnode_dfs(const int, const int, const int_t*, const int_t*, const int_t*, int*, GlobalLU_t*);
extern void ilu_cpanel_dfs(const int, const int, const int, SuperMatrix*, int*, int*, singlecomplex*, float*, int*, int*, int*, int*, int*, int_t*, GlobalLU_t*);
extern int ilu_ccolumn_dfs(const int, const int, int*, int*, int*, int*, int*, int*, int*, int_t*, GlobalLU_t*);
extern int ilu_ccopy_to_ucol(int, int, int*, int*, int*, singlecomplex*, int, milu_t, double, int, singlecomplex*, int*, GlobalLU_t*, float*);
extern int ilu_cpivotL(const int, const double, int*, int*, int, int*, int*, int*, int*, double, milu_t, singlecomplex, GlobalLU_t*, SuperLUStat_t*);
extern int ilu_cdrop_row(superlu_options_t*, int, int, double, int, int*, double*, GlobalLU_t*, float*, float*, int);

/*! \brief Driver related */

extern void cgsequ(SuperMatrix*, float*, float*, float*, float*, float*, int*);
extern void claqgs(SuperMatrix*, float*, float*, float, float, float, char*);
extern void cgscon(char*, SuperMatrix*, SuperMatrix*, float, float*, SuperLUStat_t*, int*);
extern float cPivotGrowth(int, SuperMatrix*, int*, SuperMatrix*, SuperMatrix*);
extern void cgsrfs(trans_t, SuperMatrix*, SuperMatrix*, SuperMatrix*, int*, int*, char*, float*, float*, SuperMatrix*, SuperMatrix*, float*, float*, SuperLUStat_t*, int*);

extern int sp_ctrsv(char*, char*, char*, SuperMatrix*, SuperMatrix*, singlecomplex*, SuperLUStat_t*, int*);
extern int sp_cgemv(char*, singlecomplex, SuperMatrix*, singlecomplex*, int, singlecomplex, singlecomplex*, int);

extern int sp_cgemm(char*, char*, int, int, int, singlecomplex, SuperMatrix*, singlecomplex*, int, singlecomplex, singlecomplex*, int);
extern float smach(char*); /* from C99 standard, in float.h */

/*! \brief Memory-related */
extern int_t cLUMemInit(fact_t, void*, int_t, int, int, int_t, int, float, SuperMatrix*, SuperMatrix*, GlobalLU_t*, int**, singlecomplex**);
extern void cSetRWork(int, int, singlecomplex*, singlecomplex**, singlecomplex**);
extern void cLUWorkFree(int*, singlecomplex*, GlobalLU_t*);
extern int_t cLUMemXpand(int, int_t, MemType, int_t*, GlobalLU_t*);

extern singlecomplex* singlecomplexMalloc(size_t);
extern singlecomplex* singlecomplexCalloc(size_t);
extern float* floatMalloc(size_t);
extern float* floatCalloc(size_t);
extern int_t cmemory_usage(const int_t, const int_t, const int_t, const int);
extern int cQuerySpace(SuperMatrix*, SuperMatrix*, mem_usage_t*);
extern int ilu_cQuerySpace(SuperMatrix*, SuperMatrix*, mem_usage_t*);

/*! \brief Auxiliary routines */
extern void creadhb(FILE*, int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void creadrb(int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void creadtriple(int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void creadtriple_noheader(int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void creadMM(FILE*, int*, int*, int_t*, singlecomplex**, int_t**, int_t**);
extern void cfill(singlecomplex*, int, singlecomplex);
extern void cinf_norm_error(int, SuperMatrix*, singlecomplex*);
extern float sqselect(int, float*, int);

/*! \brief Routines for debugging */
extern void cPrint_CompCol_Matrix(char*, SuperMatrix*);
extern void cPrint_SuperNode_Matrix(char*, SuperMatrix*);
extern void cPrint_Dense_Matrix(char*, SuperMatrix*);
extern void cprint_lu_col(char*, int, int, int_t*, GlobalLU_t*);
extern int print_singlecomplex_vec(const char*, int, const singlecomplex*);
extern void ccheck_tempv(int, singlecomplex*);

/*! \brief BLAS */

extern void ccopy_(int*, singlecomplex*, int*, singlecomplex*, int*);
extern void caxpy_(int*, singlecomplex*, singlecomplex*, int*, singlecomplex*, int*);
extern void cgemm_(const char*, const char*, const int*, const int*, const int*, const singlecomplex*, const singlecomplex*, const int*, const singlecomplex*, const int*, const singlecomplex*, singlecomplex*, const int*);
extern void ctrsv_(char*, char*, char*, int*, singlecomplex*, int*, singlecomplex*, int*);
extern void ctrsm_(char*, char*, char*, char*, int*, int*, singlecomplex*, singlecomplex*, int*, singlecomplex*, int*);
extern void cgemv_(char*, int*, int*, singlecomplex*, singlecomplex* a, int*, singlecomplex*, int*, singlecomplex*, singlecomplex*, int*);

extern void cusolve(int, int, singlecomplex*, singlecomplex*);
extern void clsolve(int, int, singlecomplex*, singlecomplex*);
extern void cmatvec(int, int, int, singlecomplex*, singlecomplex*, singlecomplex*);

#ifdef __cplusplus
}
#endif

#endif /* __SUPERLU_cSP_DEFS */
