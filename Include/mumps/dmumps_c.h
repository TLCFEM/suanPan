#ifndef DMUMPS_C_H
#define DMUMPS_C_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mumps_compat.h"
#include "mumps_c_types.h"

#ifndef MUMPS_VERSION
#define MUMPS_VERSION "5.2.1"
#endif
#ifndef MUMPS_VERSION_MAX_LEN
#define MUMPS_VERSION_MAX_LEN 30
#endif

typedef struct {
	MUMPS_INT sym, par, job;
	MUMPS_INT comm_fortran;
	MUMPS_INT icntl[60];
	MUMPS_INT keep[500];
	DMUMPS_REAL cntl[15];
	DMUMPS_REAL dkeep[230];
	MUMPS_INT8 keep8[150];
	MUMPS_INT n;

	MUMPS_INT nz_alloc;

	MUMPS_INT nz;
	MUMPS_INT8 nnz;
	MUMPS_INT* irn;
	MUMPS_INT* jcn;
	DMUMPS_COMPLEX* a;

	MUMPS_INT nz_loc;
	MUMPS_INT8 nnz_loc;
	MUMPS_INT* irn_loc;
	MUMPS_INT* jcn_loc;
	DMUMPS_COMPLEX* a_loc;

	MUMPS_INT nelt;
	MUMPS_INT* eltptr;
	MUMPS_INT* eltvar;
	DMUMPS_COMPLEX* a_elt;

	MUMPS_INT* perm_in;

	MUMPS_INT* sym_perm;
	MUMPS_INT* uns_perm;

	DMUMPS_REAL* colsca;
	DMUMPS_REAL* rowsca;
	MUMPS_INT colsca_from_mumps;
	MUMPS_INT rowsca_from_mumps;

	DMUMPS_COMPLEX *rhs, *redrhs, *rhs_sparse, *sol_loc, *rhs_loc;
	MUMPS_INT *irhs_sparse, *irhs_ptr, *isol_loc, *irhs_loc;
	MUMPS_INT nrhs, lrhs, lredrhs, nz_rhs, lsol_loc, nloc_rhs, lrhs_loc;
	MUMPS_INT schur_mloc, schur_nloc, schur_lld;
	MUMPS_INT mblock, nblock, nprow, npcol;
	MUMPS_INT info[80], infog[80];
	DMUMPS_REAL rinfo[40], rinfog[40];

	MUMPS_INT deficiency;
	MUMPS_INT* pivnul_list;
	MUMPS_INT* mapping;

	MUMPS_INT size_schur;
	MUMPS_INT* listvar_schur;
	DMUMPS_COMPLEX* schur;

	MUMPS_INT instance_number;
	DMUMPS_COMPLEX* wk_user;

	char version_number[MUMPS_VERSION_MAX_LEN + 1 + 1];
	char ooc_tmpdir[256];
	char ooc_prefix[64];
	char write_problem[256];
	MUMPS_INT lwk_user;
	char save_dir[256];
	char save_prefix[256];

	MUMPS_INT metis_options[40];
} DMUMPS_STRUC_C;

void MUMPS_CALL dmumps_c(DMUMPS_STRUC_C* dmumps_par);

#ifdef __cplusplus
}
#endif

#endif
