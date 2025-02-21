/*
 *
 *  This file is part of MUMPS 5.7.3, released
 *  on Mon Jul 15 11:44:21 UTC 2024
 *
 *
 *  Copyright 1991-2024 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
/* libseq is currently compiled independently of MUMPS
 * and mumps_c_types.h and default integer size relies
 * on compilation with -DINTSIZE64
 */

#include <inttypes.h>
#ifdef INTSIZE64
#define LIBSEQ_INT int64_t
#define LIBSEQ_INT8 int64_t
#else
#define LIBSEQ_INT int
#define LIBSEQ_INT8 int64_t
#endif


#if ! defined(LIBSEQ_CALL)
#if defined(_WIN32) && ! defined(__MINGW32__)
/* Choose between next lines or modify according
 * to your Windows calling conventions:
   #define LIBSEQ_CALL
   #define LIBSEQ_CALL __declspec(dllexport)
   #define LIBSEQ_CALL __declspec(dllexport) */
#define LIBSEQ_CALL
#else
#define LIBSEQ_CALL
#endif
#endif


#ifndef MUMPS_MPI_H
#define MUMPS_MPI_H

/* We define all symbols as extern "C" for users who call MUMPS with its
   libseq from a C++ driver. */
#ifdef __cplusplus
extern "C" {
#endif

/* This is the minimum to have the C interface to MUMPS work with the
 * C example provided. Other stub functions of the MPI standard may be
 * added as needed. */

typedef LIBSEQ_INT MPI_Comm; /* Simple type for MPI communicator */
static MPI_Comm MPI_COMM_WORLD=(MPI_Comm)0;

int LIBSEQ_CALL MPI_Init(int *pargc, char ***pargv);
int LIBSEQ_CALL MPI_Comm_rank(MPI_Comm  comm, int  *rank);
int LIBSEQ_CALL MPI_Finalize(void);
MPI_Comm LIBSEQ_CALL MPI_Comm_f2c(LIBSEQ_INT comm);

/* For MPI_IS_IN_PLACE tests */
void LIBSEQ_CALL MUMPS_CHECKADDREQUAL(char *a, char*b, LIBSEQ_INT8 *i);
void LIBSEQ_CALL MUMPS_CHECKADDREQUAL_(char *a, char*b, LIBSEQ_INT8 *i);
void LIBSEQ_CALL mumps_checkaddrequal_(char *a, char*b, LIBSEQ_INT8 *i);
void LIBSEQ_CALL mumps_checkaddrequal__(char *a, char*b, LIBSEQ_INT8 *i);
double LIBSEQ_CALL MPI_Wtime(void);

#ifdef __cplusplus
}
#endif

#endif /* MUMPS_MPI_H */
