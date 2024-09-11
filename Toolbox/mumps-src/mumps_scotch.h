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
#ifndef MUMPS_SCOTCH_H
#define MUMPS_SCOTCH_H
#include "mumps_common.h"
#if defined(scotch) || defined(ptscotch)
#include "scotch.h"
#if ((SCOTCH_VERSION == 6) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 7)
/* esmumpsv prototype with 64-bit integers weights of nodes in the graph are used on entry (nv) */
MUMPS_INT esmumpsv( const MUMPS_INT n, const MUMPS_INT iwlen, MUMPS_INT * const pe, const MUMPS_INT pfree,
             MUMPS_INT * const len, MUMPS_INT * const iw, MUMPS_INT * const nv, MUMPS_INT * const elen,
             MUMPS_INT * const last);
#endif
/* esmumps prototype with standard integers (weights of nodes not used on entry) */
MUMPS_INT esmumps( const MUMPS_INT n, const MUMPS_INT iwlen, MUMPS_INT * const pe, const MUMPS_INT pfree,
             MUMPS_INT * const len, MUMPS_INT * const iw, MUMPS_INT * const nv, MUMPS_INT * const elen,
             MUMPS_INT * const last);
#define  MUMPS_SCOTCH_WEIGHTUSED       \
    F_SYMBOL(scotch_weightused,SCOTCH_WEIGHTUSED)
void MUMPS_CALL
MUMPS_SCOTCH_WEIGHTUSED( MUMPS_INT  * const  weightused );
#define MUMPS_SCOTCH_ORD        \
    F_SYMBOL(scotch_ord,SCOTCH_ORD)
void MUMPS_CALL
MUMPS_SCOTCH_ORD( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
              MUMPS_INT * const        petab,
              const MUMPS_INT * const  pfree,
              MUMPS_INT * const        lentab,
              MUMPS_INT * const        iwtab,
              MUMPS_INT * const        nvtab,
              MUMPS_INT * const        elentab,
              MUMPS_INT * const        lasttab,
              MUMPS_INT * const        ncmpa, 
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
              SCOTCH_Context * const contextptr,
#endif
              MUMPS_INT * const        weightused,
              MUMPS_INT * const        weightrequested );
#define MUMPS_SCOTCH        \
    F_SYMBOL(scotch,SCOTCH)
void MUMPS_CALL
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
              MUMPS_INT * const        petab,
              const MUMPS_INT * const  pfree,
              MUMPS_INT * const        lentab,
              MUMPS_INT * const        iwtab,
              MUMPS_INT * const        nvtab,
              MUMPS_INT * const        elentab,
              MUMPS_INT * const        lasttab,
              MUMPS_INT * const        ncmpa, 
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
              SCOTCH_Context * const contextptr,
#endif
              MUMPS_INT * const        weightused,
              MUMPS_INT * const        weightrequested );
#define MUMPS_SCOTCH_VERSION        \
    F_SYMBOL(scotch_version,SCOTCH_VERSION)
void MUMPS_CALL
MUMPS_SCOTCH_VERSION(MUMPS_INT *version);
#define MUMPS_SCOTCH_GET_PTHREAD_NUMBER        \
    F_SYMBOL(scotch_get_pthread_number,SCOTCH_GET_PTHREAD_NUMBER)
void MUMPS_CALL
MUMPS_SCOTCH_GET_PTHREAD_NUMBER (MUMPS_INT *PTHREAD_NUMBER);
#define MUMPS_SCOTCH_SET_PTHREAD_NUMBER        \
    F_SYMBOL(scotch_set_pthread_number,SCOTCH_SET_PTHREAD_NUMBER)
void MUMPS_CALL
MUMPS_SCOTCH_SET_PTHREAD_NUMBER (MUMPS_INT *PTHREAD_NUMBER);
#endif /*scotch or ptscotch*/
#if defined(ptscotch)
#include "mpi.h"
#include "ptscotch.h"
#define MUMPS_DGRAPHINIT \
  F_SYMBOL(dgraphinit,DGRAPHINIT)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr);
#endif /*ptscotch*/
#endif
