/*
 *
 *  This file is part of MUMPS 5.2.1, released
 *  on Fri Jun 14 14:46:05 UTC 2019
 *
 *
 *  Copyright 1991-2019 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */
#ifndef MUMPS_SCOTCH_H
#define MUMPS_SCOTCH_H
#include "mumps_common.h"
#if defined(scotch) || defined(ptscotch)
/* esmumps prototype with standard integers */
MUMPS_INT esmumps( const MUMPS_INT n, const MUMPS_INT iwlen, MUMPS_INT * const pe, const MUMPS_INT pfree,
             MUMPS_INT * const len, MUMPS_INT * const iw, MUMPS_INT * const nv, MUMPS_INT * const elen,
             MUMPS_INT * const last);
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
              MUMPS_INT * const        ncmpa );
#endif /*scotch or ptscotch*/
#if defined(ptscotch)
#include "mpi.h"
#include <stdio.h>
#include "ptscotch.h"
#define MUMPS_DGRAPHINIT \
  F_SYMBOL(dgraphinit,DGRAPHINIT)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr);
#endif /*ptscotch*/
#endif
