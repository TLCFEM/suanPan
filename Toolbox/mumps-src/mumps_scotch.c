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
/* Interfacing with SCOTCH and pt-SCOTCH */
#include <stdio.h>
#include "mumps_scotch.h"
#if defined(scotch) || defined(ptscotch)
void MUMPS_CALL
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
                    MUMPS_INT * const  petab,
              const MUMPS_INT * const  pfree,
                    MUMPS_INT * const  lentab,
                    MUMPS_INT * const  iwtab,
                    MUMPS_INT * const  nvtab,
                    MUMPS_INT * const  elentab,
                    MUMPS_INT * const  lasttab,
                    MUMPS_INT * const  ncmpa )
{
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
}
#endif /* scotch */
#if defined(ptscotch)
void MUMPS_CALL
MUMPS_DGRAPHINIT(SCOTCH_Dgraph *graphptr, MPI_Fint *comm, MPI_Fint *ierr)
{
  MPI_Comm  int_comm;
  int_comm = MPI_Comm_f2c(*comm);
  *ierr = SCOTCH_dgraphInit(graphptr, int_comm);
  return;
}
#endif
