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
/* Interfacing with 64-bit SCOTCH and pt-SCOTCH */
#include "mumps_scotch64.h"
#if defined(scotch) || defined(ptscotch)
void MUMPS_CALL
MUMPS_SCOTCH_64( const MUMPS_INT8 * const  n,        /* in    */
                 const MUMPS_INT8 * const  iwlen,    /* in    */
                       MUMPS_INT8 * const  petab,    /* inout */
                 const MUMPS_INT8 * const  pfree,    /* in    */
                       MUMPS_INT8 * const  lentab,   /* in (modified in ANA_H) */
                       MUMPS_INT8 * const  iwtab,    /* in (modified in ANA_H) */
                       MUMPS_INT8 * const  nvtab,    /* out   */
                       MUMPS_INT8 * const  elentab,  /* out   */
                       MUMPS_INT8 * const  lasttab,  /* out   */
                       MUMPS_INT  * const  ncmpa )   /* out   */
{
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
}
#endif
