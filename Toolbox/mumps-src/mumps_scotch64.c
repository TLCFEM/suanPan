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
/* Interfacing with 64-bit SCOTCH and pt-SCOTCH */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mumps_scotch64.h"
#if defined(scotch) || defined(ptscotch)
void MUMPS_CALL
MUMPS_SCOTCH_ORD_64( const MUMPS_INT8 * const  n,        /* in    */
              const MUMPS_INT8 * const  iwlen,    /* in    */
                    MUMPS_INT8 * const  petab,    /* inout */
              const MUMPS_INT8 * const  pfree,    /* in    */
                    MUMPS_INT8 * const  lentab,   /* in (modified in ANA_H) */
                    MUMPS_INT8 * const  iwtab,    /* in (modified in ANA_H) */
                    MUMPS_INT8 * const  nvtab,    /* out or inout if weight used on entry   */
                    MUMPS_INT8 * const  elentab,  /* out permutation on output  */
                    MUMPS_INT8 * const  lasttab,  /* out   */
                    MUMPS_INT  * const  ncmpa,    /* out   */
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
                    SCOTCH_Context * const contextptr,
#endif
                    MUMPS_INT * const  weightused,         /* out   */
                    MUMPS_INT * const  weightrequested )   /* in   */
{
/* weightused(out) = weightrequested since it is always used to build graph 
   FIXME it is not exploited on output and could be suppressed from interface 
                   = 0 otherwise
*/
/* weightused(out) = weightrequested since it is always used to build graph 
   FIXME it is not exploited on output and could be suppressed from interface 
                   = 0 otherwise
*/
  MUMPS_INT8 * vendtab ;                    /* Vertex end array */
  SCOTCH_Graph grafdat;                    /* Graph            */
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
  SCOTCH_Graph grafdat_with_context;
#endif
  SCOTCH_Strat stratdat;
  MUMPS_INT8 vertnum;
  int ierr;
  *weightused = *weightrequested;
  vendtab=malloc(*n * sizeof(MUMPS_INT8));
  for (vertnum = 0; vertnum < *n; vertnum ++)
    vendtab[vertnum] = petab[vertnum] + lentab[vertnum];
  ierr=SCOTCH_graphInit        (&grafdat);
  if ( *weightrequested == 1 )
  {
    ierr=SCOTCH_graphBuild (&grafdat, 1, *n, (SCOTCH_Num *) petab, (SCOTCH_Num *) vendtab, (SCOTCH_Num *) nvtab, NULL, *iwlen, (SCOTCH_Num *) iwtab, NULL); /* Assume Fortran-based indexing */
  }
  else
  {
    ierr=SCOTCH_graphBuild (&grafdat, 1, *n, (SCOTCH_Num *) petab, (SCOTCH_Num *) vendtab, NULL, NULL, *iwlen, (SCOTCH_Num *) iwtab, NULL); /* Assume Fortran-based indexing */
  }
  ierr=SCOTCH_stratInit(&stratdat);
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
  /* Initialize and bind grafdat_with_context */
  ierr=SCOTCH_graphInit        (&grafdat_with_context);
  ierr=SCOTCH_contextBindGraph(contextptr, &grafdat, &grafdat_with_context);
  *ncmpa=SCOTCH_graphOrder(&grafdat_with_context, &stratdat, (SCOTCH_Num *) elentab, (SCOTCH_Num *) lasttab, NULL, NULL, NULL);
#else
  /* order grafdat without threads context */
  *ncmpa=SCOTCH_graphOrder(&grafdat, &stratdat, (SCOTCH_Num *) elentab, (SCOTCH_Num *) lasttab, NULL, NULL, NULL);
#endif
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
  SCOTCH_graphExit(&grafdat_with_context);
#endif
  SCOTCH_stratExit(&stratdat);
  SCOTCH_graphExit(&grafdat);
  free(vendtab);
}
void MUMPS_CALL
MUMPS_SCOTCH_64( const MUMPS_INT8 * const  n,        /* in    */
                 const MUMPS_INT8 * const  iwlen,    /* in    */
                       MUMPS_INT8 * const  petab,    /* inout */
                 const MUMPS_INT8 * const  pfree,    /* in    */
                       MUMPS_INT8 * const  lentab,   /* in (modified in ANA_H) */
                       MUMPS_INT8 * const  iwtab,    /* in (modified in ANA_H) */
                       MUMPS_INT8 * const  nvtab,    /* out or inout if weight used on entry   */
                       MUMPS_INT8 * const  elentab,  /* out   */
                       MUMPS_INT8 * const  lasttab,  /* out   */
                       MUMPS_INT  * const  ncmpa,    /* out   */
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
/* FIXME: see how to pass contextptr to esmumps/esmumpsv */
                    SCOTCH_Context * const contextptr,
#endif
                       MUMPS_INT  * const  weightused,         /* out   */
                       MUMPS_INT  * const  weightrequested )   /* in   */
{
/* weightused(out) = 1 if weight of nodes provided in nvtab are used (esmumpsv is called) 
                   = 0 otherwise
*/
#if ((SCOTCH_VERSION == 6) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 7)
/* esmumpsv prototype with 64-bit integers weights of nodes in the graph are used on entry (nvtab) */
     if ( *weightrequested == 1 )
     {
       *ncmpa = esmumpsv( *n, *iwlen, petab, *pfree,
                          lentab, iwtab, nvtab, elentab, lasttab );
       *weightused=1;
     }
     else
     {
       /* esmumps prototype with standard integers (weights of nodes not used on entry) */
       *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                         lentab, iwtab, nvtab, elentab, lasttab );
       *weightused=0;
     }
#else
     /* esmumps prototype with standard integers (weights of nodes not used on entry) */
     *ncmpa = esmumps( *n, *iwlen, petab, *pfree,
                       lentab, iwtab, nvtab, elentab, lasttab );
     *weightused=0;
#endif
}
#endif
