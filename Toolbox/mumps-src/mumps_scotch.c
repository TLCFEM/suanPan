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
/* Interfacing with SCOTCH and pt-SCOTCH */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mumps_scotch.h"
#if defined(scotch) || defined(ptscotch)
void MUMPS_CALL
MUMPS_SCOTCH_WEIGHTUSED( 
                        MUMPS_INT  * const  weightused         /* out   */
                       )
{
/* weightused(out) = 1 if weight of nodes can be used 
                   = 0 otherwise
*/
#if ((SCOTCH_VERSION == 6) && (SCOTCH_RELEASE >= 1)) || (SCOTCH_VERSION >= 7)
       *weightused=1;
#else       
       *weightused=0;
#endif
}
void MUMPS_CALL
MUMPS_SCOTCH_ORD( const MUMPS_INT * const  n,  /* NCMP or N */
              const MUMPS_INT * const  iwlen,   /* LIW8 */
                    MUMPS_INT * const  petab,   /* IPE */
              const MUMPS_INT * const  pfree,
                    MUMPS_INT * const  lentab,  /* numbers of edges for each vertex */
                    MUMPS_INT * const  iwtab,
                    MUMPS_INT * const  nvtab,   /*  weight of nodes */
                    MUMPS_INT * const  elentab, /* permutation on output (permtab) */
                    MUMPS_INT * const  lasttab,
                    MUMPS_INT * const  ncmpa,
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
  MUMPS_INT * vendtab ;                    /* Vertex end array */
  SCOTCH_Graph grafdat;                    /* Graph            */
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
  SCOTCH_Graph grafdat_with_context;
#endif
  SCOTCH_Strat stratdat;
  MUMPS_INT vertnum;
  int ierr;
  *weightused = *weightrequested;
  vendtab=malloc(*n * sizeof(MUMPS_INT));
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
MUMPS_SCOTCH( const MUMPS_INT * const  n,
              const MUMPS_INT * const  iwlen,
                    MUMPS_INT * const  petab,
              const MUMPS_INT * const  pfree,
                    MUMPS_INT * const  lentab,
                    MUMPS_INT * const  iwtab,
                    MUMPS_INT * const  nvtab,
                    MUMPS_INT * const  elentab,
                    MUMPS_INT * const  lasttab,
                    MUMPS_INT * const  ncmpa,
#if defined(MUMPS_SCOTCHIMPORTOMPTHREADS)
/* FIXME: see how to pass contextptr to esmumps/esmumpsv */
                    SCOTCH_Context * const contextptr,
#endif
                    MUMPS_INT * const  weightused,
                    MUMPS_INT * const  weightrequested )
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
void MUMPS_CALL
MUMPS_SCOTCH_VERSION(MUMPS_INT *version)
  {
  *version = SCOTCH_VERSION;
  return;
}
void MUMPS_CALL
MUMPS_SCOTCH_GET_PTHREAD_NUMBER (MUMPS_INT *PTHREAD_NUMBER)
  {
  *PTHREAD_NUMBER = -1;  /* NOT SET*/
#if (SCOTCH_VERSION>=7) 
  if (getenv("SCOTCH_PTHREAD_NUMBER"))
  {
    *PTHREAD_NUMBER = atoi(getenv("SCOTCH_PTHREAD_NUMBER"));
  }
#endif
  return;
}
void MUMPS_CALL
MUMPS_SCOTCH_SET_PTHREAD_NUMBER (MUMPS_INT *PTHREAD_NUMBER)
  {
#if (SCOTCH_VERSION>=7) 
  char param[32];
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
  int ierr;
#endif
  if (*PTHREAD_NUMBER == -1) 
  {
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
     ierr = _putenv("SCOTCH_PTHREAD_NUMBER=");
#else
     unsetenv("SCOTCH_PTHREAD_NUMBER");
#endif
  }
  else
  {
#if defined(MUMPS_WIN32) || defined(__MINGW32__)
    sprintf(param, "SCOTCH_PTHREAD_NUMBER=%d",*PTHREAD_NUMBER);
    ierr = _putenv(param);
#else
    sprintf(param, "%d", *PTHREAD_NUMBER);
    setenv("SCOTCH_PTHREAD_NUMBER",param,1);
#endif
  }
#endif
  return;
}
#endif /* scotch or ptscotch*/
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
