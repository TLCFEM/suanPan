C
C  This file is part of MUMPS 5.6.0, released
C  on Wed Apr 19 15:50:57 UTC 2023
C
C
C  Copyright 1991-2023 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      SUBROUTINE MUMPS_PRINT_IF_DEFINED(MPG)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MPG
      IF (MPG.LE.0) RETURN
      write(MPG,*) "================================================="
#if defined(ALLOC_FROM_C)
      WRITE(MPG,*) " MUMPS compiled with option -DALLOC_FROM_C"
#endif
#if defined(AFFINITY_VERBOSE)
      write(MPG, *) "MUMPS compiled with option -DAFFINITY_VERBOSE"
#endif
#if defined(MTKO_ADVANCED_VERBOSE)
#if ! defined(MTKO_ADVANCED)
      write(MPG, *) "MUMPS should not be compiled with ",
     &              "-DMTKO_ADVANCED_VERBOSE without -DMTKO_ADVANCED"
      CALL MUMPS_ABORT()
#else
      write(MPG, *) "MUMPS compiled with option -DMTKO_ADVANCED_VERBOSE"
#endif
#endif
#if defined(ZERO_TRIANGLE)
      write(MPG,*) "MUMPS compiled with option -DZERO_TRIANGLE"
#endif
#if defined(PRINT_BACKTRACE_ON_ABORT)
      write(MPG, *)
     &  "MUMPS compiled with option -DPRINT_BACKTRACE_ON_ABORT"
#endif
#if defined(GEMMT_AVAILABLE)
      write(MPG, *) "MUMPS compiled with option -DGEMMT_AVAILABLE"
#endif
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DWORKAROUNDINTELILP64MPI2INTEGER"
#endif
#if defined(WORKAROUNDILP64MPICUSTOMREDUCE)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DWORKAROUNDILP64MPICUSTOMREDUCE"
#endif
#if defined(DETERMINISTIC_PARALLEL_GRAPH)
      write(MPG,*) "MUMPS compiled with option"
     &     ," -DDETERMINISTIC_PARALLEL_GRAPH"
#endif
#if defined(metis)
      write(MPG,*) "MUMPS compiled with option -Dmetis"
#endif
#if defined(metis4)
      write(MPG,*) "MUMPS compiled with option -Dmetis4"
#endif
#if defined(MUMPS_NOF2003)
      write(MPG,*) "MUMPS compiled with option -DMUMPS_NOF2003"
#endif
#if defined(NOSCALAPACK)
      write(MPG,*) "MUMPS compiled with option -DNOSCALAPACK"
#endif
#if defined(OLD_OOC_NOPANEL)
      write(MPG,*) "MUMPS compiled with option -DOLD_OOC_NOPANEL"
#endif
#if defined(parmetis)
      write(MPG,*) "MUMPS compiled with option -Dparmetis"
#endif
#if defined(parmetis3)
      write(MPG,*) "MUMPS compiled with option -Dparmetis3"
#endif
#if defined(pord)
      write(MPG,*) "MUMPS compiled with option -Dpord"
#endif
#if defined(ptscotch)
      write(MPG,*) "MUMPS compiled with option -Dptscotch"
#endif
#if defined(scotch)
      write(MPG,*) "MUMPS compiled with option -Dscotch"
#endif
#if defined(MUMPS_USE_BLAS2)
      write(MPG,*) "MUMPS compiled with option -DMUMPS_USE_BLAS2"
#endif
#if defined(BLR_NOOPENMP)
      write(MPG,*) "MUMPS compiled with option -DBLR_NOOPENMP"
#endif
#if defined(NODYNAMICCB)
      write(MPG,*) "MUMPS compiled with option -DNODYNAMICCB"
#endif
      write(MPG,*) "================================================="
      RETURN
      END SUBROUTINE MUMPS_PRINT_IF_DEFINED
