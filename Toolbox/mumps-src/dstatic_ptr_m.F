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
      MODULE DMUMPS_STATIC_PTR_M
      PUBLIC :: DMUMPS_TMP_PTR, DMUMPS_GET_TMP_PTR
      DOUBLE PRECISION, DIMENSION(:), POINTER, SAVE :: DMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE DMUMPS_SET_STATIC_PTR(ARRAY)
      DOUBLE PRECISION, DIMENSION(:), TARGET :: ARRAY
      DMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE DMUMPS_SET_STATIC_PTR
      SUBROUTINE DMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_NOF2003)
      DOUBLE PRECISION, DIMENSION(:), POINTER :: PTR
#else
      DOUBLE PRECISION, DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#endif
      PTR => DMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE DMUMPS_GET_TMP_PTR
      END MODULE DMUMPS_STATIC_PTR_M
