C
C  This file is part of MUMPS 5.7.3, released
C  on Mon Jul 15 11:44:21 UTC 2024
C
C
C  Copyright 1991-2024 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      MODULE SMUMPS_STATIC_PTR_M
      PUBLIC :: SMUMPS_TMP_PTR, SMUMPS_GET_TMP_PTR
      REAL, DIMENSION(:), POINTER, SAVE :: SMUMPS_TMP_PTR
      CONTAINS
      SUBROUTINE SMUMPS_SET_STATIC_PTR(ARRAY)
      REAL, DIMENSION(:), TARGET :: ARRAY
      SMUMPS_TMP_PTR => ARRAY
      RETURN
      END SUBROUTINE SMUMPS_SET_STATIC_PTR
      SUBROUTINE SMUMPS_GET_TMP_PTR(PTR)
#if defined(MUMPS_NOF2003)
      REAL, DIMENSION(:), POINTER :: PTR
#else
      REAL, DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#endif
      PTR => SMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE SMUMPS_GET_TMP_PTR
      END MODULE SMUMPS_STATIC_PTR_M
