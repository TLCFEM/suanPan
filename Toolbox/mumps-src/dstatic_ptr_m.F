C
C  This file is part of MUMPS 5.2.1, released
C  on Fri Jun 14 14:46:05 UTC 2019
C
C
C  Copyright 1991-2019 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
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
#if defined(MUMPS_F2003)
      DOUBLE PRECISION, DIMENSION(:), POINTER, INTENT(OUT) :: PTR
#else
      DOUBLE PRECISION, DIMENSION(:), POINTER :: PTR
#endif
      PTR => DMUMPS_TMP_PTR
      RETURN
      END SUBROUTINE DMUMPS_GET_TMP_PTR
      END MODULE DMUMPS_STATIC_PTR_M
