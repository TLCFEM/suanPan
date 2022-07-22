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
      SUBROUTINE DMUMPS_GET_NS_OPTIONS_FACTO(N,KEEP,ICNTL,MPG)
      IMPLICIT NONE
          INTEGER N, KEEP(500), ICNTL(60), MPG
          KEEP(19)=0
          RETURN
      END SUBROUTINE DMUMPS_GET_NS_OPTIONS_FACTO
      SUBROUTINE DMUMPS_GET_NS_OPTIONS_SOLVE(ICNTL,KEEP,MPG,INFO)
      IMPLICIT NONE
      INTEGER KEEP(500), MPG, INFO(80), ICNTL(60)
      IF (KEEP(19).EQ.0.AND.KEEP(110).EQ.0) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 16  
         IF (KEEP(110).EQ.0) INFO(2) = 24 
          IF(MPG.GT.0) THEN
           WRITE( MPG,'(A)')
     &'** ERROR  : Null space computation requirement'
          WRITE( MPG,'(A)')
     &'** not consistent with factorization options'
         ENDIF
         GOTO 333
        ENDIF
      ENDIF
      IF (ICNTL(9).NE.1) THEN
        IF (KEEP(111).NE.0) THEN
         INFO(1) = -37
         INFO(2) = 9  
         IF (MPG.GT.0) THEN
          WRITE(MPG,'(A)')
     &'** ERROR  ICNTL(25) incompatible with '
          WRITE( MPG,'(A)')
     &'** option transposed system (ICNLT(9)=1) '
         ENDIF
        ENDIF
        GOTO 333
      ENDIF
 333  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_GET_NS_OPTIONS_SOLVE
      SUBROUTINE DMUMPS_RR_INIT_POINTERS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      NULLIFY(id%root%QR_TAU)
      NULLIFY(id%root%SVD_U)
      NULLIFY(id%root%SVD_VT)
      NULLIFY(id%root%SINGULAR_VALUES)
      RETURN
      END SUBROUTINE DMUMPS_RR_INIT_POINTERS
      SUBROUTINE DMUMPS_RR_FREE_POINTERS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) id
      IF (associated(id%root%QR_TAU))  THEN
        DEALLOCATE(id%root%QR_TAU)
        NULLIFY(id%root%QR_TAU)
      ENDIF
      IF (associated(id%root%SVD_U))  THEN
        DEALLOCATE(id%root%SVD_U)
        NULLIFY(id%root%SVD_U)
      ENDIF
      IF (associated(id%root%SVD_VT))  THEN
        DEALLOCATE(id%root%SVD_VT)
        NULLIFY(id%root%SVD_VT)
      ENDIF
      IF (associated(id%root%SINGULAR_VALUES))  THEN
        DEALLOCATE(id%root%SINGULAR_VALUES)
        NULLIFY(id%root%SINGULAR_VALUES)
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_RR_FREE_POINTERS
