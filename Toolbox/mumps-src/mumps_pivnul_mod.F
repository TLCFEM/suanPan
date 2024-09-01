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
      MODULE MUMPS_PIVNUL_MOD
      TYPE PIVNUL_LIST_STRUCT_T
        INTEGER :: SIZE_PIVNUL_LIST
        INTEGER, POINTER ::  PIVNUL_LIST(:) => null()
      END TYPE PIVNUL_LIST_STRUCT_T
      PUBLIC :: MUMPS_RESIZE_PIVNUL
      CONTAINS
      SUBROUTINE MUMPS_RESIZE_PIVNUL (
     &           KEEP, N, PIVNUL_LIST_STRUCT, POS_NEWENTRY, 
     &           IFLAG, IERROR
     &           )
!$    USE OMP_LIB
      IMPLICIT NONE
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      INTEGER, INTENT(IN)            :: N, POS_NEWENTRY, KEEP(500)
      INTEGER, INTENT(INOUT)         :: IFLAG, IERROR
      INTEGER, POINTER, DIMENSION(:)  :: TEMP_PTR
      INTEGER  :: NEW_SIZE, IERR, I
      INTEGER, PARAMETER :: FI=10
      IF (KEEP(405).EQ.1) THEN
!$OMP CRITICAL(critical_pivnul)
       IF (PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST.LT.POS_NEWENTRY) THEN
         NEW_SIZE =  max (POS_NEWENTRY,
     &                    PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST*FI)
         NEW_SIZE =  min(NEW_SIZE, N)
         ALLOCATE(TEMP_PTR(NEW_SIZE), stat=IERR)
         IF (IERR.GT.0) THEN
           IFLAG  = -13
           IERROR = NEW_SIZE
         ELSE
           DO I=1, PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST
             TEMP_PTR(I) =  PIVNUL_LIST_STRUCT%PIVNUL_LIST(I)
           ENDDO
           DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
           PIVNUL_LIST_STRUCT%PIVNUL_LIST      => TEMP_PTR
           PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST =  NEW_SIZE
         ENDIF
       ENDIF
!$OMP END CRITICAL(critical_pivnul)
      ELSE
         NEW_SIZE =  max (POS_NEWENTRY,
     &                    PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST*FI)
         NEW_SIZE =  min(NEW_SIZE, N)
         ALLOCATE(TEMP_PTR(NEW_SIZE), stat=IERR)
         IF (IERR.GT.0) THEN
           IFLAG  = -13
           IERROR = NEW_SIZE
         ELSE
           DO I=1, PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST
             TEMP_PTR(I) =  PIVNUL_LIST_STRUCT%PIVNUL_LIST(I)
           ENDDO
           DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
           PIVNUL_LIST_STRUCT%PIVNUL_LIST      => TEMP_PTR
           PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST =  NEW_SIZE
         ENDIF
      ENDIF
      RETURN
      END  SUBROUTINE MUMPS_RESIZE_PIVNUL
      END MODULE MUMPS_PIVNUL_MOD
