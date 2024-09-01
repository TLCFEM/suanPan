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
      MODULE SMUMPS_LR_TYPE
      IMPLICIT NONE
      TYPE LRB_TYPE
        REAL,POINTER,DIMENSION(:,:) :: Q => null()
        REAL,POINTER,DIMENSION(:,:) :: R => null()
        INTEGER :: K,M,N
        LOGICAL :: ISLR
      END TYPE LRB_TYPE
      CONTAINS
      SUBROUTINE DEALLOC_LRB(LRB_OUT, KEEP8, K34
     &                       )
        TYPE(LRB_TYPE), INTENT(INOUT) :: LRB_OUT
        INTEGER(8) :: KEEP8(150)
        INTEGER :: K34
        INTEGER :: MEM, IDUMMY, JDUMMY
        IF (LRB_OUT%M.EQ.0) RETURN
        IF (LRB_OUT%N.EQ.0) RETURN
        MEM = 0
          IF (LRB_OUT%ISLR) THEN
            IF (associated(LRB_OUT%Q)) THEN
              MEM = MEM + size(LRB_OUT%Q)
              DEALLOCATE (LRB_OUT%Q)
              NULLIFY(LRB_OUT%Q)
            ENDIF
            IF (associated(LRB_OUT%R)) THEN
              MEM = MEM + size(LRB_OUT%R)
              DEALLOCATE (LRB_OUT%R)
              NULLIFY(LRB_OUT%R)
            ENDIF
          ELSE
            IF (associated(LRB_OUT%Q)) THEN
              MEM = MEM + size(LRB_OUT%Q)
              DEALLOCATE (LRB_OUT%Q)
              NULLIFY(LRB_OUT%Q)
            ENDIF
          ENDIF
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(-int(MEM,8),
     &        .TRUE., KEEP8, 
     &        IDUMMY, JDUMMY,    
     &        .TRUE., .TRUE.)
      END SUBROUTINE DEALLOC_LRB
      SUBROUTINE DEALLOC_BLR_PANEL(BLR_PANEL, IEND, KEEP8, K34, IBEG_IN)
        INTEGER, INTENT(IN)           :: IEND
        TYPE(LRB_TYPE), INTENT(INOUT) :: BLR_PANEL(:)
        INTEGER(8) :: KEEP8(150)
        INTEGER, INTENT(IN) :: K34
        INTEGER, INTENT(IN), OPTIONAL :: IBEG_IN
        INTEGER :: I, IBEG
        IF (present(IBEG_IN)) THEN
          IBEG = IBEG_IN
        ELSE
          IBEG = 1
        ENDIF
        IF (IEND.GE.IBEG) THEN
          IF (BLR_PANEL(1)%M.NE.0) THEN
            DO I=IBEG, IEND
              CALL DEALLOC_LRB(BLR_PANEL(I), KEEP8, K34)
            ENDDO
          ENDIF
        ENDIF
      END SUBROUTINE DEALLOC_BLR_PANEL
      END MODULE SMUMPS_LR_TYPE
