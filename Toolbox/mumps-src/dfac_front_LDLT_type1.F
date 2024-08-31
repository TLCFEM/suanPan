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
      MODULE DMUMPS_FAC1_LDLT_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC1_LDLT( N, INODE, IW, LIW, A, LA,
     &                     IOLDPS, POSELT, IFLAG, IERROR,
     &                     UU, NNEGW, NNULLNEGW, NPVW, NB22T1W, NBTINYW,
     &                     DET_EXPW, DET_MANTW, DET_SIGNW,
     &                     KEEP,KEEP8,
     &                     MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP, PIVNUL_LIST_STRUCT, IWPOS
     &     , LRGROUPS
     &     , PERM
     &     )
      USE DMUMPS_FAC_FRONT_AUX_M
      USE DMUMPS_OOC
      USE DMUMPS_FAC_LR  
      USE DMUMPS_LR_TYPE
      USE DMUMPS_LR_STATS
      USE DMUMPS_ANA_LR, ONLY : GET_CUT
      USE DMUMPS_LR_DATA_M
      USE MUMPS_PIVNUL_MOD, ONLY: PIVNUL_LIST_STRUCT_T
#if ! defined(BLR_NOOPENMP)
#endif
!$    USE OMP_LIB
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, IERROR
      INTEGER, intent(inout) :: NNEGW, NNULLNEGW, NPVW, NB22T1W, NBTINYW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER MYID, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A( LA )
      INTEGER, TARGET :: IW( LIW )
      INTEGER, intent(in) :: PERM(N)
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(KEEP(280))
      INTEGER INOPV, IFINB, NFRONT, NPIV, IBEG_BLOCK, IEND_BLOCK
      INTEGER NASS, NBKJIB_ORIG, XSIZE
      INTEGER :: LDA
      DOUBLE PRECISION UUTEMP
      LOGICAL STATICMODE
      DOUBLE PRECISION SEUIL_LOC
      LOGICAL IS_MAXFROMM_AVAIL
      INTEGER PIVOT_OPTION
      INTEGER LRTRSM_OPTION
      INTEGER LAST_ROW, FIRST_ROW
      DOUBLE PRECISION MAXFROMM
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, NextPiv2beWritten, IFLAG_OOC,
     &        IDUMMY, PP_FIRST2SWAP_L, PP_LastPIVRPTRFilled
      TYPE(IO_BLOCK) :: MonBloc 
      LOGICAL LAST_CALL
      INTEGER PARPIV_T1, OFFSET
      INTEGER NFS4FATHER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: M_ARRAY
      LOGICAL LASTBL
      INTEGER CURRENT_BLR
      LOGICAL LR_ACTIVATED 
      LOGICAL COMPRESS_CB, COMPRESS_PANEL
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR, OOC_EFFECTIVE_ON_FRONT,
     &        OOC_EFF_AND_WRITE_BYPANEL
      INTEGER K473_LOC
      INTEGER INFO_TMP(2), MAXI_RANK
      INTEGER FIRST_BLOCK, LAST_BLOCK
      INTEGER HF, NPARTSASS, NPARTSCB, NB_BLR
      INTEGER MAXI_CLUSTER, LWORK, NELIM, NELIM_LOC
      TYPE(LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      INTEGER, POINTER, DIMENSION(:)          :: PTDummy
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)   :: ACC_LUA
      INTEGER, POINTER, DIMENSION(:)          :: BEGS_BLR
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)   :: BLR_L
      DOUBLE PRECISION, POINTER, DIMENSION(:)          :: DIAG
      INTEGER, POINTER, DIMENSION(:)          :: BEGS_BLR_TMP
      TYPE(LRB_TYPE), POINTER, DIMENSION(:)  :: BLR_PANEL
      INTEGER :: DIAGSIZ_STA, DIAGSIZ_DYN, DIAGPOS, I, IP, MEM, MEM_TOT
      INTEGER(8) :: POSELT_DIAG
      DOUBLE PRECISION, ALLOCATABLE :: WORK(:), TAU(:)
      INTEGER, ALLOCATABLE :: JPVT(:)
      DOUBLE PRECISION, ALLOCATABLE :: RWORK(:)
      DOUBLE PRECISION, ALLOCATABLE :: BLOCK(:,:)
      INTEGER :: allocok,J
      INTEGER :: OMP_NUM
      INTEGER :: II,JJ
      INTEGER(8) :: UPOS, LPOS, DPOS
      DOUBLE PRECISION :: ONE, MONE, ZERO
      PARAMETER (ONE = 1.0D0, MONE=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LOGICAL :: SWAP_OCCURRED
      INTEGER :: MY_NUM
      TYPE(LRB_TYPE), POINTER, DIMENSION(:) :: NEXT_BLR_L
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_STATIC
      INTEGER :: NVSCHUR, NVSCHUR_K253, IROW_L
      INCLUDE 'mumps_headers.h'
      INTEGER NBLR_ORIG, IBEG_BLR, IEND_BLR
      INTEGER Inextpiv
      INTEGER PIVSIZ,IWPOSP2
      IS_MAXFROMM_AVAIL = .FALSE.
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      SWAP_OCCURRED = .FALSE.
      INOPV = 0
      IF(KEEP(97) .EQ. 0) THEN
         STATICMODE = .FALSE.
      ELSE
         STATICMODE = .TRUE.
      ENDIF
      UUTEMP=UU
      IF (AVOID_DELAYED) THEN
        STATICMODE = .TRUE.
        SEUIL_LOC = max(SEUIL,epsilon(SEUIL))
      ELSE
        SEUIL_LOC = SEUIL
      ENDIF
      LAFAC  = -9999_8  
      XSIZE  = KEEP(IXSZ)
      NFRONT = IW(IOLDPS+XSIZE)
      LDA    = NFRONT
      NASS   = abs(IW(IOLDPS+2+XSIZE))
      IW(IOLDPS+3+XSIZE) =  -99999
      LR_ACTIVATED= .FALSE.        
      COMPRESS_PANEL = .FALSE.
      COMPRESS_CB = .FALSE.
      NULLIFY(PTDummy)
      NULLIFY(BEGS_BLR)
      NULLIFY(CB_LRB)
      NULLIFY(ACC_LUA)
      NULLIFY(BLR_L)
      NULLIFY(BEGS_BLR_TMP)
      NULLIFY(BLR_PANEL)
      NULLIFY(DIAG)
      COMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      COMPRESS_CB    = ((IW(IOLDPS+XXLR).EQ.1).OR.
     &                  (IW(IOLDPS+XXLR).EQ.3))
      LR_ACTIVATED   = (IW(IOLDPS+XXLR).GT.0)
      IF (COMPRESS_CB.AND.(.NOT.COMPRESS_PANEL)) THEN
        COMPRESS_PANEL = .TRUE.
        K473_LOC = 1
      ELSE
        K473_LOC = KEEP(473)
      ENDIF
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      OOC_EFFECTIVE_ON_FRONT= ((KEEP(201).EQ.1).AND. 
     &                         OOCWRITE_COMPATIBLE_WITH_BLR)
      CALL DMUMPS_SET_PARPIVT1 ( INODE, NFRONT, NASS, KEEP, 
     &                           LR_ACTIVATED, PARPIV_T1)
      LRTRSM_OPTION = KEEP(475)
      PIVOT_OPTION = KEEP(468)
      IF (UUTEMP.EQ.ZERO) THEN 
         PIVOT_OPTION = 0
      ELSE IF (PARPIV_T1.NE.0) THEN
         PIVOT_OPTION = min(PIVOT_OPTION,2)
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF (LRTRSM_OPTION.EQ.3) THEN
          PIVOT_OPTION = MIN(PIVOT_OPTION,1)
        ELSEIF (LRTRSM_OPTION.EQ.2) THEN
          PIVOT_OPTION = MIN(PIVOT_OPTION, 2)
        ENDIF
      ENDIF
      IF (PIVOT_OPTION.LE.1) THEN
         PARPIV_T1 = 0
      ENDIF
      IF (NASS.LT.KEEP(4)) THEN
        NBKJIB_ORIG = NASS
      ELSE IF (NASS .GT. KEEP(3)) THEN
        NBKJIB_ORIG = min( KEEP(6), NASS )
      ELSE
        NBKJIB_ORIG = min( KEEP(5), NASS )
      ENDIF
      IF (.not.LR_ACTIVATED) THEN
          NBLR_ORIG     = KEEP(420)
      ELSE
          NBLR_ORIG  = -9999 
      ENDIF
      IF ((KEEP(114).EQ.1) .AND. 
     &    (KEEP(116).GT.0) .AND. ((NFRONT-NASS-KEEP(253)).GT.0) 
     &   ) THEN
         IROW_L = IOLDPS+6+XSIZE+NASS
         CALL DMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &     N, 
     &     NFRONT-NASS-KEEP(253), 
     &     KEEP(116), 
     &     IW(IROW_L), PERM, 
     &     NVSCHUR )
      ELSE
         NVSCHUR = 0
      ENDIF
      IEND_BLOCK  = 0
      IEND_BLR    = 0
      CURRENT_BLR = 0
      LASTBL      = .FALSE.
      CALL MUMPS_GETI8(LAFAC,IW(IOLDPS+XXR))
      LIWFAC    = IW(IOLDPS+XXI)
      IF (OOC_EFFECTIVE_ON_FRONT) THEN
          IDUMMY    = -8765
          NextPiv2beWritten = 1 
          PP_FIRST2SWAP_L = NextPiv2beWritten 
          MonBloc%LastPanelWritten_L = 0 
          PP_LastPIVRPTRFilled       = 0
          MonBloc%INODE    = INODE
          MonBloc%MASTER   = .TRUE.
          MonBloc%Typenode = 1
          MonBloc%NROW     = NFRONT
          MonBloc%NCOL     = NFRONT
          MonBloc%NFS      = NASS
          MonBloc%Last     = .FALSE.   
          MonBloc%LastPiv  = -77777    
          MonBloc%INDICES  =>
     &              IW(IOLDPS+6+NFRONT+XSIZE:
     &                 IOLDPS+5+NFRONT+XSIZE+NFRONT)
      ENDIF
      IF (LR_ACTIVATED) THEN
             IF (KEEP(405) .EQ. 1) THEN
!$OMP ATOMIC UPDATE
               CNT_NODES = CNT_NODES + 1 
!$OMP END ATOMIC
             ELSE
               CNT_NODES = CNT_NODES + 1 
             ENDIF
      ELSE IF (KEEP(486).NE.0) THEN
      ENDIF
      OOC_EFF_AND_WRITE_BYPANEL  = ( (PIVOT_OPTION.GE.3) .AND.
     &                                     OOC_EFFECTIVE_ON_FRONT )
      HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
      IF (LR_ACTIVATED) THEN
         CALL GET_CUT(IW(IOLDPS+HF:IOLDPS+HF+NFRONT-1), NASS,
     &        NFRONT-NASS, LRGROUPS, NPARTSCB, 
     &        NPARTSASS, BEGS_BLR)
         CALL REGROUPING2(BEGS_BLR, NPARTSASS, NASS, NPARTSCB,
     &        NFRONT-NASS, KEEP(488), .FALSE., KEEP(472), 
     &        NFRONT, KEEP(1))     
         NB_BLR = NPARTSASS + NPARTSCB
         call MAX_CLUSTER(BEGS_BLR,NB_BLR,MAXI_CLUSTER)
         MAXI_RANK = KEEP(479)*MAXI_CLUSTER
         LWORK = MAXI_CLUSTER*MAXI_CLUSTER
         OMP_NUM = 1
#if ! defined(BLR_NOOPENMP)
!$       OMP_NUM = OMP_GET_MAX_THREADS()
#endif
         ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
         IF (allocok > 0) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           GOTO 490
         ENDIF
         ALLOCATE(ACC_LUA(OMP_NUM),stat=allocok)
         IF (allocok > 0) THEN
            IFLAG  = -13
            IERROR = OMP_NUM
            GOTO 490
         ENDIF
         IF (KEEP(480).GE.3) THEN
           DO MY_NUM=1,OMP_NUM
             CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK,
     &                      MAXI_CLUSTER, MAXI_CLUSTER, .TRUE.,
     &                      IFLAG, IERROR, KEEP8)
             IF (IFLAG.LT.0)  GOTO 500
             ACC_LUA(MY_NUM)%K = 0
           ENDDO
         ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.(KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       .OR.COMPRESS_CB
     &      )) THEN
        INFO_TMP(1) = IFLAG
        INFO_TMP(2) = IERROR
        IF (IFLAG.LT.0) GOTO 500
        CALL DMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &              .TRUE., 
     &              .FALSE., 
     &              .FALSE., 
     &              NPARTSASS, 
     &              BEGS_BLR, PTDummy, 
     &              huge(NPARTSASS),  
     &              INFO_TMP)
        IFLAG  = INFO_TMP(1) 
        IERROR = INFO_TMP(2) 
        IF (IFLAG.LT.0) GOTO 500
      ENDIF
      IF (COMPRESS_CB.AND.NPARTSCB.GT.0) THEN
        allocate(CB_LRB(NPARTSCB,NPARTSCB),stat=allocok)
        IF (allocok > 0) THEN
          IFLAG  = -13
          IERROR = NPARTSCB*NPARTSCB
          GOTO 490
        ENDIF
        DO II=1,NPARTSCB
        DO JJ=1,NPARTSCB
          CB_LRB(II,JJ)%M=0
          CB_LRB(II,JJ)%N=0
          NULLIFY(CB_LRB(II,JJ)%Q)
          NULLIFY(CB_LRB(II,JJ)%R)
          CB_LRB(II,JJ)%ISLR = .FALSE.
        ENDDO
        ENDDO
        CALL DMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
      ENDIF
      DO WHILE (IEND_BLR < NASS ) 
        CURRENT_BLR = CURRENT_BLR + 1
        IBEG_BLR = IW(IOLDPS+1+KEEP(IXSZ)) + 1 
        IF (.NOT. LR_ACTIVATED) THEN
          IEND_BLR = min(IEND_BLR + NBLR_ORIG, NASS)
        ELSE
          IEND_BLR = BEGS_BLR(CURRENT_BLR+1)-1 
          BEGS_BLR( CURRENT_BLR ) = IBEG_BLR
          IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
            MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
            LWORK = MAXI_CLUSTER*MAXI_CLUSTER
            DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
            ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              GOTO 490
            ENDIF
            IF (KEEP(480).GE.3) THEN
              DO MY_NUM=1,OMP_NUM
                CALL DEALLOC_LRB(ACC_LUA(MY_NUM), KEEP8, KEEP(34))
                CALL ALLOC_LRB(ACC_LUA(MY_NUM), MAXI_RANK,
     &                         MAXI_CLUSTER, MAXI_CLUSTER, .TRUE.,
     &                         IFLAG, IERROR, KEEP8)
                IF (IFLAG.LT.0) GOTO 500
                ACC_LUA(MY_NUM)%K = 0
              ENDDO
            ENDIF
          ENDIF
          IF (KEEP(480).GE.5) THEN
            IF (CURRENT_BLR.EQ.1) THEN
              ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR
                 GOTO 490
              ENDIF 
              IF (NB_BLR.GT.CURRENT_BLR) THEN
                BLR_L(1:NB_BLR-CURRENT_BLR)%ISLR=.FALSE.
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              ENDIF
            ELSE
              IF (NB_BLR.GT.CURRENT_BLR) THEN
                CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              ENDIF
            ENDIF
            IF (CURRENT_BLR.LT.NPARTSASS) THEN
              ALLOCATE(NEXT_BLR_L(NB_BLR-CURRENT_BLR-1),stat=allocok)
              IF (allocok > 0) THEN
                 IFLAG  = -13
                 IERROR = NB_BLR-CURRENT_BLR-1
                 GOTO 490
              ENDIF
              IF (NB_BLR.GT.CURRENT_BLR+1) THEN
                CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR+1, NEXT_BLR_L)
              ENDIF
            ENDIF
          ELSE
             ALLOCATE(BLR_L(NB_BLR-CURRENT_BLR),stat=allocok)
             IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = NB_BLR-CURRENT_BLR
                GOTO 490
             ENDIF 
          ENDIF
        ENDIF
        IF (LR_ACTIVATED) THEN
        ENDIF
        DO WHILE (IEND_BLOCK < IEND_BLR ) 
          IBEG_BLOCK = IW(IOLDPS+1+KEEP(IXSZ)) + 1
          IF (KEEP(405).EQ.0) THEN
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
          ELSE
!$OMP       ATOMIC UPDATE
            KEEP(425)=max(KEEP(425),IEND_BLOCK-IBEG_BLOCK)
!$OMP       END ATOMIC
          ENDIF
          IEND_BLOCK = min(IEND_BLOCK + NBKJIB_ORIG, IEND_BLR)
  50      CONTINUE  
            CALL DMUMPS_FAC_I_LDLT(NFRONT,NASS,N,INODE,
     &                IBEG_BLOCK, IEND_BLOCK,
     &                IW,LIW,A,LA,
     &                INOPV, NNEGW, NNULLNEGW, NB22T1W, NBTINYW,
     &                DET_EXPW, DET_MANTW, DET_SIGNW,
     &                IFLAG,IERROR,IOLDPS,POSELT,UUTEMP,
     &                SEUIL_LOC,KEEP,KEEP8,PIVSIZ,
     &      DKEEP(1), PIVNUL_LIST_STRUCT, SWAP_OCCURRED, XSIZE,
     &      PP_FIRST2SWAP_L, MonBloc%LastPanelWritten_L,
     &      PP_LastPIVRPTRFilled, MAXFROMM, IS_MAXFROMM_AVAIL,
     &      PIVOT_OPTION, IEND_BLR, Inextpiv, 
     &      OOC_EFFECTIVE_ON_FRONT,
     &      NVSCHUR, PARPIV_T1, LR_ACTIVATED
     &      )
            IF (IFLAG.LT.0) GOTO 500
          IF (INOPV.EQ.1) THEN
            IF(STATICMODE) THEN
              INOPV = -1
              GOTO 50 
            ENDIF
            LASTBL = .TRUE.
          ELSE IF ( INOPV.LE.0 ) THEN 
            INOPV = 0
            NPVW = NPVW + PIVSIZ
            NVSCHUR_K253 = 0
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_ROW = NFRONT
              NVSCHUR_K253 = NVSCHUR + KEEP(253) 
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_ROW = NASS
            ELSE
              LAST_ROW = IEND_BLR
            ENDIF
            CALL DMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &             NFRONT, NASS, IW(IOLDPS+1+XSIZE),
     &             INODE,A,LA,
     &             LDA, 
     &             POSELT,IFINB,
     &             PIVSIZ, MAXFROMM,
     &             IS_MAXFROMM_AVAIL, (UUTEMP.NE.0.0D0),
     &             PARPIV_T1, 
     &             LAST_ROW, IEND_BLR, NVSCHUR_K253,
     &             LR_ACTIVATED
     &             )
            IF(PIVSIZ .EQ. 2) THEN
              IWPOSP2 = IOLDPS+IW(IOLDPS+1+XSIZE)+6
              IW(IWPOSP2+NFRONT+XSIZE) =
     &                              -IW(IWPOSP2+NFRONT+XSIZE)
            ENDIF
            IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + PIVSIZ
            IF (IFINB.EQ.0) THEN
              GOTO 50 
            ELSE IF (IFINB.EQ.-1) THEN
              LASTBL = .TRUE.
            ENDIF
          ENDIF
          IF ( OOC_EFF_AND_WRITE_BYPANEL ) THEN
            MonBloc%Last = LASTBL
            MonBloc%LastPiv= IW(IOLDPS+1+XSIZE)
            LAST_CALL=.FALSE.
            CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE,
     &        TYPEF_L, A(POSELT),
     &        LAFAC, MonBloc, NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC,
     &        MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
            IF (IFLAG_OOC < 0 ) THEN
              IFLAG=IFLAG_OOC
              GOTO 500
            ENDIF
          ENDIF
          NPIV       =  IW(IOLDPS+1+XSIZE)
          IF ( IEND_BLR .GT. IEND_BLOCK ) THEN
            IF (PIVOT_OPTION.GE.3) THEN
              LAST_ROW = NFRONT
            ELSEIF (PIVOT_OPTION.EQ.2) THEN
              LAST_ROW = NASS
            ELSE
              LAST_ROW = IEND_BLR
            ENDIF
              CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,
     &             NPIV, NFRONT,NASS,INODE,A,LA,
     &             LDA, POSELT,
     &             KEEP, KEEP8,
     &             -6666, -6666, 
     &             IEND_BLR, LAST_ROW,
     &             .FALSE., .TRUE., LR_ACTIVATED,
     &             IW, LIW, -6666 
     &             )
          ENDIF
        END DO 
        NPIV   =  IW(IOLDPS+1+XSIZE)
        IF (.NOT. LR_ACTIVATED
     &      .OR. (.NOT. COMPRESS_PANEL)
     &     ) THEN
          IF (PIVOT_OPTION.GE.3) THEN
            LAST_ROW = NFRONT
          ELSEIF (PIVOT_OPTION.EQ.2) THEN
            LAST_ROW = NASS
          ELSE
            LAST_ROW = IEND_BLR
          ENDIF
          CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLR,IEND_BLR,NPIV,
     &             NFRONT,NASS,INODE,A,LA,
     &             LDA, POSELT,
     &             KEEP, KEEP8,
     &             IEND_BLR, NASS,
     &             NASS, LAST_ROW,
     &             (PIVOT_OPTION.LE.1), .TRUE., LR_ACTIVATED,
     &             IW, LIW, IOLDPS+6+XSIZE+NFRONT+IBEG_BLR-1) 
        ELSE
          NELIM = IEND_BLOCK - NPIV
          IF (NELIM .EQ. IEND_BLR - IBEG_BLR + 1) THEN
            IF (KEEP(480).GE.2
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &         ) THEN
              DO J=1,NB_BLR-CURRENT_BLR
                 BLR_L(J)%M=0
                 BLR_L(J)%N=0
                 BLR_L(J)%K=0
                 BLR_L(J)%ISLR=.FALSE.
                 NULLIFY(BLR_L(J)%Q)
                 NULLIFY(BLR_L(J)%R)
              ENDDO
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
              NULLIFY(BLR_L)
            IF (KEEP(480).GE.2 .AND. IEND_BLR.LT.NASS) THEN
              IF (LRTRSM_OPTION.EQ.2) THEN
                FIRST_BLOCK = NPARTSASS-CURRENT_BLR
              ELSE
                FIRST_BLOCK = 1
              ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
              CALL DMUMPS_BLR_UPD_PANEL_LEFT_LDLT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 
     &          BEGS_BLR, CURRENT_BLR, NB_BLR, NPARTSASS,
     &          NELIM,
     &          IW(HF+IOLDPS+NFRONT), BLOCK,
     &          ACC_LUA, MAXI_CLUSTER, MAXI_RANK,
     &          1, IFLAG, IERROR,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), KEEP8, FIRST_BLOCK=FIRST_BLOCK)
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
              IF (IFLAG.LT.0) GOTO 500
            ENDIF
            ENDIF
            IF (KEEP(486).EQ.3) THEN
              IF (KEEP(480).EQ.0) THEN
                DEALLOCATE(BLR_L)
                NULLIFY(BLR_L)
              ENDIF
            ENDIF
            GOTO 100
          ENDIF
          IF (PIVOT_OPTION.GE.3) THEN
            FIRST_ROW = NFRONT
          ELSEIF (PIVOT_OPTION.EQ.2) THEN
            FIRST_ROW = NASS
          ELSE
            FIRST_ROW = IEND_BLR
          ENDIF
          IF (LRTRSM_OPTION.EQ.3) THEN
            LAST_ROW = IEND_BLR
          ELSEIF (LRTRSM_OPTION.EQ.2) THEN
            LAST_ROW = NASS
          ELSE
            LAST_ROW = NFRONT
          ENDIF
          IF ((IEND_BLR.LT.NFRONT) .AND. (LAST_ROW-FIRST_ROW.GT.0)) THEN
            CALL DMUMPS_FAC_SQ_LDLT(IBEG_BLR, IEND_BLR,
     &            NPIV, NFRONT, NASS, 
     &            INODE, A, LA, LDA, POSELT, 
     &            KEEP, KEEP8, 
     &            FIRST_ROW, LAST_ROW,
     &            -6666, -6666,  
     &            .TRUE., .FALSE., LR_ACTIVATED,
     &            IW, LIW, IOLDPS+6+XSIZE+NFRONT+IBEG_BLR-1) 
          ENDIF
#if ! defined(BLR_NOOPENMP)
#endif
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL PRIVATE(UPOS,LPOS,DPOS,OFFSET) 
!$OMP&         FIRSTPRIVATE(FIRST_BLOCK,LAST_BLOCK)
#endif
          CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG, IERROR, 
     &        NFRONT,
     &        BEGS_BLR, NB_BLR, DKEEP(8), KEEP(466), KEEP(458),
     &        K473_LOC, BLR_L, 
     &        CURRENT_BLR,
     &        'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCK, MAXI_CLUSTER, NELIM,
     &        .FALSE., 0, 0,
     &        1, KEEP(483), KEEP8,
     &        K480=KEEP(480)
     &        )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif              
          IF (IFLAG.LT.0) GOTO 400
          IF (PIVOT_OPTION.LT.3) THEN
            IF (LRTRSM_OPTION.GE.2) THEN
              IF (PIVOT_OPTION.LE.1.AND.LRTRSM_OPTION.EQ.3) THEN
                FIRST_BLOCK = CURRENT_BLR+1
              ELSE
                FIRST_BLOCK = NPARTSASS+1
              ENDIF
              CALL DMUMPS_BLR_PANEL_LRTRSM(A, LA, POSELT, NFRONT,
     &              IBEG_BLR, NB_BLR, BLR_L, 
     &              CURRENT_BLR, FIRST_BLOCK, NB_BLR,
     &              1, 1, 0, 
     &              .FALSE.,
     &              IW, OFFSET_IW=IOLDPS+6+XSIZE+NFRONT+IBEG_BLR-1)
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
            ENDIF
            IF (NELIM.GT.0) THEN
              IF (PIVOT_OPTION.LE.1) THEN
                FIRST_BLOCK = CURRENT_BLR+1
              ELSE
                FIRST_BLOCK = NPARTSASS+1
              ENDIF
              LPOS = POSELT
     &         +int(BEGS_BLR(CURRENT_BLR+1)-1-NELIM,8)*int(NFRONT,8)
     &         +int(BEGS_BLR(CURRENT_BLR)-1,8)
              DPOS = POSELT 
     &         +int(BEGS_BLR(CURRENT_BLR)-1,8)*int(NFRONT,8)
     &         +int(BEGS_BLR(CURRENT_BLR)-1,8)
              OFFSET=IOLDPS+6+XSIZE+NFRONT+IBEG_BLR-1
              UPOS = POSELT+int(BEGS_BLR(CURRENT_BLR)-1,8)*int(NFRONT,8)
     &             +int(BEGS_BLR(CURRENT_BLR+1)-1-NELIM,8)
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif          
              CALL DMUMPS_FAC_LDLT_COPYSCALE_U( NELIM, 1, 
     &             KEEP(424), NFRONT, NPIV-IBEG_BLR+1, 
     &             LIW, IW, OFFSET, LA, A, POSELT, LPOS, UPOS, DPOS)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif          
              LPOS = POSELT
     &           +int(BEGS_BLR(CURRENT_BLR+1)-1,8)*int(NFRONT,8)
     &           +int(BEGS_BLR(CURRENT_BLR+1)-1-NELIM,8)
              CALL DMUMPS_BLR_UPD_NELIM_VAR_L(
     &          A, LA, UPOS, A, LA, LPOS,
     &          IFLAG, IERROR, NFRONT, NFRONT,
     &          BEGS_BLR, CURRENT_BLR, BLR_L, NB_BLR,
     &          FIRST_BLOCK, NELIM, 'N')
            ENDIF
          ENDIF
          IF (IFLAG.LT.0) GOTO 400
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif          
          IF (KEEP(480).NE.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            IF (KEEP(480).LT.5) THEN
              CALL DMUMPS_BLR_SAVE_PANEL_LORU (
     &              IW(IOLDPS+XXF),
     &              0, 
     &              CURRENT_BLR, BLR_L)
            ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif          
          IF (KEEP(480).GE.2) THEN
            IF (IEND_BLR.LT.NASS) THEN
              IF (LRTRSM_OPTION.EQ.2) THEN
                FIRST_BLOCK = NPARTSASS-CURRENT_BLR
              ELSE
                FIRST_BLOCK = 1
              ENDIF
              CALL DMUMPS_BLR_UPD_PANEL_LEFT_LDLT(A, LA, POSELT,
     &          NFRONT, IW(IOLDPS+XXF), 
     &          BEGS_BLR, CURRENT_BLR, NB_BLR, NPARTSASS,
     &          NELIM,
     &          IW(HF+IOLDPS+NFRONT), BLOCK,
     &          ACC_LUA, MAXI_CLUSTER, MAXI_RANK,
     &          1, IFLAG, IERROR,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(483), KEEP8, FIRST_BLOCK=FIRST_BLOCK)
            ENDIF
          ELSE
            CALL DMUMPS_BLR_UPDATE_TRAILING_LDLT(A, LA, POSELT, 
     &        IFLAG, IERROR, NFRONT,
     &        BEGS_BLR, NB_BLR, CURRENT_BLR, BLR_L, NELIM,
     &        IW(HF+IOLDPS+NFRONT+IBEG_BLR-1), BLOCK,
     &        MAXI_CLUSTER, NPIV,
     &        1, 
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 400
          IF (LRTRSM_OPTION.GE.2) THEN
            IF (LRTRSM_OPTION.EQ.2) THEN
              FIRST_BLOCK = NPARTSASS+1
            ELSE
              FIRST_BLOCK = CURRENT_BLR+1
            ENDIF
            IF (KEEP(486).NE.2) THEN
              LAST_BLOCK = NB_BLR
            ELSEIF(UU.GT.0) THEN
              LAST_BLOCK = NPARTSASS
            ELSE
              LAST_BLOCK = CURRENT_BLR
            ENDIF
            CALL DMUMPS_DECOMPRESS_PANEL(A, LA, POSELT, NFRONT, NFRONT,
     &        .TRUE.,   
     &        BEGS_BLR(CURRENT_BLR),
     &        BEGS_BLR(CURRENT_BLR+1), NB_BLR, BLR_L, CURRENT_BLR, 'V',
     &        1,
     &        BEG_I_IN=FIRST_BLOCK, END_I_IN=LAST_BLOCK)
          ENDIF
 400      CONTINUE         
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 500
          IF (KEEP(486).EQ.3) THEN
            IF (KEEP(480).EQ.0) THEN
              CALL DEALLOC_BLR_PANEL(BLR_L, NB_BLR-CURRENT_BLR, KEEP8,
     &                               KEEP(34))
              DEALLOCATE(BLR_L)
            ELSE
              NULLIFY(NEXT_BLR_L)
            ENDIF
          ENDIF
          NULLIFY(BLR_L)
        ENDIF
        IF ( OOC_EFF_AND_WRITE_BYPANEL ) THEN
             MonBloc%Last = LASTBL
             MonBloc%LastPiv= NPIV
             LAST_CALL=.FALSE.
             CALL DMUMPS_OOC_IO_LU_PANEL(
     &          STRAT_TRY_WRITE,
     &          TYPEF_L, A(POSELT),
     &          LAFAC, MonBloc, NextPiv2beWritten, IDUMMY, IW(IOLDPS),
     &          LIWFAC, MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
                IF (IFLAG_OOC .LT. 0 ) THEN
                  IFLAG = IFLAG_OOC
                  GOTO 500
                ENDIF
        ENDIF
 100    CONTINUE
      END DO 
      IF (LR_ACTIVATED) THEN
        IBEG_BLR = IW(IOLDPS+1+XSIZE) + 1 
        BEGS_BLR( CURRENT_BLR + 1 ) = IBEG_BLR
        IF ( (KEEP(486).EQ.2) 
     &       ) THEN
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW(IOLDPS+XXF),
     &                      BEGS_BLR_STATIC)
          IF (UU.GT.0 .AND. SWAP_OCCURRED) THEN
            allocate(BEGS_BLR_TMP(NB_BLR+1),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = NB_BLR+1
              GOTO 500
            ENDIF
            DO J=1,NB_BLR+1
               BEGS_BLR_TMP(J) = BEGS_BLR_STATIC(J)
            ENDDO
          ENDIF
        ENDIF
        MEM_TOT = 0
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
!$OMP& PRIVATE(IP, NELIM_LOC, BLR_PANEL)
#endif
        IF (
     &         (KEEP(486).EQ.2) 
     &       ) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP DO PRIVATE(DIAG, DIAGSIZ_STA, DIAGSIZ_DYN, DIAGPOS, POSELT_DIAG,
!$OMP&           MEM, allocok)
!$OMP&   REDUCTION(+:MEM_TOT)
#endif
          DO IP=1,NPARTSASS
            IF (IFLAG.LT.0) CYCLE
            DIAGSIZ_DYN = BEGS_BLR(IP+1)-BEGS_BLR(IP)
            DIAGSIZ_STA = BEGS_BLR_STATIC(IP+1)-BEGS_BLR(IP)
            MEM = DIAGSIZ_DYN*DIAGSIZ_STA
            MEM_TOT = MEM_TOT + MEM
            ALLOCATE(DIAG(MEM),stat=allocok)
            IF (allocok > 0) THEN
              IFLAG  = -13
              IERROR = MEM
              CYCLE
            ENDIF 
            DIAGPOS = 1
            POSELT_DIAG = POSELT + int(BEGS_BLR(IP)-1,8)*int(NFRONT,8)
     &                             + int(BEGS_BLR(IP)-1,8)
            DO I=1,DIAGSIZ_STA
              DIAG(DIAGPOS:DIAGPOS+DIAGSIZ_DYN-1) =
     &                 A(POSELT_DIAG:POSELT_DIAG+int(DIAGSIZ_DYN-1,8))
              DIAGPOS = DIAGPOS + DIAGSIZ_DYN
              POSELT_DIAG = POSELT_DIAG + int(NFRONT,8)
            ENDDO
            CALL DMUMPS_BLR_SAVE_DIAG_BLOCK(
     &            IW(IOLDPS+XXF),
     &            IP, DIAG)
          ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP ENDDO
!$OMP SINGLE
#endif
          CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(int(MEM_TOT,8),
     &        (KEEP(405).NE.0), KEEP8, IFLAG, IERROR, .TRUE., .TRUE.)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
          IF (IFLAG.LT.0) GOTO 447
          IF (UU.GT.0 .AND. SWAP_OCCURRED) THEN
            DO IP=1,NPARTSASS
              NELIM_LOC = BEGS_BLR_TMP(IP+1)-BEGS_BLR(IP+1)
              CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &         IW(IOLDPS+XXF), 0, IP, BLR_PANEL)
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif
              CALL DEALLOC_BLR_PANEL(BLR_PANEL, NPARTSASS-IP, KEEP8,
     &                                 KEEP(34))
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
              CALL DMUMPS_COMPRESS_PANEL(A, LA, POSELT, IFLAG,
     &            IERROR, NFRONT, BEGS_BLR_TMP,
     &            NB_BLR, DKEEP(8), KEEP(466), KEEP(458), K473_LOC,
     &            BLR_PANEL, IP,
     &            'V', WORK, TAU, JPVT, LWORK, RWORK,
     &            BLOCK, MAXI_CLUSTER, NELIM_LOC,
     &            .FALSE., 0, 0,
     &            1, KEEP(483), KEEP8,
     &          END_I_IN=NPARTSASS, FRSWAP=.TRUE.
     &          )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
              IF (IFLAG.LT.0) GOTO 445
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif
              BEGS_BLR_TMP(IP+1) = BEGS_BLR(IP+1)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
            ENDDO
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
 445        CONTINUE
          ENDIF 
 447      CONTINUE
        ENDIF 
        IF (IFLAG .LT. 0) GOTO 450
        IF (KEEP(480) .GE. 2) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP SINGLE
#endif
          CALL DMUMPS_BLR_RETRIEVE_BEGSBLR_STA(IW(IOLDPS+XXF),
     &                        BEGS_BLR_STATIC)
#if ! defined(BLR_NOOPENMP)
!$OMP END SINGLE
#endif
          CALL DMUMPS_BLR_UPD_CB_LEFT_LDLT(A, LA, POSELT, NFRONT,
     &          BEGS_BLR_STATIC, BEGS_BLR, NPARTSCB, NPARTSASS, NASS,
     &          IW(IOLDPS+XXF),
     &          IW(HF+IOLDPS+NFRONT), BLOCK,
     &          ACC_LUA, MAXI_CLUSTER, MAXI_RANK,
     &          1, IFLAG, IERROR,
     &          KEEP(481), DKEEP(11), KEEP(466), KEEP(477), 
     &          KEEP(480), KEEP(479), KEEP(478), KEEP(476), 
     &          KEEP(484), KEEP8)
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
        ENDIF
        IF (IFLAG.LT.0) GOTO 450
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
        IF (COMPRESS_CB
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &       ) THEN
            CALL DMUMPS_BLR_SAVE_BEGS_BLR_DYN(IW(IOLDPS+XXF),
     &        BEGS_BLR)
        ENDIF
        IF (COMPRESS_CB) THEN
            IEND_BLR = BEGS_BLR(CURRENT_BLR+2)
            IF ( IEND_BLR - IBEG_BLR + 1 .GT. MAXI_CLUSTER ) THEN
              MAXI_CLUSTER = IEND_BLR - IBEG_BLR + 1
              LWORK = MAXI_CLUSTER*MAXI_CLUSTER
              DEALLOCATE(BLOCK, WORK, RWORK, TAU, JPVT)
              ALLOCATE(BLOCK(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &             RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &             TAU(MAXI_CLUSTER*OMP_NUM),
     &             JPVT(MAXI_CLUSTER*OMP_NUM), 
     &             WORK(LWORK*OMP_NUM),stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
              ENDIF
            ENDIF
        ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
        IF (IFLAG.LT.0) GOTO 450
        IF (COMPRESS_CB) THEN
#if ! defined(BLR_NOOPENMP)
!$OMP MASTER
#endif
            NFS4FATHER = -9999
            IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2) ) THEN
             CALL DMUMPS_BLR_RETRIEVE_NFS4FATHER ( IW(IOLDPS+XXF),
     &             NFS4FATHER )
             IF (NFS4FATHER.GE.0) NFS4FATHER = NFS4FATHER + NELIM
            ENDIF
            ALLOCATE(M_ARRAY(max(NFS4FATHER,1)), stat=allocok)
            IF ( allocok.GT.0 ) THEN
                  IFLAG = -13
                  IERROR = max(NFS4FATHER,1)
            ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
            IF (IFLAG.LT.0) GOTO 448
            CALL DMUMPS_COMPRESS_CB(A, LA, POSELT, NFRONT,
     &      BEGS_BLR, BEGS_BLR, NPARTSCB, NPARTSCB, NPARTSASS,
     &      NFRONT-NASS, NFRONT-NASS, INODE,
     &      IW(IOLDPS+XXF), 2, 1, IFLAG, IERROR,
     &      DKEEP(12), KEEP(466), KEEP(484), KEEP(489), CB_LRB,
     &      WORK, TAU, JPVT, LWORK, RWORK, BLOCK,
     &      MAXI_CLUSTER, KEEP8,
     &      NFS4FATHER, NPIV, NVSCHUR+KEEP(253), KEEP(1), 
     &      M_ARRAY=M_ARRAY,
     &      NELIM=NELIM )
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
            IF (IFLAG.LT.0) GOTO 448
#if ! defined(BLR_NOOPENMP)
!$OMP  MASTER
#endif
            IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
                 INFO_TMP(1) = IFLAG
                 INFO_TMP(2) = IERROR
                 CALL DMUMPS_BLR_SAVE_M_ARRAY( IW(IOLDPS+XXF),
     &            M_ARRAY, INFO_TMP)
                 IFLAG  = INFO_TMP(1) 
                 IERROR = INFO_TMP(2) 
            ENDIF
            DEALLOCATE(M_ARRAY)
#if ! defined(BLR_NOOPENMP)
!$OMP END MASTER
!$OMP BARRIER
#endif
 448        CONTINUE         
        ENDIF
 450    CONTINUE          
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
        IF (
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND.UU.GT.0 .AND. SWAP_OCCURRED
     &       ) THEN
            deallocate(BEGS_BLR_TMP)
        ENDIF
        IF (IFLAG.LT.0) GOTO 500
        CALL UPD_MRY_LU_FR(NASS, NFRONT-NASS, 1, NASS-NPIV)
        DO IP=1,NPARTSASS
          CALL DMUMPS_BLR_RETRIEVE_PANEL_LORU(
     &             IW(IOLDPS+XXF), 0, IP, BLR_PANEL)
          CALL UPD_MRY_LU_LRGAIN(BLR_PANEL, NB_BLR-IP
     &             )
        ENDDO
        CALL UPD_FLOP_FACTO_FR(NFRONT, NASS, NPIV, 2, 1)
      ENDIF
      IF (.NOT. COMPRESS_PANEL)  THEN
        CALL DMUMPS_FAC_T_LDLT(NFRONT,NASS,IW,LIW,A,LA,
     &         LDA, IOLDPS,POSELT, KEEP,KEEP8,
     &         (PIVOT_OPTION.NE.3), ETATASS,
     &         TYPEF_L, LAFAC, MonBloc, NextPiv2beWritten,
     &         LIWFAC, MYID, IFLAG, IOLDPS+6+XSIZE+NFRONT, INODE )
      ENDIF
      IF (KEEP(486).NE.0) THEN
        IF (.NOT.LR_ACTIVATED) THEN
          CALL UPD_FLOP_FRFRONTS(NFRONT, NPIV, NASS, 1, 1)
        ENDIF
      ENDIF
      IF (OOC_EFFECTIVE_ON_FRONT) THEN
          STRAT            = STRAT_WRITE_MAX   
          MonBloc%Last     = .TRUE.
          MonBloc%LastPiv  = IW(IOLDPS+1+XSIZE)
          LAST_CALL    = .TRUE.
          CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEF_L, 
     &           A(POSELT), LAFAC, MonBloc,
     &           NextPiv2beWritten, IDUMMY,
     &           IW(IOLDPS), LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
          IF (IFLAG_OOC .LT. 0 ) THEN
            IFLAG = IFLAG_OOC
            GOTO 500
          ENDIF
          CALL DMUMPS_OOC_PP_TRYRELEASE_SPACE (IWPOS, 
     &      IOLDPS, IW, LIW, MonBloc , NFRONT, KEEP)
      ENDIF
      GOTO 600
 490  CONTINUE
 500  CONTINUE
 600  CONTINUE
      IF (LR_ACTIVATED) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK))  DEALLOCATE(WORK)
        IF (allocated(TAU))   deallocate(TAU)
        IF (allocated(JPVT))  deallocate(JPVT)
        IF (allocated(BLOCK)) deallocate(BLOCK)
        IF (associated(ACC_LUA)) THEN
          IF (KEEP(480).GE.3) THEN
            DO MY_NUM=1,OMP_NUM
              CALL DEALLOC_LRB(ACC_LUA(MY_NUM), KEEP8, KEEP(34))
            ENDDO
          ENDIF
          DEALLOCATE(ACC_LUA)
          NULLIFY(ACC_LUA)
        ENDIF
        IF (associated(BEGS_BLR)) THEN
          DEALLOCATE(BEGS_BLR)
          NULLIFY(BEGS_BLR)
        ENDIF
      ENDIF
      IF (LR_ACTIVATED.AND.KEEP(480).NE.0) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &     ) THEN
          CALL DMUMPS_BLR_FREE_ALL_PANELS(IW(IOLDPS+XXF), 0, 
     &                        KEEP8, KEEP(34))
        ENDIF
      ENDIF
      IF (LR_ACTIVATED) THEN
        IF (.NOT.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &    .AND. .NOT.COMPRESS_CB) THEN
          CALL DMUMPS_BLR_END_FRONT(IW(IOLDPS+XXF), IFLAG, KEEP8,
     &                              KEEP(34), MTK405=KEEP(405))
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC1_LDLT
      END MODULE DMUMPS_FAC1_LDLT_M
      SUBROUTINE DMUMPS_FAC1_LDLT_I( N, INODE, IW, LIW, A, LA,
     &                     IOLDPS, POSELT, IFLAG, IERROR,
     &                     UU, NNEGW, NNULLNEGW, NPVW, NB22T1W, NBTINYW,
     &                     DET_EXPW, DET_MANTW, DET_SIGNW,
     &                     KEEP,KEEP8,
     &                     MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP, PIVNUL_LIST_STRUCT, IWPOS
     &     , LRGROUPS
     &     , PERM
     &     )
      USE MUMPS_PIVNUL_MOD, ONLY: PIVNUL_LIST_STRUCT_T
      USE DMUMPS_FAC1_LDLT_M, ONLY: DMUMPS_FAC1_LDLT
      IMPLICIT NONE
      INTEGER(8) :: LA, POSELT
      INTEGER N, INODE, LIW, IFLAG, IERROR
      INTEGER, intent(inout) :: NNEGW, NNULLNEGW, NPVW, NB22T1W, NBTINYW
      INTEGER, intent(inout) :: DET_EXPW, DET_SIGNW
      DOUBLE PRECISION, intent(inout) :: DET_MANTW
      INTEGER MYID, IOLDPS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A( LA )
      INTEGER IW( LIW )
      INTEGER, intent(in) :: PERM(N)
      LOGICAL AVOID_DELAYED
      INTEGER ETATASS, IWPOS
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION DKEEP(230)
      INTEGER :: LRGROUPS(KEEP(280))
      CALL DMUMPS_FAC1_LDLT( N, INODE, IW, LIW, A, LA,
     &                     IOLDPS, POSELT, IFLAG, IERROR,
     &                     UU, NNEGW, NNULLNEGW, NPVW, NB22T1W, NBTINYW,
     &                     DET_EXPW, DET_MANTW, DET_SIGNW,
     &                     KEEP,KEEP8,
     &                     MYID, SEUIL, AVOID_DELAYED, ETATASS,
     &     DKEEP, PIVNUL_LIST_STRUCT, IWPOS
     &     , LRGROUPS
     &     , PERM
     &     )
      RETURN
      END SUBROUTINE DMUMPS_FAC1_LDLT_I
