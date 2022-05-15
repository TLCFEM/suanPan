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
      MODULE DMUMPS_LR_STATS
      USE DMUMPS_LR_TYPE
      IMPLICIT NONE
      DOUBLE PRECISION :: ACC_MRY_CB_GAIN, 
     &                    ACC_MRY_CB_FR,
     &                    FRONT_L11_BLR_SAVINGS,  
     &                    FRONT_U11_BLR_SAVINGS, 
     &                    FRONT_L21_BLR_SAVINGS,  
     &                    FRONT_U12_BLR_SAVINGS,
     &                    ACC_FR_MRY,
     &                    GLOBAL_BLR_SAVINGS,
     &                    GLOBAL_MRY_LPRO_COMPR,
     &                    GLOBAL_MRY_LTOT_COMPR
      INTEGER :: CNT_NODES
      DOUBLE PRECISION :: FLOP_FR_UPDT,
     &                    FLOP_LR_UPDT,
     &                    FLOP_LR_UPDT_OUT,
     &                    FLOP_RMB,
     &                    FLOP_FR_TRSM,
     &                    FLOP_LR_TRSM,
     &                    FLOP_PANEL,
     &                    FLOP_TRSM,
     &                    FLOP_DEC_ACC,
     &                    FLOP_REC_ACC,
     &                    FLOP_FR_SWAP,
     &                    FLOP_DEMOTE,
     &                    FLOP_CB_DEMOTE,
     &                    FLOP_CB_PROMOTE,
     &                    FLOP_PROMOTE,
     &                    LR_FLOP_GAIN 
      DOUBLE PRECISION :: ACC_LR_FLOP_GAIN
      DOUBLE PRECISION :: ACC_FLOP_FR_FACTO,
     &                    ACC_FLOP_LR_FACTO,
     &                    ACC_FLOP_FR_TRSM,
     &                    ACC_FLOP_LR_TRSM,
     &                    ACC_FLOP_FR_UPDT,
     &                    ACC_FLOP_LR_UPDT,
     &                    ACC_FLOP_LR_UPDT_OUT,
     &                    ACC_FLOP_RMB,
     &                    ACC_FLOP_DEMOTE,
     &                    ACC_FLOP_CB_DEMOTE,
     &                    ACC_FLOP_CB_PROMOTE,
     &                    ACC_FLOP_PROMOTE,
     &                    ACC_FLOP_TRSM,
     &                    ACC_FLOP_DEC_ACC,
     &                    ACC_FLOP_REC_ACC,
     &                    ACC_FLOP_FR_SWAP,
     &                    ACC_FLOP_PANEL,
     &                    ACC_FLOP_FRFRONTS,
     &                    ACC_FLOP_FR_SOLVE,
     &                    ACC_FLOP_LR_SOLVE
      DOUBLE PRECISION :: FACTOR_PROCESSED_FRACTION
      INTEGER(KIND=8)  :: FACTOR_SIZE
      DOUBLE PRECISION :: TOTAL_FLOP
      DOUBLE PRECISION :: BLR_TIME_LRGROUPING
      DOUBLE PRECISION :: BLR_TIME_SEPGROUPING
      DOUBLE PRECISION :: BLR_TIME_GETHALO
      DOUBLE PRECISION :: BLR_TIME_KWAY
      DOUBLE PRECISION :: BLR_TIME_GNEW
      DOUBLE PRECISION :: ACC_UPDT_TIME
      DOUBLE PRECISION :: ACC_RMB_TIME
      DOUBLE PRECISION :: ACC_UPDT_TIME_OUT
      DOUBLE PRECISION :: ACC_PROMOTING_TIME
      DOUBLE PRECISION :: ACC_PROMOTING_UCFS_TIME
      DOUBLE PRECISION :: ACC_DEMOTING_TIME
      DOUBLE PRECISION :: ACC_CB_DEMOTING_TIME
      DOUBLE PRECISION :: ACC_LR_MODULE_TIME
      DOUBLE PRECISION :: ACC_TRSM_TIME
      DOUBLE PRECISION :: ACC_FRPANELS_TIME
      DOUBLE PRECISION :: ACC_FAC_I_TIME
      DOUBLE PRECISION :: ACC_FAC_MQ_TIME
      DOUBLE PRECISION :: ACC_FAC_SQ_TIME
      DOUBLE PRECISION :: ACC_FRFRONTS_TIME
      DOUBLE PRECISION :: ACC_DIAGCOPY_TIME
      DOUBLE PRECISION :: ACC_FRSWAP_TIME
      DOUBLE PRECISION :: ACC_DEC_ASM1_TIME
      DOUBLE PRECISION :: ACC_DEC_LOCASM2_TIME
      DOUBLE PRECISION :: ACC_DEC_MAPLIG1_TIME
      DOUBLE PRECISION :: ACC_DEC_ASMS2S_TIME
      DOUBLE PRECISION :: ACC_DEC_ASMS2M_TIME
      DOUBLE PRECISION :: AVG_ACC_FLOP_LR_FACTO
      DOUBLE PRECISION :: MIN_ACC_FLOP_LR_FACTO
      DOUBLE PRECISION :: MAX_ACC_FLOP_LR_FACTO
      INTEGER :: TOTAL_NBLOCKS_ASS, TOTAL_NBLOCKS_CB
      INTEGER :: MIN_BLOCKSIZE_ASS, MAX_BLOCKSIZE_ASS
      INTEGER :: MIN_BLOCKSIZE_CB, MAX_BLOCKSIZE_CB
      DOUBLE PRECISION :: AVG_BLOCKSIZE_ASS, AVG_BLOCKSIZE_CB
      INTEGER, POINTER :: STEP_STATS(:)
      CONTAINS
      SUBROUTINE COLLECT_BLOCKSIZES(CUT,NPARTSASS,NPARTSCB)
        INTEGER, INTENT(IN) :: NPARTSASS, NPARTSCB
        INTEGER, POINTER, DIMENSION(:) :: CUT
        INTEGER :: LOC_MIN_ASS, LOC_MIN_CB, LOC_MAX_ASS, LOC_MAX_CB,
     &             LOC_TOT_ASS, LOC_TOT_CB
        DOUBLE PRECISION :: LOC_AVG_ASS, LOC_AVG_CB 
        INTEGER :: I
        LOC_TOT_ASS = 0
        LOC_TOT_CB = 0
        LOC_AVG_ASS = 0.D0
        LOC_AVG_CB = 0.D0
        LOC_MIN_ASS = 100000
        LOC_MIN_CB = 100000
        LOC_MAX_ASS = 0
        LOC_MAX_CB = 0
        DO I = 1,NPARTSASS
          LOC_AVG_ASS = ( LOC_TOT_ASS * LOC_AVG_ASS
     &                        + CUT(I+1) - CUT(I) )
     &                        / (LOC_TOT_ASS + 1)
          LOC_TOT_ASS = LOC_TOT_ASS + 1
          IF (CUT(I+1) - CUT(I) .LE. LOC_MIN_ASS) THEN
            LOC_MIN_ASS = CUT(I+1) - CUT(I)
          END IF
          IF (CUT(I+1) - CUT(I) .GE. LOC_MAX_ASS) THEN
            LOC_MAX_ASS = CUT(I+1) - CUT(I)
          END IF
        END DO
        DO I = NPARTSASS+1,NPARTSASS+NPARTSCB
          LOC_AVG_CB = ( LOC_TOT_CB * LOC_AVG_CB
     &                        + CUT(I+1) - CUT(I) )
     &                        / (LOC_TOT_CB + 1)
          LOC_TOT_CB = LOC_TOT_CB + 1
          IF (CUT(I+1) - CUT(I) .LE. LOC_MIN_CB) THEN
            LOC_MIN_CB = CUT(I+1) - CUT(I)
          END IF
          IF (CUT(I+1) - CUT(I) .GE. LOC_MAX_CB) THEN
            LOC_MAX_CB = CUT(I+1) - CUT(I)
          END IF
        END DO
        AVG_BLOCKSIZE_ASS = (TOTAL_NBLOCKS_ASS*AVG_BLOCKSIZE_ASS
     &     + LOC_TOT_ASS*LOC_AVG_ASS) / (TOTAL_NBLOCKS_ASS+LOC_TOT_ASS)
        AVG_BLOCKSIZE_CB = (TOTAL_NBLOCKS_CB*AVG_BLOCKSIZE_CB
     &     + LOC_TOT_CB*LOC_AVG_CB) / (TOTAL_NBLOCKS_CB+LOC_TOT_CB)
        TOTAL_NBLOCKS_ASS = TOTAL_NBLOCKS_ASS + LOC_TOT_ASS
        TOTAL_NBLOCKS_CB = TOTAL_NBLOCKS_CB + LOC_TOT_CB
        MIN_BLOCKSIZE_ASS = min(MIN_BLOCKSIZE_ASS,LOC_MIN_ASS)
        MIN_BLOCKSIZE_CB = min(MIN_BLOCKSIZE_CB,LOC_MIN_CB)
        MAX_BLOCKSIZE_ASS = max(MAX_BLOCKSIZE_ASS,LOC_MAX_ASS)
        MAX_BLOCKSIZE_CB = max(MAX_BLOCKSIZE_CB,LOC_MAX_CB)
      END SUBROUTINE COLLECT_BLOCKSIZES
      SUBROUTINE UPDATE_ALL_TIMES(INODE, LOC_FACTO_TIME,
     &      LOC_PROMOTING_TIME, LOC_DEMOTING_TIME, LOC_CB_DEMOTING_TIME,
     &      LOC_FRPANELS_TIME, LOC_FRFRONTS_TIME,
     &      LOC_TRSM_TIME, LOC_LR_MODULE_TIME,
     &      LOC_FAC_I_TIME, LOC_FAC_MQ_TIME, LOC_FAC_SQ_TIME)
          INTEGER, INTENT(IN) :: INODE
          DOUBLE PRECISION, INTENT(IN) :: LOC_FACTO_TIME, 
     &            LOC_PROMOTING_TIME, LOC_DEMOTING_TIME,
     &            LOC_CB_DEMOTING_TIME, LOC_FRPANELS_TIME,
     &            LOC_FRFRONTS_TIME, LOC_TRSM_TIME, LOC_LR_MODULE_TIME,
     &            LOC_FAC_I_TIME, LOC_FAC_MQ_TIME, LOC_FAC_SQ_TIME
        ACC_UPDT_TIME = ACC_UPDT_TIME + LOC_FACTO_TIME
        ACC_PROMOTING_TIME = ACC_PROMOTING_TIME + LOC_PROMOTING_TIME
        ACC_DEMOTING_TIME = ACC_DEMOTING_TIME + LOC_DEMOTING_TIME
        ACC_CB_DEMOTING_TIME = ACC_CB_DEMOTING_TIME + 
     &      LOC_CB_DEMOTING_TIME
        ACC_FRPANELS_TIME = ACC_FRPANELS_TIME + LOC_FRPANELS_TIME
        ACC_FAC_I_TIME = ACC_FAC_I_TIME + LOC_FAC_I_TIME
        ACC_FAC_MQ_TIME = ACC_FAC_MQ_TIME + LOC_FAC_MQ_TIME
        ACC_FAC_SQ_TIME = ACC_FAC_SQ_TIME + LOC_FAC_SQ_TIME
        ACC_FRFRONTS_TIME = ACC_FRFRONTS_TIME + LOC_FRFRONTS_TIME
        ACC_TRSM_TIME = ACC_TRSM_TIME + LOC_TRSM_TIME
        ACC_LR_MODULE_TIME = ACC_LR_MODULE_TIME + LOC_LR_MODULE_TIME
      END SUBROUTINE UPDATE_ALL_TIMES
      SUBROUTINE UPDATE_CB_DEMOTING_TIME(INODE, LOC_CB_DEMOTING_TIME)
          INTEGER, INTENT(IN) :: INODE
          DOUBLE PRECISION, INTENT(IN) :: LOC_CB_DEMOTING_TIME
          ACC_CB_DEMOTING_TIME = ACC_CB_DEMOTING_TIME + 
     &      LOC_CB_DEMOTING_TIME
      END SUBROUTINE UPDATE_CB_DEMOTING_TIME
      SUBROUTINE UPDATE_UPDT_TIME(INODE, LOC_UPDT_TIME)
          INTEGER, INTENT(IN) :: INODE
          DOUBLE PRECISION, INTENT(IN) :: LOC_UPDT_TIME
          ACC_UPDT_TIME = ACC_UPDT_TIME + LOC_UPDT_TIME
      END SUBROUTINE UPDATE_UPDT_TIME
      SUBROUTINE UPDATE_UPDT_TIME_OUT(LOC_UPDT_TIME_OUT)
          DOUBLE PRECISION, INTENT(IN) :: LOC_UPDT_TIME_OUT
          ACC_UPDT_TIME_OUT = ACC_UPDT_TIME_OUT + LOC_UPDT_TIME_OUT
      END SUBROUTINE UPDATE_UPDT_TIME_OUT
      SUBROUTINE UPDATE_RMB_TIME(LOC_RMB_TIME)
          DOUBLE PRECISION, INTENT(IN) :: LOC_RMB_TIME
          ACC_RMB_TIME = ACC_RMB_TIME + LOC_RMB_TIME
      END SUBROUTINE UPDATE_RMB_TIME
      SUBROUTINE UPDATE_PROMOTING_TIME(INODE, LOC_PROMOTING_TIME)
          INTEGER, INTENT(IN) :: INODE
          DOUBLE PRECISION, INTENT(IN) :: LOC_PROMOTING_TIME
          ACC_PROMOTING_TIME = ACC_PROMOTING_TIME + 
     &      LOC_PROMOTING_TIME
      END SUBROUTINE UPDATE_PROMOTING_TIME
      SUBROUTINE UPDATE_FLOP_STATS_CB_PROMOTE(COST, NIV)
        DOUBLE PRECISION :: COST
        INTEGER :: NIV
        IF (NIV.EQ.1) THEN
!$OMP CRITICAL(cb_flop_cost_pro_cri)
                  FLOP_CB_PROMOTE = FLOP_CB_PROMOTE + COST
!$OMP END CRITICAL(cb_flop_cost_pro_cri)
        ELSE
!$OMP CRITICAL(acc_cb_flop_cost_pro_cri)
          ACC_FLOP_CB_PROMOTE = ACC_FLOP_CB_PROMOTE + COST
!$OMP END CRITICAL(acc_cb_flop_cost_pro_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_CB_PROMOTE
      SUBROUTINE UPDATE_FLOP_STATS_PROMOTE(COST, NIV)
        DOUBLE PRECISION :: COST
        INTEGER :: NIV
        IF (NIV.EQ.1) THEN
!$OMP CRITICAL(flop_cost_pro_cri)
                  FLOP_PROMOTE = FLOP_PROMOTE + COST
!$OMP END CRITICAL(flop_cost_pro_cri)
        ELSE
!$OMP CRITICAL(acc_flop_cost_pro_cri)
          ACC_FLOP_PROMOTE = ACC_FLOP_PROMOTE + COST
!$OMP END CRITICAL(acc_flop_cost_pro_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_PROMOTE
      SUBROUTINE UPDATE_FLOP_STATS_CB_DEMOTE(COST, NIV)
        DOUBLE PRECISION :: COST
        INTEGER :: NIV
        IF (NIV.EQ.1) THEN
!$OMP CRITICAL(cb_flop_cost_dem_cri)
                  FLOP_CB_DEMOTE = FLOP_CB_DEMOTE + COST
                  FLOP_DEMOTE    = FLOP_DEMOTE    + COST
!$OMP END CRITICAL(cb_flop_cost_dem_cri)
        ELSE
!$OMP CRITICAL(acc_cb_flop_cost_dem_cri)
          ACC_FLOP_CB_DEMOTE = ACC_FLOP_CB_DEMOTE + COST
          ACC_FLOP_DEMOTE    = ACC_FLOP_DEMOTE    + COST
!$OMP END CRITICAL(acc_cb_flop_cost_dem_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_CB_DEMOTE
      SUBROUTINE UPDATE_FLOP_STATS_DEMOTE(LR_B, NIV, REC_ACC, 
     &         CB_COMPRESS, FR_SWAP)
        TYPE(LRB_TYPE),INTENT(IN) :: LR_B
        INTEGER(8) :: M,N,K
        INTEGER :: NIV
        DOUBLE PRECISION :: HR_COST,BUILDQ_COST
        LOGICAL, OPTIONAL :: REC_ACC, CB_COMPRESS, FR_SWAP
        M = int(LR_B%M,8)
        N = int(LR_B%N,8)
        K = int(LR_B%K,8)
        HR_COST =  dble(4_8*K*K*K/3_8 + 4_8*K*M*N - 2_8*(M+N)*K*K)
        IF (LR_B%ISLR) THEN 
          BUILDQ_COST = dble(4_8*K*K*M - K*K*K)
        ELSE 
          BUILDQ_COST = 0.0d0
        END IF
        IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
                  FLOP_DEMOTE = FLOP_DEMOTE + HR_COST + BUILDQ_COST
                  IF (present(REC_ACC)) THEN
                    IF (REC_ACC) THEN
                      FLOP_REC_ACC = FLOP_REC_ACC + HR_COST+BUILDQ_COST
                    ENDIF
                  ENDIF
                  IF (present(CB_COMPRESS)) THEN
                    IF (CB_COMPRESS) THEN
                      FLOP_CB_DEMOTE = FLOP_CB_DEMOTE + HR_COST 
     &                                                + BUILDQ_COST
                    ENDIF
                  ENDIF
                  IF (present(FR_SWAP)) THEN
                    IF (FR_SWAP) THEN
                      FLOP_FR_SWAP = FLOP_FR_SWAP + HR_COST+BUILDQ_COST
                    ENDIF
                  ENDIF
!$OMP END CRITICAL(lr_flop_gain_cri)
        ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
           ACC_FLOP_DEMOTE = ACC_FLOP_DEMOTE + (HR_COST + BUILDQ_COST) 
           IF (present(REC_ACC)) THEN
             IF (REC_ACC) THEN
               ACC_FLOP_REC_ACC = ACC_FLOP_REC_ACC +HR_COST+BUILDQ_COST
             ENDIF
           ENDIF
           IF (present(CB_COMPRESS)) THEN
             IF (CB_COMPRESS) THEN
               ACC_FLOP_CB_DEMOTE = ACC_FLOP_CB_DEMOTE + HR_COST 
     &                                                 + BUILDQ_COST
             ENDIF
           ENDIF
           IF (present(FR_SWAP)) THEN
             IF (FR_SWAP) THEN
               ACC_FLOP_FR_SWAP = ACC_FLOP_FR_SWAP +HR_COST+BUILDQ_COST
             ENDIF
           ENDIF
!$OMP END CRITICAL(lr_flop_gain_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_DEMOTE
      SUBROUTINE UPDATE_FLOP_STATS_REC_ACC(LR_B, NIV, K1, K2, BUILDQ1)
        TYPE(LRB_TYPE),INTENT(IN) :: LR_B
        INTEGER,INTENT(IN) :: NIV, K1, K2
        LOGICAL,INTENT(IN) :: BUILDQ1
        INTEGER(8) :: M,N,K
        DOUBLE PRECISION :: HR_COST, BUILDQ_COST, GS_COST, UPDT_COST,
     &                      TOT_COST
        M = int(LR_B%M,8)
        N = int(LR_B%N,8)
        K = int(LR_B%K - K1,8)
        GS_COST = dble((4_8*(K1)+1_8)*M*K2) 
        HR_COST = dble(4_8*K*K*K/3_8 + 4_8*K*M*K2 - 2_8*(M+K2)*K*K)
        IF (BUILDQ1) THEN
          BUILDQ_COST = dble(4_8*K*K*M - K*K*K)
          UPDT_COST = dble(2_8*K*K2*N)
        ELSE
          BUILDQ_COST = 0.0d0
          UPDT_COST = 0.0d0
        ENDIF
        TOT_COST = BUILDQ_COST + HR_COST + GS_COST + UPDT_COST
        IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
                  FLOP_DEMOTE  = FLOP_DEMOTE  + TOT_COST
                  FLOP_REC_ACC = FLOP_REC_ACC + TOT_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
        ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
           ACC_FLOP_DEMOTE  = ACC_FLOP_DEMOTE  + TOT_COST
           ACC_FLOP_REC_ACC = ACC_FLOP_REC_ACC + TOT_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_REC_ACC
      SUBROUTINE UPDATE_FLOP_STATS_PANEL(NFRONT, NPIV, NIV, SYM)
          INTEGER :: NFRONT, NPIV, NIV, SYM
          DOUBLE PRECISION :: COST_PANEL, COST_TRSM
          IF (SYM.EQ.0) THEN
            COST_TRSM  = dble(2 * NPIV-1) * dble(NPIV) 
     &                                    * dble(NFRONT-NPIV)
            COST_PANEL = dble(NPIV) * dble(NPIV - 1) 
     &                              * dble(4 * NPIV + 1)/dble(6)
          ELSE
            COST_TRSM  = dble(NPIV) * dble(NPIV) * dble(NFRONT-NPIV)
            COST_PANEL = dble(NPIV) * dble(NPIV - 1) 
     &                              * dble(2 * NPIV + 1)/dble(6)
          ENDIF
          IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
             FLOP_PANEL = FLOP_PANEL + COST_PANEL 
             FLOP_TRSM  = FLOP_TRSM  + COST_TRSM
!$OMP END CRITICAL(lr_flop_gain_cri)
          ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
             ACC_FLOP_PANEL = ACC_FLOP_PANEL + COST_PANEL
             ACC_FLOP_TRSM  = ACC_FLOP_TRSM  + COST_TRSM
!$OMP END CRITICAL(lr_flop_gain_cri)
          ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_PANEL
      SUBROUTINE UPDATE_FLOP_STATS_TRSM(LRB, NIV, LorU)
          TYPE(LRB_TYPE),INTENT(IN) :: LRB
          INTEGER,INTENT(IN) :: NIV, LorU
          DOUBLE PRECISION :: LR_FLOP_COST, FR_FLOP_COST
          IF (LorU.EQ.0) THEN 
            FR_FLOP_COST = dble(LRB%M)*dble(LRB%N)*dble(LRB%N)
            IF (LRB%ISLR) THEN
              LR_FLOP_COST = dble(LRB%K)*dble(LRB%N)*dble(LRB%N)
            ELSE
              LR_FLOP_COST = FR_FLOP_COST
            ENDIF
          ELSE 
            FR_FLOP_COST = dble(LRB%M-1)*dble(LRB%N)*dble(LRB%N)
            IF (LRB%ISLR) THEN
              LR_FLOP_COST = dble(LRB%N-1)*dble(LRB%N)*dble(LRB%K)
            ELSE
              LR_FLOP_COST = FR_FLOP_COST
            ENDIF
          ENDIF
          IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
                  FLOP_FR_TRSM = FLOP_FR_TRSM + FR_FLOP_COST
                  FLOP_LR_TRSM = FLOP_LR_TRSM + LR_FLOP_COST
                  LR_FLOP_GAIN = LR_FLOP_GAIN + FR_FLOP_COST
     &                                    - LR_FLOP_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
          ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
              ACC_FLOP_FR_TRSM = ACC_FLOP_FR_TRSM + FR_FLOP_COST
              ACC_FLOP_LR_TRSM = ACC_FLOP_LR_TRSM + LR_FLOP_COST
              ACC_LR_FLOP_GAIN = ACC_LR_FLOP_GAIN + FR_FLOP_COST 
     &                                            - LR_FLOP_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
          END IF
      END SUBROUTINE UPDATE_FLOP_STATS_TRSM
      SUBROUTINE UPDATE_FLOP_STATS_LRB_PRODUCT(LRB1, LRB2, TRANSB1,
     &      TRANSB2, NIV, COMPRESS_MID_PRODUCT, RANK_IN, BUILDQ,
     &      IS_DIAG, K480, REC_ACC_IN)
!$    USE OMP_LIB
        TYPE(LRB_TYPE),INTENT(IN) :: LRB1,LRB2
        CHARACTER(len=1), INTENT(IN) :: TRANSB1, TRANSB2
        LOGICAL, INTENT(IN), OPTIONAL :: BUILDQ, IS_DIAG, REC_ACC_IN
        INTEGER, INTENT(IN), OPTIONAL :: NIV, RANK_IN,
     &                    COMPRESS_MID_PRODUCT, K480
        LOGICAL :: REC_ACC
        DOUBLE PRECISION :: LR_FLOP_COST, LR_FLOP_COST_OUT, FR_FLOP_COST
        DOUBLE PRECISION :: HR_COST, BUILDQ_COST
        DOUBLE PRECISION :: M1,N1,K1,M2,N2,K2,RANK
        CHARACTER(len=2) :: PROD, TRANS
         IF(present(K480).AND.present(REC_ACC_IN)) THEN
           IF (K480.GE.4) THEN
             REC_ACC = REC_ACC_IN
           ELSE
             REC_ACC = .FALSE.
           ENDIF
         ELSE
           REC_ACC = .FALSE.
         ENDIF
        M1 = dble(LRB1%M)
        N1 = dble(LRB1%N)
        K1 = dble(LRB1%K)
        M2 = dble(LRB2%M)
        N2 = dble(LRB2%N)
        K2 = dble(LRB2%K)
        RANK = dble(RANK_IN)
        IF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==0)) THEN 
          PROD = '00'
        ELSE IF ((LRB1%LRFORM==1).AND.(LRB2%LRFORM==0)) THEN
          PROD = '10'
        ELSE IF ((LRB1%LRFORM==0).AND.(LRB2%LRFORM==1)) THEN
          PROD = '01'
        ELSE 
          PROD = '11'
        END IF
        IF ((TRANSB1=='N').AND.(TRANSB2=='N')) THEN
          TRANS = 'NN'
        ELSE IF ((TRANSB1=='T').AND.(TRANSB2=='N')) THEN
          TRANS = 'TN'
        ELSE IF ((TRANSB1=='N').AND.(TRANSB2=='T')) THEN
          TRANS = 'NT'
        ELSE 
          TRANS = 'TT'
        END IF
        LR_FLOP_COST_OUT = 0.0D0
        HR_COST = 0.0D0
        BUILDQ_COST = 0.0D0
        SELECT CASE (PROD)  
        CASE('00') 
          SELECT CASE (TRANS)
          CASE('NN')
        FR_FLOP_COST = 2.0D0*M1*N2*N1
        LR_FLOP_COST = 2.0D0*M1*N2*N1
          CASE('TN')
        FR_FLOP_COST = 2.0D0*N1*N2*M1
        LR_FLOP_COST = 2.0D0*M1*N2*N1
          CASE('NT')
        FR_FLOP_COST = 2.0D0*M1*M2*N1
        LR_FLOP_COST = 2.0D0*M1*M2*N1
          CASE('TT')
        FR_FLOP_COST = 2.0D0*N1*M2*M1
        LR_FLOP_COST = 2.0D0*N1*M2*M1
          END SELECT
        CASE('10')
          SELECT CASE (TRANS)
          CASE('NN')
        FR_FLOP_COST = 2.0D0*M1*N2*N1
        LR_FLOP_COST = 2.0D0*K1*N2*N1 + 2.0D0*M1*N2*K1
        LR_FLOP_COST_OUT = 2.0D0*M1*N2*K1
          CASE('TN')
        FR_FLOP_COST = 2.0D0*N1*N2*M1
        LR_FLOP_COST = 2.0D0*K1*N2*M1 + 2.0D0*N1*N2*K1
        LR_FLOP_COST_OUT = 2.0D0*N1*N2*K1
          CASE('NT')
        FR_FLOP_COST = 2.0D0*M1*M2*N1
        LR_FLOP_COST = 2.0D0*K1*M2*N1 + 2.0D0*M1*M2*K1
        LR_FLOP_COST_OUT = 2.0D0*M1*M2*K1
          CASE('TT')
        FR_FLOP_COST = 2.0D0*N1*M2*M1
        LR_FLOP_COST = 2.0D0*K1*M2*M1 + 2.0D0*N1*M2*K1
        LR_FLOP_COST_OUT = 2.0D0*N1*M2*K1
          END SELECT
        CASE('01')
          SELECT CASE (TRANS)
          CASE('NN')
        FR_FLOP_COST = 2.0D0*M1*N2*N1
        LR_FLOP_COST = 2.0D0*M1*K2*N1 + 2.0D0*M1*N2*K2
        LR_FLOP_COST_OUT = 2.0D0*M1*N2*K2
          CASE('TN')
        FR_FLOP_COST = 2.0D0*N1*N2*M1
        LR_FLOP_COST = 2.0D0*N1*K2*M1 + 2.0D0*N1*N2*K2
        LR_FLOP_COST_OUT = 2.0D0*N1*N2*K2
          CASE('NT')
        FR_FLOP_COST = 2.0D0*M1*M2*N1
        LR_FLOP_COST = 2.0D0*M1*K2*N1 + 2.0D0*M1*M2*K2
        LR_FLOP_COST_OUT = 2.0D0*M1*M2*K2
          CASE('TT')
        FR_FLOP_COST = 2*N1*M2*M1
        LR_FLOP_COST = 2.0D0*N1*K2*M1 + 2.0D0*N1*M2*K2
        LR_FLOP_COST_OUT = 2.0D0*N1*M2*K2
          END SELECT
        CASE('11')
        IF (COMPRESS_MID_PRODUCT.GE.1) THEN
            HR_COST =  4.0D0*RANK*RANK*RANK/3.0D0 + 
     &                 4.0D0*RANK*K1*K2 - 
     &                 2.0D0*(K1+K2)*RANK*RANK    
            IF (BUILDQ) THEN
                BUILDQ_COST = 4.0D0*RANK*RANK*K1 - RANK*RANK*RANK
            ENDIF
        ENDIF
          SELECT CASE (TRANS)
          CASE('NN')
        FR_FLOP_COST = 2.0D0*M1*N2*N1
        IF ((COMPRESS_MID_PRODUCT.GE.1).AND.BUILDQ) THEN
            LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                     2.0D0*K1*M1*RANK + 2.0D0*K2*N2*RANK +
     &                     2.0D0*M1*N2*RANK
            LR_FLOP_COST_OUT = 2.0D0*M1*N2*RANK
        ELSE
            IF (K1 .GE. K2) THEN
                LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                         2.0D0*K1*M1*K2 + 2.0D0*M1*N2*K2
                LR_FLOP_COST_OUT = 2.0D0*M1*N2*K2
            ELSE
                LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                         2.0D0*K1*N2*K2 + 2.0D0*M1*N2*K1
                LR_FLOP_COST_OUT = 2.0D0*M1*N2*K1
            ENDIF
        ENDIF
          CASE('TN')
        FR_FLOP_COST = 2.0D0*N1*N2*M1
        IF ((COMPRESS_MID_PRODUCT.GE.1).AND.BUILDQ) THEN
            LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                     2.0D0*K1*N1*RANK + 2.0D0*K2*N2*RANK +
     &                     2.0D0*N1*N2*RANK
            LR_FLOP_COST_OUT = 2.0D0*N1*N2*RANK
        ELSE
            IF (K1 .GE. K2) THEN
                LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                         2.0D0*K1*N1*K2 + 2.0D0*N1*N2*K2
                LR_FLOP_COST_OUT = 2.0D0*N1*N2*K2
            ELSE
                LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                         2.0D0*K1*N2*K2 + 2.0D0*N1*N2*K1
                LR_FLOP_COST_OUT = 2.0D0*N1*N2*K1
            ENDIF
        ENDIF
          CASE('NT')
        FR_FLOP_COST = 2.0D0*M1*M2*N1
        IF ((COMPRESS_MID_PRODUCT.GE.1).AND.BUILDQ) THEN
            LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                     2.0D0*K1*M1*RANK + 2.0D0*K2*M2*RANK +
     &                     2.0D0*M1*M2*RANK
            LR_FLOP_COST_OUT = 2.0D0*M1*M2*RANK
        ELSE
            IF (K1 .GE. K2) THEN
                LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                         2.0D0*K1*M1*K2 + 2.0D0*M1*M2*K2
                LR_FLOP_COST_OUT = 2.0D0*M1*M2*K2
            ELSE
                LR_FLOP_COST = 2.0D0*K1*K2*N1 + 
     &                         2.0D0*K1*M2*K2 + 2.0D0*M1*M2*K1
                LR_FLOP_COST_OUT = 2.0D0*M1*M2*K1
            ENDIF
        ENDIF
          CASE('TT')
        FR_FLOP_COST = 2.0D0*N1*M2*M1
        IF ((COMPRESS_MID_PRODUCT.GE.1).AND.BUILDQ) THEN
            LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                     2.0D0*K1*N1*RANK + 2.0D0*K2*M2*RANK +
     &                     2.0D0*N1*M2*RANK
            LR_FLOP_COST_OUT = 2.0D0*N1*M2*RANK
        ELSE
            IF (K1 .GE. K2) THEN
                LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                         2.0D0*K1*N1*K2 + 2.0D0*N1*M2*K2
                LR_FLOP_COST_OUT = 2.0D0*N1*M2*K2
            ELSE
                LR_FLOP_COST = 2.0D0*K1*K2*M1 + 
     &                         2.0D0*K1*M2*K2 + 2.0D0*N1*M2*K1
                LR_FLOP_COST_OUT = 2.0D0*N1*M2*K1
            ENDIF
        ENDIF
          END SELECT
        END SELECT
        IF (present(IS_DIAG)) THEN
          IF (IS_DIAG) THEN
            FR_FLOP_COST = FR_FLOP_COST/2.0D0
            LR_FLOP_COST = LR_FLOP_COST/2.0D0
          ENDIF
        ENDIF
        IF (present(K480)) THEN
          IF (K480.GE.3) THEN
            LR_FLOP_COST = LR_FLOP_COST - LR_FLOP_COST_OUT
            LR_FLOP_COST_OUT = 0.0D0
            IF (REC_ACC) THEN
              IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
                  FLOP_REC_ACC = FLOP_REC_ACC + LR_FLOP_COST
     &                           + HR_COST + BUILDQ_COST
                  FLOP_DEMOTE  = FLOP_DEMOTE + LR_FLOP_COST
     &                           + HR_COST + BUILDQ_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
              ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
                ACC_FLOP_REC_ACC = ACC_FLOP_REC_ACC + LR_FLOP_COST
     &                           + HR_COST + BUILDQ_COST
                ACC_FLOP_DEMOTE  = ACC_FLOP_DEMOTE + LR_FLOP_COST
     &                           + HR_COST + BUILDQ_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF (.NOT.REC_ACC) THEN
          IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
              LR_FLOP_GAIN = LR_FLOP_GAIN + FR_FLOP_COST - LR_FLOP_COST
              FLOP_FR_UPDT = FLOP_FR_UPDT + FR_FLOP_COST
              FLOP_LR_UPDT = FLOP_LR_UPDT + LR_FLOP_COST
              FLOP_LR_UPDT_OUT = FLOP_LR_UPDT_OUT + LR_FLOP_COST_OUT
              FLOP_DEMOTE  = FLOP_DEMOTE  + HR_COST + BUILDQ_COST
              FLOP_RMB     = FLOP_RMB     + HR_COST + BUILDQ_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
          ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
            ACC_LR_FLOP_GAIN = ACC_LR_FLOP_GAIN + 
     &                     FR_FLOP_COST - LR_FLOP_COST 
            ACC_FLOP_FR_UPDT = ACC_FLOP_FR_UPDT + FR_FLOP_COST
            ACC_FLOP_LR_UPDT = ACC_FLOP_LR_UPDT + LR_FLOP_COST
            ACC_FLOP_LR_UPDT_OUT = ACC_FLOP_LR_UPDT_OUT + 
     &                         LR_FLOP_COST_OUT
            ACC_FLOP_DEMOTE  = ACC_FLOP_DEMOTE  + HR_COST + BUILDQ_COST
            ACC_FLOP_RMB     = ACC_FLOP_RMB     + HR_COST + BUILDQ_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
          ENDIF
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_LRB_PRODUCT
      SUBROUTINE UPDATE_FLOP_STATS_DEC_ACC(LRB, NIV)
        TYPE(LRB_TYPE),INTENT(IN) :: LRB
        INTEGER,INTENT(IN) :: NIV
        DOUBLE PRECISION :: FLOP_COST
        FLOP_COST = 2.0D0*dble(LRB%M)*dble(LRB%N)*dble(LRB%K)
        IF (NIV .EQ. 1) THEN
!$OMP CRITICAL(lr_flop_gain_cri)
                  LR_FLOP_GAIN = LR_FLOP_GAIN - FLOP_COST
                  FLOP_LR_UPDT = FLOP_LR_UPDT + FLOP_COST
                  FLOP_LR_UPDT_OUT = FLOP_LR_UPDT_OUT + FLOP_COST
                  FLOP_DEC_ACC = FLOP_DEC_ACC + FLOP_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
        ELSE
!$OMP CRITICAL(lr_flop_gain_cri)
           ACC_LR_FLOP_GAIN = ACC_LR_FLOP_GAIN - FLOP_COST
           ACC_FLOP_LR_UPDT = ACC_FLOP_LR_UPDT + FLOP_COST
           ACC_FLOP_LR_UPDT_OUT = ACC_FLOP_LR_UPDT_OUT + 
     &                            FLOP_COST
           ACC_FLOP_DEC_ACC = ACC_FLOP_DEC_ACC + FLOP_COST
!$OMP END CRITICAL(lr_flop_gain_cri)
        ENDIF
      END SUBROUTINE UPDATE_FLOP_STATS_DEC_ACC
      SUBROUTINE UPDATE_FLOPS_STATS_ROOT(KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID)
        INTEGER, intent(in) :: KEEP50, NFRONT, NPIV,
     &           NPROW, NPCOL, MYID
        DOUBLE PRECISION :: COST, COST_PER_PROC
        INTEGER, PARAMETER :: LEVEL3 = 3
        CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NFRONT, KEEP50, LEVEL3, 
     &                          COST)
        COST_PER_PROC = dble(int( COST,8) / int(NPROW * NPCOL,8))
        ACC_FLOP_FRFRONTS = ACC_FLOP_FRFRONTS + COST_PER_PROC
        RETURN
      END SUBROUTINE UPDATE_FLOPS_STATS_ROOT
      SUBROUTINE INIT_STATS_FRONT(NFRONT,INODE,NASS,NCB) 
        INTEGER,INTENT(IN) :: NFRONT,INODE,NASS,NCB
        FRONT_L11_BLR_SAVINGS = 0.D0 
        FRONT_U11_BLR_SAVINGS = 0.D0 
        FRONT_L21_BLR_SAVINGS = 0.D0 
        FRONT_U12_BLR_SAVINGS = 0.D0 
                  LR_FLOP_GAIN = 0.D0
                  FLOP_CB_DEMOTE = 0.D0
                  FLOP_CB_PROMOTE = 0.D0
                  FLOP_PROMOTE = 0.D0
                  FLOP_FR_UPDT = 0.D0 
                  FLOP_LR_UPDT = 0.D0
                  FLOP_LR_UPDT_OUT = 0.D0
                  FLOP_RMB = 0.D0
                  FLOP_FR_TRSM = 0.D0 
                  FLOP_LR_TRSM = 0.D0
                  FLOP_DEMOTE = 0.D0
                  FLOP_DEC_ACC = 0.D0
                  FLOP_REC_ACC = 0.D0
                  FLOP_FR_SWAP = 0.D0
                FLOP_PANEL = 0.D0
                FLOP_TRSM = 0.D0
      END SUBROUTINE INIT_STATS_FRONT
      SUBROUTINE INIT_STATS_GLOBAL(id)
        use DMUMPS_STRUC_DEF
        TYPE (DMUMPS_STRUC), TARGET :: id
        ACC_MRY_CB_GAIN = 0.D0
        ACC_MRY_CB_FR = 0.D0
        ACC_FLOP_CB_DEMOTE = 0.D0
        ACC_FLOP_CB_PROMOTE = 0.D0
        ACC_FLOP_PROMOTE = 0.D0
        ACC_FLOP_FR_FACTO = 0.D0 
        ACC_FLOP_LR_FACTO = 0.D0
        ACC_FLOP_FR_UPDT = 0.D0 
        ACC_FLOP_LR_UPDT = 0.D0
        ACC_FLOP_LR_UPDT_OUT = 0.D0
        ACC_FLOP_RMB = 0.D0
        ACC_FLOP_FR_TRSM = 0.D0 
        ACC_FLOP_LR_TRSM = 0.D0
        ACC_FLOP_DEMOTE = 0.D0
        ACC_FLOP_TRSM = 0.D0
        ACC_FLOP_DEC_ACC = 0.D0
        ACC_FLOP_REC_ACC = 0.D0
        ACC_FLOP_FR_SWAP = 0.D0
        ACC_FLOP_PANEL = 0.D0
        ACC_FLOP_FRFRONTS = 0.D0
        ACC_FLOP_FR_SOLVE = 0.D0
        ACC_FLOP_LR_SOLVE = 0.D0
        ACC_LR_FLOP_GAIN = 0.D0
        TOTAL_NBLOCKS_ASS = 0
        TOTAL_NBLOCKS_CB = 0
        AVG_BLOCKSIZE_ASS = 0.D0
        AVG_BLOCKSIZE_CB = 0.D0
        MIN_BLOCKSIZE_ASS = huge(1)
        MAX_BLOCKSIZE_ASS = 0
        MIN_BLOCKSIZE_CB = huge(1)
        MAX_BLOCKSIZE_CB = 0
        ACC_FR_MRY = 0.D0
        GLOBAL_BLR_SAVINGS = 0.D0 
        ACC_UPDT_TIME = 0.D0
        ACC_UPDT_TIME_OUT = 0.D0
        ACC_RMB_TIME = 0.D0
        ACC_PROMOTING_TIME = 0.D0
        ACC_PROMOTING_UCFS_TIME = 0.D0
        ACC_DEMOTING_TIME = 0.D0
        ACC_CB_DEMOTING_TIME = 0.D0
        ACC_FRPANELS_TIME = 0.0D0
        ACC_FAC_I_TIME = 0.0D0
        ACC_FAC_MQ_TIME = 0.0D0
        ACC_FAC_SQ_TIME = 0.0D0
        ACC_FRFRONTS_TIME = 0.0D0
        ACC_DIAGCOPY_TIME = 0.0D0
        ACC_FRSWAP_TIME = 0.0D0
        ACC_TRSM_TIME = 0.D0
        ACC_LR_MODULE_TIME = 0.D0
        ACC_DEC_ASM1_TIME = 0.0D0
        ACC_DEC_MAPLIG1_TIME = 0.0D0
        ACC_DEC_LOCASM2_TIME = 0.0D0
        ACC_DEC_ASMS2S_TIME = 0.0D0
        ACC_DEC_ASMS2M_TIME = 0.0D0
        CNT_NODES = 0
        STEP_STATS => id%STEP
      END SUBROUTINE INIT_STATS_GLOBAL
      SUBROUTINE STATS_COMPUTE_MRY_FRONT_TYPE1(NASS, NCB,
     &                       SYM, INODE, NELIM)
        INTEGER,INTENT(IN) :: NASS, NCB, SYM, INODE, NELIM
        DOUBLE PRECISION :: FRONT_BLR_SAVINGS, FRONT_FR_MRY
        IF (SYM .GT. 0) THEN
           FRONT_BLR_SAVINGS  = FRONT_L11_BLR_SAVINGS 
     &                        + FRONT_L21_BLR_SAVINGS
           FRONT_FR_MRY   = dble(NASS-NELIM) *
     &                            (dble(NASS-NELIM)+1.D0)/2.D0
     &                    + dble(NASS-NELIM) * dble(NCB+NELIM)
        ELSE
           FRONT_BLR_SAVINGS  = FRONT_L11_BLR_SAVINGS 
     &                        + FRONT_L21_BLR_SAVINGS
     &                        + FRONT_U11_BLR_SAVINGS
     &                        + FRONT_U12_BLR_SAVINGS
           FRONT_FR_MRY   = dble(NASS-NELIM) * dble(NASS-NELIM) 
     &                      + 2.0D0 * dble(NASS-NELIM) * dble(NCB+NELIM)
        END IF
           ACC_FR_MRY  = ACC_FR_MRY + FRONT_FR_MRY
           GLOBAL_BLR_SAVINGS = GLOBAL_BLR_SAVINGS + FRONT_BLR_SAVINGS
      END SUBROUTINE STATS_COMPUTE_MRY_FRONT_TYPE1
      SUBROUTINE STATS_COMPUTE_MRY_FRONT_TYPE2(NASS, NFRONT,
     &                       SYM, INODE, NELIM)
        INTEGER,INTENT(IN) :: NASS, NFRONT, SYM, INODE, NELIM
        IF (SYM .GT. 0) THEN
          ACC_FR_MRY   = ACC_FR_MRY +
     &                       dble(NASS-NELIM) *
     &                      (dble(NASS-NELIM)+1.D0)/2.D0
     &                    + dble(NASS-NELIM) * dble(NFRONT-NASS+NELIM)
        ELSE
          ACC_FR_MRY = ACC_FR_MRY +
     &                   dble(NASS-NELIM) * dble(NASS-NELIM) 
     &         + 2.0D0 * dble(NASS-NELIM) * dble(NFRONT-NASS+NELIM)
        ENDIF
      END SUBROUTINE STATS_COMPUTE_MRY_FRONT_TYPE2
      SUBROUTINE STATS_COMPUTE_MRY_FRONT_CB(NROWS, NCOLS,
     &                       SYM, NIV, INODE,
     &                       FRONT_CB_BLR_SAVINGS)
        INTEGER,INTENT(IN) :: NROWS, NCOLS, SYM, NIV, INODE,
     &                        FRONT_CB_BLR_SAVINGS
        DOUBLE PRECISION :: MRY_CB_FR
        IF (SYM.EQ.0) THEN
          MRY_CB_FR = dble(NCOLS)*dble(NROWS)
        ELSE
          MRY_CB_FR = dble(NCOLS-NROWS)*dble(NROWS) +
     &                dble(NROWS)*dble(NROWS+1)/2.D0
        ENDIF
!$OMP CRITICAL(lr_flop_gain_cri)
        ACC_MRY_CB_FR = ACC_MRY_CB_FR + MRY_CB_FR
        ACC_MRY_CB_GAIN = ACC_MRY_CB_GAIN + dble(FRONT_CB_BLR_SAVINGS)
!$OMP END CRITICAL(lr_flop_gain_cri)
      END SUBROUTINE STATS_COMPUTE_MRY_FRONT_CB
      SUBROUTINE STATS_STORE_BLR_PANEL_MRY(BLR_PANEL, NB_INASM,
     &                NB_INCB, DIR, NIV)
        INTEGER,INTENT(IN) :: NB_INASM, NB_INCB, NIV
        TYPE(LRB_TYPE),  INTENT(IN) :: BLR_PANEL(:)
        CHARACTER(len=1) :: DIR
        INTEGER :: I
        IF (NB_INASM.GT.0.AND.DIR .EQ.'V') THEN 
          ACC_FLOP_FR_SOLVE = ACC_FLOP_FR_SOLVE + 
     &      dble(BLR_PANEL(1)%N)*dble(BLR_PANEL(1)%N-1)
          ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + 
     &      dble(BLR_PANEL(1)%N)*dble(BLR_PANEL(1)%N-1)
        ENDIF
        DO I = 1 , NB_INASM
            IF (DIR .EQ. 'V') THEN 
              ACC_FLOP_FR_SOLVE = ACC_FLOP_FR_SOLVE + 
     &            2.0D0*dble(BLR_PANEL(I)%M)*dble(BLR_PANEL(I)%N) 
            ENDIF
            IF (BLR_PANEL(I)%ISLR) THEN
              IF (DIR .EQ. 'H') THEN 
                 IF (NIV .EQ. 1) THEN
                    FRONT_U11_BLR_SAVINGS   = 
     &                FRONT_U11_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K )    *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                 ELSE
                    GLOBAL_BLR_SAVINGS   = 
     &                GLOBAL_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K )    *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                 ENDIF
              ELSE
                ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + 
     &             2.0D0*(dble(BLR_PANEL(I)%M)+dble(BLR_PANEL(I)%N))*
     &             dble(BLR_PANEL(I)%K)
                 IF (NIV .EQ. 1) THEN
                    FRONT_L11_BLR_SAVINGS = 
     &                FRONT_L11_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K )    *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                 ELSE
                    GLOBAL_BLR_SAVINGS = 
     &                GLOBAL_BLR_SAVINGS + 
     &            dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &            dble( BLR_PANEL(I)%K )   *
     &            dble( BLR_PANEL(I)%M) + dble(BLR_PANEL(I)%N ) 
                 ENDIF
              ENDIF
            ELSE
              IF (DIR .EQ. 'V') THEN 
                ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + 
     &          2.0D0*dble(BLR_PANEL(I)%M)*dble(BLR_PANEL(I)%N) 
              ENDIF
            ENDIF
        END DO
        DO I = NB_INASM + 1 , NB_INASM + NB_INCB
            IF (DIR .EQ. 'V') THEN 
              ACC_FLOP_FR_SOLVE = ACC_FLOP_FR_SOLVE + 
     &            2.0D0*dble(BLR_PANEL(I)%M)*dble(BLR_PANEL(I)%N) 
            ENDIF
            IF (BLR_PANEL(I)%ISLR) THEN
              IF (DIR .EQ. 'H') THEN 
                 IF (NIV .EQ. 1) THEN
                    FRONT_U12_BLR_SAVINGS = 
     &                FRONT_U12_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K )   *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                 ELSE
                    GLOBAL_BLR_SAVINGS = 
     &                GLOBAL_BLR_SAVINGS + 
     &             dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &             dble( BLR_PANEL(I)%K )    *
     &             dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                 ENDIF
              ELSE
                ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + 
     &             2.0D0*(dble(BLR_PANEL(I)%M)+dble(BLR_PANEL(I)%N))*
     &             dble(BLR_PANEL(I)%K)
                IF (NIV .EQ. 1) THEN
                  FRONT_L21_BLR_SAVINGS = 
     &                FRONT_L21_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K ) *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                ELSE
                  GLOBAL_BLR_SAVINGS = 
     &                GLOBAL_BLR_SAVINGS + 
     &              dble( BLR_PANEL(I)%M ) * dble ( BLR_PANEL(I)%N ) - 
     &              dble( BLR_PANEL(I)%K ) *
     &              dble( BLR_PANEL(I)%M + BLR_PANEL(I)%N ) 
                ENDIF
              ENDIF
            ELSE
              IF (DIR .EQ. 'V') THEN 
                ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + 
     &          2.0D0*dble(BLR_PANEL(I)%M)*dble(BLR_PANEL(I)%N) 
              ENDIF
            ENDIF
          END DO
      END SUBROUTINE STATS_STORE_BLR_PANEL_MRY
      SUBROUTINE STATS_COMPUTE_FLOP_FRONT_TYPE1( NFRONT, NASS, NPIV, 
     &                  KEEP50, INODE)
          INTEGER,INTENT(IN) :: NFRONT, KEEP50, NASS, NPIV, INODE
          DOUBLE PRECISION   :: FLOP_FR_FACTO
          CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NASS, 
     &                              KEEP50, 1, FLOP_FR_FACTO)
          ACC_FLOP_FR_FACTO   = ACC_FLOP_FR_FACTO   + FLOP_FR_FACTO
          ACC_LR_FLOP_GAIN    = ACC_LR_FLOP_GAIN    + LR_FLOP_GAIN
          ACC_FLOP_FR_UPDT    = ACC_FLOP_FR_UPDT    + FLOP_FR_UPDT
          ACC_FLOP_LR_UPDT    = ACC_FLOP_LR_UPDT    + FLOP_LR_UPDT
          ACC_FLOP_LR_UPDT_OUT= ACC_FLOP_LR_UPDT_OUT+ FLOP_LR_UPDT_OUT
          ACC_FLOP_RMB        = ACC_FLOP_RMB        + FLOP_RMB
          ACC_FLOP_FR_TRSM    = ACC_FLOP_FR_TRSM    + FLOP_FR_TRSM
          ACC_FLOP_LR_TRSM    = ACC_FLOP_LR_TRSM    + FLOP_LR_TRSM
          ACC_FLOP_DEMOTE     = ACC_FLOP_DEMOTE     + FLOP_DEMOTE
          ACC_FLOP_CB_DEMOTE  = ACC_FLOP_CB_DEMOTE  + FLOP_CB_DEMOTE
          ACC_FLOP_CB_PROMOTE = ACC_FLOP_CB_PROMOTE + FLOP_CB_PROMOTE
          ACC_FLOP_PROMOTE    = ACC_FLOP_PROMOTE    + FLOP_PROMOTE
          ACC_FLOP_DEC_ACC    = ACC_FLOP_DEC_ACC    + FLOP_DEC_ACC
          ACC_FLOP_REC_ACC    = ACC_FLOP_REC_ACC    + FLOP_REC_ACC
          ACC_FLOP_FR_SWAP    = ACC_FLOP_FR_SWAP    + FLOP_FR_SWAP
          ACC_FLOP_TRSM       = ACC_FLOP_TRSM       + FLOP_TRSM
          ACC_FLOP_PANEL      = ACC_FLOP_PANEL      + FLOP_PANEL
      END SUBROUTINE STATS_COMPUTE_FLOP_FRONT_TYPE1
      SUBROUTINE STATS_COMPUTE_FLOP_FRONT_TYPE2( NFRONT, NASS,
     &                KEEP50, INODE, NELIM)
          INTEGER,INTENT(IN) :: NFRONT, KEEP50, NASS, INODE, NELIM
          DOUBLE PRECISION   :: FLOP_FR_FACTO
          CALL MUMPS_GET_FLOPS_COST(NFRONT, NASS-NELIM, NASS, 
     &                              KEEP50, 2, FLOP_FR_FACTO)
          ACC_FLOP_FR_FACTO = ACC_FLOP_FR_FACTO + FLOP_FR_FACTO
      END SUBROUTINE STATS_COMPUTE_FLOP_FRONT_TYPE2
      SUBROUTINE STATS_COMPUTE_FLOP_SLAVE_TYPE2( NROW1, NCOL1,
     &                NASS1, KEEP50, INODE)
          INTEGER,INTENT(IN) :: NROW1, NCOL1, KEEP50, NASS1, INODE
          DOUBLE PRECISION   :: NROW2, NCOL2, NASS2
          DOUBLE PRECISION   :: FLOP_FR_FACTO
          NROW2 = dble(NROW1)
          NCOL2 = dble(NCOL1)
          NASS2 = dble(NASS1)
          IF (KEEP50.EQ.0) THEN
            FLOP_FR_FACTO = NROW2*NASS2*NASS2         
     &              + 2.0D0*NROW2*NASS2*(NCOL2-NASS2) 
          ELSE
             FLOP_FR_FACTO =
     &            NROW2*NASS2*NASS2
     &            +       NROW2*NASS2*NROW2
     &            + 2.0D0*NROW2*NASS2*(NCOL2-NASS2-NROW2)
          ENDIF
          ACC_FLOP_FR_FACTO = ACC_FLOP_FR_FACTO + FLOP_FR_FACTO
      END SUBROUTINE STATS_COMPUTE_FLOP_SLAVE_TYPE2
      SUBROUTINE UPDATE_FLOP_STATS_FRFRONTS(NFRONT, NPIV, NASS, SYM, 
     &                                        NIV)
          INTEGER, INTENT(IN) :: NFRONT, NPIV, NASS, SYM, NIV
          DOUBLE PRECISION    :: FLOP_FRFRONTS, FLOP_SOLVE
          CALL MUMPS_GET_FLOPS_COST(NFRONT, NPIV, NASS, 
     &                              SYM, NIV, FLOP_FRFRONTS)
          ACC_FLOP_FRFRONTS = ACC_FLOP_FRFRONTS + FLOP_FRFRONTS
          FLOP_SOLVE = dble(NPIV)*dble(NPIV-1) +
     &                 2.0D0*dble(NFRONT-NPIV)*dble(NPIV)
          ACC_FLOP_FR_SOLVE = ACC_FLOP_FR_SOLVE + FLOP_SOLVE
          ACC_FLOP_LR_SOLVE = ACC_FLOP_LR_SOLVE + FLOP_SOLVE
      END SUBROUTINE UPDATE_FLOP_STATS_FRFRONTS
      SUBROUTINE UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &                                        KEEP50, INODE)
          INTEGER,INTENT(IN) :: NROW1, NCOL1, KEEP50, NASS1, INODE
          DOUBLE PRECISION   :: NROW2, NCOL2, NASS2
          DOUBLE PRECISION   :: FLOP_FRFRONTS
          NROW2 = dble(NROW1)
          NCOL2 = dble(NCOL1)
          NASS2 = dble(NASS1)
          IF (KEEP50.EQ.0) THEN
            FLOP_FRFRONTS = NROW2*NASS2*NASS2         
     &              + 2.0D0*NROW2*NASS2*(NCOL2-NASS2) 
          ELSE
             FLOP_FRFRONTS =
     &            NROW2*NASS2*NASS2
     &            +       NROW2*NASS2*NROW2
     &            + 2.0D0*NROW2*NASS2*(NCOL2-NASS2-NROW2)
          ENDIF
          ACC_FLOP_FRFRONTS = ACC_FLOP_FRFRONTS + FLOP_FRFRONTS
      END SUBROUTINE UPD_FLOP_FRFRONT_SLAVE
      SUBROUTINE COMPUTE_GLOBAL_GAINS(NB_ENTRIES_FACTOR, 
     &                FLOP_NUMBER, NB_ENTRIES_FACTOR_withLR,
     &                PROKG, MPG)
        INTEGER(8), INTENT(IN) :: NB_ENTRIES_FACTOR   
        INTEGER, INTENT(IN)    :: MPG
        LOGICAL, INTENT(IN)    :: PROKG
        DOUBLE PRECISION, INTENT(IN)        :: FLOP_NUMBER   
        INTEGER(8), INTENT(OUT) :: 
     &                  NB_ENTRIES_FACTOR_withLR 
        IF (NB_ENTRIES_FACTOR < 0) THEN
         IF (PROKG.AND.MPG.GT.0) THEN
          WRITE(MPG,*) "NEGATIVE NUMBER OF ENTRIES IN FACTOR"
          WRITE(MPG,*) "===> OVERFLOW ?"
         END IF
        END IF
        IF (ACC_FR_MRY .EQ. 0) THEN
           GLOBAL_MRY_LPRO_COMPR = 100.0D0
        ELSE
           GLOBAL_MRY_LPRO_COMPR = 100.0D0 *
     &                             GLOBAL_BLR_SAVINGS/ACC_FR_MRY
        ENDIF
        IF (ACC_MRY_CB_FR .EQ. 0) THEN 
          ACC_MRY_CB_FR = 100.0D0
        END IF
        NB_ENTRIES_FACTOR_withLR = NB_ENTRIES_FACTOR -
     &                             int(GLOBAL_BLR_SAVINGS,8)
        IF (NB_ENTRIES_FACTOR.EQ.0) THEN
          FACTOR_PROCESSED_FRACTION = 100.0D0
          GLOBAL_MRY_LTOT_COMPR = 100.0D0
        ELSE
          FACTOR_PROCESSED_FRACTION = 100.0D0 *
     &                            ACC_FR_MRY/dble(NB_ENTRIES_FACTOR)
          GLOBAL_MRY_LTOT_COMPR = 
     &            100.0D0*GLOBAL_BLR_SAVINGS/dble(NB_ENTRIES_FACTOR)
        ENDIF
        TOTAL_FLOP = FLOP_NUMBER
        ACC_FLOP_LR_FACTO    = ACC_FLOP_FR_FACTO - ACC_LR_FLOP_GAIN
     &                             + ACC_FLOP_DEMOTE + ACC_FLOP_PROMOTE
        RETURN
      END SUBROUTINE COMPUTE_GLOBAL_GAINS
      SUBROUTINE SAVEandWRITE_GAINS(LOCAL, K489, DKEEP, N, 
     &         ICNTL36,
     &         DEPTH, BCKSZ, NASSMIN, NFRONTMIN, SYM, K486,
     &         K472, K475, K478, K480, K481, K483, K484, 
     &         K8110, K849, 
     &         NBTREENODES, NPROCS, MPG, PROKG)
        INTEGER, INTENT(IN) :: LOCAL,K489,DEPTH, N,
     &      ICNTL36, BCKSZ,NASSMIN,
     &      NFRONTMIN, K486, NBTREENODES, MPG, 
     &      K472, K475, K478, K480, K481, K483, K484, 
     &      SYM, NPROCS
        INTEGER(8), INTENT(IN) :: K8110, K849
        LOGICAL, INTENT(IN) :: PROKG
        DOUBLE PRECISION :: DKEEP(230)
        LOGICAL PROK
        PROK = (PROKG.AND.(MPG.GE.0))
        IF (PROK) THEN
        WRITE(MPG,'(/A,A)') 
     & '-------------- Beginning of BLR statistics -------------------',
     & '--------------'
        WRITE(MPG,'(A,I2)') 
     & ' ICNTL(36) BLR variant                            = ', ICNTL36
        WRITE(MPG,'(A,ES8.1)')  
     & ' CNTL(7)   Dropping parameter controlling accuray = ',
     &                          DKEEP(8)
        WRITE(MPG,'(A)') 
     &          ' Statistics after BLR factorization :'
        WRITE(MPG,'(A,I8)') 
     &    '     Number of BLR fronts                     = ',
     &                          CNT_NODES
        ENDIF  
        IF (PROK) WRITE(MPG,'(A,F8.1,A)')
     &    '     Fraction of factors in BLR fronts        =',
     &                FACTOR_PROCESSED_FRACTION,'% ' 
        IF (PROK) THEN
          WRITE(MPG,'(A)') 
     &  '     Statistics on the number of entries in factors :'
          WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     INFOG(29) Theoretical nb of entries in factors      ='
     &     ,dble(K8110),' (100.0%)'
          WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     INFOG(35) Effective nb of entries  (% of INFOG(29)) ='
     &     ,dble(K849),' ('
     &     ,dble(100)*(dble(K849)/dble(max(K8110,1_8)))
     &     ,'%)'
        ENDIF
        IF (PROK) WRITE(MPG,'(A)')
     &  '     Statistics on operation counts (OPC):'
        TOTAL_FLOP = MAX(TOTAL_FLOP,EPSILON(1.0D0))
        DKEEP(55)=dble(TOTAL_FLOP)
        DKEEP(60)=dble(100)
        DKEEP(56)=dble(ACC_FLOP_LR_FACTO+ACC_FLOP_FRFRONTS)
        DKEEP(61)=dble(100*(ACC_FLOP_LR_FACTO+
     &       ACC_FLOP_FRFRONTS) /TOTAL_FLOP)
        IF (PROK) THEN 
        WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     RINFOG(3) Total theoretical operations counts       ='
     &     ,TOTAL_FLOP,' (',100*TOTAL_FLOP/TOTAL_FLOP,'%)'
        WRITE(MPG,'(A,ES10.3,A,F5.1,A)') 
     &  '     RINFOG(14) Total effective OPC     (% of RINFOG(3)) ='
     &     ,ACC_FLOP_LR_FACTO+ACC_FLOP_FRFRONTS,' ('
     &,100*(ACC_FLOP_LR_FACTO+ACC_FLOP_FRFRONTS)/TOTAL_FLOP
     &,'%)'
        ENDIF
      IF (PROK) WRITE(MPG,'(A,A)') 
     & '-------------- End of BLR statistics -------------------------',
     & '--------------'
      RETURN
      END SUBROUTINE SAVEandWRITE_GAINS
      END MODULE DMUMPS_LR_STATS
