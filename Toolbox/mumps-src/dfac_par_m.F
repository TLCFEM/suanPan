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
      MODULE DMUMPS_FAC_PAR_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_PAR(N, IW, LIW, A, LA, NSTK_STEPS,
     & ND,FILS,STEP, FRERE, DAD, CAND, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     & NSTEPSDONE, OPASS, OPELI, NELVA, COMP, MAXFRT, NMAXNPIV, NTOTPV,
     & NOFFNEGPV, NULLNEGPV, NB22T1, NB22T2, NBTINY,
     & DET_EXP, DET_MANT, DET_SIGN,
     & PTRIST, PTRAST, PIMASTER, PAMASTER, PTRARW, PTRAIW, PTR8ARR,
     & NINCOLARR, NINROWARR, PTRDEBARR, ITLOC, RHS_MUMPS, IPOOL, LPOOL,
     & L0_OMP_MAPPING, LL0_OMP_MAPPING,
     & MUMPS_TPS_ARR, DMUMPS_TPS_ARR, LTPS_ARR,
     & RINFO, POSFAC ,IWPOS, LRLU, IPTRLU, LRLUS, LEAF, NBROOT, NBRTOT,
     & NBROOT_UNDER_L0, 
     & UU, ICNTL, PTLUST, PTRFAC, INFO, KEEP,KEEP8, PROCNODE_STEPS,
     & SLAVEF,MYID, COMM_NODES, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES,
     & INTARR, DBLARR, root, PERM, NELT, FRTPTR, FRTELT, LPTRAR,
     & COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2, MEM_DISTRIB, NE,
     & DKEEP, PIVNUL_LIST_STRUCT, LRGROUPS )
!$    USE OMP_LIB
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY : 
     &                                    DMUMPS_DM_FREEALLDYNAMICCB
      USE DMUMPS_BUF, ONLY : DMUMPS_BUF_TEST
      USE DMUMPS_LOAD
      USE DMUMPS_OOC
      USE DMUMPS_FAC_ASM_MASTER_M
      USE DMUMPS_FAC_ASM_MASTER_ELT_M
      USE DMUMPS_FAC1_LDLT_M
      USE DMUMPS_FAC2_LDLT_M
      USE DMUMPS_FAC1_LU_M
      USE DMUMPS_FAC2_LU_M
      USE OMP_LIB
      USE MUMPS_TPS_M
      USE DMUMPS_TPS_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE MUMPS_PIVNUL_MOD
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER N, LIW, LPTRAR, NSTEPSDONE, INFO(80)
      DOUBLE PRECISION, INTENT(INOUT) :: OPASS, OPELI
      INTEGER, INTENT(INOUT) :: NELVA, COMP
      INTEGER, INTENT(INOUT) :: MAXFRT, NTOTPV, NMAXNPIV, NOFFNEGPV,
     &                          NULLNEGPV
      INTEGER, INTENT(INOUT) :: NB22T1, NB22T2, NBTINY
      INTEGER, INTENT(INOUT) :: DET_SIGN, DET_EXP
      DOUBLE PRECISION, INTENT(INOUT) :: DET_MANT
      INTEGER(8) :: LA
      DOUBLE PRECISION, TARGET :: A(LA)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      INTEGER LPOOL
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N+KEEP(253))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP8(85))
      INTEGER IW(LIW), NSTK_STEPS(KEEP(28))
      INTEGER(8), INTENT(IN) :: PTRARW(LPTRAR), PTRAIW(LPTRAR)
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER PTLUST(KEEP(28)), PERM(N)
      INTEGER CAND(SLAVEF+1,max(1,KEEP(56)))
      INTEGER   ISTEP_TO_INIV2(KEEP(71)),
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IPOOL(LPOOL)
      INTEGER NE(KEEP(28))
      DOUBLE PRECISION RINFO(40)
      INTEGER(8) :: PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER(8) :: POSFAC, LRLU, LRLUS, IPTRLU
      INTEGER IWPOS, LEAF, NBROOT, NBRTOT
      INTEGER, INTENT(in) :: NBROOT_UNDER_L0
      INTEGER COMM_LOAD, ASS_IRECV
      DOUBLE PRECISION UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER INTARR( KEEP8(27) )
      LOGICAL IS_ISOLATED_NODE
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION DKEEP(230)
      INTEGER LRGROUPS(KEEP(280))
      INTEGER, INTENT( IN ) :: LTPS_ARR
      TYPE (MUMPS_TPS_T), TARGET  :: MUMPS_TPS_ARR( LTPS_ARR )
      TYPE (DMUMPS_TPS_T), TARGET :: DMUMPS_TPS_ARR( LTPS_ARR )
      INTEGER, INTENT( IN ) :: LL0_OMP_MAPPING
      INTEGER, INTENT( IN ) :: L0_OMP_MAPPING( LL0_OMP_MAPPING )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0, DONE = 1.0D0
      INTEGER INODE
      INTEGER IWPOSCB
      INTEGER FPERE, TYPEF
      INTEGER MP, LP, DUMMY(1)
      INTEGER NBFIN, NBROOT_TRAITEES
      INTEGER NFRONT, IOLDPS
      INTEGER(8) :: NFRONT8
      INTEGER(8) :: POSELT
      INTEGER IPOSROOT, IPOSROOTROWINDICES
      INTEGER GLOBK109
      INTEGER(8) :: LBUFRX
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: BUFRX
      LOGICAL :: IS_BUFRX_ALLOCATED
      DOUBLE PRECISION FLOP1
      INTEGER TYPE
      LOGICAL SON_LEVEL2, SET_IRECV, BLOCKING,
     &        MESSAGE_RECEIVED
      LOGICAL AVOID_DELAYED
      LOGICAL LAST_CALL
      INTEGER MASTER_ROOT
      INTEGER LOCAL_M, LOCAL_N
      INTEGER LRHS_CNTR_MASTER_ROOT, FWD_LOCAL_N_RHS
      LOGICAL ROOT_OWNER
!$    INTEGER :: NOMP
!$    INTEGER(8) :: CHUNK8  
      INTEGER  LIWK_RR, PHASE, MBLOCK, NBLOCK
      INTEGER(8) :: LWK_RR
      INTEGER(8) :: I8
      INTEGER  I, K, KEEP17_LU
      INTEGER  NOFFNEGPV_ROOT, NTOTPV_ROOT, NB22T1_ROOT, NBTINY_ROOT,
     &         NULLNEGPV_ROOT,
     &         DET_EXP_ROOT, DET_SIGN_ROOT, 
     &         LRecord, Header_ROOT(5)
      DOUBLE PRECISION  DET_MANT_ROOT  
      DOUBLE PRECISION DKEEP_SAVE(230)
      DOUBLE PRECISION, DIMENSION(:), POINTER     :: A_ROOT_SAVE
      LOGICAL                            :: IS_A_ROOT_SAVE_ALLOCATED
      INTEGER, DIMENSION(:), ALLOCATABLE :: RECORD_ROOT
      INTEGER KEEP_SAVE(500)             
      INTEGER(8) KEEP8_SAVE(150)
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      LOGICAL MUMPS_INSSARBR
      EXTERNAL MUMPS_INSSARBR
      LOGICAL DMUMPS_POOL_EMPTY
      EXTERNAL DMUMPS_POOL_EMPTY, DMUMPS_EXTRACT_POOL
      LOGICAL STACK_RIGHT_AUTHORIZED
      INTEGER, EXTERNAL :: MUMPS_NUMROC
      INTEGER JOBASS, ETATASS
      INTEGER(8) :: LAFAC
      INTEGER LIWFAC, STRAT, TYPEFile, NextPiv2beWritten,
     &        IDUMMY
      INTEGER(8) :: ITMP8
      TYPE(IO_BLOCK) :: MonBloc
      INCLUDE 'mumps_headers.h'
      INTEGER MPA
      DOUBLE PRECISION  OPLAST_PRINTED
      DOUBLE PRECISION :: ROOTTIME
      INTEGER:: ITH
      DOUBLE PRECISION :: DUMMY_FLOP_ESTIM_ACC
      DUMMY_FLOP_ESTIM_ACC = 0.0d0
      ITLOC(1:N+KEEP(253)) =0
      ASS_IRECV = MPI_REQUEST_NULL
      MP = ICNTL(2)
      LP = ICNTL(1)
      IWPOSCB = LIW
      NULLIFY(BUFRX)
      IS_BUFRX_ALLOCATED = .FALSE.
      KEEP(143) = -1
      KEEP17_LU = -1
      NULLIFY(A_ROOT_SAVE)
      IS_A_ROOT_SAVE_ALLOCATED = .FALSE.
      IF ( INFO(1) .LT. 0 ) THEN
        GOTO 640
      ENDIF
      OPLAST_PRINTED = DONE
      MPA            = ICNTL(2)
      IF (ICNTL(4).LT.2) MPA=0
      IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
      CALL DMUMPS_UPDATE_PROGRESS( OPELI, KEEP8 )
      STACK_RIGHT_AUTHORIZED = .TRUE.
      CALL DMUMPS_ALLOC_CB( .FALSE., 0_8,
     &     .FALSE., .FALSE., MYID_NODES, N, KEEP, KEEP8, DKEEP,
     &     IW, LIW, A, LA, LRLU, IPTRLU, IWPOS, IWPOSCB,
     &     SLAVEF, PROCNODE_STEPS, DAD,
     &     PTRIST, PTRAST, STEP, PIMASTER,
     &     PAMASTER, KEEP(IXSZ), 0_8, -444, -444, .true.,
     &     COMP, LRLUS, KEEP8(67),
     &     INFO(1), INFO(2)
     &     )
      JOBASS  = 0
      ETATASS = 0
      NBFIN = NBRTOT
      NBROOT_TRAITEES = 0
      KEEP(121)=0
      IF ( KEEP(38).NE.0 ) THEN
        IF (root%yes) THEN
            CALL DMUMPS_ROOT_ALLOC_STATIC(
     &        root, KEEP(38), N, IW, LIW,
     &        A, LA,
     &        FILS, DAD, MYID_NODES, SLAVEF, PROCNODE_STEPS,
     &        LPTRAR, NELT, FRTPTR, FRTELT,
     &        PTRAIW, PTRARW,
     &        PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &        INTARR, DBLARR,
     &        LRLU, IPTRLU,
     &        IWPOS, IWPOSCB, PTRIST, PTRAST,
     &        STEP, PIMASTER, PAMASTER, ITLOC, RHS_MUMPS,
     &        COMP, LRLUS, INFO(1), KEEP,KEEP8, DKEEP, INFO(2) )
        ENDIF
        IF ( INFO(1) .LT. 0 ) GOTO 635
      END IF
      IF (KEEP(400).GT.0) THEN
        NBROOT_TRAITEES = NBROOT_UNDER_L0
        IF (NBROOT_TRAITEES .GT.0) THEN
          IF (NBROOT_TRAITEES.EQ.NBROOT) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF .GT. 1) THEN
               CALL DMUMPS_MCAST2( NBROOT, 1, MPI_INTEGER,
     &         MYID_NODES, COMM_NODES, RACINE, SLAVEF, KEEP )
            ENDIF
          ENDIF
        ENDIF
        IF (NBFIN .EQ. 0) GOTO 640
      ENDIF
      KEEP(429)=0
 20   CONTINUE
      CALL DMUMPS_UPDATE_PROGRESS( OPELI, KEEP8 )
      CALL MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL, INFO, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 635 
      NIV1_FLAG=0
      SET_IRECV = .TRUE.
      BLOCKING = .FALSE.
      MESSAGE_RECEIVED = .FALSE.
      IF ( SLAVEF .GT. 1 ) THEN 
        CALL DMUMPS_TRY_RECVTREAT(
     &      COMM_LOAD, ASS_IRECV, BLOCKING, SET_IRECV,
     &      MESSAGE_RECEIVED,
     &      MPI_ANY_SOURCE, MPI_ANY_TAG,
     &      STATUS, BUFR, LBUFR,
     &      LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA,
     &      PTRIST, PTLUST, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS,
     &      COMP, INFO(1), INFO(2), COMM_NODES, PERM,
     &      IPOOL, LPOOL, LEAF, NBFIN, MYID_NODES, SLAVEF,
     &      root, OPASS, OPELI, ITLOC, RHS_MUMPS, FILS, DAD,
     &      PTRARW, PTRAIW,
     &      PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &      INTARR, DBLARR, ICNTL, KEEP,KEEP8,DKEEP, ND, FRERE,
     &      LPTRAR, NELT, FRTPTR, FRTELT,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &      STACK_RIGHT_AUTHORIZED
     &               , LRGROUPS
     &       )
        CALL DMUMPS_LOAD_RECV_MSGS(COMM_LOAD)
      ENDIF
      IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
      IF (MESSAGE_RECEIVED) THEN
          IF ( INFO(1) .LT. 0 ) GO TO 640
          IF ( NBFIN .eq. 0 ) GOTO 640
      ELSE
          IF ( .NOT. DMUMPS_POOL_EMPTY( IPOOL, LPOOL) )THEN
            CALL DMUMPS_EXTRACT_POOL( N, IPOOL, LPOOL,
     &      PROCNODE_STEPS,
     &      SLAVEF, STEP, INODE, KEEP,KEEP8, MYID_NODES, ND,
     &      (.NOT. STACK_RIGHT_AUTHORIZED) )
            STACK_RIGHT_AUTHORIZED = .TRUE.
            IF (KEEP(47) .GE. 3) THEN
              CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &              IPOOL, LPOOL,
     &              PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &              MYID_NODES, STEP, N, ND, FILS )
            ENDIF
            IF (KEEP(47).EQ.4) THEN
               IF(INODE.GT.0.AND.INODE.LE.N)THEN
                  IF((NE(STEP(INODE)).EQ.0).AND.
     &                 (FRERE(STEP(INODE)).EQ.0))THEN
                     IS_ISOLATED_NODE=.TRUE.
                  ELSE
                     IS_ISOLATED_NODE=.FALSE.
                  ENDIF
               ENDIF
               CALL DMUMPS_LOAD_SBTR_UPD_NEW_POOL(
     &              IS_ISOLATED_NODE,INODE,IPOOL,LPOOL,
     &              MYID_NODES,SLAVEF,COMM_LOAD,KEEP,KEEP8)
            ENDIF
            IF ((( KEEP(80) == 2 .OR. KEEP(80)==3 ) .AND.
     &           ( KEEP(47) == 4 )).OR.
     &           (KEEP(80) == 1 .AND. KEEP(47) .GE. 1)) THEN
               CALL DMUMPS_UPPER_PREDICT(INODE,STEP,KEEP(28),
     &         PROCNODE_STEPS,FRERE,ND,COMM_LOAD,SLAVEF,
     &         MYID_NODES,KEEP,KEEP8,N)
            END IF
            GOTO 30
          ELSE
            CALL DMUMPS_BUF_TEST()
          ENDIF
      ENDIF
      GO TO 20
 30   CONTINUE
      IF ( INODE .LT. 0 ) THEN
        INODE = -INODE
        FPERE = DAD(STEP(INODE))
        GOTO 130
      ELSE IF (INODE.GT.N) THEN
       INODE = INODE - N
       IF (INODE.EQ.KEEP(38)) THEN
         NBROOT_TRAITEES = NBROOT_TRAITEES + 1
         IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
            NBFIN = NBFIN - NBROOT
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID_NODES,
     &          COMM_NODES, RACINE, SLAVEF, KEEP )
            END IF
         ENDIF
         IF (NBFIN.EQ.0) GOTO 640
         GOTO 20
       ENDIF
       TYPE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
       IF (TYPE.EQ.1) GOTO 100
       FPERE = DAD(STEP(INODE))
       AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
       IF ( KEEP(50) .eq. 0 ) THEN
         CALL  DMUMPS_FAC2_LU( COMM_LOAD, ASS_IRECV,
     &        N, INODE, FPERE, IW, LIW, A, LA, UU,
     &        NOFFNEGPV, NTOTPV, NBTINY,
     &        DET_EXP, DET_MANT, DET_SIGN,
     &        COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &        NBFIN,LEAF, INFO(1), INFO(2), IPOOL,LPOOL,
     &        SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &        LRLUS, COMP, PTRIST, PTRAST, PTLUST, PTRFAC,
     &        STEP, PIMASTER, PAMASTER,
     &        NSTK_STEPS,PERM,PROCNODE_STEPS,
     &        root, OPASS, OPELI, ITLOC, RHS_MUMPS,
     &        FILS, DAD, PTRARW, PTRAIW,
     &        PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &        INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &        LPTRAR, NELT, FRTPTR, FRTELT, SEUIL,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &        DKEEP(1), PIVNUL_LIST_STRUCT
     &        , LRGROUPS
     &             )
        IF ( INFO(1) .LT. 0 ) GOTO 640
        IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
       ELSE
         CALL DMUMPS_FAC2_LDLT( COMM_LOAD, ASS_IRECV,
     &             N, INODE, FPERE, IW, LIW, A, LA, UU,
     &             NOFFNEGPV, NULLNEGPV, NTOTPV,
     &             NB22T2, NBTINY, DET_EXP, DET_MANT, DET_SIGN,
     &             COMM_NODES, MYID_NODES, BUFR, LBUFR,LBUFR_BYTES,
     &             NBFIN,LEAF, INFO(1), INFO(2), IPOOL,LPOOL,
     &             SLAVEF, POSFAC, IWPOS, IWPOSCB, IPTRLU, LRLU,
     &             LRLUS, COMP, PTRIST, PTRAST, PTLUST, PTRFAC,
     &             STEP, PIMASTER, PAMASTER,
     &             NSTK_STEPS,PERM,PROCNODE_STEPS,
     &             root, OPASS, OPELI, ITLOC, RHS_MUMPS,
     &             FILS, DAD, PTRARW, PTRAIW, 
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             INTARR, DBLARR, ICNTL, KEEP,KEEP8, ND, FRERE,
     &             LPTRAR, NELT, FRTPTR, FRTELT, SEUIL_LDLT_NIV2,
     &             ISTEP_TO_INIV2, TAB_POS_IN_PERE, AVOID_DELAYED,
     &              DKEEP(1), PIVNUL_LIST_STRUCT
     &           , LRGROUPS
     &             )
        IF ( INFO(1) .LT. 0 ) GOTO 640
        IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
        IF ( IW( PTLUST(STEP(INODE)) + KEEP(IXSZ) + 5 ) .GT. 1 ) THEN
             GOTO 20
        END IF
       END IF
       GOTO 130
      ENDIF
      IF (INODE.EQ.KEEP(38)) THEN
         CALL  DMUMPS_LAST_RTNELIND( COMM_LOAD, ASS_IRECV,
     &    root, FRERE,
     &    INODE,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_STEPS, COMP,
     &    INFO(1), INFO(2), COMM_NODES,
     &    PERM,
     &    IPOOL, LPOOL, LEAF,
     &    NBFIN, MYID_NODES, SLAVEF,
     &
     &    OPASS, OPELI, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP, ND,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &      )
         IF ( INFO(1) .LT. 0 ) GOTO 640
         IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
        GOTO 20
      ENDIF
      TYPE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
      IF (TYPE.EQ.1) THEN
        IF (KEEP(55).NE.0) THEN
         CALL DMUMPS_FAC_ASM_NIV1_ELT( COMM_LOAD, ASS_IRECV,
     &        NELT, FRTPTR, FRTELT,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRT,root,OPASS, OPELI,
     &        PTRIST,PTLUST,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,
     &        PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &        ITLOC, RHS_MUMPS, NSTEPSDONE, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS, KEEP8(67),
     &        ICNTL, KEEP,KEEP8,DKEEP,
     &        INTARR,KEEP8(27),DBLARR,KEEP8(26),
     &    NSTK_STEPS,PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    , MUMPS_TPS_ARR, DMUMPS_TPS_ARR,
     &    L0_OMP_MAPPING
     &    )
        ELSE
         JOBASS = 0
         CALL DMUMPS_FAC_ASM_NIV1(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRT,root,OPASS, OPELI,
     &        PTRIST,PTLUST,PTRFAC,PTRAST,STEP, PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW, PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &        ITLOC, RHS_MUMPS, NSTEPSDONE, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS, KEEP8(67),
     &        ICNTL, KEEP,KEEP8,DKEEP, INTARR,KEEP8(27),
     &        DBLARR,KEEP8(26),
     &    NSTK_STEPS,PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &    PERM,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, JOBASS,ETATASS
     &               , LRGROUPS
     &    , MUMPS_TPS_ARR, DMUMPS_TPS_ARR,
     &    L0_OMP_MAPPING
     &    )
        ENDIF
        IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
       IF ( INFO(1) .LT. 0 ) GOTO 640
        IF ((IW(PTLUST(STEP(INODE))+XXNBPR).GT.0).OR.(SON_LEVEL2)) THEN
          GOTO 20
        ENDIF
      ELSE
        IF ( KEEP(55) .eq. 0 ) THEN
          CALL DMUMPS_FAC_ASM_NIV2(COMM_LOAD, ASS_IRECV,
     &    N, INODE, IW, LIW, A, LA,
     &    INFO(1),
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRT,
     &    root, OPASS, OPELI, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW, PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    ITLOC, RHS_MUMPS, NSTEPSDONE,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,KEEP8(27),DBLARR,KEEP8(26),
     &    PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0)
     &               , LRGROUPS
     &    )
        ELSE
          CALL DMUMPS_FAC_ASM_NIV2_ELT( COMM_LOAD, ASS_IRECV,
     &    NELT, FRTPTR, FRTELT,
     &    N, INODE, IW, LIW, A, LA, INFO(1),
     &    ND, FILS, FRERE, DAD, CAND,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    MAXFRT,
     &    root, OPASS, OPELI, PTRIST, PTLUST, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, PTRARW, NSTK_STEPS,
     &    PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    ITLOC, RHS_MUMPS, NSTEPSDONE,
     &    COMP, LRLU, IPTRLU, IWPOS, IWPOSCB, POSFAC, LRLUS,
     &    ICNTL, KEEP,KEEP8,DKEEP,INTARR,KEEP8(27),DBLARR,KEEP8(26),
     &    PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    NBFIN, LEAF, IPOOL, LPOOL, PERM,
     &    MEM_DISTRIB(0)
     &               , LRGROUPS
     &     )
        END IF
        IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
        IF (INFO(1).LT.0) GOTO 640
        GOTO 20
      ENDIF
 100  CONTINUE
       FPERE = DAD(STEP(INODE))
      IF ( INODE .eq. KEEP(20) ) THEN
        POSELT = PTRAST(STEP(INODE))
        IF (PTRFAC(STEP(INODE)).NE.POSELT) THEN
          WRITE(*,*) "ERROR 2 in DMUMPS_FAC_PAR", POSELT
          GOTO 635
        ENDIF
        IF (KEEP(118).GE.40) THEN
          IOLDPS = PTLUST(STEP(INODE))
          LRecord = IW(IOLDPS+XXI)
          IF (allocated(RECORD_ROOT)) DEALLOCATE(RECORD_ROOT)
          ALLOCATE(RECORD_ROOT(LRecord), stat=IERR)
          IF (IERR.GT.0) THEN
            INFO(1)= -13
            INFO(2)= LRecord
            IF (LP > 0)
     &      write(LP,*) "ERROR allocate RECORD_ROOT"
            GOTO 635
          ENDIF
          RECORD_ROOT(1:LRecord) = IW(IOLDPS:IOLDPS+LRecord-1)
        ENDIF
        CALL DMUMPS_CHANGE_HEADER
     &       ( IW(PTLUST(STEP(INODE))+KEEP(IXSZ)), KEEP(253) )
        IF (KEEP(118).GE.40) THEN
          Header_ROOT(1:5) = IW(PTLUST(STEP(INODE))+KEEP(IXSZ):
     &                        PTLUST(STEP(INODE))+KEEP(IXSZ)+4)
        ENDIF
        GOTO 200
      END IF
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST(STEP(INODE))
      AVOID_DELAYED = ( (FPERE .eq. KEEP(20) .OR. FPERE .eq. KEEP(38))
     &                   .AND. KEEP(60).ne.0 )
      IF (KEEP(50).EQ.0) THEN
           CALL DMUMPS_FAC1_LU ( 
     &               N, INODE, IW, LIW, A, LA,
     &               IOLDPS, POSELT,
     &               INFO(1), INFO(2), UU, NOFFNEGPV, NTOTPV, NBTINY,
     &               DET_EXP, DET_MANT, DET_SIGN,
     &               KEEP,KEEP8,
     &               STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     &               SEUIL, AVOID_DELAYED, ETATASS,
     &              DKEEP(1), PIVNUL_LIST_STRUCT, IWPOS 
     &           , LRGROUPS
     &           , PERM
     &           )
           IF (INFO(1).LT.0) GOTO 635
      ELSE  
            IW( IOLDPS+4+KEEP(IXSZ) ) = 1
              CALL DMUMPS_FAC1_LDLT( N, INODE,
     &           IW, LIW, A, LA,
     &           IOLDPS, POSELT,
     &           INFO(1), INFO(2), UU, NOFFNEGPV, NULLNEGPV, NTOTPV,
     &           NB22T1, NBTINY, DET_EXP, DET_MANT, DET_SIGN,
     &           KEEP,KEEP8, MYID_NODES, SEUIL, AVOID_DELAYED,
     &           ETATASS,
     &           DKEEP(1), PIVNUL_LIST_STRUCT, IWPOS
     &           , LRGROUPS
     &           , PERM
     &           )
             IF (INFO(1).LT.0) GOTO 635
            IW( IOLDPS+4+KEEP(IXSZ) ) = STEP(INODE)
          ENDIF 
          JOBASS = ETATASS          
          IF (JOBASS.EQ.1) THEN
              CALL DMUMPS_FAC_ASM_NIV1(COMM_LOAD, ASS_IRECV,
     &        N,INODE,IW,LIW,A,LA,
     &        INFO(1),ND,
     &        FILS,FRERE,DAD,MAXFRT,root,OPASS, OPELI,
     &        PTRIST,PTLUST,PTRFAC,PTRAST,STEP,PIMASTER,PAMASTER,
     &        PTRARW,PTRAIW,PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &        ITLOC, RHS_MUMPS, NSTEPSDONE, SON_LEVEL2,
     &        COMP, LRLU, IPTRLU,
     &        IWPOS,IWPOSCB, POSFAC, LRLUS, KEEP8(67),
     &        ICNTL, KEEP,KEEP8,DKEEP,INTARR,KEEP8(27),DBLARR,KEEP8(26),
     &        NSTK_STEPS, PROCNODE_STEPS, SLAVEF,
     &        COMM_NODES, MYID_NODES,
     &        BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL, LPOOL, LEAF,
     &        PERM,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &        JOBASS,ETATASS
     &               , LRGROUPS
     &         )
          ENDIF
          IF (MPA.GT.0)
     &    CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
      IF (INFO(1).LT.0) GOTO 635
 130  CONTINUE
      TYPE  = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
      IF ( FPERE .NE. 0 ) THEN
        TYPEF = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(FPERE)),KEEP(199))
      ELSE
        TYPEF = -9999
      END IF
      CALL DMUMPS_FAC_STACK( COMM_LOAD, ASS_IRECV,
     &       N,INODE,TYPE,TYPEF,LA,IW,LIW,A,
     &       INFO(1),INFO(2),OPELI,NELVA,NMAXNPIV,
     &       PTRIST,PTLUST,PTRFAC,
     &       PTRAST, STEP, PIMASTER, PAMASTER,
     &       NE, POSFAC,LRLU, LRLUS,KEEP8(67),
     &       IPTRLU,ICNTL,KEEP,KEEP8,DKEEP,COMP,IWPOS,IWPOSCB,
     &       PROCNODE_STEPS,SLAVEF,FPERE,COMM_NODES,MYID_NODES,
     &       IPOOL, LPOOL, LEAF,
     &       NSTK_STEPS, PERM, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     &       root, OPASS, ITLOC, RHS_MUMPS, FILS, DAD, PTRARW, PTRAIW,
     &       PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &       INTARR, DBLARR,
     &       ND, FRERE, LPTRAR, NELT, FRTPTR, FRTELT,
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &       ,DUMMY_FLOP_ESTIM_ACC
     &       )
      IF (MPA.GT.0)
     &   CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                        dble(DKEEP(17)), 
     &                        OPELI,
     &                        OPLAST_PRINTED, MPA)
      IF (INFO(1).LT.0) GOTO 640
 200  CONTINUE
      IF ( INODE .eq. KEEP(38) ) THEN
        WRITE(*,*) 'Error .. in DMUMPS_FAC_PAR: ',
     &             ' INODE == KEEP(38)'
        CALL MUMPS_ABORT()
      END IF
      IF ( FPERE.EQ.0 ) THEN
        NBROOT_TRAITEES = NBROOT_TRAITEES + 1
        IF ( NBROOT_TRAITEES .EQ. NBROOT ) THEN
           IF (KEEP(201).EQ.1) THEN 
              CALL DMUMPS_OOC_FORCE_WRT_BUF_PANEL(IERR)
           ELSE IF ( KEEP(201).EQ.2) THEN 
              CALL DMUMPS_FORCE_WRITE_BUF(IERR)
           ENDIF
            NBFIN = NBFIN - NBROOT
            IF ( NBFIN .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_FAC_PAR: ',
     &                   ' NBFIN=', NBFIN
              CALL MUMPS_ABORT()
            END IF
            IF ( NBROOT .LT. 0 ) THEN
              WRITE(*,*) ' ERROR 1 in DMUMPS_FAC_PAR: ',
     &                   ' NBROOT=', NBROOT
              CALL MUMPS_ABORT()
            END IF
            IF (SLAVEF.GT.1) THEN
                DUMMY(1) = NBROOT
                CALL DMUMPS_MCAST2( DUMMY(1), 1, MPI_INTEGER,
     &          MYID_NODES, COMM_NODES, RACINE, SLAVEF, KEEP )
            END IF
        ENDIF
        IF (NBFIN.EQ.0)THEN
           GOTO 640
        ENDIF
      ELSEIF ( FPERE.NE.KEEP(38) .AND.
     &         MUMPS_PROCNODE(PROCNODE_STEPS(STEP(FPERE)),
     &         KEEP(199)) .EQ.  MYID_NODES ) THEN
        NSTK_STEPS(STEP(FPERE)) = NSTK_STEPS(STEP(FPERE))-1
        IF ( NSTK_STEPS( STEP( FPERE )).EQ.0) THEN
          IF (KEEP(234).NE.0 .AND.
     &      MUMPS_INSSARBR(PROCNODE_STEPS(STEP(INODE)),KEEP(199)))
     &      THEN
            STACK_RIGHT_AUTHORIZED = .FALSE.
          ENDIF
          CALL DMUMPS_INSERT_POOL_N(N, IPOOL, LPOOL,
     &         PROCNODE_STEPS, SLAVEF, KEEP(199), KEEP(28), KEEP(76),
     &         KEEP(80), KEEP(47), STEP, FPERE )
          IF (KEEP(47) .GE. 3) THEN
             CALL DMUMPS_LOAD_POOL_UPD_NEW_POOL(
     &            IPOOL, LPOOL,
     &            PROCNODE_STEPS, KEEP,KEEP8, SLAVEF, COMM_LOAD,
     &            MYID_NODES, STEP, N, ND, FILS )
          ENDIF
          CALL MUMPS_ESTIM_FLOPS( FPERE, N, PROCNODE_STEPS,KEEP(199),
     &           ND, FILS, FRERE, STEP, PIMASTER, KEEP(28),
     &           KEEP(50), KEEP(253), FLOP1,
     &           IW, LIW, KEEP(IXSZ) )
          IF (FPERE.NE.KEEP(20))
     &    CALL DMUMPS_LOAD_UPDATE(1,.FALSE.,FLOP1,KEEP,KEEP8)
        ENDIF
      ENDIF
      GO TO 20
 635  CONTINUE
      IF (allocated(RECORD_ROOT)) DEALLOCATE(RECORD_ROOT)
      CALL DMUMPS_BDC_ERROR( MYID_NODES, SLAVEF, COMM_NODES, KEEP )
 640  CONTINUE
        CALL DMUMPS_CANCEL_IRECV( INFO(1),
     &       KEEP,
     &       ASS_IRECV, BUFR, LBUFR,
     &       LBUFR_BYTES,
     &       COMM_NODES,
     &       MYID_NODES, SLAVEF)
       CALL DMUMPS_CLEAN_PENDING( INFO(1), KEEP,
     &      BUFR, LBUFR,
     &      LBUFR_BYTES,
     &      COMM_NODES, COMM_LOAD, SLAVEF,
     &      .TRUE.,
     &      .TRUE.)
      CALL MPI_BARRIER( COMM_NODES, IERR )
      IF (INFO(1) .LT. 0) THEN
        CALL DMUMPS_DM_FREEALLDYNAMICCB( MYID, N, SLAVEF, KEEP, KEEP8,
     &  IW, LIW, IWPOSCB, IWPOS,
     &  STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &  .FALSE. ) 
        IF ( KEEP(400) .GT. 0
     &     ) THEN
!$OMP PARALLEL DO SCHEDULE(STATIC,1)
            DO ITH = 1, KEEP(400)
              IF (associated(MUMPS_TPS_ARR(ITH)%IW)) THEN
                CALL DMUMPS_DM_FREEALLDYNAMICCB_I( MYID, N, SLAVEF,
     &          KEEP, KEEP8,
     &          MUMPS_TPS_ARR(ITH)%IW(1), MUMPS_TPS_ARR(ITH)%LIW,
     &          MUMPS_TPS_ARR(ITH)%IWPOSCB, MUMPS_TPS_ARR(ITH)%IWPOS,
     &          STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &          .TRUE. ) 
              ENDIF
            ENDDO
!$OMP END PARALLEL DO
        ENDIF
      ENDIF
      IF ( INFO(1) .GE. 0 ) THEN
          IF( KEEP(38) .NE. 0 .OR. KEEP(20).NE.0) THEN
            CALL MUMPS_SECDEB(ROOTTIME)
            MASTER_ROOT = MUMPS_PROCNODE(
     &                  PROCNODE_STEPS(STEP(max(KEEP(38),KEEP(20)))),
     &                  KEEP(199))
            ROOT_OWNER  = (MASTER_ROOT .EQ. MYID_NODES)
            IF ( KEEP(38) .NE. 0 ) THEN
               IF (KEEP(60).EQ.0) THEN
                 IOLDPS  = PTLUST(STEP(KEEP(38)))
                 LOCAL_M = IW(IOLDPS+2+KEEP(IXSZ))
                 LOCAL_N = IW(IOLDPS+1+KEEP(IXSZ))
               ELSE
                 IOLDPS  = -999
                 LOCAL_M = root%SCHUR_MLOC
                 LOCAL_N = root%SCHUR_NLOC
               ENDIF
               ITMP8   = int(LOCAL_M,8)*int(LOCAL_N,8)
               LBUFRX = min(int(root%MBLOCK,8)*int(root%NBLOCK,8),
     &            int(root%TOT_ROOT_SIZE,8)*int(root%TOT_ROOT_SIZE,8) )
               IS_BUFRX_ALLOCATED = .FALSE.
               IF ( LRLU .GT. LBUFRX ) THEN
                   BUFRX => A(POSFAC:POSFAC+LRLU-1_8)
                   LBUFRX=LRLU
               ELSE
                   ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                   IF (IERR.gt.0) THEN
                         INFO(1) = -13
                         CALL MUMPS_SET_IERROR(LBUFRX, INFO(2) )
                         IF (LP > 0 )
     &                   write(LP,*) ' Error allocating, real array ',
     &                   'of size before DMUMPS_FACTO_ROOT',  LBUFRX
                   ELSE
                     IS_BUFRX_ALLOCATED = .TRUE.
                   ENDIF
               ENDIF
               CALL MUMPS_PROPINFO( ICNTL, INFO,
     &              COMM_NODES, MYID_NODES )
               IF (INFO(1).GE.0) THEN
                 CALL DMUMPS_FACTO_ROOT( 
     &               MPA, MYID_NODES, MASTER_ROOT,
     &               root, N, KEEP(38),
     &               COMM_NODES, IW, LIW, IWPOS + 1,
     &               A, LA, PTRAST, PTLUST, PTRFAC, STEP,
     &               INFO(1), KEEP(50), KEEP(19),
     &               BUFRX(1), LBUFRX, KEEP,KEEP8, DKEEP,
     &               OPELI, DET_EXP, DET_MANT, DET_SIGN )
                 CALL DMUMPS_UPDATE_PROGRESS( OPELI, KEEP8 )
                 IF (IS_BUFRX_ALLOCATED) DEALLOCATE ( BUFRX )
                 NULLIFY(BUFRX)
                 IS_BUFRX_ALLOCATED = .FALSE.
                 CALL MUMPS_PROPINFO( ICNTL, INFO,
     &                COMM_NODES, MYID_NODES )
               ENDIF
                IF ( MYID_NODES .eq. 
     &               MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(38))),
     &                              KEEP(199))
     &             ) THEN
                   IF ( INFO(1) .EQ. -10 .OR. INFO(1) .EQ. -40 ) THEN
                      NTOTPV = NTOTPV + INFO(2)
                   ELSE IF ( INFO(1) .GE. 0 ) THEN
                      NTOTPV = NTOTPV + root%TOT_ROOT_SIZE
                      NMAXNPIV = max(NMAXNPIV,root%TOT_ROOT_SIZE)
                   END IF
                END IF
                IF (INFO(1).GE.0.AND.KEEP(60).EQ.0) THEN
                 IF (root%yes) THEN
                  IF (KEEP(252).EQ.0) THEN
                  IF (KEEP(201).EQ.1) THEN 
                    CALL MUMPS_GETI8(LAFAC, IW(IOLDPS+XXR))
                    LIWFAC    = IW(IOLDPS+XXI)
                    TYPEFile  = TYPEF_L
                    NextPiv2beWritten = 1 
                    MonBloc%INODE    = KEEP(38)   
                    MonBloc%MASTER   = .TRUE.
                    MonBloc%Typenode = 3
                    MonBloc%NROW     = LOCAL_M
                    MonBloc%NCOL     = LOCAL_N
                    MonBloc%NFS      = MonBloc%NCOL
                    MonBloc%Last     = .TRUE.   
                    MonBloc%LastPiv  =  MonBloc%NCOL
                    MonBloc%LastPanelWritten_L=-9999 
                    MonBloc%LastPanelWritten_U=-9999 
                    NULLIFY(MonBloc%INDICES)
                    STRAT        = STRAT_WRITE_MAX
                    MonBloc%Last = .TRUE.
                    LAST_CALL = .TRUE.
                    CALL DMUMPS_OOC_IO_LU_PANEL
     &                                 ( STRAT, TYPEFile,
     &                                  A(PTRFAC(STEP(KEEP(38)))),
     &                                  LAFAC, MonBloc,
     &                                  NextPiv2beWritten, IDUMMY,
     &                                  IW(IOLDPS), LIWFAC,
     &                                  MYID, KEEP8(31), IERR,LAST_CALL)
                    IF (IERR .LT.0) THEN
                      INFO(1) = IERR
                      IF (LP > 0 ) THEN
                        WRITE(LP,*)MYID,
     &                  ': Error in DMUMPS_OOC_IO_LU_PANEL',IERR
                      ENDIF
                    ENDIF
                  ELSE IF (KEEP(201).EQ.2) THEN
                    KEEP8(31)=KEEP8(31)+ ITMP8
                    CALL DMUMPS_NEW_FACTOR(KEEP(38),PTRFAC,
     &              KEEP,KEEP8,A,LA, ITMP8, IERR)
                    IF(IERR.LT.0)THEN
                      INFO(1)=IERR
                      IF (LP > 0 ) THEN
                        WRITE(LP,*)MYID,
     &                  ': Error in DMUMPS_NEW_FACTOR',IERR
                      ENDIF
                    ENDIF
                  ENDIF 
                  ENDIF 
                  IF (KEEP(201).NE.0 .OR. KEEP(252).NE.0) THEN
                     LRLUS = LRLUS + ITMP8 
                     KEEP8(69) = KEEP8(69) - ITMP8 
                     IF (KEEP(252).NE.0) THEN
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,0_8,-ITMP8,
     &                 KEEP,KEEP8,LRLUS)
                     ELSE         
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLUS)
                     ENDIF
                     IF (PTRFAC(STEP(KEEP(38))).EQ.POSFAC-ITMP8) THEN
                       POSFAC = POSFAC  - ITMP8
                       LRLU   = LRLU    + ITMP8
                     ENDIF
                  ELSE
                       CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                 LA-LRLUS
     &                 ,ITMP8,    
     &                 0_8,
     &                 KEEP,KEEP8,LRLUS)
                  ENDIF
                 ENDIF  
                 CALL MUMPS_PROPINFO( ICNTL, INFO,
     &               COMM_NODES, MYID_NODES )
                ENDIF  
                IF ( INFO(1).GE.0 .AND. KEEP(252) .NE. 0 .AND.
     &              (KEEP(60).EQ.0 .OR. KEEP(221).EQ.1)) THEN
                  IF (root%yes) THEN
                    IF (MYID_NODES .EQ. MASTER_ROOT) THEN
                      LRHS_CNTR_MASTER_ROOT = root%TOT_ROOT_SIZE*
     &                                        KEEP(253)
                    ELSE
                      LRHS_CNTR_MASTER_ROOT = 1
                    ENDIF
                    ALLOCATE(root%RHS_CNTR_MASTER_ROOT( 
     &                       LRHS_CNTR_MASTER_ROOT), stat=IERR )
                    IF (IERR.gt.0) THEN
                      INFO(1) = -13
                      INFO(2) = LRHS_CNTR_MASTER_ROOT
                      IF (LP > 0 )
     &                write(LP,*) ' Error allocating, real array ',
     &                'CNTR_MASTER_ROOT of size',
     &                LRHS_CNTR_MASTER_ROOT
                    ENDIF
                  ENDIF 
                  CALL MUMPS_PROPINFO( ICNTL, INFO, COMM_NODES,
     &            MYID_NODES)
                  IF (root%yes .AND. INFO(1).GE.0) THEN
                    FWD_LOCAL_N_RHS = MUMPS_NUMROC(KEEP(253),
     &              root%NBLOCK, root%MYCOL, 0, root%NPCOL)
                    FWD_LOCAL_N_RHS = max(1,FWD_LOCAL_N_RHS)
                    CALL DMUMPS_GATHER_ROOT( MYID_NODES,
     &              root%TOT_ROOT_SIZE, KEEP(253),
     &              root%RHS_CNTR_MASTER_ROOT(1), LOCAL_M,
     &              FWD_LOCAL_N_RHS, root%MBLOCK, root%NBLOCK,
     &              root%RHS_ROOT(1,1), MASTER_ROOT,
     &              root%NPROW, root%NPCOL, COMM_NODES )
                  ENDIF
                ENDIF
            ELSE
                IF (KEEP(19).NE.0) THEN
                  CALL MPI_REDUCE(KEEP(109), GLOBK109, 1,
     &                 MPI_INTEGER, MPI_SUM,
     &                 MASTER_ROOT,
     &                 COMM_NODES, IERR)
                ENDIF
                IF (ROOT_OWNER) THEN
                   IPOSROOT = PTLUST(STEP(KEEP(20)))
                   NFRONT   = IW(IPOSROOT+KEEP(IXSZ)+3)   
                   NFRONT8  = int(NFRONT,8)
                   IPOSROOTROWINDICES=IPOSROOT+6+KEEP(IXSZ)+ 
     &                             IW(IPOSROOT+5+KEEP(IXSZ)) 
                   NTOTPV = NTOTPV + NFRONT 
                   NMAXNPIV = max(NMAXNPIV,NFRONT)
                END IF
                IF ( KEEP(60) .EQ. 0 ) THEN
                 IF ( ROOT_OWNER ) THEN
                  IF (KEEP(118).GE.40) THEN
                   NOFFNEGPV_ROOT = 0        
                   NULLNEGPV_ROOT = 0        
                   NTOTPV_ROOT = 0
                   NB22T1_ROOT = 0
                   NBTINY_ROOT = 0
                   DET_SIGN_ROOT = 1
                   DET_EXP_ROOT = 0
                   DET_MANT_ROOT = 1.0D0
                   DKEEP_SAVE(:) = DKEEP(:)
                   KEEP_SAVE(:)  = KEEP(:)
                   KEEP8_SAVE(:) = KEEP8(:)
                   KEEP_SAVE(201) = 0
                   IF (KEEP(110).EQ.0) THEN
                     KEEP_SAVE(110)= 1        
                     IF (KEEP(118).EQ.40) THEN
                      IF ((DKEEP(10).LE.0).OR.(DKEEP(10).GT.1)) THEN
                             DKEEP_SAVE(1) = DKEEP(9)*1D-1
                      ELSE
                             DKEEP_SAVE(1) = DKEEP(9)*DKEEP(10) 
                      ENDIF    
                     ELSE IF (KEEP(118).EQ.41) THEN
                        DKEEP_SAVE(1) = DKEEP(9)
                     ELSE IF (KEEP(118).EQ.42) THEN
                      IF (DKEEP(13).LT.1)  THEN
                          DKEEP_SAVE(1) = DKEEP(9)*10 
                      ELSE
                         DKEEP_SAVE(1) = DKEEP(9)*DKEEP(13)       
                      ENDIF
                     ENDIF  
                   ELSE 
                     DKEEP_SAVE(1) = DKEEP(9)
                   ENDIF 
                   IS_A_ROOT_SAVE_ALLOCATED = .FALSE.
                   IF (LRLU.GT.NFRONT8*NFRONT8) THEN
                     A_ROOT_SAVE  => A(POSFAC:POSFAC+LRLU-1_8)
                   ELSE
                    IF (associated(A_ROOT_SAVE)) 
     &                               DEALLOCATE(A_ROOT_SAVE)
                    ALLOCATE(A_ROOT_SAVE(NFRONT8*NFRONT8),stat=IERR)
                    IF (IERR.GT.0) THEN
                     INFO(1) = -13
                     CALL MUMPS_SET_IERROR(NFRONT8*NFRONT8, INFO(2) )
                     IF (LP > 0 )
     &                 write(LP,*) "ERROR allocating A_ROOT_SAVE ",
     &                 " of size ", NFRONT*NFRONT
                     GOTO 735
                    ENDIF
                    IS_A_ROOT_SAVE_ALLOCATED = .TRUE.
                   ENDIF
!$        NOMP = OMP_GET_MAX_THREADS()
!$        CHUNK8=int(KEEP(361),8)
!$OMP     PARALLEL DO PRIVATE(I8) SCHEDULE(STATIC, CHUNK8)
!$OMP&    IF ( NFRONT8*NFRONT8 > int(KEEP(361),8) .AND. NOMP .GT. 1)
                   DO I8 =1_8, NFRONT8*NFRONT8
                       A_ROOT_SAVE(I8) = 
     &                   A(PTRAST(STEP(KEEP(20)))+I8-1_8)
                   ENDDO
                   IW(PTLUST(STEP(INODE))+KEEP(IXSZ): 
     &                    PTLUST(STEP(INODE))+KEEP(IXSZ)+4)
     &                 = RECORD_ROOT(KEEP(IXSZ)+1:KEEP(IXSZ)+5)
                   IW(PTLUST(STEP(INODE))+XXLR) = 0
                   IF (KEEP(50).EQ.0) THEN
                    CALL DMUMPS_FAC1_LU_I ( 
     &               N, INODE, IW, LIW, A_ROOT_SAVE(1),
     &               NFRONT8*NFRONT8, IPOSROOT, 1_8,
     &               INFO(1), INFO(2), UU, NOFFNEGPV_ROOT, NTOTPV_ROOT,
     &               NBTINY_ROOT,
     &               DET_EXP_ROOT, DET_MANT_ROOT, DET_SIGN_ROOT,
     &               KEEP_SAVE,KEEP8_SAVE,
     &               STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     &               SEUIL, AVOID_DELAYED, ETATASS,
     &               DKEEP_SAVE(1), PIVNUL_LIST_STRUCT, IWPOS 
     &               , LRGROUPS
     &               , PERM
     &               )
                    IF (INFO(1).LT.0) THEN
                     IF (LP.GT.0) 
     &               write(LP,*) "ERROR after DMUMPS_FAC1_LU ",
     &               "on the root INFO(1)= ", INFO(1)
                     GOTO 735
                    ENDIF 
                   ELSE
                    CALL DMUMPS_FAC1_LDLT_I (N,KEEP_SAVE(20),
     &               IW, LIW, A_ROOT_SAVE(1), NFRONT8*NFRONT8,
     &               IPOSROOT, 1_8,
     &               INFO(1), INFO(2), UU,
     &               NOFFNEGPV_ROOT, NULLNEGPV_ROOT, NTOTPV_ROOT,
     &               NB22T1_ROOT, NBTINY_ROOT, 
     &               DET_EXP_ROOT, DET_MANT_ROOT, DET_SIGN_ROOT,
     &               KEEP_SAVE,KEEP8_SAVE, MYID_NODES, SEUIL, 
     &               AVOID_DELAYED, ETATASS, DKEEP_SAVE(1),
     &               PIVNUL_LIST_STRUCT, IWPOS
     &               , LRGROUPS
     &               , PERM
     &               )
                     IF (INFO(1).LT.0) THEN
                      IF (LP.GT.0) 
     &                write(LP,*) "ERROR after DMUMPS_FAC1_LDLT ",
     &                "on the root INFO(1)= ", INFO(1)
                      GOTO 735
                     ENDIF 
                   ENDIF 
                   LRecord = IW(IOLDPS+XXI)
                   IW(PTLUST(STEP(INODE)):
     &                PTLUST(STEP(INODE))+LRecord-1) = 
     &                                    RECORD_ROOT(1:LRecord) 
                   IW(PTLUST(STEP(INODE))+KEEP(IXSZ): 
     &                PTLUST(STEP(INODE))+KEEP(IXSZ)+4) =
     &                                          Header_ROOT(1:5)
                   KEEP17_LU = KEEP_SAVE(109)-KEEP(109)
                   IF (KEEP_SAVE(109).GT.KEEP(109)) THEN
                     K = 1
                     DO I = KEEP(109)+1, KEEP(109)+KEEP17_LU
                        RECORD_ROOT(K) = 
     &                          PIVNUL_LIST_STRUCT%PIVNUL_LIST(I)
                        K = K+1
                     ENDDO
                   ENDIF
                   IF (IS_A_ROOT_SAVE_ALLOCATED) DEALLOCATE(A_ROOT_SAVE)
                   NULLIFY(A_ROOT_SAVE)
                   IS_A_ROOT_SAVE_ALLOCATED = .FALSE.
                   DET_SIGN = DET_SIGN * DET_SIGN_ROOT
                   DET_EXP  = DET_EXP + DET_EXP_ROOT
                   CALL DMUMPS_UPDATEDETER ( DET_MANT_ROOT, 
     &                                       DET_MANT, DET_EXP)
                   NOFFNEGPV = NOFFNEGPV + NOFFNEGPV_ROOT
                   NULLNEGPV = NULLNEGPV + NULLNEGPV_ROOT
                  ENDIF 
                  LOCAL_M = 0
                  LOCAL_N = 0
                  MBLOCK  = 0
                  NBLOCK  = 0
                  PHASE   = 1  
                  CALL DMUMPS_SVD_QR_ESTIM_WK( PHASE,
     &             MBLOCK, NBLOCK, NFRONT, LOCAL_M, LOCAL_N,
     &             ROOT_OWNER, KEEP,KEEP8,
     &             LIWK_RR, LWK_RR )    
                  LBUFRX             = LWK_RR
                  IS_BUFRX_ALLOCATED = .FALSE.
                  IF ( LRLU .GT. LBUFRX ) THEN
                      BUFRX => A(POSFAC:POSFAC+LBUFRX-1_8)
                  ELSE
                      ALLOCATE( BUFRX( LBUFRX ), stat = IERR )
                      IF (IERR.gt.0) THEN
                        INFO(1) = -13
                        CALL MUMPS_SET_IERROR(LBUFRX, INFO(2))
                        IF (LP.GT.0)
     &                  write(LP,*) ' Error allocating, real
     &                        array ','of size ', LBUFRX,
     &                        ' before DMUMPS_SEQ_FACTO_ROOT_SVD_QR'
                        GOTO 735
                      ENDIF
                      IS_BUFRX_ALLOCATED = .TRUE.
                  ENDIF
                 IF (PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST .LT. 
     &               KEEP(109)+NFRONT) THEN
                  CALL MUMPS_RESIZE_PIVNUL(KEEP, N, PIVNUL_LIST_STRUCT, 
     &                             KEEP(109)+NFRONT, INFO(1), INFO(2) ) 
                  IF (INFO(1).LT.0) GOTO 735
                 ENDIF
                  CALL DMUMPS_SEQ_FACTO_ROOT_SVD_QR(
     &                 NFRONT,A(PTRAST(STEP(KEEP(20)))),
     &                 root,
     &                 BUFRX(1), int(LBUFRX),
     &                 KEEP,KEEP8, INFO, LP, DKEEP,
     &                 GLOBK109, OPELI,
     &                 PIVNUL_LIST_STRUCT%PIVNUL_LIST(KEEP(109)+1),
     &                 PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST- KEEP(109), 
     &                 IW(IPOSROOTROWINDICES))
                  IF (IS_BUFRX_ALLOCATED) DEALLOCATE ( BUFRX )
                  NULLIFY(BUFRX)
                  IS_BUFRX_ALLOCATED = .FALSE.
                  IF (INFO(1).LT.0) GOTO 735
                  IF (MPA.GT.0)
     &            CALL MUMPS_PRINT_STILL_ACTIVE (MYID_NODES, KEEP, 
     &                     dble(DKEEP(17)), 
     &                     OPELI,
     &                     OPLAST_PRINTED, MPA)
                  CALL DMUMPS_UPDATE_PROGRESS( OPELI, KEEP8 )
                  KEEP(143) = KEEP17_LU
                  IF (KEEP(118).GE.40) THEN
                    K = 1 
                    IF (KEEP(17).GT.0) THEN
                     DO I = KEEP(109)+1, KEEP(109)+KEEP(17)
                      IF ( K .GT. KEEP17_LU ) THEN
                        PIVNUL_LIST_STRUCT%PIVNUL_LIST(I) = -1
                      ELSE
                        PIVNUL_LIST_STRUCT%PIVNUL_LIST(I) = 
     &                                 RECORD_ROOT(K)
                      ENDIF
                      K = K+1
                     ENDDO
                    ENDIF
                  ENDIF
                  IF (allocated(RECORD_ROOT)) DEALLOCATE(RECORD_ROOT)
                      IF (KEEP(201).EQ.1) THEN 
                        CALL MUMPS_GETI8(LAFAC, IW(IPOSROOT+XXR))
                        LIWFAC    = IW(IPOSROOT+XXI)
                        TYPEFile  = TYPEF_L
                        NextPiv2beWritten = 1 
                        MonBloc%INODE    = KEEP(20)   
                        MonBloc%MASTER   = .TRUE.
                        MonBloc%Typenode = 3
                        MonBloc%NCOL     = NFRONT
                        MonBloc%NROW     = NFRONT
                        MonBloc%NFS      = NFRONT
                        MonBloc%Last     = .TRUE.   
                        MonBloc%LastPiv  = MonBloc%NCOL
                        MonBloc%LastPanelWritten_L=-9999 
                        MonBloc%LastPanelWritten_U=-9999 
                        NULLIFY(MonBloc%INDICES)
                        STRAT        = STRAT_WRITE_MAX
                        MonBloc%Last = .TRUE.
                        LAST_CALL    = .TRUE.
                        CALL DMUMPS_OOC_IO_LU_PANEL
     &                                 ( STRAT, TYPEFile,
     &                                  A(PTRAST(STEP(KEEP(20)))),
     &                                  LAFAC, MonBloc,
     &                                  NextPiv2beWritten, IDUMMY,
     &                                  IW(IPOSROOT), LIWFAC,
     &                              MYID, KEEP8(31), IERR, LAST_CALL)
                        IF(IERR.LT.0)THEN
                        IF (LP > 0)
     &                  WRITE(LP,*)MYID,
     &                  ': Error raised in DMUMPS_OOC_IO_LU_PANEL',
     &                  IERR
                        INFO(1)=IERR
                        ENDIF
                      ELSE IF (KEEP(201).EQ.2) THEN 
                        KEEP8(31)=KEEP8(31)+NFRONT8*NFRONT8
                        CALL DMUMPS_NEW_FACTOR(KEEP(20),PTRFAC,
     &                  KEEP,KEEP8,A,LA, NFRONT8*NFRONT8, IERR)
                        IF(IERR.LT.0)THEN
                        WRITE(*,*)MYID,
     &                  ': Internal error in DMUMPS_NEW_FACTOR',
     &                  IERR
                        GOTO 735
                        ENDIF
                      ENDIF
                      ITMP8  =  NFRONT8*NFRONT8
                      IF(KEEP(201).NE.0)THEN 
                        IF (PTRFAC(STEP(KEEP(20))).EQ.
     &                      POSFAC-ITMP8) THEN
                          POSFAC = POSFAC  - ITMP8
                          LRLU   = LRLU    + ITMP8
                          LRLUS  = LRLUS   + ITMP8
                          KEEP8(69) = KEEP8(69) - ITMP8
                        ELSE
                          IF (LP.GT.0) 
     &                    WRITE(LP,*) "Internal error",
     &                    POSFAC,NFRONT8,
     &                    "root KEEP(20) not on top in OOC"
                          GOTO 735
                         ENDIF
                      ENDIF
                      CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                LA-LRLUS,ITMP8,0_8,KEEP,KEEP8,LRLUS)
                ENDIF
 735         CONTINUE
             CALL MUMPS_PROPINFO( ICNTL, INFO, COMM_NODES, MYID_NODES )
             IF (INFO(1).LT.0) GOTO 745
                 CALL MPI_BCAST( KEEP(17), 1, MPI_INTEGER,
     &                MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(20))),
     &                               KEEP(199)),
     &                COMM_NODES, IERR )
                 CALL MPI_BCAST( KEEP(143), 1, MPI_INTEGER,
     &                MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(20))),
     &                               KEEP(199)),
     &                COMM_NODES, IERR )
                END IF
                IF (ROOT_OWNER.AND.KEEP(60).NE.0) THEN 
                  ITMP8 = NFRONT8*NFRONT8
                  IF ( PTRFAC(STEP(KEEP(20))) .EQ. POSFAC -
     &                 ITMP8 ) THEN
                    POSFAC = POSFAC - ITMP8
                    LRLUS  = LRLUS  + ITMP8
                    LRLU   = LRLUS  + ITMP8
                    KEEP8(69) = KEEP8(69) - ITMP8
                    CALL DMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &              LA-LRLUS,0_8,-ITMP8,KEEP,KEEP8,LRLUS)
                  ENDIF
                ENDIF
            END IF
            GOTO 750
 745  CONTINUE
            IF (allocated(RECORD_ROOT)) DEALLOCATE(RECORD_ROOT)
            IF (IS_A_ROOT_SAVE_ALLOCATED) DEALLOCATE(A_ROOT_SAVE)
            NULLIFY(A_ROOT_SAVE)
 750  CONTINUE
            IF (INFO(1).LT.0) GOTO 500
            CALL MUMPS_SECFIN(ROOTTIME)
            DKEEP(99)=ROOTTIME
          END IF  
       END IF 
 500   CONTINUE
       IF ( KEEP(38) .NE. 0 ) THEN
         IF (MYID_NODES.EQ. 
     &        MUMPS_PROCNODE(PROCNODE_STEPS(STEP(KEEP(38))),KEEP(199))
     &      ) THEN
           MAXFRT = max ( MAXFRT, root%TOT_ROOT_SIZE)
         END IF
       END IF
      RETURN
      END SUBROUTINE DMUMPS_FAC_PAR
      SUBROUTINE DMUMPS_CHANGE_HEADER( HEADER, KEEP253 )
        INTEGER HEADER( 6 ), KEEP253
        INTEGER NFRONT, NASS
        NFRONT = HEADER(1)
        IF ( HEADER(2) .ne. 0 ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 1 :',HEADER(2)
          CALL MUMPS_ABORT()
        END IF
        NASS   = abs( HEADER( 3 ) )
        IF ( NASS .NE. abs( HEADER( 4 ) ) ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 2 :',HEADER(3:4)
          CALL MUMPS_ABORT()
        END IF
        IF ( NASS+KEEP253 .NE. NFRONT ) THEN
          WRITE(*,*) ' *** CHG_HEADER ERROR 3 : not root',
     &    NASS, KEEP253, NFRONT
          CALL MUMPS_ABORT()
        END IF
        HEADER( 1 ) = KEEP253 
        HEADER( 2 ) = 0
        HEADER( 3 ) = NFRONT 
        HEADER( 4 ) = NFRONT-KEEP253    
        RETURN
      END SUBROUTINE DMUMPS_CHANGE_HEADER
      END MODULE DMUMPS_FAC_PAR_M
