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
      RECURSIVE SUBROUTINE SMUMPS_PROCESS_SYM_BLOCFACTO( 
     &   COMM_LOAD, ASS_IRECV,
     &   BUFR, LBUFR,
     &   LBUFR_BYTES, PROCNODE_STEPS, MSGSOU,
     &   SLAVEF, IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, N, IW, LIW,
     &   A, LA, PTRIST, PTRAST, NSTK_S, PERM,
     &   COMP, STEP, PIMASTER, PAMASTER, POSFAC,
     &   MYID, COMM, IFLAG, IERROR, NBFIN,
     &
     &    PTLUST_S, PTRFAC, root, OPASSW, OPELIW,
     &    ITLOC, RHS_MUMPS, FILS, DAD,
     &    PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,
     &    ICNTL,KEEP,KEEP8,DKEEP,IPOOL, LPOOL, LEAF, ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &    )
      USE SMUMPS_OOC, ONLY : IO_BLOCK
      USE MUMPS_OOC_COMMON, ONLY : TYPEF_L, STRAT_TRY_WRITE
      USE SMUMPS_LOAD
      USE SMUMPS_BUF
      USE SMUMPS_LR_CORE
      USE SMUMPS_LR_TYPE
      USE SMUMPS_LR_STATS
      USE SMUMPS_FAC_LR
      USE SMUMPS_ANA_LR, ONLY : GET_CUT
      USE SMUMPS_LR_DATA_M
      USE SMUMPS_STRUC_DEF, ONLY : SMUMPS_ROOT_STRUC
      USE SMUMPS_DYNAMIC_MEMORY_M, ONLY : SMUMPS_DM_SET_DYNPTR
      USE SMUMPS_FAC_FRONT_AUX_M, ONLY : SMUMPS_GET_SIZE_SCHUR_IN_FRONT
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      TYPE (SMUMPS_ROOT_STRUC) :: root
      INTEGER ICNTL( 60 ), KEEP( 500 )
      INTEGER(8) KEEP8(150)
      REAL    DKEEP(230)
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER N, SLAVEF, IWPOS, IWPOSCB, LIW
      INTEGER(8) IPTRLU, LRLU, LRLUS, LA, POSFAC
      INTEGER COMP
      INTEGER IFLAG, IERROR, NBFIN, MSGSOU
      INTEGER PROCNODE_STEPS(KEEP(28)), PTRIST(KEEP(28)),
     &        NSTK_S(KEEP(28))
      INTEGER(8) PTRAST(KEEP(28)), PTRFAC(KEEP(28)), PAMASTER(KEEP(28))
      INTEGER PERM(N), STEP(N), 
     & PIMASTER(KEEP(28))
      INTEGER IW( LIW )
      REAL A( LA )
      INTEGER, intent(in) :: LRGROUPS(KEEP(280))
      INTEGER LPTRAR, NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER COMM, MYID
      INTEGER PTLUST_S(KEEP(28)),
     &        ITLOC(N+KEEP(253)), FILS(N), DAD(KEEP(28)), ND(KEEP(28))
      REAL :: RHS_MUMPS(KEEP8(85))
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER FRERE_STEPS(KEEP(28))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION FLOP1
      INTEGER INTARR( KEEP8(27) )
      REAL DBLARR( KEEP8(26) )
      INTEGER LEAF, LPOOL 
      INTEGER IPOOL( LPOOL )
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER PIVI
      INTEGER (8) POSPV1,POSPV2,OFFDAG,LPOS1
      INTEGER J2
      REAL MULT1,MULT2, A11, DETPIV, A22, A12
      INTEGER :: NFS4FATHER, NVSCHUR_K253, NSLAVES_L, IROW_L
      REAL, ALLOCATABLE, DIMENSION(:) :: M_ARRAY
      INTEGER NBROWSinF
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER LP
      INTEGER INODE, POSITION, NPIV, IERR
      INTEGER NCOL
      INTEGER :: IFLAG_OOC
      INTEGER(8) :: POSBLOCFACTO
      INTEGER :: LD_BLOCFACTO 
      INTEGER(8) :: LA_BLOCFACTO 
      INTEGER(8) :: LA_PTR 
      INTEGER(8) :: POSELT
      REAL, DIMENSION(:), POINTER :: A_PTR
      INTEGER IOLDPS, LCONT1, NASS1, NROW1, NCOL1, NPIV1
      INTEGER NSLAV1, HS, ISW, DEST
      INTEGER ICT11
      INTEGER(8) LPOS, LPOS2, DPOS, UPOS
      INTEGER (8) IPOS, KPOS
      INTEGER I, IPIV, FPERE, NSLAVES_TOT,
     &        NSLAVES_FOLLOW, NB_BLOC_FAC
      INTEGER IPOSK, JPOSK, NPIVSENT, Block, IROW, BLSIZE
      INTEGER allocok, TO_UPDATE_CPT_END
      REAL, DIMENSION(:), ALLOCATABLE :: UIP21K
      REAL, DIMENSION(:), ALLOCATABLE :: DYN_BLOCFACTO
      INTEGER, DIMENSION(:), ALLOCATABLE :: LIST_SLAVES_FOLLOW
      INTEGER, DIMENSION(:), ALLOCATABLE :: DYN_PIVINFO
      LOGICAL LASTBL
      LOGICAL BLOCKING, SET_IRECV, MESSAGE_RECEIVED
      REAL ONE,ALPHA
      PARAMETER (ONE = 1.0E0, ALPHA=-1.0E0)
      INTEGER LIWFAC, STRAT, NextPivDummy
      TYPE(IO_BLOCK) :: MonBloc
      LOGICAL LAST_CALL
      INTEGER LRELAY_INFO
      LOGICAL COUNTER_WAS_HUGE
      INTEGER TO_UPDATE_CPT_RECUR
      INTEGER :: LR_ACTIVATED_INT
      LOGICAL :: LR_ACTIVATED, COMPRESS_CB, COMPRESS_PANEL
      LOGICAL :: DYNAMIC_ALLOC
      LOGICAL OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER :: CURRENT_BLR, NSLAVES_PREC, INFO_TMP(2)
      INTEGER :: NELIM, NB_BLR_LM, NB_BLR_LS,  
     &           MAXI_CLUSTER_LM, MAXI_CLUSTER_LS, MAXI_CLUSTER, 
     &           NPARTSASS, NPARTSCB, NPARTSCB_COL, NPARTSASS_COL, 
     &           NB_BLR_COL, MAXI_CLUSTER_COL
       INTEGER :: NPARTSASS_MASTER, IPANEL, NB_ACCESSES_INIT
      TYPE (LRB_TYPE), DIMENSION(:), ALLOCATABLE :: BLR_LM 
      TYPE (LRB_TYPE), DIMENSION(:), POINTER     :: BLR_LS
      TYPE(LRB_TYPE), POINTER, DIMENSION(:,:) :: CB_LRB
      INTEGER, POINTER, DIMENSION(:) :: BEGS_BLR_LM, BEGS_BLR_LS, 
     &                                  BEGS_BLR_COL, BEGS_BLR_COL_TMP
      LOGICAL KEEP_BEGS_BLR_LS, KEEP_BEGS_BLR_COL, KEEP_BLR_LS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK, TAU
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JPVT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BLOCKLR
      REAL,ALLOCATABLE,DIMENSION(:) :: RWORK
      INTEGER :: OMP_NUM, LWORK
      INTEGER :: II,JJ
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      LP = ICNTL(1)
      IF (ICNTL(4) .LE. 0) LP = -1
      POSITION = 0
      TO_UPDATE_CPT_END = -654321
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, INODE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NPIV, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, FPERE, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LASTBL = (NPIV.LE.0)
      IF (LASTBL) THEN 
         NPIV = -NPIV
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
         CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NB_BLOC_FAC, 1,
     &                 MPI_INTEGER, COMM, IERR )
      ENDIF
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NCOL, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NELIM, 1,
     &                 MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 NPARTSASS_MASTER, 1,
     &                 MPI_INTEGER, COMM, IERR )
       NPARTSASS_COL = NPARTSASS_MASTER
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, IPANEL,
     &                 1, MPI_INTEGER, COMM, IERR )
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, LR_ACTIVATED_INT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      LR_ACTIVATED    = (LR_ACTIVATED_INT.EQ.1)
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, NSLAVES_TOT, 1,
     &                 MPI_INTEGER, COMM, IERR )
      KEEP_BEGS_BLR_LS  =.FALSE. 
      NULLIFY(BEGS_BLR_LS)
      KEEP_BEGS_BLR_COL =.FALSE.
      NULLIFY(BEGS_BLR_COL)
      KEEP_BLR_LS       =.FALSE.
      NULLIFY(BLR_LS)
      NULLIFY(BEGS_BLR_LM)
      IF ( LR_ACTIVATED ) THEN
        LA_BLOCFACTO = int(NPIV,8) * int(NPIV+NELIM,8)
        LD_BLOCFACTO = max(NPIV+NELIM,1)
      ELSE
        LA_BLOCFACTO = int(NPIV,8) * int(NCOL,8)
        LD_BLOCFACTO = max(NCOL,1)
      ENDIF
      IF (LR_ACTIVATED) THEN
        DYNAMIC_ALLOC = .TRUE.
      ELSE
        DYNAMIC_ALLOC = .FALSE.
      ENDIF
      IF ( .NOT. DYNAMIC_ALLOC ) THEN
        IF ( NPIV .EQ. 0 ) THEN
          IPIV = 1  
          POSBLOCFACTO = 1_8 
        ELSE 
          CALL SMUMPS_GET_SIZE_NEEDED(
     &       NPIV, LA_BLOCFACTO, .FALSE.,
     &       KEEP(1), KEEP8(1), 
     &       N, IW, LIW, A, LA,
     &       LRLU, IPTRLU,
     &       IWPOS, IWPOSCB, PTRIST, PTRAST,
     &       STEP, PIMASTER, PAMASTER, LRLUS,
     &       KEEP(IXSZ),COMP,DKEEP(97),
     &       MYID, SLAVEF, PROCNODE_STEPS, DAD,
     &       IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
          LRLU  = LRLU - LA_BLOCFACTO
          LRLUS = LRLUS - LA_BLOCFACTO
          KEEP8(69) = KEEP8(69) + LA_BLOCFACTO
          KEEP8(67) = min(LRLUS, KEEP8(67))
          KEEP8(68) = max(KEEP8(69), KEEP8(68))
          POSBLOCFACTO = POSFAC
          POSFAC = POSFAC + LA_BLOCFACTO
          IPIV = IWPOS
          IWPOS = IWPOS + NPIV
          CALL SMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &               LA-LRLUS,0_8,LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        ENDIF
      ELSE 
        ALLOCATE(DYN_PIVINFO(max(1,NPIV)),
     &           DYN_BLOCFACTO(max(1_8,LA_BLOCFACTO)),
     &           stat=allocok)
        IF (allocok.GT.0) THEN
           IF (LP > 0 ) WRITE(LP,*) MYID,
     &     ": ALLOCATION FAILURE FOR DYN_PIVINFO and DYN_BLOCFACTO IN ",
     &     "SMUMPS_PROCESS_SYM_BLOCFACTO"
           IFLAG = -13
           CALL MUMPS_SET_IERROR(max(1_8,LA_BLOCFACTO), IERROR)
           GOTO 700
        ENDIF
        POSBLOCFACTO = 1_8
        IPIV = 1
      ENDIF
      IF (NPIV.GT.0) THEN
        IF (DYNAMIC_ALLOC) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 DYN_PIVINFO, NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 IW( IPIV ), NPIV,
     &                 MPI_INTEGER, COMM, IERR )
        ENDIF
        IF (DYNAMIC_ALLOC) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                DYN_BLOCFACTO, int(LA_BLOCFACTO),
     &                 MPI_REAL,
     &                 COMM, IERR )
        ELSE
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 A(POSBLOCFACTO), int(LA_BLOCFACTO),
     &                 MPI_REAL,
     &                 COMM, IERR )
        ENDIF
        IF ( LR_ACTIVATED ) THEN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                 NB_BLR_LM, 1, MPI_INTEGER,
     &                 COMM, IERR )
          ALLOCATE(BLR_LM(max(NB_BLR_LM,1)), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
             IF (LP > 0 ) WRITE(LP,*) MYID,
     &            ": ALLOCATION FAILURE FOR BLR_LM IN ",
     &            "SMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -13
             IERROR = max(NB_BLR_LM,1)
             GOTO 700
          END IF
          ALLOCATE(BEGS_BLR_LM(NB_BLR_LM+2), stat=allocok)
          IF ( allocok .GT. 0 ) THEN
             IF (LP > 0 ) WRITE(LP,*) MYID,
     &            ": ALLOCATION FAILURE FOR BEGS_BLR_LM IN ",
     &            "SMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -13
             IERROR = NB_BLR_LM+2
             GOTO 700
          END IF
          CALL SMUMPS_MPI_UNPACK_LR(
     &          BUFR, LBUFR, LBUFR_BYTES, POSITION, NPIV, NELIM, 
     &          'V', BLR_LM, NB_BLR_LM,
     &          BEGS_BLR_LM(1), KEEP8, COMM, IERR, IFLAG, IERROR)
          IF (IFLAG.LT.0) GOTO 700
        ENDIF 
      ENDIF 
      CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION, 
     &                 LRELAY_INFO, 1,
     &                 MPI_INTEGER, COMM, IERR )
      IF (PTRIST(STEP( INODE )) .EQ. 0) THEN
          CALL SMUMPS_TREAT_DESCBAND( INODE, COMM_LOAD, ASS_IRECV,
     &      BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &      IWPOS, IWPOSCB, IPTRLU,
     &      LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &      PTLUST_S, PTRFAC,
     &      PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &      IFLAG, IERROR, COMM,
     &      PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &      root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &      FILS, DAD, PTRARW, PTRAIW,
     &      PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &      INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &      LPTRAR, NELT, FRTPTR, FRTELT, 
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE.
     &               , LRGROUPS
     &        )
          IF ( IFLAG .LT. 0 ) GOTO 550
      ENDIF
      IF ( IW( PTRIST(STEP(INODE)) + 3 + KEEP(IXSZ)) .EQ. 0 ) THEN
       DO WHILE ( IW(PTRIST(STEP(INODE)) + XXNBPR) .NE. 0)
        BLOCKING  = .TRUE.
        SET_IRECV = .FALSE.
        MESSAGE_RECEIVED = .FALSE.
        CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD,
     &    ASS_IRECV, BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, CONTRIB_TYPE2,
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR, ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &      )
        IF ( IFLAG .LT. 0 ) GOTO 550
      END  DO
      ENDIF
        SET_IRECV = .TRUE.
        BLOCKING  = .FALSE.
        MESSAGE_RECEIVED = .TRUE.
        CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &    BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &    MPI_ANY_SOURCE, MPI_ANY_TAG, 
     &    STATUS,
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU,
     &    LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &    PTLUST_S, PTRFAC,
     &    PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &    IFLAG, IERROR, COMM,
     &    PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &
     &    root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &    FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT,
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &    )
      IOLDPS = PTRIST(STEP(INODE))
      CALL SMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      LCONT1 = IW( IOLDPS + KEEP(IXSZ))
      NASS1  = IW( IOLDPS + 1 + KEEP(IXSZ))
      COMPRESS_PANEL = (IW(IOLDPS+XXLR).GE.2)
      OOCWRITE_COMPATIBLE_WITH_BLR = 
     &          ( .NOT.LR_ACTIVATED.OR.  (.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(486).NE.2) 
     &          )
      IF ( NASS1 < 0 ) THEN
        NASS1 = -NASS1
        IW( IOLDPS + 1 + KEEP(IXSZ)) = NASS1
        IF (KEEP(55) .EQ. 0) THEN 
          CALL SMUMPS_ASM_SLAVE_ARROWHEADS(INODE, STEP(INODE),
     &       N, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS,
     &       PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &       INTARR, DBLARR, KEEP8(27), KEEP8(26), RHS_MUMPS, 
     &       LRGROUPS)
        ELSE
          CALL SMUMPS_ASM_SLAVE_ELEMENTS(INODE, N, NELT, IW, LIW,
     &       IOLDPS, A_PTR(POSELT), LA_PTR, 1_8, KEEP, KEEP8, ITLOC,
     &       FILS, PTRAIW,
     &       PTRARW, INTARR, DBLARR, KEEP8(27), KEEP8(26),
     &       FRTPTR, FRTELT, RHS_MUMPS, LRGROUPS)
        ENDIF
      ENDIF
      NROW1  = IW( IOLDPS + 2 +KEEP(IXSZ))
      NPIV1  = IW( IOLDPS + 3 +KEEP(IXSZ))
      NSLAV1 = IW( IOLDPS + 5 + KEEP(IXSZ))
      NSLAVES_FOLLOW = NSLAV1 - XTRA_SLAVES_SYM
      HS     = 6 + NSLAV1 + KEEP(IXSZ)
      NCOL1  = LCONT1 + NPIV1
      IF ( LASTBL ) THEN
        TO_UPDATE_CPT_END = ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) * 
     &                       NB_BLOC_FAC
      END IF
      IF (NPIV.GT.0) THEN
        ICT11 = IOLDPS+HS+NROW1+NPIV1 - 1
        DO I = 1, NPIV
          IF (DYNAMIC_ALLOC) THEN
            PIVI = abs(DYN_PIVINFO(I))
          ELSE
            PIVI = abs(IW(IPIV+I-1))
          ENDIF
          IF (PIVI.EQ.I) CYCLE
          ISW = IW(ICT11+I)
          IW(ICT11+I) = IW(ICT11+PIVI)
          IW(ICT11+PIVI) = ISW
          IPOS = POSELT + int(NPIV1 + I - 1,8)
          KPOS = POSELT + int(NPIV1 + PIVI - 1,8)
          CALL sswap(NROW1, A_PTR(IPOS), NCOL1, A_PTR(KPOS), NCOL1)
        ENDDO
        IF (.NOT.LR_ACTIVATED) THEN
          ALLOCATE( UIP21K( NPIV * NROW1 ), stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN SMUMPS_PROCESS_SYM_BLOCFACTO"
            IFLAG = -13
            IERROR = NPIV * NROW1
            GOTO 700
          END IF
        ELSE
         ALLOCATE( UIP21K( 1 ), stat = allocok )
         IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR UIP21K IN SMUMPS_PROCESS_SYM_BLOCFACTO"
          IFLAG = -13
          IERROR = NPIV * 1
          GOTO 700
        END IF
        ENDIF
        IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
          ALLOCATE( LIST_SLAVES_FOLLOW ( NSLAVES_FOLLOW ),
     &            stat = allocok )
          IF ( allocok .GT. 0 ) THEN
            IF (LP > 0 ) WRITE(LP,*) MYID,
     &": ALLOCATION FAILURE FOR LIST_SLAVES_FOLLOW
     & IN SMUMPS_PROCESS_SYM_BLOCFACTO"
            IFLAG = -13
            IERROR = NSLAVES_FOLLOW
            GOTO 700
          END IF
          LIST_SLAVES_FOLLOW(1:NSLAVES_FOLLOW)=
     &    IW(IOLDPS+6+XTRA_SLAVES_SYM+KEEP(IXSZ):
     &     IOLDPS+5+XTRA_SLAVES_SYM+KEEP(IXSZ)+NSLAVES_FOLLOW)
        END IF
        IF ((.NOT. LR_ACTIVATED).OR.KEEP(475).EQ.0) THEN
            IF (DYNAMIC_ALLOC) THEN
              CALL strsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &             DYN_BLOCFACTO, LD_BLOCFACTO,
     &             A_PTR(POSELT+int(NPIV1,8)), NCOL1)
            ELSE
              CALL strsm( 'L', 'U', 'T', 'U', NPIV, NROW1, ONE,
     &             A( POSBLOCFACTO ), LD_BLOCFACTO,
     &             A_PTR(POSELT+int(NPIV1,8)), NCOL1)
            ENDIF
        ENDIF
        IF (.NOT.LR_ACTIVATED) THEN
         LPOS = POSELT + int(NPIV1,8)
         UPOS = 1_8
         DO I = 1, NROW1
          UIP21K( UPOS: UPOS + int(NPIV-1,8) ) = 
     &                       A_PTR(LPOS: LPOS+int(NPIV-1,8))
          LPOS = LPOS + int(NCOL1,8)
          UPOS = UPOS + int(NPIV,8)
         END DO
        ENDIF
        IF ((.NOT. LR_ACTIVATED).OR.KEEP(475).EQ.0) THEN
          LPOS = POSELT + int(NPIV1,8)
          IF (DYNAMIC_ALLOC) THEN
            DPOS = 1_8
          ELSE
            DPOS = POSBLOCFACTO
          ENDIF
          I = 1
          DO
            IF(I .GT. NPIV) EXIT
            IF (DYNAMIC_ALLOC) THEN
              PIVI = DYN_PIVINFO(I)
            ELSE
              PIVI = IW(IPIV+I-1)
            ENDIF
            IF(PIVI .GT. 0) THEN
              IF (DYNAMIC_ALLOC) THEN
                A11 = ONE/DYN_BLOCFACTO(DPOS)
              ELSE
                A11 = ONE/A(DPOS)
              ENDIF
              CALL sscal( NROW1, A11, A_PTR(LPOS), NCOL1 )
              LPOS = LPOS + 1_8
              DPOS = DPOS + int(LD_BLOCFACTO + 1,8)
              I = I+1
            ELSE
              POSPV1 = DPOS
              POSPV2 = DPOS+ int(LD_BLOCFACTO + 1,8)
              OFFDAG = POSPV1+1_8
              IF (DYNAMIC_ALLOC) THEN
                A11 = DYN_BLOCFACTO(POSPV1)
                A22 = DYN_BLOCFACTO(POSPV2)
                A12 = DYN_BLOCFACTO(OFFDAG)
                DETPIV = A11*A22 - A12**2
                A22 = A11/DETPIV
                A11 = DYN_BLOCFACTO(POSPV2)/DETPIV
                A12 = -A12/DETPIV
              ELSE
                A11 = A(POSPV1)
                A22 = A(POSPV2)
                A12 = A(OFFDAG)
                DETPIV = A11*A22 - A12**2
                A22 = A11/DETPIV
                A11 = A(POSPV2)/DETPIV
                A12 = -A12/DETPIV
              ENDIF
              LPOS1 = LPOS
              DO J2 = 1,NROW1
                 MULT1 = A11*A_PTR(LPOS1)+A12*A_PTR(LPOS1+1_8)
                 MULT2 = A12*A_PTR(LPOS1)+A22*A_PTR(LPOS1+1_8)
                 A_PTR(LPOS1) = MULT1
                 A_PTR(LPOS1+1_8) = MULT2
                 LPOS1 = LPOS1 + int(NCOL1,8)
              ENDDO
              LPOS = LPOS + 2_8
              DPOS = POSPV2 + int(LD_BLOCFACTO + 1,8)
              I = I+2
            ENDIF
          ENDDO
        ENDIF
      ENDIF 
      COMPRESS_CB = .FALSE.
      IF ( LR_ACTIVATED) THEN
        NSLAVES_PREC = NSLAVES_TOT - NSLAVES_FOLLOW -1
        COMPRESS_CB    = ((IW(IOLDPS+XXLR).EQ.1).OR.
     &                    (IW(IOLDPS+XXLR).EQ.3))
      ENDIF
      IF (COMPRESS_CB.AND.NPIV.EQ.0) THEN
           COMPRESS_CB = .FALSE.
           IW(IOLDPS+XXLR) = IW(IOLDPS+XXLR) -1
      ENDIF
      IF (NPIV.GT.0) THEN
         IF (NROW1.LE.0) THEN
            CALL MUMPS_ABORT()  
         ENDIF
       IF (LR_ACTIVATED) THEN
        IF (NPIV1.NE.0) THEN
           CALL SMUMPS_BLR_RETRIEVE_BEGS_BLR_L (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_LS)
           KEEP_BEGS_BLR_LS = .TRUE.  
           NB_BLR_LS = size(BEGS_BLR_LS) - 2
           NPARTSCB  = NB_BLR_LS
        ELSE
          CALL GET_CUT(IW(IOLDPS+HS:IOLDPS+HS+NROW1-1), 0,
     &                    NROW1, LRGROUPS, NPARTSCB, 
     &                    NPARTSASS, BEGS_BLR_LS)
              CALL REGROUPING2(BEGS_BLR_LS, NPARTSASS, NASS1, NPARTSCB,
     &                        NROW1-0, KEEP(488), .TRUE., KEEP(472),
     &                        NCOL1, KEEP(1))
             NB_BLR_LS = NPARTSCB
        ENDIF
        call MAX_CLUSTER(BEGS_BLR_LM,NB_BLR_LM+1,MAXI_CLUSTER_LM)
        call MAX_CLUSTER(BEGS_BLR_LS,NB_BLR_LS+1,MAXI_CLUSTER_LS)
        MAXI_CLUSTER=max(MAXI_CLUSTER_LS,MAXI_CLUSTER_LM,NPIV)
        IF (COMPRESS_CB) THEN
         IF (NPIV1.EQ.0) THEN
          CALL GET_CUT(IW(IOLDPS+HS+NROW1:IOLDPS+HS+NROW1+NCOL1-1), 
     &                    NASS1,
     &                    NCOL1-NASS1, LRGROUPS, NPARTSCB_COL, 
     &                    NPARTSASS_COL, BEGS_BLR_COL)
          CALL REGROUPING2(BEGS_BLR_COL, NPARTSASS_COL, NASS1, 
     &                     NPARTSCB_COL,
     &                     NCOL1-NASS1, KEEP(488), .FALSE., KEEP(472),
     &                     NCOL1, KEEP(1))
          NB_BLR_COL = NPARTSCB_COL + NPARTSASS_COL
          IF (NPARTSASS_MASTER.NE.NPARTSASS_COL) THEN
             ALLOCATE(BEGS_BLR_COL_TMP(
     &         size(BEGS_BLR_COL)-NPARTSASS_COL+NPARTSASS_MASTER),
     &         stat=allocok)
             IF ( allocok .GT. 0 ) THEN
                IF (LP > 0 ) WRITE(LP,*) MYID,
     &            ": ALLOCATION FAILURE FOR BEGS_BLR_COL_TMP in",
     &            "SMUMPS_PROCESS_SYM_BLOCFACTO"
                IFLAG = -13
                IERROR = size(BEGS_BLR_COL)
     &                   -NPARTSASS_COL+NPARTSASS_MASTER
                GOTO 700
             END IF
             IF ( size(BEGS_BLR_COL).GT. NPARTSASS_COL) THEN
              DO II=1, size(BEGS_BLR_COL) - NPARTSASS_COL
                BEGS_BLR_COL_TMP (II+NPARTSASS_MASTER) =
     &                     BEGS_BLR_COL(II+NPARTSASS_COL)
              ENDDO
             ENDIF
             DO II= 1, NPARTSASS_MASTER
                BEGS_BLR_COL_TMP (II) = 
     &                 BEGS_BLR_COL(max(NPARTSASS_COL,1)+1)
             ENDDO
             DEALLOCATE(BEGS_BLR_COL)
             BEGS_BLR_COL => BEGS_BLR_COL_TMP
             NPARTSASS_COL = NPARTSASS_MASTER
             NB_BLR_COL = NPARTSCB_COL + NPARTSASS_COL
          ENDIF
         ELSE
            CALL SMUMPS_BLR_RETRIEVE_BEGS_BLR_C (IW(IOLDPS+XXF), 
     &                  BEGS_BLR_COL, NPARTSASS_COL )
            KEEP_BEGS_BLR_COL = .TRUE.  
            NB_BLR_COL   = size(BEGS_BLR_COL) - 1
            NPARTSCB_COL = NB_BLR_COL - NPARTSASS_COL 
         ENDIF
        ELSE
         NULLIFY(BEGS_BLR_COL)
        ENDIF
        IF (NPIV1.EQ.0)  THEN
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          NB_ACCESSES_INIT=0
            IF (NSLAVES_PREC.GT.0) THEN
              NB_ACCESSES_INIT=NSLAVES_PREC+1
            ENDIF
          IF ( (KEEP(486).EQ.2) 
     &       ) THEN
            NB_ACCESSES_INIT = huge(NPARTSASS_MASTER)
          END IF
          INFO_TMP(1) = IFLAG
          INFO_TMP(2) = IERROR
          IF (IFLAG.LT.0) GOTO 700
          CALL SMUMPS_BLR_SAVE_INIT(IW(IOLDPS+XXF), 
     &              .TRUE., .TRUE., .TRUE., NPARTSASS_COL, 
     &              BEGS_BLR_LS, BEGS_BLR_COL, NB_ACCESSES_INIT, 
     &              INFO_TMP)
          IFLAG  = INFO_TMP(1) 
          IERROR = INFO_TMP(2) 
          IF (IFLAG.LT.0) GOTO 700
        ENDIF
        LWORK = MAXI_CLUSTER*MAXI_CLUSTER
        OMP_NUM = 1
#if ! defined(BLR_NOOPENMP)
!$      OMP_NUM = OMP_GET_MAX_THREADS()
#endif
        ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &            RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &            TAU(MAXI_CLUSTER*OMP_NUM),
     &            JPVT(MAXI_CLUSTER*OMP_NUM), 
     &            WORK(LWORK*OMP_NUM),
     &            stat=allocok)
        IF (allocok > 0 ) THEN
           IFLAG  = -13
           IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
           GOTO 700
        ENDIF
          CURRENT_BLR = 1
          ALLOCATE(BLR_LS(NB_BLR_LS), stat=allocok)
          IF (allocok > 0 ) THEN
             IFLAG  = -13
             IERROR = NB_BLR_LS
             GOTO 700
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
          CALL SMUMPS_COMPRESS_PANEL_I_NOOPT
     &        (A_PTR(POSELT), LA_PTR, 1_8,
     &        IFLAG, IERROR, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        DKEEP(8), KEEP(466), 0, 
     &        KEEP(473), BLR_LS(1), 
     &        CURRENT_BLR, 'V', WORK, TAU, JPVT, LWORK, RWORK,
     &        BLOCKLR, MAXI_CLUSTER, NELIM, 
     &        .TRUE.,  
     &        NPIV, NPIV1,
     &        2, KEEP(483), KEEP8, 
     &        OMP_NUM)
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif
          IF (IFLAG.LT.0) GOTO 300
          IF (KEEP(475).GE.1) THEN
            IF (DYNAMIC_ALLOC) THEN
              CALL SMUMPS_BLR_PANEL_LRTRSM(
     &              DYN_BLOCFACTO, LA_BLOCFACTO, 1_8,
     &              LD_BLOCFACTO, -6666, 
     &              NB_BLR_LS+1, 
     &              BLR_LS, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_LS+1, 
     &              2, 1, 0, 
     &              .TRUE., 
     &              DYN_PIVINFO, OFFSET_IW=1)
            ELSE
            CALL SMUMPS_BLR_PANEL_LRTRSM(A, LA, POSBLOCFACTO, 
     &              LD_BLOCFACTO, -6666, 
     &              NB_BLR_LS+1, 
     &              BLR_LS, CURRENT_BLR, CURRENT_BLR+1, NB_BLR_LS+1, 
     &              2, 1, 0, 
     &              .TRUE., 
     &              IW, OFFSET_IW=IPIV)
            ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP BARRIER
#endif          
            IF (KEEP(486).NE.2) THEN
              CALL SMUMPS_DECOMPRESS_PANEL_I_NOOPT(
     &        A_PTR(POSELT), LA_PTR, 1_8,
     &        NCOL1, NCOL1,
     &        .TRUE.,    
     &        NPIV1+1,   
     &        1,         
     &        NB_BLR_LS+1, BLR_LS(1), CURRENT_BLR, 'V', 1)
            ENDIF
          ENDIF
 300      CONTINUE         
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
        ENDIF
      ENDIF
      IF ( (KEEP(201).eq.1) .AND. 
     &    (OOCWRITE_COMPATIBLE_WITH_BLR .OR. NPIV.EQ.0) ) THEN
        MonBloc%INODE = INODE
        MonBloc%MASTER = .FALSE.
        MonBloc%Typenode = 2
        MonBloc%NROW = NROW1  
        MonBloc%NCOL = NCOL1  
        MonBloc%NFS  = NASS1
        MonBloc%LastPiv = NPIV1 + NPIV 
        MonBloc%LastPanelWritten_L = -9999 
        MonBloc%LastPanelWritten_U = -9999 
        NULLIFY(MonBloc%INDICES)
        MonBloc%Last = LASTBL
        STRAT = STRAT_TRY_WRITE 
        NextPivDummy      = -8888 
        LIWFAC = IW(IOLDPS+XXI)
        LAST_CALL = .FALSE.
        CALL SMUMPS_OOC_IO_LU_PANEL_I( STRAT, TYPEF_L,
     &       A_PTR(POSELT),
     &       LA_PTR, MonBloc, NextPivDummy, NextPivDummy,
     &       IW(IOLDPS), LIWFAC, MYID, KEEP8(31), IFLAG_OOC, LAST_CALL)
        IF ( IFLAG_OOC .LT. 0 )THEN
          IFLAG = IFLAG_OOC
          GOTO 700
        ENDIF
      ENDIF
      IF (NPIV.GT.0) THEN
        IF (LR_ACTIVATED) THEN
          IF (NELIM.GT.0) THEN
            LPOS2 = POSELT + int(NPIV1,8)
            UPOS = 1_8+int(NPIV,8)
            LPOS  = LPOS2 + int(NPIV,8)
            IF (DYNAMIC_ALLOC) THEN
              CALL SMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        DYN_BLOCFACTO, LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        CURRENT_BLR, BLR_LS(1), NB_BLR_LS+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ELSE
              CALL SMUMPS_BLR_UPD_NELIM_VAR_L_I(
     &        A(POSBLOCFACTO), LA_BLOCFACTO, UPOS,
     &        A_PTR(POSELT), LA_PTR, LPOS-POSELT+1_8,
     &        IFLAG, IERROR, LD_BLOCFACTO, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        CURRENT_BLR, BLR_LS(1), NB_BLR_LS+1, 
     &        CURRENT_BLR+1, NELIM, 'N')
            ENDIF
          ENDIF
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif          
          IF (DYNAMIC_ALLOC) THEN
            CALL SMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1, NROW1,
     &        DYN_BLOCFACTO, LA_BLOCFACTO,
     &        LD_BLOCFACTO, 
     &        BEGS_BLR_LM(1), size(BEGS_BLR_LM), NB_BLR_LM+1,
     &        BLR_LM(1), NPIV1, 
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        BLR_LS(1), 0, 
     &        CURRENT_BLR, CURRENT_BLR,   
     &        DYN_PIVINFO, 
     &        BLOCKLR,
     &        MAXI_CLUSTER, OMP_NUM,
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ELSE
            CALL SMUMPS_BLR_SLV_UPD_TRAIL_LDLT_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, 
     &        IFLAG, IERROR, NCOL1, NROW1,
     &        A(POSBLOCFACTO), LA_BLOCFACTO,
     &        LD_BLOCFACTO, 
     &        BEGS_BLR_LM(1), size(BEGS_BLR_LM), NB_BLR_LM+1,
     &        BLR_LM(1), NPIV1, 
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS), NB_BLR_LS+1,
     &        BLR_LS(1), 0, 
     &        CURRENT_BLR, CURRENT_BLR,   
     &        IW(IPIV), 
     &        BLOCKLR,
     &        MAXI_CLUSTER, OMP_NUM,
     &        KEEP(481), DKEEP(11), KEEP(466), KEEP(477) 
     &        )
          ENDIF
          IF (IFLAG.LT.0) GOTO 400
 400      CONTINUE          
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif          
          IF (IFLAG.LT.0) GOTO 700
          CALL UPD_MRY_LU_LRGAIN(BLR_LS, NPARTSCB
     &               )
          CALL DEALLOC_BLR_PANEL(BLR_LM, NB_BLR_LM, KEEP8, KEEP(34))
          DEALLOCATE(BLR_LM)
          IF (NSLAVES_PREC.GT.0
     &       .OR.
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &     ) THEN
            CALL SMUMPS_BLR_SAVE_PANEL_LORU(
     &          IW(IOLDPS+XXF),
     &          0,   
     &          IPANEL, BLR_LS)
             KEEP_BLR_LS = .TRUE.
          ENDIF
        ELSE 
          IF (NPIV .GT. 0 .AND. NCOL-NPIV.GT.0)THEN
          LPOS2 = POSELT + int(NPIV1,8)
          LPOS  = LPOS2 + int(NPIV,8)
            IF (DYNAMIC_ALLOC) THEN
              UPOS = int(NPIV+1,8)
              CALL sgemm('N','N', NCOL-NPIV, NROW1, NPIV,
     &             ALPHA, DYN_BLOCFACTO(UPOS), NCOL,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ELSE
              UPOS = POSBLOCFACTO+int(NPIV,8)
              CALL sgemm('N','N', NCOL-NPIV, NROW1, NPIV,
     &             ALPHA,A(UPOS), NCOL,
     &             A_PTR(LPOS2), NCOL1, ONE, A_PTR(LPOS), NCOL1)
            ENDIF
          ENDIF
          DPOS = POSELT + int(NCOL1 - NROW1,8)
#if defined(GEMMT_AVAILABLE)
            IF ( KEEP(421).EQ. -1 .OR.
     &        ( KEEP(421) .GT. 0 .AND. NROW1 .GT. KEEP(421) ) ) THEN
              LPOS2 = POSELT + int(NPIV1,8)
              UPOS  = 1_8
              CALL sgemmt( 'U', 'T', 'N', NROW1, NPIV, ALPHA,
     &         UIP21K( UPOS ), NPIV,
     &         A_PTR( LPOS2 ), NCOL1, ONE,
     &         A_PTR( DPOS ), NCOL1 )
            ELSE
#endif
              IF ( NROW1 .GT. KEEP(7) ) THEN
                BLSIZE = KEEP(8)
              ELSE
                BLSIZE = NROW1
              ENDIF
              IF ( NROW1 .GT. 0 ) THEN
                DO IROW = 1, NROW1, BLSIZE
                  Block = min( BLSIZE, NROW1 - IROW + 1 )
                  DPOS  = POSELT + int(NCOL1 - NROW1,8)
     &                + int( IROW - 1, 8 ) * int( NCOL1 + 1, 8 )
                  LPOS2 = POSELT + int(NPIV1,8)
     &                + int( IROW - 1, 8 ) * int( NCOL1, 8 )
                  UPOS  = int( IROW - 1, 8 ) * int(NPIV, 8) + 1_8
                  DO I = 1, Block
                  CALL sgemv( 'T', NPIV, Block-I+1, ALPHA,
     &             A_PTR( LPOS2 + int(I - 1,8) * int(NCOL1,8) ), NCOL1,
     &             UIP21K( UPOS + int(NPIV,8) * int( I - 1, 8 ) ),
     &             1, ONE, A_PTR(DPOS+int(NCOL1+1,8)*int(I-1,8)),NCOL1 )
                  END DO
                  IF ( NROW1-IROW+1-Block .ne. 0 )
     &            CALL sgemm( 'T', 'N', Block, NROW1-IROW+1-Block,
     &            NPIV, ALPHA,
     &            UIP21K( UPOS ), NPIV,
     &            A_PTR( LPOS2 + int(Block,8) * int(NCOL1,8) ), NCOL1,
     &            ONE,
     &            A_PTR( DPOS + int(Block,8) * int(NCOL1,8) ), NCOL1 )
                ENDDO
              ENDIF
#if defined(GEMMT_AVAILABLE)
            ENDIF 
#endif
        ENDIF 
        FLOP1 = dble(NROW1) * dble(NPIV) *
     &           dble( 2 * NCOL  - NPIV + NROW1 +1 )
        FLOP1 = -FLOP1
        CALL SMUMPS_LOAD_UPDATE( 1, .FALSE., FLOP1, KEEP,KEEP8 )
      ENDIF 
      IW(IOLDPS+KEEP(IXSZ)) = IW(IOLDPS+KEEP(IXSZ)) - NPIV
      IW(IOLDPS+3+KEEP(IXSZ)) = IW(IOLDPS+3+KEEP(IXSZ)) + NPIV
      IF (LASTBL) IW(IOLDPS+1+KEEP(IXSZ)) = IW(IOLDPS + 3+KEEP(IXSZ))
      IF ( .NOT. LR_ACTIVATED ) THEN
      IF (DYNAMIC_ALLOC) THEN
        IF (allocated(DYN_PIVINFO)  ) DEALLOCATE(DYN_PIVINFO)
        IF (allocated(DYN_BLOCFACTO)) THEN
           DEALLOCATE(DYN_BLOCFACTO)
        ENDIF
      ELSE
      LRLU  = LRLU + LA_BLOCFACTO
      LRLUS = LRLUS + LA_BLOCFACTO
      KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
      POSFAC = POSFAC - LA_BLOCFACTO
      IWPOS = IWPOS - NPIV
      CALL SMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &                LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
      ENDIF 
      ENDIF 
      IF (LR_ACTIVATED) THEN
        IF (allocated(RWORK))  DEALLOCATE(RWORK)
        IF (allocated(WORK)) DEALLOCATE(WORK)
        IF (allocated(TAU)) DEALLOCATE(TAU)
        IF (allocated(JPVT)) DEALLOCATE(JPVT)
        IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
      ENDIF
      IF ( NSLAVES_FOLLOW .NE. 0 .and. NPIV .NE. 0 ) THEN
         IPOSK = NPIV1 + 1
         JPOSK = NCOL1 - NROW1 + 1
           NPIVSENT = NPIV
           IERR = -1
           DO WHILE ( IERR .eq. -1 )
            IF (DYNAMIC_ALLOC) THEN
              CALL SMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, KEEP,
     &             LR_ACTIVATED, BLR_LS, IPANEL, 
     &             DYN_BLOCFACTO, LA_BLOCFACTO,
     &             1_8, LD_BLOCFACTO,
     &             DYN_PIVINFO, MAXI_CLUSTER,
     &                    IERR )
            ELSE
              CALL SMUMPS_BUF_SEND_BLFAC_SLAVE(
     &                    INODE, NPIVSENT, FPERE,
     &                    IPOSK, JPOSK,
     &                    UIP21K, NROW1,
     &                    NSLAVES_FOLLOW,
     &                    LIST_SLAVES_FOLLOW(1),
     &                    COMM, KEEP,
     &             LR_ACTIVATED, BLR_LS, IPANEL, 
     &             A, LA, 
     &             POSBLOCFACTO, LD_BLOCFACTO,
     &             IW(IPIV), MAXI_CLUSTER,
     &                    IERR )
            ENDIF
            IF (IERR .EQ. -1 ) THEN
              IOLDPS = PTRIST(STEP(INODE))
              IF ( IW(IOLDPS+6+KEEP(IXSZ)) .eq.
     &              huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
                    COUNTER_WAS_HUGE=.TRUE.
                    IW(IOLDPS+6+KEEP(IXSZ)) = 1
              ELSE
                    COUNTER_WAS_HUGE=.FALSE.
              ENDIF
              TO_UPDATE_CPT_RECUR =
     &                      ( NSLAVES_TOT - NSLAVES_FOLLOW - 1 ) *
     &                       (2*NASS1/KEEP(6))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) - TO_UPDATE_CPT_RECUR - 10
              BLOCKING = .FALSE.
              SET_IRECV= .TRUE.
              MESSAGE_RECEIVED = .FALSE.
              CALL SMUMPS_TRY_RECVTREAT( COMM_LOAD, ASS_IRECV,
     &         BLOCKING, SET_IRECV, MESSAGE_RECEIVED,
     &         MPI_ANY_SOURCE, MPI_ANY_TAG,
     &         STATUS, 
     &         BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &         IWPOS, IWPOSCB, IPTRLU,
     &         LRLU, LRLUS, N, IW, LIW, A, LA, PTRIST,
     &         PTLUST_S, PTRFAC,
     &         PTRAST, STEP, PIMASTER, PAMASTER, NSTK_S, COMP,
     &         IFLAG, IERROR, COMM,
     &         PERM, IPOOL, LPOOL, LEAF, NBFIN, MYID, SLAVEF,
     &         root, OPASSW, OPELIW, ITLOC, RHS_MUMPS,
     &         FILS, DAD, PTRARW, PTRAIW,
     &         PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &         INTARR, DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND, FRERE_STEPS,
     &         LPTRAR, NELT, FRTPTR, FRTELT, 
     &         ISTEP_TO_INIV2, TAB_POS_IN_PERE, .TRUE. 
     &               , LRGROUPS
     &           )
              IOLDPS = PTRIST(STEP(INODE))
              IW(IOLDPS+6+KEEP(IXSZ)) =
     &             IW(IOLDPS+6+KEEP(IXSZ)) + TO_UPDATE_CPT_RECUR + 10
              IF ( COUNTER_WAS_HUGE .AND.
     &             IW(IOLDPS+6+KEEP(IXSZ)).EQ.1 ) THEN
                IW(IOLDPS+6+KEEP(IXSZ)) = huge(IW(IOLDPS+6+KEEP(IXSZ)))
              ENDIF
              IF ( IFLAG .LT. 0 ) GOTO 550
            END IF
           END DO
           IF ( IERR .eq. -2 ) THEN
             IF (LP > 0 ) THEN 
               WRITE(LP,*) MYID,
     &": FAILURE, SEND BUFFER TOO SMALL DURING
     & SMUMPS_PROCESS_SYM_BLOCFACTO"
             ENDIF
             IFLAG = -17
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           IF ( IERR .eq. -3 ) THEN
              IF (LP > 0 ) WRITE(LP,*) MYID,
     &": FAILURE, RECV BUFFER TOO SMALL DURING
     & SMUMPS_PROCESS_SYM_BLOCFACTO"
             IFLAG = -20
             IERROR = 5 * KEEP(34) + NPIV * NROW1 * KEEP(35)
             GOTO 700
           END IF
           DEALLOCATE(LIST_SLAVES_FOLLOW)
      END IF
      IF ( LR_ACTIVATED ) THEN
        IF (NPIV.GT.0 .AND. NSLAVES_PREC.GT.0
     &     .AND. KEEP(486).EQ.3
     &    ) THEN
          IOLDPS = PTRIST(STEP(INODE))
          CALL SMUMPS_BLR_DEC_AND_TRYFREE_L(IW(IOLDPS+XXF),IPANEL,
     &                       KEEP8, KEEP(34))
        ENDIF 
        IF (DYNAMIC_ALLOC) THEN
          IF (allocated(DYN_PIVINFO))   DEALLOCATE(DYN_PIVINFO)
          IF (allocated(DYN_BLOCFACTO)) THEN
            DEALLOCATE(DYN_BLOCFACTO)
          ENDIF
        ELSE IF (NPIV .GT. 0) THEN
          LRLU  = LRLU + LA_BLOCFACTO
          LRLUS = LRLUS + LA_BLOCFACTO
          KEEP8(69) = KEEP8(69) - LA_BLOCFACTO
          POSFAC = POSFAC - LA_BLOCFACTO
          IWPOS = IWPOS - NPIV
      CALL SMUMPS_LOAD_MEM_UPDATE(.FALSE.,.FALSE.,
     &             LA-LRLUS,0_8,-LA_BLOCFACTO,KEEP,KEEP8,LRLUS)
        ENDIF 
      ENDIF 
      IF ( NPIV .NE. 0 )  THEN
        IF (allocated(UIP21K)) THEN
          DEALLOCATE( UIP21K )
        ENDIF
      ENDIF
      IOLDPS = PTRIST(STEP(INODE))
      CALL SMUMPS_DM_SET_DYNPTR( IW(IOLDPS+XXS), A, LA,
     &     PTRAST(STEP(INODE)), IW(IOLDPS+XXD), IW(IOLDPS+XXR),
     &     A_PTR, POSELT, LA_PTR )
      IF (LASTBL) THEN
         IF ( KEEP(486) .NE. 0) THEN
           IF (LR_ACTIVATED) THEN
             CALL STATS_COMPUTE_FLOP_SLAVE_TYPE2(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
           ELSE
             CALL UPD_FLOP_FRFRONT_SLAVE(NROW1, NCOL1, NASS1,
     &             KEEP(50), INODE)
           ENDIF
         ENDIF
         IF ( IW(IOLDPS+6+KEEP(IXSZ)).EQ.
     &     huge(IW(IOLDPS+6+KEEP(IXSZ))) ) THEN
           IW(IOLDPS+6+KEEP(IXSZ)) =  1
         ENDIF
         IW(IOLDPS+6+KEEP(IXSZ)) = IW(IOLDPS+6+KEEP(IXSZ))
     &                           - TO_UPDATE_CPT_END 
     &                           - 1 
         IF ( IW(IOLDPS+6+KEEP(IXSZ) ) .eq. 0
     &       .and. KEEP(50) .ne. 0 .and. NSLAVES_FOLLOW .eq. 0
     &       .and. NSLAVES_TOT.NE.1 ) THEN
          DEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(INODE)),
     &                           KEEP(199) )
          CALL SMUMPS_BUF_SEND_1INT( INODE, DEST, END_NIV2_LDLT,
     &                              COMM, KEEP, IERR )
          IF ( IERR .LT. 0 ) THEN
            write(*,*) ' Internal error in PROCESS_SYM_BLOCFACTO.'
            IFLAG = -99
            GOTO 700
          END IF
        ENDIF
      END IF
        IF (IW(IOLDPS+6+KEEP(IXSZ)) .eq. 0 ) THEN 
          NELIM =  IW( IOLDPS + 4 + KEEP(IXSZ)) - 
     &              IW( IOLDPS + 3  + KEEP(IXSZ))
          IF (LR_ACTIVATED) THEN
            IF (COMPRESS_CB) THEN
              allocate(CB_LRB(NB_BLR_LS,NB_BLR_COL-NPARTSASS_COL),
     &                 stat=allocok)
              IF (allocok > 0) THEN
                IFLAG  = -13
                IERROR = NB_BLR_LS*(NB_BLR_COL-NPARTSASS_COL)
                GOTO 700
              ENDIF
              DO II=1,NB_BLR_LS
              DO JJ=1,NB_BLR_COL-NPARTSASS_COL
                CB_LRB(II,JJ)%M=0
                CB_LRB(II,JJ)%N=0
                NULLIFY(CB_LRB(II,JJ)%Q)
                NULLIFY(CB_LRB(II,JJ)%R)
                CB_LRB(II,JJ)%ISLR = .FALSE.
              ENDDO
              ENDDO
              CALL SMUMPS_BLR_SAVE_CB_LRB(IW(IOLDPS+XXF),CB_LRB)
            ENDIF
            IF (COMPRESS_CB) THEN
              NFS4FATHER = -9999
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2) ) THEN
               CALL SMUMPS_BLR_RETRIEVE_NFS4FATHER ( IW(IOLDPS+XXF),
     &             NFS4FATHER )
                NFS4FATHER = max(NFS4FATHER,0) + NELIM
              ENDIF
              ALLOCATE(M_ARRAY(max(1,NFS4FATHER)), stat=allocok)
              IF ( allocok .GT. 0 ) THEN
                IF (LP > 0 ) WRITE(LP,*) MYID,
     &          ": ALLOCATION FAILURE FOR M_ARRAY ",
     &          "SMUMPS_PROCESS_SYM_BLOCFACTO"
                IFLAG = -13
                IERROR = max(1,NFS4FATHER)
              ENDIF
              BEGS_BLR_COL(1+NPARTSASS_COL) = 
     &               BEGS_BLR_COL(1+NPARTSASS_COL) - NELIM
              CALL MAX_CLUSTER(BEGS_BLR_COL,NB_BLR_COL,MAXI_CLUSTER_COL)
              MAXI_CLUSTER=max(MAXI_CLUSTER_LS, MAXI_CLUSTER_COL)
              ALLOCATE(BLOCKLR(MAXI_CLUSTER, OMP_NUM*MAXI_CLUSTER),
     &                  RWORK(2*MAXI_CLUSTER*OMP_NUM), 
     &                  TAU(MAXI_CLUSTER*OMP_NUM),
     &                  JPVT(MAXI_CLUSTER*OMP_NUM), 
     &                  WORK(LWORK*OMP_NUM),
     &                  stat=allocok)
              IF (allocok > 0 ) THEN
                IFLAG  = -13
                IERROR = OMP_NUM*(LWORK + MAXI_CLUSTER*(MAXI_CLUSTER+4))
                GOTO 700
              ENDIF
              NBROWSinF    = 0
              NVSCHUR_K253 =  0
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
               CALL SMUMPS_COMPUTE_NBROWSinF (
     &                N, INODE, FPERE, KEEP, 
     &                IOLDPS, HS, 
     &                IW, LIW, 
     &                NROW1, NCOL1, NPIV+NPIV1,
     &                NELIM, NFS4FATHER,
     &                NBROWSinF
     &                )
               IF ((KEEP(114).EQ.1) .AND. (KEEP(116).GT.0) ) THEN
                  NSLAVES_L = IW(PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
                  IROW_L    = PTRIST(STEP(INODE)) + 6 + NSLAVES_L + 
     &                        KEEP(IXSZ)
                  CALL SMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &                 N, 
     &                 NROW1,
     &                 KEEP(116), 
     &                 IW(IROW_L),
     &                 PERM, NVSCHUR_K253 )
               ELSE IF (KEEP(253).NE.0) THEN
                  NSLAVES_L = IW(PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ))
                  IROW_L    = PTRIST(STEP(INODE)) + 6 + NSLAVES_L + 
     &                        KEEP(IXSZ)
                  CALL SMUMPS_GET_SIZE_SCHUR_IN_FRONT ( 
     &                 N, 
     &                 NROW1,
     &                 0,  
     &                 IW(IROW_L),
     &                 PERM, NVSCHUR_K253 )
               ENDIF
              ENDIF
              IF (IFLAG.LT.0) GOTO 700
#if ! defined(BLR_NOOPENMP)
!$OMP PARALLEL
#endif
              CALL SMUMPS_COMPRESS_CB_I(
     &        A_PTR(POSELT), LA_PTR, 1_8, NCOL1,
     &        BEGS_BLR_LS(1), size(BEGS_BLR_LS),
     &        BEGS_BLR_COL(1), size(BEGS_BLR_COL),
     &        NB_BLR_LS, NB_BLR_COL-NPARTSASS_COL,
     &        NPARTSASS_COL, 
     &        NROW1, NCOL1-NPIV1-NPIV, INODE,
     &        IW(IOLDPS+XXF), 1, 2, IFLAG, IERROR,
     &        DKEEP(12), KEEP(466), KEEP(484), KEEP(489),
     &        CB_LRB(1,1),
     &        WORK, TAU, JPVT, LWORK, RWORK, BLOCKLR,
     &        MAXI_CLUSTER, KEEP8, OMP_NUM,
     &        NFS4FATHER, NPIV1+NPIV, NVSCHUR_K253, KEEP(1), 
     &        M_ARRAY
     &        , NELIM, NBROWSinF
     &        )
#if ! defined(BLR_NOOPENMP)
!$OMP END PARALLEL
#endif
              IF (IFLAG.LT.0) GOTO 650
              IF ( (KEEP(219).NE.0).AND.(KEEP(50).EQ.2).AND.
     &             NFS4FATHER.GT.0  ) THEN
                 INFO_TMP(1) = IFLAG
                 INFO_TMP(2) = IERROR
                 CALL SMUMPS_BLR_SAVE_M_ARRAY( IW(IOLDPS+XXF),
     &            M_ARRAY, INFO_TMP)
                 IFLAG  = INFO_TMP(1) 
                 IERROR = INFO_TMP(2) 
              ENDIF
              DEALLOCATE(M_ARRAY)
 650          CONTINUE         
            ENDIF 
            IF (IFLAG.LT.0) GOTO 700
          ENDIF
          CALL SMUMPS_END_FACTO_SLAVE( COMM_LOAD, ASS_IRECV, 
     &    N, INODE, FPERE, 
     &    root,
     &    MYID, COMM,
     &    
     &    BUFR, LBUFR, LBUFR_BYTES, PROCNODE_STEPS, POSFAC,
     &    IWPOS, IWPOSCB, IPTRLU, LRLU, LRLUS, IW, LIW, A, LA,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER,
     &    PAMASTER,
     &    NSTK_S, COMP, IFLAG, IERROR, PERM,
     &    IPOOL, LPOOL, LEAF, NBFIN, SLAVEF,
     &    OPASSW, OPELIW, ITLOC, RHS_MUMPS, FILS, DAD, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    INTARR,DBLARR,ICNTL,KEEP,KEEP8,DKEEP,ND,FRERE_STEPS,
     &    LPTRAR, NELT, FRTPTR, FRTELT, 
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &               , LRGROUPS
     &      )
        ENDIF
        GOTO 550
 700    CONTINUE
        CALL SMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 550    CONTINUE
        IF (LR_ACTIVATED) THEN
          IF (allocated(RWORK))  DEALLOCATE(RWORK)
          IF (allocated(WORK)) DEALLOCATE(WORK)
          IF (allocated(TAU)) DEALLOCATE(TAU)
          IF (allocated(JPVT)) DEALLOCATE(JPVT)
          IF (allocated(BLOCKLR)) DEALLOCATE(BLOCKLR)
          IF (NPIV.GT.0) THEN
            IF (.NOT.KEEP_BEGS_BLR_LS) THEN
              IF (associated(BEGS_BLR_LS)) DEALLOCATE(BEGS_BLR_LS)
            ENDIF
            IF (.NOT.KEEP_BLR_LS) THEN
             IF (associated(BLR_LS)) THEN
              CALL DEALLOC_BLR_PANEL(BLR_LS, NB_BLR_LS, KEEP8, KEEP(34))
              DEALLOCATE(BLR_LS)
             ENDIF
            ENDIF
            IF (associated(BEGS_BLR_LM)) DEALLOCATE(BEGS_BLR_LM)
            IF (.NOT.KEEP_BEGS_BLR_COL) THEN
              IF (COMPRESS_CB) THEN
                IF (associated(BEGS_BLR_COL)) THEN 
                  DEALLOCATE( BEGS_BLR_COL)
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      RETURN
      END SUBROUTINE SMUMPS_PROCESS_SYM_BLOCFACTO
