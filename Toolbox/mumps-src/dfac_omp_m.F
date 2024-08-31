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
      MODULE DMUMPS_FAC_OMP_M
      INTEGER(8), PARAMETER :: UnderL0 = -20_8
      INTEGER(8), PARAMETER :: CopyNotStarted = -19_8
      INTEGER(8), PARAMETER :: WaitMem = -18_8
      INTEGER(8), PARAMETER :: CopyFactorsFinished = -17_8
      INTEGER(8), PARAMETER :: AllocateViderCBEnCours = -16_8
      INTEGER(8), PARAMETER :: Finished = -15_8
      CONTAINS
      SUBROUTINE DMUMPS_FAC_L0_OMP(N,LIW, NSTK_STEPS, ND,
     & FILS,STEP, FRERE, DAD, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     & PTRIST, PTRAST, PIMASTER, PAMASTER, PTRARW, PTRAIW, PTR8ARR,
     & NINCOLARR, NINROWARR, PTRDEBARR, ITLOC,
     & RHS_MUMPS, RINFO, NBROOT, NBRTOT, NBROOT_UNDER_L0, UU, ICNTL,
     & PTLUST_S, PTRFAC, INFO, KEEP, KEEP8, PROCNODE_STEPS,SLAVEF,
     & COMM_NODES, MYID, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES,
     & INTARR, DBLARR, root, PERM, NELT, FRTPTR, FRTELT, LPTRAR,
     & COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     & MEM_DISTRIB, NE, DKEEP, PIVNUL_LIST_STRUCT,
     &  LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &  L_VIRT_L0_OMP, VIRT_L0_OMP, VIRT_L0_OMP_MAPPING,
     &  L_PHYS_L0_OMP, PHYS_L0_OMP, PERM_L0_OMP, PTR_LEAFS_L0_OMP,
     &  L0_OMP_MAPPING, LL0_OMP_MAPPING, THREAD_LA,
     &  MUMPS_TPS_ARR, DMUMPS_TPS_ARR,
     &  NSTEPSW, OPASSW, OPELIW, NELVAW, COMP,
     &  MAXFRW, NMAXNPIVW, NPVW, NOFFNEGW, NULLNEGW,
     &  NB22T1W, NBTINYW, DET_EXPW, DET_MANTW, DET_SIGNW,
     &  LRGROUPS, L0_OMP_FACTORS, LL0_OMP_FACTORS,
     &  I4_L0_OMP, NBSTATS_I4, NBCOLS_I4,
     &  I8_L0_OMP, NBSTATS_I8, NBCOLS_I8 )
      USE DMUMPS_LOAD
!$    USE OMP_LIB
      USE MUMPS_TPS_M, ONLY : MUMPS_TPS_T
      USE DMUMPS_TPS_M, ONLY : DMUMPS_TPS_T
      USE DMUMPS_LR_STATS
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC,
     &                             DMUMPS_L0OMPFAC_T
      USE DMUMPS_DYNAMIC_MEMORY_M, ONLY :
     &                             DMUMPS_DM_FAC_ALLOC_ALLOWED
      USE MUMPS_PIVNUL_MOD, ONLY : PIVNUL_LIST_STRUCT_T
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER N,LIW, LPTRAR, 
     &        NSTEPSW, INFO(80)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER, DIMENSION(0: SLAVEF - 1) :: MEM_DISTRIB
      INTEGER KEEP(500), ICNTL(60)
      INTEGER(8) KEEP8(150)
      INTEGER(8), INTENT(IN) :: THREAD_LA
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER ITLOC(N+KEEP(253))
      DOUBLE PRECISION :: RHS_MUMPS(KEEP8(85))
      INTEGER NSTK_STEPS(KEEP(28))
      INTEGER(8), INTENT(IN) :: PTRARW(LPTRAR), PTRAIW(LPTRAR)
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER ND(KEEP(28))
      INTEGER FILS(N),PTRIST(KEEP(28))
      INTEGER STEP(N), FRERE(KEEP(28)), DAD(KEEP(28))
      INTEGER PIMASTER(KEEP(28))
      INTEGER PTLUST_S(KEEP(28)), PERM(N)
      INTEGER   ISTEP_TO_INIV2(KEEP(71)),
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NE(KEEP(28))
      DOUBLE PRECISION RINFO(40)
      INTEGER(8) :: PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NBROOT 
      INTEGER NBRTOT 
      INTEGER, intent(out) :: NBROOT_UNDER_L0 
      INTEGER COMM_LOAD, ASS_IRECV
      DOUBLE PRECISION UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER NELT
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      DOUBLE PRECISION DBLARR( KEEP8(26) )
      INTEGER INTARR( KEEP8(27) )
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION DKEEP(230)
      DOUBLE PRECISION :: OPASSW, OPELIW
      INTEGER LRGROUPS(KEEP(280))
      INTEGER, INTENT ( IN )  :: LPOOL_B_L0_OMP
      INTEGER, INTENT ( IN )  :: IPOOL_B_L0_OMP
     &                           ( LPOOL_B_L0_OMP )
      INTEGER, INTENT ( IN )  :: L_PHYS_L0_OMP
      INTEGER, INTENT ( IN )  :: PHYS_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT ( IN )  :: L_VIRT_L0_OMP
      INTEGER, INTENT ( IN )  :: VIRT_L0_OMP( L_VIRT_L0_OMP )
      INTEGER, INTENT ( IN )  :: VIRT_L0_OMP_MAPPING( L_VIRT_L0_OMP )
      INTEGER, INTENT ( IN )  :: PERM_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT ( IN )  :: PTR_LEAFS_L0_OMP( L_PHYS_L0_OMP + 1)
      INTEGER, INTENT ( IN )  :: LL0_OMP_MAPPING
      INTEGER, INTENT ( OUT ) :: L0_OMP_MAPPING( LL0_OMP_MAPPING )
      TYPE (MUMPS_TPS_T), DIMENSION(:) :: MUMPS_TPS_ARR
      TYPE (DMUMPS_TPS_T), DIMENSION(:) :: DMUMPS_TPS_ARR
      INTEGER, INTENT ( IN )  :: LL0_OMP_FACTORS
      TYPE (DMUMPS_L0OMPFAC_T), INTENT(INOUT) ::
     &                        L0_OMP_FACTORS(LL0_OMP_FACTORS)
      INTEGER, INTENT (IN)    :: NBSTATS_I4, NBSTATS_I8
      INTEGER, INTENT (IN)    :: NBCOLS_I4, NBCOLS_I8
      INTEGER, INTENT (IN)    :: I4_L0_OMP(NBSTATS_I4, NBCOLS_I4)
      INTEGER(8), INTENT (IN) :: I8_L0_OMP(NBSTATS_I8, NBCOLS_I8)
      LOGICAL DMUMPS_POOL_EMPTY
      EXTERNAL DMUMPS_POOL_EMPTY, DMUMPS_EXTRACT_POOL
      INTEGER :: MYTHREAD_ID, ITH  
      INTEGER :: THREAD_ID_P       
      DOUBLE PRECISION, PARAMETER :: DZERO = 0.0D0, DONE = 1.0D0
      INTEGER INODE, LEAF
      INTEGER TYPEF
      INTEGER NBFIN
      INTEGER TYPE
      INTEGER NBROOT_PROCESSED
      INTEGER MAXFRW, NPVW, NMAXNPIVW, NOFFNEGW, NULLNEGW, NELVAW, COMP
      INTEGER :: NB22T1W, NBTINYW, DET_EXPW, DET_SIGNW
      DOUBLE PRECISION :: DET_MANTW
      DOUBLE PRECISION FLOP_ESTIM_ACC
      INTEGER :: LPOOL_P
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPOOL_P
      INTEGER(8) :: TO_ALLOCATE
      INTEGER, DIMENSION(:), ALLOCATABLE    :: ID
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: VAL
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: STATE, SIZE_COPIED
      INTEGER ::  NbWaitMem, NbFinished, NbOnGoingCopies, NbUnderL0
      INTEGER(8) :: KEEP8_77_SAVE
      DOUBLE PRECISION :: GTIME
      INTEGER(8) :: MEMDISPO_UNDERL0, MEMDISPO_PERTHREAD
      INTEGER    :: BLR_STRAT
      INTEGER :: VIRTUAL_TASK, PHYSICAL_TASK
      INTEGER :: IFATH
      INTEGER :: I, INFO_P(2), allocok
      INTEGER(8) :: I8
!$    INTEGER(8) :: CHUNK8 
!$    LOGICAL    :: OMP_FLAG
!$    INTEGER    :: NOMP_SAVE
      CALL MUMPS_LOAD_DISABLE()
      GTIME = MPI_WTIME()
      L0_OMP_MAPPING = 0
      NBROOT_PROCESSED  = 0
      NSTEPSW           = 0
      OPASSW         = DZERO
      OPELIW         = DZERO
      NELVAW            = 0
      COMP              = 0
      MAXFRW            = 0 
      NMAXNPIVW         = 0 
      NOFFNEGW          = 0 
      NULLNEGW          = 0 
      FLOP_ESTIM_ACC = DZERO
      NPVW              = 0 
      NB22T1W           = 0 
      NBTINYW           = 0 
      DET_EXPW          = 0 
      DET_MANTW         = cmplx(1.0D0,0.0D0, kind=kind(1.0D0))
      DET_SIGNW         = 1
      DO ITH = 1, KEEP(400)
        NULLIFY(MUMPS_TPS_ARR(ITH)%IW)
        NULLIFY(MUMPS_TPS_ARR(ITH)%ITLOC)
        NULLIFY(DMUMPS_TPS_ARR(ITH)%A)
        CALL DMUMPS_SET_MAXS_MAXIS_THREAD(
     &     MUMPS_TPS_ARR(ITH)%LA,
     &     MUMPS_TPS_ARR(ITH)%LIW, BLR_STRAT,
     &     KEEP,
     &     I4_L0_OMP(1,ITH), NBSTATS_I4,
     &     I8_L0_OMP(1,ITH), NBSTATS_I8)
      ENDDO
      IF (KEEP8(4) .NE. 0_8) THEN
        CALL DMUMPS_MA_EFF_MEM_DISPO ( 
     &       MUMPS_TPS_ARR, KEEP(400),KEEP8, KEEP, 
     &       N, BLR_STRAT, LPOOL_B_L0_OMP, 
     &       I8_L0_OMP, NBSTATS_I8,
     &       MEMDISPO_UNDERL0   
     &      )
        IF (KEEP(486).EQ.2) THEN
         MEMDISPO_UNDERL0 = max(MEMDISPO_UNDERL0/20_8,0_8)
        ELSE
         MEMDISPO_UNDERL0 = max(MEMDISPO_UNDERL0/4_8,0_8)
        ENDIF
        KEEP8(77) = KEEP8(77) + MEMDISPO_UNDERL0
        MEMDISPO_PERTHREAD = 0_8
        IF (MEMDISPO_UNDERL0.GT.0) THEN
         MEMDISPO_PERTHREAD = MEMDISPO_UNDERL0/(int(KEEP(400),8))
        ENDIF
        DO ITH = 1, KEEP(400)
          MUMPS_TPS_ARR(ITH)%LA =  MUMPS_TPS_ARR(ITH)%LA + 
     &                             MEMDISPO_PERTHREAD
        ENDDO
      ENDIF
      DO ITH = 1, KEEP(400)
        MUMPS_TPS_ARR(ITH)%LRLU  = MUMPS_TPS_ARR(ITH)%LA
        MUMPS_TPS_ARR(ITH)%LRLUS = MUMPS_TPS_ARR(ITH)%LA
        MUMPS_TPS_ARR(ITH)%LRLUSM  = MUMPS_TPS_ARR(ITH)%LA
        MUMPS_TPS_ARR(ITH)%IPTRLU  = MUMPS_TPS_ARR(ITH)%LA
        MUMPS_TPS_ARR(ITH)%POSFAC  = 1_8
        MUMPS_TPS_ARR(ITH)%IWPOS   = 1
        MUMPS_TPS_ARR(ITH)%IWPOSCB = MUMPS_TPS_ARR(ITH)%LIW
      ENDDO
      IF (KEEP(406) .EQ. 2 ) THEN
        ALLOCATE(STATE(KEEP(400)), SIZE_COPIED(KEEP(400)), stat=allocok)
        IF (allocok .GT. 0 ) THEN
          WRITE(*,*) "Problem allocating STATE/SIZE_COPIED", KEEP(400)
          CALL MUMPS_ABORT()
        ENDIF
        CALL DMUMPS_PERFORM_COPIES_INIT( STATE, KEEP8_77_SAVE,
     &  NbWaitMem, NbFinished, NbOnGoingCopies, NbUnderL0, 
     &  KEEP, KEEP8 )
      ENDIF
!$OMP PARALLEL
!$OMP& SHARED ( IPOOL_B_L0_OMP, LPOOL_B_L0_OMP )
!$OMP& PRIVATE ( VIRTUAL_TASK, PHYSICAL_TASK,
!$OMP&     IPOOL_P, LPOOL_P, LEAF, INODE, IFATH, INFO_P, I, I8,
!$OMP&     TO_ALLOCATE, THREAD_ID_P,
!$OMP&     TYPE, TYPEF, NOMP_SAVE, allocok )
!$OMP& REDUCTION ( + : NPVW, OPASSW, OPELIW, NOFFNEGW, NELVAW, COMP,
!$OMP&                 NB22T1W, NBTINYW, DET_EXPW, NULLNEGW,
!$OMP&                 FLOP_ESTIM_ACC, NBROOT_PROCESSED, NSTEPSW )
!$OMP& REDUCTION ( * : DET_MANTW, DET_SIGNW )
!$OMP& REDUCTION ( max : MAXFRW, NMAXNPIVW )
      THREAD_ID_P = 1
!$    THREAD_ID_P = OMP_GET_THREAD_NUM () + 1
!$OMP BARRIER
!$    NOMP_SAVE = omp_get_max_threads()
!$    CALL omp_set_num_threads(1)
      LPOOL_P           = LPOOL_B_L0_OMP
      LEAF              = 1
      INFO_P            = 0
      VIRTUAL_TASK      = 0
!$    IF ( omp_get_num_threads() .NE. KEEP(400) ) THEN
!$      INFO_P(1)=-58
!$      INFO_P(2)=-100-omp_get_num_threads()
!$      GOTO 700
!$    ENDIF
      CALL DMUMPS_DM_FAC_ALLOC_ALLOWED( MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &                                  KEEP8,
     &                                  INFO_P(1), INFO_P(2) )
      IF (INFO_P(1) .LT. 0) GOTO 700
      ALLOCATE ( DMUMPS_TPS_ARR(THREAD_ID_P)%A(
     &           max(1_8,MUMPS_TPS_ARR(THREAD_ID_P)%LA) ),
     &           stat=allocok)
      IF (allocok.GT.0) THEN
          INFO_P(1) = -13
          CALL MUMPS_SETI8TOI4( MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &                          INFO_P(2))
          GOTO 700
      ELSE
          CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &       MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &       KEEP(405).EQ.1,
     &       KEEP8, INFO_P(1), INFO_P(2),
     &       .FALSE., .FALSE. ) 
          IF (INFO_P(1) .LT. 0) THEN
            GOTO 700
          ENDIF
      ENDIF
      TO_ALLOCATE = 
     & ((int(MUMPS_TPS_ARR(THREAD_ID_P)%LIW,8) * int(KEEP(34),8 )) /
     &   int(KEEP(35),8 ))+
     & ((int(LPOOL_P,8) * int(KEEP(34),8 )) / int(KEEP(35),8 ))+
     & ((int(N+KEEP(253),8) * int(KEEP(34),8 )) / int(KEEP(35),8 ))
      CALL DMUMPS_DM_FAC_ALLOC_ALLOWED( TO_ALLOCATE,
     &            KEEP8, INFO_P(1), INFO_P(2) )
      IF ( INFO_P(1) .LT. 0 ) GOTO 700
      ALLOCATE ( MUMPS_TPS_ARR(THREAD_ID_P)%IW(
     &                          MUMPS_TPS_ARR(THREAD_ID_P)%LIW ),
     &           IPOOL_P ( LPOOL_P ),
     &           MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC ( N + KEEP(253) ),
     &           stat = allocok )
      IF ( allocok .GT. 0 ) THEN
          INFO_P(1) = -13
          INFO_P(2) = MUMPS_TPS_ARR(THREAD_ID_P)%LIW +
     &                LPOOL_P + N+KEEP(253)
          GOTO 700
      ELSE
          CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS( TO_ALLOCATE,
     &       KEEP(405).EQ.1, KEEP8,
     &       INFO_P(1), INFO_P(2), .TRUE., .FALSE. )
          IF (INFO_P(1) .LT. 0) THEN
            GOTO 700
          ENDIF
      ENDIF
      CALL DMUMPS_ALLOC_CB( .FALSE., 0_8,
     &     .FALSE., .FALSE., MYID_NODES, N, KEEP, KEEP8, DKEEP,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%IW(1),
     &     MUMPS_TPS_ARR(THREAD_ID_P)%LIW,
     &     DMUMPS_TPS_ARR(THREAD_ID_P)%A(1),
     &     MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%LRLU,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%IPTRLU,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%IWPOS,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%IWPOSCB,
     &     SLAVEF, PROCNODE_STEPS, DAD,
     &     PTRIST, PTRAST, STEP, PIMASTER,
     &     PAMASTER, KEEP(IXSZ), 0_8, -444, -444, .true.,
     &     COMP, MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS,
     &     MUMPS_TPS_ARR(THREAD_ID_P)%LRLUSM,
     &     INFO_P(1), INFO_P(2)
     &     )
      CALL DMUMPS_INIT_POOL_LAST3( IPOOL_P(1), LPOOL_P,
     &           LEAF )
      MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC = 0
600   CONTINUE
      VIRTUAL_TASK = VIRTUAL_TASK + 1
      IF ( VIRTUAL_TASK .LT. L_VIRT_L0_OMP ) THEN
        IF ( VIRT_L0_OMP_MAPPING( VIRTUAL_TASK ) .EQ. THREAD_ID_P ) THEN
       DO PHYSICAL_TASK =
     &       VIRT_L0_OMP ( VIRTUAL_TASK ),
     &       VIRT_L0_OMP ( VIRTUAL_TASK + 1 ) - 1
        DO I = PTR_LEAFS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK )+1 ) + 1,
     &         PTR_LEAFS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK ) )
          IF (IPOOL_B_L0_OMP(I) .GT. 0) THEN 
            CALL DMUMPS_INSERT_POOL_N( N, IPOOL_P(1),
     &      LPOOL_P,
     &      PROCNODE_STEPS, SLAVEF, KEEP(199), KEEP(28), 3, 0, 1, STEP,
     &      IPOOL_B_L0_OMP(I) )
          END IF
        END DO
        DO WHILE (
     &    .NOT. DMUMPS_POOL_EMPTY( IPOOL_P(1), LPOOL_P )
     &    .AND. INFO_P(1) .GE. 0 ) 
          CALL DMUMPS_EXTRACT_POOL( N, IPOOL_P(1), LPOOL_P,
     &    PROCNODE_STEPS, SLAVEF, STEP, INODE, KEEP, KEEP8, MYID_NODES,
     &    ND, .FALSE. ) 
 10       CONTINUE
          L0_OMP_MAPPING ( STEP ( INODE ) ) = THREAD_ID_P
          IFATH = DAD ( STEP ( INODE ) )
          TYPE = 1
          IF ( IFATH .NE. 0 ) THEN
              TYPEF = 1
          ELSE
              TYPEF = -9999
          ENDIF
          CALL MUMPS_STOP_ON_USER_REQUEST(KEEP, KEEP8, ICNTL,
     &         INFO_P, MYID)
          IF (INFO_P(1) .LT. 0) THEN
            GOTO 700
          ENDIF
          IF (THREAD_ID_P .EQ. KEEP(400)-1) THEN
            CALL DMUMPS_UPDATE_PROGRESS( OPELIW*KEEP(400), KEEP8 )
          ENDIF
          CALL DMUMPS_PROCESS_FRONT_NIV1(COMM_LOAD, ASS_IRECV, N, INODE,
     &  TYPE, TYPEF, MUMPS_TPS_ARR(THREAD_ID_P)%LA, MUMPS_TPS_ARR(THREAD
     &  _ID_P)%IW(1), MUMPS_TPS_ARR(THREAD_ID_P)%LIW, DMUMPS_TPS_ARR(
     &  THREAD_ID_P)%A(1), MAXFRW, NOFFNEGW, NULLNEGW, NPVW,
     &  NB22T1W, NBTINYW, DET_EXPW, DET_MANTW, DET_SIGNW, INFO_P, UU,
     &  SEUIL, SEUIL_LDLT_NIV2, OPELIW, NELVAW, NMAXNPIVW, NSTEPSW,
     &  PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &  NE, MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC,
     &  MUMPS_TPS_ARR(THREAD_ID_P)%LRLU,
     &  MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS, MUMPS_TPS_ARR(THREAD_ID_P)%
     %  LRLUSM, MUMPS_TPS_ARR(THREAD_ID_P)%IPTRLU, ICNTL, KEEP, KEEP8,
     &  DKEEP, PIVNUL_LIST_STRUCT, COMP, MUMPS_TPS_ARR(THREAD_ID_P)%
     &  IWPOS, MUMPS_TPS_ARR(THREAD_ID_P)%IWPOSCB, PROCNODE_STEPS,
     &  SLAVEF, IFATH, COMM_NODES, MYID_NODES, IPOOL_P, LPOOL_P, LEAF,
     &  PERM, NSTK_STEPS, BUFR, LBUFR, LBUFR_BYTES,
     &  NBFIN, root, OPASSW, MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC(1),
     &  RHS_MUMPS, FILS, PTRARW, PTRAIW, PTR8ARR, NINCOLARR, NINROWARR,
     &  PTRDEBARR, INTARR, DBLARR, ND, FRERE, DAD,
     &  LPTRAR, NELT, FRTPTR, FRTELT, ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     &  LRGROUPS, FLOP_ESTIM_ACC )
          IF (INFO_P(1) .LT. 0) THEN
            GOTO 700
          ENDIF
          IF ( IFATH .NE. 0 ) THEN
              IF ( PHYS_L0_OMP ( PERM_L0_OMP ( PHYSICAL_TASK ) )
     &             .NE. INODE ) THEN
                  NSTK_STEPS ( STEP ( IFATH ) ) =
     &             NSTK_STEPS ( STEP ( IFATH ) ) - 1
                  IF ( NSTK_STEPS ( STEP ( IFATH ) ) .EQ. 0 ) THEN
                       INODE = IFATH
                       GOTO 10
                  ENDIF
              ELSE
!$OMP ATOMIC UPDATE
                  NSTK_STEPS ( STEP ( IFATH ) ) =
     &             NSTK_STEPS ( STEP ( IFATH ) ) - 1
!$OMP END ATOMIC
              END IF
          ELSE
            NBROOT_PROCESSED = NBROOT_PROCESSED + 1
          END IF
        END DO
       END DO
      ENDIF
      GOTO 600
      ENDIF
 700  CONTINUE
      IF (associated(MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC)) THEN
        DEALLOCATE(MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC)
        NULLIFY(MUMPS_TPS_ARR(THREAD_ID_P)%ITLOC)
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &      -(int(N+KEEP(253),8) * int(KEEP(34),8 )) / int(KEEP(35),8),
     &      KEEP(405).EQ.1, KEEP8,
     &      INFO_P(1), INFO_P(2), .TRUE., .FALSE. )
      ENDIF
      IF (allocated(IPOOL_P)) THEN
        DEALLOCATE(IPOOL_P);
        CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &      -(int(LPOOL_P,8) * int(KEEP(34),8 )) / int(KEEP(35),8),
     &      KEEP(405).EQ.1, KEEP8,
     &      INFO_P(1), INFO_P(2), .TRUE., .FALSE. )
      ENDIF
      IF ( KEEP(406) .EQ. 2) THEN
        CALL DMUMPS_PERFORM_COPIES( THREAD_ID_P,
     &  MUMPS_TPS_ARR, DMUMPS_TPS_ARR,
     &  L0_OMP_FACTORS, LL0_OMP_FACTORS,
     &  STATE, SIZE_COPIED,
     &  NbWaitMem, NbFinished, NbOngoingCopies, NbUnderL0,
     &  MYID_NODES, N, SLAVEF,
     &  STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &  KEEP, KEEP8, INFO_P
     &  )
      ELSE
       IF ((KEEP(407) .EQ. 1) .OR. (KEEP(406) .EQ.1) ) THEN
        IF (INFO_P(1) .GE. 0) THEN
          CALL DMUMPS_DM_CBSTATIC2DYNAMIC_I  
     &    (2,   
     &     0_8, 
     &    .FALSE., 
     &    MYID_NODES, N, SLAVEF, KEEP, KEEP8,
     &    MUMPS_TPS_ARR(THREAD_ID_P)%IW(1),
     &    MUMPS_TPS_ARR(THREAD_ID_P)%LIW,
     &    MUMPS_TPS_ARR(THREAD_ID_P)%IWPOSCB, 
     &    MUMPS_TPS_ARR(THREAD_ID_P)%IWPOS,
     &   DMUMPS_TPS_ARR(THREAD_ID_P)%A(1),
     &    MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &    MUMPS_TPS_ARR(THREAD_ID_P)%LRLU,
     &    MUMPS_TPS_ARR(THREAD_ID_P)%IPTRLU,
     &    MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS,
     &    STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &    INFO_P(1), INFO_P(2) )
        ENDIF
       ENDIF
       IF (KEEP(406) .EQ.1) THEN
        IF (INFO_P(1) .GE.0 )THEN 
          TO_ALLOCATE = max(MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1,1_8)
          CALL DMUMPS_DM_FAC_ALLOC_ALLOWED( TO_ALLOCATE,
     &             KEEP8, INFO_P(1), INFO_P(2) )
        ENDIF
        IF (INFO_P(1) .GE.0 )THEN 
          ALLOCATE(L0_OMP_FACTORS(THREAD_ID_P)%A(TO_ALLOCATE),
     &           stat=allocok)
          IF (allocok .GT. 0) THEN
            INFO_P(1) = -13
            CALL MUMPS_SETI8TOI4(TO_ALLOCATE, INFO_P(2))
            L0_OMP_FACTORS(THREAD_ID_P)%LA = 0_8
          ELSE
            L0_OMP_FACTORS(THREAD_ID_P)%LA =
     &                           MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1_8
            CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &         L0_OMP_FACTORS(THREAD_ID_P)%LA, KEEP(405).EQ.1, KEEP8,
     &         INFO_P(1), INFO_P(2), .TRUE., .FALSE. )
          ENDIF
        ENDIF
        IF (INFO_P(1) .GE.0 ) THEN
          DO I8 = 1_8, MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1_8
            L0_OMP_FACTORS(THREAD_ID_P)%A(I8) =
     &                             DMUMPS_TPS_ARR(THREAD_ID_P)%A(I8)
          ENDDO
        ENDIF
        IF ( associated(DMUMPS_TPS_ARR(THREAD_ID_P)%A)) THEN
          DEALLOCATE(DMUMPS_TPS_ARR(THREAD_ID_P)%A)
          NULLIFY(DMUMPS_TPS_ARR(THREAD_ID_P)%A)
          CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &                  -MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &                  KEEP(405).EQ.1, 
     &                  KEEP8,
     &                  INFO_P(1), INFO_P(2),
     &                  .FALSE., .FALSE. ) 
          IF (INFO_P(1) .GE. 0) THEN
!$OMP       ATOMIC UPDATE
            KEEP8(69) = KEEP8(69) - L0_OMP_FACTORS(THREAD_ID_P)%LA
!$OMP       END ATOMIC 
          ENDIF
        ENDIF
      ENDIF  
      ENDIF
      IF (INFO_P(1) .LT.0) THEN
!$OMP CRITICAL(critical_info)
        INFO(1) = INFO_P(1)
        INFO(2) = INFO_P(2)
!$OMP END CRITICAL(critical_info)
      ELSE IF (INFO_P(1) .GE. 0) THEN
!$OMP CRITICAL(critical_info)
        IF (INFO(1) .EQ. 0) THEN
          INFO(1) = INFO_P(1)
          INFO(2) = INFO_P(2)
        ENDIF
!$OMP END CRITICAL(critical_info)
      ENDIF
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$    CALL omp_set_num_threads(int(NOMP_SAVE,4))
#else
!$    CALL omp_set_num_threads(NOMP_SAVE)
#endif
!$OMP END PARALLEL
      IF (INFO(1) .LT. 0) THEN
        IF (ICNTL(1) .GT. 0 .AND. ICNTL(4) .GE.1 ) THEN
          WRITE(ICNTL(1),'(A,I6,I16,A,I5,A)')
     &    "** ERROR DURING L0_OMP: INFO(1:2)=",
     &    INFO(1), INFO(2), " (MPI worker ", MYID_NODES,")"
        ENDIF
      ENDIF
      IF ( KEEP(406) .EQ. 0 ) THEN
       ALLOCATE(ID(KEEP(400)), VAL(KEEP(400)),
     &          stat = allocok)
       IF ( allocok .GT. 0 ) THEN
          INFO(1) = -13
          INFO(2) = KEEP(400)
          GOTO 800
       ENDIF
       DO MYTHREAD_ID = 1, KEEP(400)
         VAL (MYTHREAD_ID) =  MUMPS_TPS_ARR( MYTHREAD_ID )%POSFAC-1_8
         ID (MYTHREAD_ID)  =  MYTHREAD_ID
       ENDDO
       CALL MUMPS_SORT_INT8(KEEP(400), VAL, ID) 
       DO ITH=1, KEEP(400)
        MYTHREAD_ID = ID(ITH)
        IF ((KEEP(407).NE.1) .AND. (KEEP(406).EQ.0)) THEN
          IF (INFO(1) .GE. 0) THEN
          CALL DMUMPS_DM_CBSTATIC2DYNAMIC_I  
     &    (2,   
     &     0_8, 
     &    .FALSE., 
     &    MYID_NODES, N, SLAVEF, KEEP, KEEP8,
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%IW(1),
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%LIW,
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%IWPOSCB, 
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%IWPOS,
     &   DMUMPS_TPS_ARR(MYTHREAD_ID)%A(1),
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%LA,
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%LRLU,
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%IPTRLU,
     &    MUMPS_TPS_ARR(MYTHREAD_ID)%LRLUS,
     &    STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &    INFO(1), INFO(2) )
          ENDIF
        ENDIF
        IF (KEEP(406).EQ.0) THEN
         IF (INFO(1) .GE. 0 )THEN
          TO_ALLOCATE = max(MUMPS_TPS_ARR(MYTHREAD_ID)%POSFAC-1,1_8)
          CALL DMUMPS_DM_FAC_ALLOC_ALLOWED( TO_ALLOCATE,
     &             KEEP8, INFO(1), INFO(2) )
         ENDIF
         IF (INFO(1) .GE.0 ) THEN 
           ALLOCATE(L0_OMP_FACTORS(MYTHREAD_ID)%A(TO_ALLOCATE),
     &            stat=allocok)
           IF (allocok .GT. 0) THEN
            INFO(1) = -13
            CALL MUMPS_SETI8TOI4(TO_ALLOCATE, INFO(2))
            L0_OMP_FACTORS(MYTHREAD_ID)%LA = 0_8
           ELSE
            L0_OMP_FACTORS(MYTHREAD_ID)%LA =
     &                          MUMPS_TPS_ARR(MYTHREAD_ID)%POSFAC-1_8
            CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &        L0_OMP_FACTORS(MYTHREAD_ID)%LA,
     &        KEEP(405).EQ.1, KEEP8,
     &        INFO(1), INFO(2), .TRUE., .FALSE. )
           ENDIF
         ENDIF
         IF (INFO(1) .GE. 0) THEN
!$         CHUNK8 = max( int(KEEP(361),8),
!$   &     (MUMPS_TPS_ARR(MYTHREAD_ID)%POSFAC+KEEP(400)-2_8) /
!$   &     KEEP(400) )
!$         OMP_FLAG = ( (MUMPS_TPS_ARR(MYTHREAD_ID)%POSFAC-1_8 >
!$   &                                               int(KEEP(361),8))
!$   &                  .AND. (KEEP(400).GT.1)
!$   &                 )
!$OMP      PARALLEL DO PRIVATE(I8) SCHEDULE(STATIC, CHUNK8)
!$OMP&     IF (OMP_FLAG) 
           DO I8 = 1_8, MUMPS_TPS_ARR(MYTHREAD_ID)%POSFAC-1_8
             L0_OMP_FACTORS(MYTHREAD_ID)%A(I8) =
     &                               DMUMPS_TPS_ARR(MYTHREAD_ID)%A(I8)
           ENDDO
!$OMP      END PARALLEL DO 
         ENDIF
         IF ( associated(DMUMPS_TPS_ARR(MYTHREAD_ID)%A)) THEN
           DEALLOCATE(DMUMPS_TPS_ARR(MYTHREAD_ID)%A)
           NULLIFY(DMUMPS_TPS_ARR(MYTHREAD_ID)%A)
           CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS( 
     &                     -MUMPS_TPS_ARR(MYTHREAD_ID)%LA,
     &                     KEEP(405).EQ.1, KEEP8,
     &                     INFO(1), INFO(2),
     &                     .FALSE., .FALSE. ) 
           IF (INFO(1).GE.0) THEN
             KEEP8(69) = KEEP8(69) - L0_OMP_FACTORS(MYTHREAD_ID)%LA
           ENDIF
         ENDIF
        ENDIF
       ENDDO
       IF (ALLOCATED(ID)) DEALLOCATE(ID)
       IF (ALLOCATED(VAL)) DEALLOCATE(VAL)
      ENDIF
  800 CONTINUE
      DO ITH = 1, KEEP(400)
        IF ( associated(DMUMPS_TPS_ARR(ITH)%A)) THEN
           DEALLOCATE(DMUMPS_TPS_ARR(ITH)%A)
           NULLIFY(DMUMPS_TPS_ARR(ITH)%A)
           CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS( 
     &                     -MUMPS_TPS_ARR(ITH)%LA,
     &                     KEEP(405).EQ.1, KEEP8,
     &                     INFO(1), INFO(2),
     &                     .FALSE., .FALSE. ) 
           IF (INFO(1).GE.0) THEN
             KEEP8(69) = KEEP8(69) - L0_OMP_FACTORS(ITH)%LA
           ENDIF
        ENDIF
      ENDDO
      KEEP8(64) = 0_8
      DO I = 1, KEEP(400)
        KEEP8(64) = KEEP8(64) + MUMPS_TPS_ARR(I)%POSFAC - 1_8
      ENDDO
      KEEP8(62) = 0_8
      DO I = 1, KEEP(400)
        KEEP8(62) = KEEP8(62) + MUMPS_TPS_ARR(I)%LRLUSM
      ENDDO
      NBROOT_UNDER_L0 = NBROOT_PROCESSED
      DKEEP(95) = MPI_WTIME() - GTIME 
      IF (KEEP(486) .NE. 0) THEN
      TIME_UPDATE          = TIME_UPDATE/dble(KEEP(400))
      TIME_COMPRESS        = TIME_COMPRESS/dble(KEEP(400))
      TIME_FRSWAP_COMPRESS = TIME_FRSWAP_COMPRESS/dble(KEEP(400))
      TIME_CB_COMPRESS     = TIME_CB_COMPRESS/dble(KEEP(400))
      TIME_PANEL           = TIME_PANEL/dble(KEEP(400))
      TIME_FAC_I           = TIME_FAC_I/dble(KEEP(400))
      TIME_FAC_MQ          = TIME_FAC_MQ/dble(KEEP(400))
      TIME_FAC_SQ          = TIME_FAC_SQ/dble(KEEP(400))
      TIME_FRFRONTS        = TIME_FRFRONTS/dble(KEEP(400))
      TIME_LRTRSM          = TIME_LRTRSM/dble(KEEP(400))
      TIME_FRTRSM          = TIME_FRTRSM/dble(KEEP(400))
      TIME_LR_MODULE       = TIME_LR_MODULE/dble(KEEP(400))
      TIME_DECOMP          = TIME_DECOMP/dble(KEEP(400))
      TIME_DIAGCOPY        = TIME_DIAGCOPY/dble(KEEP(400))
      TIME_DECOMP_UCFS     = TIME_DECOMP_UCFS/dble(KEEP(400))
      TIME_LRASM_NIV1      = TIME_LRASM_NIV1/dble(KEEP(400))
      TIME_LRASM_LOCASM2   = TIME_LRASM_LOCASM2/dble(KEEP(400))
      TIME_LRASM_MAPLIG1   = TIME_LRASM_MAPLIG1/dble(KEEP(400))
      TIME_LRASM_CONTRIB2  = TIME_LRASM_CONTRIB2/dble(KEEP(400))
      TIME_FRASM_LOCASM2   = TIME_FRASM_LOCASM2/dble(KEEP(400))
      TIME_FRASM_MAPLIG1   = TIME_FRASM_MAPLIG1/dble(KEEP(400))
      TIME_FRASM_CONTRIB2  = TIME_FRASM_CONTRIB2/dble(KEEP(400))
      ENDIF
      DKEEP(97)            = DKEEP(97) / dble(KEEP(400))
      CALL MUMPS_LOAD_ENABLE()
      CALL DMUMPS_LOAD_UPDATE(0,.FALSE., FLOP_ESTIM_ACC,KEEP,KEEP8)
      RETURN
      END SUBROUTINE DMUMPS_FAC_L0_OMP
      SUBROUTINE DMUMPS_SET_MAXS_MAXIS_THREAD(MAXS_BASE_RELAXED8TH,
     &     MAXIS_BASE_RELAXEDTH, BLR_STRAT,
     &     KEEP,
     &     I4_L0_OMPTH, NBSTATS_I4,
     &     I8_L0_OMPTH, NBSTATS_I8)
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: KEEP(500), NBSTATS_I4, NBSTATS_I8
      INTEGER, INTENT(IN)     :: I4_L0_OMPTH(NBSTATS_I4)
      INTEGER(8), INTENT(IN)  :: I8_L0_OMPTH(NBSTATS_I8)
      INTEGER(8), INTENT(OUT) :: MAXS_BASE_RELAXED8TH
      INTEGER, INTENT(OUT)    :: MAXIS_BASE_RELAXEDTH
      INTEGER, INTENT(OUT)    :: BLR_STRAT 
      INTEGER    :: PERLU
      INTEGER(8) :: MAXS_BASE8TH
      INTEGER(8) :: MAXIS_BASE_RELAXEDTH8
      PERLU = KEEP(12)
      CALL DMUMPS_SET_BLRSTRAT_AND_MAXS ( MAXS_BASE8TH,
     &     MAXS_BASE_RELAXED8TH, BLR_STRAT, KEEP(1),
     &     I8_L0_OMPTH(2), I8_L0_OMPTH(3), I8_L0_OMPTH(5),
     &     I8_L0_OMPTH(6), I8_L0_OMPTH(7), I8_L0_OMPTH(8) )
      IF ( KEEP(201) .EQ. 0 ) THEN
        MAXIS_BASE_RELAXEDTH8 = int(I4_L0_OMPTH(2),8)  
      ELSE
        MAXIS_BASE_RELAXEDTH8 = int(I4_L0_OMPTH(4),8)  
      ENDIF
      MAXIS_BASE_RELAXEDTH8 = max( 1_8,
     &            MAXIS_BASE_RELAXEDTH8 + 3 * max(PERLU,10) *
     &            ( MAXIS_BASE_RELAXEDTH8 / 100 + 1 )
     &  )
      MAXIS_BASE_RELAXEDTH8 = min(MAXIS_BASE_RELAXEDTH8,
     &                            int( huge( MAXIS_BASE_RELAXEDTH ) ,8)
     &                           )
      MAXIS_BASE_RELAXEDTH  = int( MAXIS_BASE_RELAXEDTH8 )
      RETURN
      END SUBROUTINE DMUMPS_SET_MAXS_MAXIS_THREAD
      SUBROUTINE DMUMPS_MA_EFF_MEM_DISPO(
     &     MUMPS_TPS_ARR, NBTHREADS, KEEP8, KEEP,
     &     N, BLR_STRAT, LPOOL_P,
     &     I8_L0_OMP, NBSTATS_I8,
     &     MEMDISPO_UNDERL0)
      USE MUMPS_TPS_M, ONLY : MUMPS_TPS_T
      IMPLICIT NONE
      INTEGER, INTENT(IN)     :: N, BLR_STRAT, KEEP(500)
      INTEGER, INTENT(IN)     :: NBSTATS_I8, NBTHREADS, LPOOL_P
      INTEGER(8), INTENT(IN)  :: KEEP8(150)
      INTEGER(8), INTENT(IN)  :: I8_L0_OMP(NBSTATS_I8,NBTHREADS)
      INTEGER(8), INTENT(OUT) :: MEMDISPO_UNDERL0
      TYPE (MUMPS_TPS_T), INTENT(IN) :: MUMPS_TPS_ARR(:)
      INTEGER    :: PERLU, ITH, ITHMIN, ITHMIN_if_LRLU, OOC_STRAT
      INTEGER(8) :: TO_ALLOCATE, BLR_RELATED, COPY_RELATED
      INTEGER(8) :: MIN_NRLADU_underL0, MIN_NRLADU_if_LR_LU_underL0
      PERLU     = KEEP(12)
      OOC_STRAT = KEEP(201)
      TO_ALLOCATE = 0_8
      DO ITH = 1, NBTHREADS
         TO_ALLOCATE = TO_ALLOCATE +
     &   ((int(MUMPS_TPS_ARR(ITH)%LIW,8) * int(KEEP(34),8 )) /
     &     int(KEEP(35),8 ))
     &   + MUMPS_TPS_ARR(ITH)%LA
      ENDDO
      TO_ALLOCATE = TO_ALLOCATE + int(NBTHREADS,8)* (
     & ((int(LPOOL_P,8) * int(KEEP(34),8 )) / int(KEEP(35),8 )) +
     & ((int(N+KEEP(253),8) * int(KEEP(34),8 )) / int(KEEP(35),8 ))
     & )
      BLR_RELATED = 0_8
      DO ITH = 1, NBTHREADS
         IF (BLR_STRAT.EQ.1) THEN
           BLR_RELATED = BLR_RELATED +
     &        I8_L0_OMP(10,ITH) + 
     &        int(PERLU,8) * ( I8_L0_OMP(10,ITH) / 100_8 + 1_8)
         ELSE  IF (BLR_STRAT.EQ.2) THEN
           BLR_RELATED = BLR_RELATED +
     &        I8_L0_OMP(13,ITH) + 
     &        int(PERLU,8) * ( I8_L0_OMP(13,ITH) / 100_8 + 1_8)
         ELSE IF (BLR_STRAT.EQ.3) THEN
           BLR_RELATED = BLR_RELATED +
     &        I8_L0_OMP(8,ITH) + 
     &        int(PERLU,8) * ( I8_L0_OMP(8,ITH) / 100_8 + 1_8)
         ENDIF
      ENDDO
      COPY_RELATED = 0_8
      ITHMIN = 1
      ITHMIN_if_LRLU = 1 
      MIN_NRLADU_if_LR_LU_underL0 = I8_L0_OMP(4,1)
      MIN_NRLADU_underL0          = I8_L0_OMP(1,1)
      DO ITH = 1, NBTHREADS
        IF (I8_L0_OMP(1,ITH).LT.MIN_NRLADU_underL0)
     &      THEN
             MIN_NRLADU_underL0 = I8_L0_OMP(1,ITH)
             ITHMIN = ITH
        ENDIF
        IF (I8_L0_OMP(4,ITH).LT.MIN_NRLADU_if_LR_LU_underL0) 
     &      THEN
            MIN_NRLADU_if_LR_LU_underL0 = I8_L0_OMP(4,ITH)
            ITHMIN_if_LRLU = ITH
        ENDIF
      ENDDO
      IF (BLR_STRAT.EQ.0) THEN
        IF (OOC_STRAT .LE. 0 .AND. OOC_STRAT .NE. -1 ) THEN
         COPY_RELATED = COPY_RELATED +     
     &         I8_L0_OMP(1,ITHMIN) +  
     &         I8_L0_OMP(23, ITHMIN)  
        ELSE
         COPY_RELATED = COPY_RELATED +     
     &         I8_L0_OMP(23, ITHMIN)  
        ENDIF
      ELSE
        IF (OOC_STRAT .LE. 0 .AND. OOC_STRAT .NE. -1 ) THEN
         COPY_RELATED = COPY_RELATED +     
     &    I8_L0_OMP(4,ITHMIN_if_LRLU) + 
     &    I8_L0_OMP(23,ITHMIN_if_LRLU )  
        ELSE
         COPY_RELATED = COPY_RELATED +     
     &    I8_L0_OMP(23, ITHMIN_if_LRLU)  
        ENDIF
      ENDIF
      COPY_RELATED = COPY_RELATED +
     &   int(PERLU,8)*(COPY_RELATED / 100_8 + 1_8 )
      TO_ALLOCATE       = TO_ALLOCATE + COPY_RELATED + BLR_RELATED
      MEMDISPO_UNDERL0  = KEEP8(75) - TO_ALLOCATE
      RETURN
      END SUBROUTINE DMUMPS_MA_EFF_MEM_DISPO
      SUBROUTINE DMUMPS_L0OMP_COPY_IW( IW, LIW, IWPOS,
     &                                 MUMPS_TPS_ARR, KEEP,
     &                                 PTLUST, ICNTL, INFO )
      USE MUMPS_TPS_M, ONLY : MUMPS_TPS_T
      IMPLICIT NONE
      INTEGER                     :: KEEP(500)
      INTEGER, INTENT( IN )       :: LIW
      INTEGER, INTENT( INOUT )    :: IW(:) 
      INTEGER, INTENT( INOUT )    :: IWPOS
      INTEGER, INTENT( INOUT )    :: PTLUST(KEEP(28))
      INTEGER, INTENT( IN )       :: ICNTL(60)
      INTEGER, INTENT( INOUT )    :: INFO(80)
      TYPE (MUMPS_TPS_T), TARGET  :: MUMPS_TPS_ARR(:)
      INTEGER :: ITHREAD, JTHREAD
      INTEGER :: REQUESTED_SIZE
      INTEGER :: IWPOS_TO_COPY 
      INTEGER :: LOC_IPOS
      INTEGER :: LOC_SIZE, LOC_ISTEP
      TYPE (MUMPS_TPS_T), POINTER :: MUMPS_TPS
      INCLUDE 'mumps_headers.h'
      REQUESTED_SIZE = 0
      DO ITHREAD = 1, size(MUMPS_TPS_ARR)
        MUMPS_TPS => MUMPS_TPS_ARR(ITHREAD)
        REQUESTED_SIZE = REQUESTED_SIZE + MUMPS_TPS%IWPOS - 1
      ENDDO
      IF ( LIW - IWPOS + 1 .LT. REQUESTED_SIZE ) THEN
        WRITE(*,*) " LIW too small in DMUMPS_L0OMP_COPY_IW !!", LIW,
     &  REQUESTED_SIZE
        INFO(1) = -8
        INFO(2) = REQUESTED_SIZE-LIW+IWPOS-1
        IF ( ICNTL(1) .GT. 0 .AND. ICNTL(4).GE.1) THEN
          WRITE(ICNTL(1),*) " ** ERROR IN DMUMPS_L0OMP_COPY_IW: ",
     &    "LIW TOO SMALL TO COPY LOCAL FACTOR INFORMATION",
     &    INFO(2)
        ENDIF
        GOTO 500
      ENDIF
      DO ITHREAD = 1, size(MUMPS_TPS_ARR)
        MUMPS_TPS => MUMPS_TPS_ARR(ITHREAD)
        IWPOS_TO_COPY = IWPOS
        DO JTHREAD=1, ITHREAD - 1
          IWPOS_TO_COPY = IWPOS_TO_COPY+MUMPS_TPS_ARR(JTHREAD)%IWPOS-1
        ENDDO
        IW(IWPOS_TO_COPY: IWPOS_TO_COPY+MUMPS_TPS%IWPOS - 2) =
     &                              MUMPS_TPS%IW(1:MUMPS_TPS%IWPOS-1)
        LOC_IPOS = 1
        DO WHILE ( LOC_IPOS .NE. MUMPS_TPS%IWPOS )
          LOC_SIZE  = MUMPS_TPS%IW(LOC_IPOS+XXI)
          LOC_ISTEP = MUMPS_TPS%IW(LOC_IPOS+KEEP(IXSZ)+4)
          PTLUST(LOC_ISTEP) = IWPOS_TO_COPY+LOC_IPOS-1
          LOC_IPOS = LOC_IPOS + LOC_SIZE
        ENDDO
      ENDDO
      IWPOS = IWPOS + REQUESTED_SIZE
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_L0OMP_COPY_IW
      SUBROUTINE DMUMPS_PERFORM_COPIES_INIT( STATE, KEEP8_77_SAVE,
     & NbWaitMem, NbFinished, NbOnGoingCopies, NbUnderL0,
     & KEEP, KEEP8 )
      INTEGER, INTENT(IN)    :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(OUT)   :: NbWaitMem,
     &                          NbFinished,
     &                          NbOnGoingCopies,
     &                          NbUnderL0
      INTEGER(8), INTENT(OUT) :: STATE(KEEP(400)), KEEP8_77_SAVE
      INTEGER :: ITH
      NbWaitMem       = 0
      NbFinished      = 0
      NbOnGoingCopies = 0
      NbUnderL0       = KEEP(400)
      DO ITH=1, KEEP(400)
        STATE(ITH) = UnderL0
      ENDDO
      KEEP8_77_SAVE = KEEP8(77)
      RETURN
      END SUBROUTINE DMUMPS_PERFORM_COPIES_INIT
      SUBROUTINE DMUMPS_PERFORM_COPIES( THREAD_ID_P, 
     &           MUMPS_TPS_ARR, DMUMPS_TPS_ARR,
     &           L0_OMP_FACTORS, LL0_OMP_FACTORS,
     &           STATE, SIZE_COPIED,
     &           NbWaitMem, NbFinished, NbOngoingCopies, NbUnderL0,
     &           MYID_NODES, N, SLAVEF,
     &           STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &           KEEP, KEEP8, INFO_P
     &           )
      USE MUMPS_TPS_M, ONLY : MUMPS_TPS_T
      USE DMUMPS_TPS_M, ONLY : DMUMPS_TPS_T
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_L0OMPFAC_T
      INTEGER :: KEEP(500)
      INTEGER(8) :: KEEP8(150)
      INTEGER, INTENT(IN) :: THREAD_ID_P
      INTEGER, INTENT(INOUT) :: INFO_P(2)
      INTEGER, INTENT(IN) :: MYID_NODES, N, SLAVEF
      INTEGER, INTENT(IN) ::  STEP(N), DAD(KEEP(28))
      INTEGER(8), INTENT(IN) :: PAMASTER(KEEP(28)), PTRAST(KEEP(28))
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT(INOUT) :: NbWaitMem,
     &                          NbFinished,
     &                          NbOnGoingCopies,
     &                          NbUnderL0
      INTEGER(8), INTENT(INOUT) :: STATE( KEEP(400) )
      INTEGER(8), INTENT(INOUT) :: SIZE_COPIED(KEEP(400) )
      TYPE (MUMPS_TPS_T), DIMENSION(:) :: MUMPS_TPS_ARR
      TYPE (DMUMPS_TPS_T), DIMENSION(:) :: DMUMPS_TPS_ARR
      INTEGER, INTENT ( IN )  :: LL0_OMP_FACTORS
      TYPE (DMUMPS_L0OMPFAC_T), INTENT(INOUT) ::
     &                        L0_OMP_FACTORS(LL0_OMP_FACTORS)
      INTEGER :: NbFinishedPrivateCopy 
      INTEGER :: LOCAL_ACTION
      INTEGER, PARAMETER :: NOTHING            = 0
      INTEGER, PARAMETER :: FREE_WORK_MYID     = 1
      INTEGER, PARAMETER :: COPY_FACTORS       = 2
      INTEGER, PARAMETER :: AllocateViderCB    = 3
      INTEGER, PARAMETER :: DORMIR             = 4
      INTEGER(8)       :: COPY_START, CHUNK8, I8, TO_ALLOCATE
      INTEGER          :: ITH, K
      INTEGER          :: allocok
      INTEGER(8) :: PeakAuthorized_P
      INTEGER(8) :: MemNeeded_P, MemNeededForCB_P, MemDispo_P,
     &        CBCopiedToDynamic_P, LRLUS_SAVE_P
      INTEGER(8) :: KEEP8_71, KEEP8_73
!$OMP CRITICAL(L0_COPIES)
      STATE(THREAD_ID_P) = CopyNotStarted
      IF ( INFO_P(1) .LT. 0 ) THEN
        NbFinished = NbFinished + 1
        STATE(THREAD_ID_P) = Finished
      ENDIF
      DO ITH = 1, KEEP(400)
       IF ( STATE(ITH) .EQ. WaitMem ) THEN
         STATE(ITH)=CopyNotStarted
       ENDIF
      ENDDO
      NbWaitMem       = 0
      NbFinishedPrivateCopy = NbFinished
!$OMP END CRITICAL(L0_COPIES)
      DO WHILE ( NbFinishedPrivateCopy .NE. KEEP(400) )
        LOCAL_ACTION = DORMIR
!$OMP CRITICAL(L0_COPIES)
        NbFinishedPrivateCopy = NbFinished
        IF ( NbFinished.EQ. KEEP(400)) THEN
         LOCAL_ACTION = NOTHING 
        ELSE IF ( (NbFinished+NbWaitMem) .EQ. KEEP(400) ) THEN
!$OMP ATOMIC READ
            KEEP8_73 = KEEP8(73)
!$OMP END ATOMIC
!$OMP ATOMIC READ
            KEEP8_71 = KEEP8(71)
!$OMP END ATOMIC
            MemDispo_P  = KEEP8(77)  - (KEEP8_73 -KEEP8_71)
            MemDispo_P  = min(MemDispo_P, KEEP8(75)-KEEP8_73)
            MemNeeded_P = huge(MemNeeded_P)
            DO ITH = 1, KEEP(400)
             IF (STATE(ITH).EQ.WaitMem) THEN
               MemNeeded_P = min( MemNeeded_P,
     &                     MUMPS_TPS_ARR(ITH)%LA -
     &                     MUMPS_TPS_ARR(ITH)%LRLUS )
             ENDIF
            ENDDO
            IF ((KEEP8(75)-KEEP8_73).LT.MemNeeded_P) THEN
              INFO_P(1)  = -19
              CALL MUMPS_SET_IERROR (
     &             MemNeeded_P-(KEEP8(75)-KEEP8_73), INFO_P(2))
              DO ITH = 1, KEEP(400)
               STATE(ITH) = Finished
              ENDDO
              NbFinished = KEEP(400)
            ELSE
              KEEP8(77) = MemNeeded_P + (KEEP8_73 -KEEP8_71)
              DO ITH = 1, KEEP(400)
                IF ( STATE(ITH) .EQ. WaitMem ) THEN
                  STATE(ITH)=CopyNotStarted
                ENDIF
              ENDDO
              NbWaitMem = 0
            ENDIF
            LOCAL_ACTION = NOTHING
        ELSE
         SELECT CASE (STATE(THREAD_ID_P))
          CASE ( CopyFactorsFinished )
           LOCAL_ACTION = FREE_WORK_MYID
          CASE ( CopyNotStarted )
!$OMP ATOMIC READ
           KEEP8_73 = KEEP8(73)
!$OMP END ATOMIC
!$OMP ATOMIC READ
           KEEP8_71 = KEEP8(71)
!$OMP END ATOMIC
           PeakAuthorized_P = KEEP8(77)
           MemDispo_P  = PeakAuthorized_P - (KEEP8_73 -KEEP8_71)
           MemDispo_P  = min(MemDispo_P, KEEP8(75)-KEEP8_73)
           MemNeeded_P = MUMPS_TPS_ARR(THREAD_ID_P)%LA -
     &                  MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS
           MemNeededForCB_P = MemNeeded_P - 
     &                   ( MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC - 1_8 )
           IF ( MemDispo_P .GE. MemNeeded_P ) THEN
!$OMP ATOMIC CAPTURE
             KEEP8(73) = KEEP8(73) + MemNeeded_P
             KEEP8_73  = KEEP8(73)
!$OMP END ATOMIC
!$OMP ATOMIC READ
             KEEP8_71 = KEEP8(71)
!$OMP END ATOMIC
             MemDispo_P  = PeakAuthorized_P - (KEEP8_73 -KEEP8_71)
             MemDispo_P  = min(MemDispo_P, KEEP8(75)-KEEP8_73)
             IF ( MemDispo_P .LT. 0 ) THEN
!$OMP ATOMIC UPDATE
               KEEP8(73) = KEEP8(73) - MemNeeded_P
!$OMP END ATOMIC
               IF ( STATE(THREAD_ID_P) .NE. WaitMem ) THEN
                 STATE( THREAD_ID_P ) = WaitMem
                 NbWaitMem = NbWaitMem + 1
               ENDIF
             ELSE
!$OMP ATOMIC UPDATE
              KEEP8(74) = max(KEEP8(74), KEEP8_73 )
!$OMP END ATOMIC
              IF ( STATE( THREAD_ID_P ) .EQ. WaitMem ) THEN
               NbWaitMem = NbWaitMem - 1
              ENDIF
              STATE( THREAD_ID_P ) = AllocateViderCBEnCours
              LOCAL_ACTION    = AllocateViderCB
              NbOngoingCopies = NbOnGoingCopies + 1
             ENDIF
           ELSE
            IF ( STATE(THREAD_ID_P) .NE. WaitMem ) THEN
              STATE( THREAD_ID_P ) = WaitMem
              NbWaitMem = NbWaitMem + 1
            ENDIF
           ENDIF
          CASE DEFAULT
           ITH = -1
           DO K = THREAD_ID_P, THREAD_ID_P + KEEP(400) - 1
            IF ( K > KEEP(400) ) THEN
              ITH = K - KEEP(400)
            ELSE 
              ITH = K
            ENDIF
            IF ( STATE(ITH) .GE. 0 .AND. 
     &           STATE(ITH) .LT. MUMPS_TPS_ARR(ITH)%POSFAC - 1_8 ) THEN
              EXIT
            ELSE
              ITH = -1
            ENDIF
           ENDDO
           IF ( ITH .GT. 0 ) THEN
            LOCAL_ACTION = COPY_FACTORS
            COPY_START = STATE(ITH) + 1
            CHUNK8 = max(
     &
     &      int(KEEP(361),8),
     &
     &       (MUMPS_TPS_ARR(ITH)%POSFAC+KEEP(400)-2_8) / 
     &        (int(KEEP(400)*2,8))
     &
     &      )
            IF (KEEP(72) .EQ. 1) THEN
              CHUNK8 = 4_8
            ENDIF
            CHUNK8 = min( CHUNK8,
     &      MUMPS_TPS_ARR(ITH)%POSFAC - 1_8 - COPY_START + 1_8
     &      )
            STATE(ITH) = STATE(ITH) + CHUNK8
           ENDIF
         END SELECT
        ENDIF
!$OMP END CRITICAL(L0_COPIES)
        SELECT CASE ( LOCAL_ACTION )
          CASE ( FREE_WORK_MYID )
            IF ( associated(DMUMPS_TPS_ARR(THREAD_ID_P)%A) ) THEN
              DEALLOCATE(DMUMPS_TPS_ARR(THREAD_ID_P)%A)
              NULLIFY(DMUMPS_TPS_ARR(THREAD_ID_P)%A)
              CALL MUMPS_DM_FAC_UPD_DYN_MEMCNTS(
     &                  -MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &                  KEEP(405).EQ.1, 
     &                  KEEP8,
     &                  INFO_P(1), INFO_P(2),
     &                  .FALSE., .FALSE. ) 
              IF (INFO_P(1) .GE. 0) THEN
!$OMP           ATOMIC UPDATE
                KEEP8(69) = KEEP8(69) - L0_OMP_FACTORS(THREAD_ID_P)%LA
!$OMP           END ATOMIC
!$OMP CRITICAL(L0_COPIES)
                DO ITH = 1, KEEP(400)
                  IF ( STATE(ITH) .EQ. WaitMem ) THEN
                    STATE(ITH)=CopyNotStarted
                  ENDIF
                ENDDO
                NbWaitMem       = 0
                NbFinished      = NbFinished + 1
                STATE( THREAD_ID_P ) = Finished
                NbOnGoingCopies = NbOnGoingCopies -1
!$OMP END CRITICAL(L0_COPIES)
              ENDIF
            ENDIF
          CASE ( AllocateViderCB )
          TO_ALLOCATE = max(MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1_8,1_8)
          ALLOCATE( L0_OMP_FACTORS(THREAD_ID_P)%A(TO_ALLOCATE),
     &              stat=allocok )
          IF ( allocok .GT. 0 ) THEN
            INFO_P(1) = -13
            CALL MUMPS_SETI8TOI4(TO_ALLOCATE, INFO_P(2))
            L0_OMP_FACTORS(THREAD_ID_P)%LA = 0_8
!$OMP CRITICAL(L0_COPIES)
!$OMP ATOMIC UPDATE
            KEEP8(73) = KEEP8(73) - MemNeeded_P
!$OMP END ATOMIC
            STATE(THREAD_ID_P) = Finished
            NbFinished         = NbFinished + 1
!$OMP END CRITICAL(L0_COPIES)
          ELSE
            L0_OMP_FACTORS(THREAD_ID_P)%LA =
     &                          MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1_8
!$OMP ATOMIC UPDATE
            KEEP8(69) = KEEP8(69) + L0_OMP_FACTORS(THREAD_ID_P)%LA
!$OMP END ATOMIC
!$OMP CRITICAL(L0_COPIES)
            IF ( MUMPS_TPS_ARR(THREAD_ID_P)%POSFAC-1_8 == 0_8 ) THEN
              STATE(THREAD_ID_P) = CopyFactorsFinished
            ELSE
              STATE      ( THREAD_ID_P ) = 0
              SIZE_COPIED( THREAD_ID_P ) = 0
            ENDIF
!$OMP END CRITICAL(L0_COPIES)
            LRLUS_SAVE_P = MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS
            CALL DMUMPS_DM_CBSTATIC2DYNAMIC_I  
     &      (3,   
     &       0_8, 
     &      .FALSE., 
     &      MYID_NODES, N, SLAVEF, KEEP, KEEP8,
     &      MUMPS_TPS_ARR(THREAD_ID_P)%IW(1),
     &      MUMPS_TPS_ARR(THREAD_ID_P)%LIW,
     &      MUMPS_TPS_ARR(THREAD_ID_P)%IWPOSCB, 
     &      MUMPS_TPS_ARR(THREAD_ID_P)%IWPOS,
     &     DMUMPS_TPS_ARR(THREAD_ID_P)%A(1),
     &      MUMPS_TPS_ARR(THREAD_ID_P)%LA,
     &      MUMPS_TPS_ARR(THREAD_ID_P)%LRLU,
     &      MUMPS_TPS_ARR(THREAD_ID_P)%IPTRLU,
     &      MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS,
     &      STEP, PTRAST, PAMASTER, PROCNODE_STEPS, DAD,
     &      INFO_P(1), INFO_P(2) )
            CBCopiedToDynamic_P = 
     &      MUMPS_TPS_ARR(THREAD_ID_P)%LRLUS - LRLUS_SAVE_P
            IF (INFO_P(1) .LT. 0 ) THEN
!$OMP CRITICAL(L0_COPIES)
!$OMP ATOMIC UPDATE
              KEEP8(73) = KEEP8(73) -
     &                      ( MemNeededForCB_P - CBCopiedToDynamic_P )
!$OMP END ATOMIC
              STATE(THREAD_ID_P) = Finished
              NbFinished         = NbFinished+1
!$OMP END CRITICAL(L0_COPIES)
            ELSE
            ENDIF
          ENDIF
          CASE ( COPY_FACTORS )
          DO I8 = COPY_START, COPY_START + CHUNK8 - 1
            L0_OMP_FACTORS(ITH)%A(I8) = DMUMPS_TPS_ARR(ITH)%A(I8)
          ENDDO
!$OMP CRITICAL(L0_COPIES)
          SIZE_COPIED(ITH) = SIZE_COPIED(ITH) + CHUNK8
          IF ( SIZE_COPIED(ITH) .EQ. L0_OMP_FACTORS(ITH)%LA ) THEN
            STATE(ITH) = CopyFactorsFinished
          ENDIF
!$OMP END CRITICAL(L0_COPIES)
          CASE ( NOTHING )
          CASE ( DORMIR )
          CALL MUMPS_USLEEP(1000)
          CASE DEFAULT
          WRITE(*,*) " Internal error in DMUMPS_PERFORM_COPIES",
     &    LOCAL_ACTION
        END SELECT
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_PERFORM_COPIES
      END MODULE DMUMPS_FAC_OMP_M
      RECURSIVE SUBROUTINE DMUMPS_PROCESS_FRONT_NIV1( COMM_LOAD,
     & ASS_IRECV, N, INODE, TYPE, TYPEF, LA, IW, LIW, A,
     & MAXFRW, NOFFNEGW, NULLNEGW, NPVW, NB22T1W, NBTINYW,
     & DET_EXPW, DET_MANTW, DET_SIGNW,
     & INFO_P, UU, SEUIL, SEUIL_LDLT_NIV2,
     & OPELIW, NELVAW, NMAXNPIVW, NSTEPSW, PTRIST, PTLUST_S,
     & PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER, NE, POSFAC, 
     & LRLU, LRLUS, LRLUSM, IPTRLU, ICNTL, KEEP, KEEP8, DKEEP,
     & PIVNUL_LIST_STRUCT, COMP, IWPOS, IWPOSCB, PROCNODE_STEPS,
     & SLAVEF, IFATH, COMM_NODES, MYID_NODES, IPOOL_P,
     & LPOOL_P, LEAF, PERM, NSTK_STEPS, BUFR, LBUFR,
     & LBUFR_BYTES, NBFIN, root, OPASSW, ITLOC, RHS_MUMPS, FILS,
     & PTRARW, PTRAIW, PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     & INTARR, DBLARR, ND, FRERE, DAD, LPTRAR, NELT,
     & FRTPTR, FRTELT, ISTEP_TO_INIV2, TAB_POS_IN_PERE, LRGROUPS,
     & FLOP_ESTIM_ACC )
      USE DMUMPS_FAC_ASM_MASTER_M
      USE DMUMPS_FAC_ASM_MASTER_ELT_M
      USE DMUMPS_FAC1_LU_M
      USE DMUMPS_FAC1_LDLT_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      USE MUMPS_PIVNUL_MOD, ONLY : PIVNUL_LIST_STRUCT_T
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER COMM_LOAD, ASS_IRECV
      INTEGER COMM_NODES, MYID_NODES, TYPE, TYPEF
      INTEGER N, LIW, INODE,INFO_P(2)
      INTEGER ICNTL(60), KEEP(500)
      DOUBLE PRECISION    DKEEP(230)
      DOUBLE PRECISION UU, SEUIL, SEUIL_LDLT_NIV2
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: LA, POSFAC, LRLU, LRLUS, LRLUSM, IPTRLU
      INTEGER IWPOSCB, IWPOS,
     &        IFATH, SLAVEF, NELVAW, NMAXNPIVW, NSTEPSW
      INTEGER IW(LIW),PROCNODE_STEPS(KEEP(28))
      INTEGER(8) :: PTRAST  (KEEP(28))
      INTEGER(8) :: PTRFAC  (KEEP(28))
      INTEGER(8) :: PAMASTER(KEEP(28))
      INTEGER PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER STEP(N), PIMASTER(KEEP(28)), NE(KEEP(28))
      DOUBLE PRECISION    A(LA)
      INTEGER :: MAXFRW, NOFFNEGW, NULLNEGW, NPVW, NBTINYW
      INTEGER, intent(in) :: LRGROUPS(KEEP(280))
      DOUBLE PRECISION OPASSW, OPELIW
      DOUBLE PRECISION DBLARR(KEEP8(26))
      INTEGER INTARR(KEEP8(27)) 
      INTEGER ITLOC( N + KEEP(253) ), FILS( N ),
     &        ND( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER DAD( KEEP(28) )
      DOUBLE PRECISION :: RHS_MUMPS(KEEP8(85))
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER NELT, LPTRAR
      INTEGER FRTPTR( N+1 ), FRTELT( NELT )
      INTEGER(8), INTENT(IN) :: PTRAIW( LPTRAR ), PTRARW( LPTRAR )
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER LEAF, COMP
      INTEGER :: NB22T1W, DET_EXPW, DET_SIGNW
      DOUBLE PRECISION :: DET_MANTW
      INTEGER PERM( N )
      INTEGER NSTK_STEPS( KEEP(28) )
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER NBFIN
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
      DOUBLE PRECISION FLOP_ESTIM_ACC
      INTEGER, INTENT(IN) :: LPOOL_P
      INTEGER, INTENT(IN) :: IPOOL_P(LPOOL_P)
      INTEGER :: IOLDPS, JOBASS, ETATASS
      INTEGER(8) :: POSELT
      LOGICAL :: AVOID_DELAYED, SON_LEVEL2
      JOBASS  = 0
      ETATASS = 0
      IF ( KEEP(55) .EQ. 0 ) THEN
          JOBASS = 0
          CALL DMUMPS_FAC_ASM_NIV1( COMM_LOAD, ASS_IRECV, N, INODE,
     &    IW, LIW, A, LA,
     &    INFO_P, ND, FILS, FRERE, DAD, MAXFRW, root, OPASSW, OPELIW,
     &    PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER,
     &    PTRARW, PTRAIW, PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    ITLOC, RHS_MUMPS, NSTEPSW,
     &    SON_LEVEL2,COMP, LRLU, IPTRLU,
     &    IWPOS, IWPOSCB, POSFAC,
     &    LRLUS, LRLUSM, ICNTL, KEEP, KEEP8, DKEEP, 
     &    INTARR, KEEP8(27), DBLARR, KEEP8(26),
     &    NSTK_STEPS, PROCNODE_STEPS, SLAVEF, COMM_NODES,
     &    MYID_NODES, BUFR, LBUFR, LBUFR_BYTES, NBFIN, IPOOL_P,
     &    LPOOL_P, LEAF, PERM, ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &    JOBASS, ETATASS
     &    , LRGROUPS
     &    )
      ELSE
          CALL DMUMPS_FAC_ASM_NIV1_ELT(COMM_LOAD,ASS_IRECV,NELT,FRTPTR,
     &    FRTELT, N, INODE, IW, LIW, A,
     &    LA, INFO_P, ND, FILS, FRERE, DAD, MAXFRW,
     &    root, OPASSW, OPELIW, PTRIST, PTLUST_S, PTRFAC, PTRAST, STEP,
     &    PIMASTER, PAMASTER, PTRARW, PTRAIW,
     &    PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &    ITLOC, RHS_MUMPS, NSTEPSW, SON_LEVEL2, COMP, LRLU,
     &    IPTRLU, IWPOS, IWPOSCB,
     &    POSFAC, LRLUS, LRLUSM, ICNTL, KEEP, KEEP8, DKEEP,
     &    INTARR, KEEP8(27), DBLARR, KEEP8(26),
     &    NSTK_STEPS, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES, NBFIN,
     &    IPOOL_P, LPOOL_P, LEAF, PERM, ISTEP_TO_INIV2,
     &    TAB_POS_IN_PERE
     &    , LRGROUPS
     &    )
      ENDIF
      IF (INFO_P(1) .LT. 0) THEN
        RETURN
      ENDIF
      AVOID_DELAYED = ( ( IFATH .EQ. KEEP(20)
     &                  .OR.
     &                    IFATH .EQ. KEEP(38) )
     &                .AND.
     &                  ( KEEP(60) .NE. 0 ) )
      POSELT = PTRAST(STEP(INODE))
      IOLDPS = PTLUST_S(STEP(INODE))
      IF ( KEEP(50) .EQ. 0 ) THEN
          CALL DMUMPS_FAC1_LU( N, INODE,
     &    IW, LIW,
     &    A, LA, IOLDPS,
     &    POSELT,
     &    INFO_P(1), INFO_P(2), UU, NOFFNEGW, NPVW, NBTINYW,
     &    DET_EXPW, DET_MANTW, DET_SIGNW,
     &    KEEP, KEEP8,
     &    STEP, PROCNODE_STEPS, MYID_NODES, SLAVEF,
     &    SEUIL, AVOID_DELAYED, ETATASS,
     &    DKEEP(1), PIVNUL_LIST_STRUCT, IWPOS 
     &    , LRGROUPS
     &    , PERM
     &    )
      ELSE
          IW( IOLDPS + 4 + KEEP(IXSZ) ) = 1
          CALL DMUMPS_FAC1_LDLT( N, INODE,
     &    IW, LIW, A,
     &    LA,
     &    IOLDPS, POSELT,
     &    INFO_P(1), INFO_P(2), UU, NOFFNEGW, NULLNEGW, NPVW,
     &    NB22T1W, NBTINYW, DET_EXPW, DET_MANTW, DET_SIGNW,
     &    KEEP, KEEP8, MYID_NODES, SEUIL,
     &    AVOID_DELAYED,
     &    ETATASS,
     &    DKEEP(1), PIVNUL_LIST_STRUCT, IWPOS
     &    , LRGROUPS
     &    , PERM
     &    )
          IW(IOLDPS + 4 + KEEP(IXSZ)) = STEP(INODE)
      ENDIF
      IF (INFO_P(1) .LT. 0) THEN
        RETURN
      ENDIF
      CALL DMUMPS_FAC_STACK( COMM_LOAD, ASS_IRECV, N, INODE, TYPE,
     &TYPEF, LA, IW, LIW, A,
     &INFO_P(1), INFO_P(2), OPELIW, NELVAW, NMAXNPIVW, PTRIST, PTLUST_S,
     &PTRFAC, PTRAST, STEP, PIMASTER, PAMASTER, NE, POSFAC, 
     &LRLU, LRLUS, LRLUSM, IPTRLU, ICNTL, KEEP,
     &KEEP8, DKEEP,
     &COMP,IWPOS, IWPOSCB, PROCNODE_STEPS,
     &SLAVEF, IFATH, COMM_NODES, MYID_NODES, IPOOL_P,
     &LPOOL_P, LEAF, NSTK_STEPS, PERM, BUFR, LBUFR,
     &LBUFR_BYTES, NBFIN, root, OPASSW, ITLOC, RHS_MUMPS,
     &FILS, DAD, PTRARW, PTRAIW,
     &PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &INTARR, DBLARR, ND, FRERE,
     &LPTRAR, NELT, FRTPTR, FRTELT, ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     & LRGROUPS,
     & FLOP_ESTIM_ACC
     &)
      RETURN
      END SUBROUTINE DMUMPS_PROCESS_FRONT_NIV1
