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
      MODULE SMUMPS_SOL_L0OMP_M
!$    USE OMP_LIB, ONLY: OMP_LOCK_KIND
      INTEGER, PARAMETER :: NB_LOCK_MAX = 18
!$    INTEGER(OMP_LOCK_KIND),
!$   &ALLOCATABLE, DIMENSION(:), SAVE :: LOCK_FOR_SCATTER
      CONTAINS
      SUBROUTINE SMUMPS_SOL_L0OMP_LI( K400 )
!$    USE OMP_LIB, ONLY: OMP_INIT_LOCK
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K400
!$    INTEGER :: I
!$    IF (K400 .GT. 0) THEN
!$      ALLOCATE(LOCK_FOR_SCATTER(min(NB_LOCK_MAX,K400)))
!$      DO I = 1, min(NB_LOCK_MAX,K400)
!$        CALL OMP_INIT_LOCK(LOCK_FOR_SCATTER(I))
!$      ENDDO
!$    ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SOL_L0OMP_LI
      SUBROUTINE SMUMPS_SOL_L0OMP_LD( K400 )
!$    USE OMP_LIB, ONLY : OMP_DESTROY_LOCK
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: K400
!$    INTEGER :: I
!$    IF (allocated(LOCK_FOR_SCATTER)) THEN
!$      IF (K400 .GT. 0) THEN
!$        DO I = 1, min(NB_LOCK_MAX,K400)
!$          CALL OMP_DESTROY_LOCK(LOCK_FOR_SCATTER(I))
!$        ENDDO
!$        DEALLOCATE(LOCK_FOR_SCATTER)
!$      ENDIF
!$    ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SOL_L0OMP_LD
      SUBROUTINE SMUMPS_SOL_L0OMP_R(N, MTYPE,
     &  NRHS, LIW, IW, PTRICB, RHSINTR, LRHSINTR, POSINRHSINTR_FWD,
     &  STEP, FRERE, DAD, FILS, NSTK, PTRIST, PTRFAC, INFO,
     &  KEEP, KEEP8, DKEEP, PROCNODE_STEPS, SLAVEF,
     &  COMM, MYID,              
     &  BUFR, LBUFR, LBUFR_BYTES,
     &  RHS_ROOT, LRHS_ROOT,             
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     &  RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE,
     &  FROM_PP, 
     &  NBROOT_UNDER_L0, LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &  L_VIRT_L0_OMP, VIRT_L0_OMP, L_PHYS_L0_OMP, PHYS_L0_OMP,
     &  PERM_L0_OMP, PTR_LEAFS_L0_OMP,
     &  L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &  L0_OMP_FACTORS, LL0_OMP_FACTORS,
     &  DO_PRUN, TO_PROCESS )
      USE SMUMPS_STRUC_DEF, ONLY : SMUMPS_L0OMPFAC_T
!$      USE OMP_LIB
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: N, MTYPE, NRHS, SLAVEF, LIW
      INTEGER, INTENT( in ) :: IW(LIW)
      INTEGER               :: INFO( 80 ), KEEP(500)
      INTEGER(8)            :: KEEP8(150)
      REAL                  :: DKEEP(230)
      INTEGER, INTENT( in ) :: PROCNODE_STEPS( KEEP(28) )
      INTEGER               :: PTRICB( KEEP(28) ) 
      INTEGER, INTENT( in ) :: POSINRHSINTR_FWD(N), LRHSINTR 
      REAL, INTENT(inout):: RHSINTR(LRHSINTR,NRHS)
      INTEGER, INTENT( in ) :: STEP(N), FRERE( KEEP(28) ), FILS( N ),
     &                         DAD( KEEP(28) )
      INTEGER, INTENT( inout ) :: NSTK(KEEP(28))
      INTEGER, INTENT( in )    :: PTRIST(KEEP(28))
      INTEGER(8), INTENT( in ) :: PTRFAC(KEEP(28))
      INTEGER,    INTENT( IN ) :: COMM, MYID
      INTEGER,    INTENT( IN ) :: LBUFR, LBUFR_BYTES
      INTEGER                  :: BUFR(LBUFR)
      INTEGER(8), INTENT(IN)   :: LRHS_ROOT
      REAL                  :: RHS_ROOT(LRHS_ROOT)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      LOGICAL, INTENT( in ) :: DO_NBSPARSE
      INTEGER, INTENT( in ) :: LRHS_BOUNDS
      INTEGER, INTENT( in ) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, INTENT( in ) :: FROM_PP
      INTEGER, INTENT( out ):: NBROOT_UNDER_L0
      INTEGER, INTENT( in ) :: LPOOL_B_L0_OMP
      INTEGER, INTENT( in ) :: IPOOL_B_L0_OMP
     &                         ( LPOOL_B_L0_OMP )
      INTEGER, INTENT( in ) :: L_PHYS_L0_OMP
      INTEGER, INTENT( in ) :: PHYS_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT( in ) :: L_VIRT_L0_OMP
      INTEGER, INTENT( in ) :: VIRT_L0_OMP( L_VIRT_L0_OMP )
      INTEGER, INTENT( in ) :: PERM_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT( in ) :: PTR_LEAFS_L0_OMP( L_PHYS_L0_OMP + 1)
      INTEGER, INTENT( in ) :: LL0_OMP_MAPPING
      INTEGER, INTENT( in ) :: L0_OMP_MAPPING( LL0_OMP_MAPPING )
      INTEGER, INTENT( in ) :: LL0_OMP_FACTORS
      LOGICAL, INTENT( in ) :: DO_PRUN
      LOGICAL, INTENT( in ) :: TO_PROCESS( KEEP(28) )
      TYPE (SMUMPS_L0OMPFAC_T), INTENT(IN) ::
     &                        L0_OMP_FACTORS(LL0_OMP_FACTORS)
      INTEGER :: LASTFSSBTRSTA_P, LASTFSSBTRDYN_P
      INTEGER :: THREAD_ID, IL0OMPFAC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOOL_P
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWCB_P
      REAL, ALLOCATABLE, DIMENSION(:) :: WCB_P
      INTEGER :: LPOOL_P, LEAF_P, LIWCB_P
      INTEGER(8) :: LWCB_P
      INTEGER(8) :: POSWCB_P, PLEFTWCB_P
      INTEGER    :: POSIWCB_P
      LOGICAL :: IS_INODE_PROCESSED_P
      LOGICAL :: ERROR_WAS_BROADCASTED_P
      INTEGER :: INFO_P(2), allocok
      INTEGER :: I, VIRTUAL_TASK, PHYSICAL_TASK
      INTEGER :: INODE, IFATH, IROOT_SBTR
      INTEGER :: NBROOT_PROCESSED 
      INTEGER :: NEXT_TASK_DYN 
!$    INTEGER :: NOMP_SAVE
      INTEGER :: NBFIN_DUMMY
      NBFIN_DUMMY = huge(NBFIN_DUMMY)
      NBROOT_PROCESSED = 0
      PTRICB = 0
      NEXT_TASK_DYN = KEEP(400)+1
!$OMP PARALLEL
!$OMP& SHARED ( NEXT_TASK_DYN, IPOOL_B_L0_OMP,
!$OMP&          LPOOL_B_L0_OMP, NBFIN_DUMMY )
!$OMP& PRIVATE ( THREAD_ID, IL0OMPFAC, VIRTUAL_TASK, PHYSICAL_TASK,
!$OMP&           IPOOL_P, LPOOL_P, LEAF_P,
!$OMP&           LIWCB_P, LWCB_P, IWCB_P, WCB_P,
!$OMP&           PLEFTWCB_P, POSWCB_P, POSIWCB_P,
!$OMP&           LASTFSSBTRSTA_P, LASTFSSBTRDYN_P,
!$OMP&           INODE, IROOT_SBTR, IFATH,
!$OMP&           IS_INODE_PROCESSED_P,
!$OMP&           INFO_P, ERROR_WAS_BROADCASTED_P, NOMP_SAVE, allocok )
!$OMP& REDUCTION( + : NBROOT_PROCESSED )
!$    NOMP_SAVE = omp_get_max_threads()
      THREAD_ID = 1
!$    THREAD_ID = OMP_GET_THREAD_NUM() + 1
!$OMP BARRIER
!$    CALL omp_set_num_threads(1)
      LPOOL_P = LPOOL_B_L0_OMP
      INFO_P(1:2) = 0
      LWCB_P  = int(KEEP(133),8)*int(NRHS,8)
      LIWCB_P = KEEP(133)
      PLEFTWCB_P = 1_8
      POSWCB_P = LWCB_P
      POSIWCB_P = LIWCB_P
      ALLOCATE(IPOOL_P(LPOOL_P), IWCB_P(LIWCB_P), WCB_P( LWCB_P),
     &         stat=allocok)
      IF ( allocok > 0 ) THEN
        INFO_P(1) = -13
        CALL MUMPS_SETI8TOI4(LPOOL_P + LIWCB_P + LWCB_P,
     &                       INFO(2))
!$OMP CRITICAL(critical_info)
        INFO(1) = -13
        INFO(2) = INFO_P(2)
!$OMP END CRITICAL(critical_info)
      ENDIF
!$OMP BARRIER
      IF (INFO(1) .LT. 0) THEN
        GOTO 50
      ENDIF
      VIRTUAL_TASK = THREAD_ID
 600  CONTINUE
      IF (VIRTUAL_TASK .LT. L_VIRT_L0_OMP) THEN
        DO PHYSICAL_TASK = VIRT_L0_OMP( VIRTUAL_TASK ),
     &       VIRT_L0_OMP ( VIRTUAL_TASK + 1 ) - 1
          LEAF_P = 1
          DO I = PTR_LEAFS_L0_OMP( PERM_L0_OMP( PHYSICAL_TASK )+1 )+1,
     &      PTR_LEAFS_L0_OMP( PERM_L0_OMP( PHYSICAL_TASK ) )
            IF ( IPOOL_B_L0_OMP(I) .GT. 0 ) THEN
              IPOOL_P(LEAF_P) = IPOOL_B_L0_OMP(I)
              LEAF_P = LEAF_P + 1
            ENDIF
          ENDDO
          IF ( LEAF_P .EQ. 1 ) THEN
                WRITE(*,*) " Internal error 1 in SMUMPS_SOL_L0OMP_R",
     &             LEAF_P
          ENDIF
          IROOT_SBTR = PHYS_L0_OMP( PERM_L0_OMP( PHYSICAL_TASK ))
          IF (DO_PRUN) THEN
            IF (.NOT. TO_PROCESS(STEP(IROOT_SBTR))) THEN
              CYCLE
            ENDIF
          ENDIF
          INODE = IROOT_SBTR
          DO WHILE (INODE .GT. 0)
            LASTFSSBTRSTA_P = INODE
            INODE=FILS(INODE)
          ENDDO
          CALL MUMPS_COMPUTE_LASTFS_DYN( IROOT_SBTR, LASTFSSBTRDYN_P,
     &    MTYPE, KEEP, IW, LIW, N, STEP, PTRIST, FILS, FRERE )
          DO WHILE (LEAF_P .NE.1 .AND. INFO_P(1) .GE. 0)
            LEAF_P = LEAF_P - 1
            INODE = IPOOL_P(LEAF_P)
            IFATH = DAD(STEP(INODE) )
            IL0OMPFAC = L0_OMP_MAPPING(STEP(INODE))
            IF (IL0OMPFAC .NE. THREAD_ID) THEN
            ENDIF
            IF (DO_PRUN) THEN
              IS_INODE_PROCESSED_P = TO_PROCESS(STEP(INODE))
            ELSE
              IS_INODE_PROCESSED_P = .TRUE.
            ENDIF
            IF ( IS_INODE_PROCESSED_P ) THEN
              CALL SMUMPS_SOLVE_NODE_FWD( INODE,
     &        LASTFSSBTRSTA_P, LASTFSSBTRDYN_P,
     &        BUFR, LBUFR, LBUFR_BYTES, MYID, SLAVEF, COMM,   
     &        N, IPOOL_P, LPOOL_P, LEAF_P,  NBFIN_DUMMY, NSTK,
     &        IWCB_P, LIWCB_P, WCB_P, LWCB_P,
     &        L0_OMP_FACTORS(IL0OMPFAC)%A(1),
     &        L0_OMP_FACTORS(IL0OMPFAC)%LA,
     &        IW, LIW, 
     &        NRHS, POSWCB_P, PLEFTWCB_P, POSIWCB_P, 
     &        PTRICB, PTRIST, PTRFAC, PROCNODE_STEPS,
     &        FILS, STEP, FRERE, DAD, INFO_P, KEEP, KEEP8, DKEEP,
     &        RHS_ROOT, LRHS_ROOT, MTYPE,
     &        RHSINTR, LRHSINTR, POSINRHSINTR_FWD,
     &        ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     &        RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &        , ERROR_WAS_BROADCASTED_P )
              IF (INFO_P(1) .LT. 0) THEN
!$OMP CRITICAL(critical_info)
                INFO(1) = INFO_P(1)
                INFO(2) = INFO_P(2)
!$OMP END CRITICAL(critical_info)
              ENDIF
              IF ( INFO(1) .LT. 0 ) GOTO 50
              IF (ERROR_WAS_BROADCASTED_P) THEN
                WRITE(*,*) " Internal error 2 in SMUMPS_SOL_L0OMP_R",
     &          ERROR_WAS_BROADCASTED_P
              ENDIF
            ENDIF
            IF ( IFATH .EQ. 0 ) THEN
              IF ( IS_INODE_PROCESSED_P ) THEN
                NBROOT_PROCESSED = NBROOT_PROCESSED + 1
              ENDIF
            ELSE
              PTRICB(STEP(INODE)) = 0
              IF (IFATH .NE. 0) THEN
                IF ( INODE .NE. IROOT_SBTR ) THEN
                  IF ( IS_INODE_PROCESSED_P ) THEN
                    NSTK(STEP(IFATH)) = NSTK(STEP(IFATH)) - 1
                  ENDIF
                  IF (NSTK(STEP(IFATH)) .EQ. 0 .OR.
     &                NSTK(STEP(IFATH)) .EQ. -1 ) THEN
                    IPOOL_P( LEAF_P ) = IFATH
                    LEAF_P = LEAF_P + 1
                    IF (DO_PRUN) THEN
                      NSTK(STEP(IFATH)) = huge(NSTK(STEP(IFATH)))
                    ENDIF
                  ENDIF
                ELSE
                  IF ( IS_INODE_PROCESSED_P ) THEN
!$OMP ATOMIC UPDATE
                    NSTK(STEP(IFATH)) = NSTK(STEP(IFATH)) - 1
!$OMP END ATOMIC
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!$OMP ATOMIC CAPTURE
        VIRTUAL_TASK = NEXT_TASK_DYN
        NEXT_TASK_DYN = NEXT_TASK_DYN + 1
!$OMP END ATOMIC
        GOTO 600
      ENDIF
  50  CONTINUE
      IF (allocated(IPOOL_P)) DEALLOCATE(IPOOL_P)
      IF (allocated(IWCB_P)) DEALLOCATE(IWCB_P)
      IF (allocated(WCB_P)) DEALLOCATE(WCB_P)
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$    CALL omp_set_num_threads(int(NOMP_SAVE,4))
#else
!$    CALL omp_set_num_threads(NOMP_SAVE)
#endif
!$OMP END PARALLEL
      NBROOT_UNDER_L0 = NBROOT_PROCESSED
      RETURN
      END SUBROUTINE SMUMPS_SOL_L0OMP_R
      SUBROUTINE SMUMPS_SOL_L0OMP_S(N, MTYPE, NRHS, LIW, IW,
     &  PTRICB, PTRACB, RHSINTR, LRHSINTR, POSINRHSINTR_BWD,
     &  STEP, FRERE, FILS, NE_STEPS, PTRIST, PTRFAC, INFO,
     &  KEEP, KEEP8, DKEEP, PROCNODE_STEPS, SLAVEF,
     &  COMM, MYID, BUFR, LBUFR, LBUFR_BYTES, RHS_ROOT, LRHS_ROOT, 
     &  ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS,
     &  PRUN_BELOW_BWD, TO_PROCESS, SIZE_TO_PROCESS,
     &  RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP, LPOOL_B_L0_OMP,
     &  L_VIRT_L0_OMP, VIRT_L0_OMP, L_PHYS_L0_OMP, PHYS_L0_OMP,
     &  PERM_L0_OMP, PTR_LEAFS_L0_OMP, L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &  L0_OMP_FACTORS, LL0_OMP_FACTORS )
      USE SMUMPS_STRUC_DEF, ONLY : SMUMPS_L0OMPFAC_T
      USE OMP_LIB
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: N, MTYPE, NRHS, SLAVEF, LIW
      INTEGER, INTENT( in ) :: IW(LIW)
      INTEGER               :: INFO( 80 ), KEEP(500)
      INTEGER(8)            :: KEEP8(150)
      REAL                  :: DKEEP(230)
      INTEGER, INTENT( in ) :: PROCNODE_STEPS( KEEP(28) )
      INTEGER               :: PTRICB( KEEP(28) ) 
      INTEGER(8)            :: PTRACB( KEEP(28) ) 
      INTEGER, INTENT( in ) :: POSINRHSINTR_BWD(N), LRHSINTR 
      REAL, INTENT(inout):: RHSINTR(LRHSINTR,NRHS)
      INTEGER, INTENT( in ) :: STEP(N), FRERE( KEEP(28) ), FILS( N )
      INTEGER, INTENT( inout ) :: NE_STEPS(KEEP(28))
      INTEGER, INTENT( in )    :: PTRIST(KEEP(28))
      INTEGER(8), INTENT( in ) :: PTRFAC(KEEP(28))
      INTEGER,    INTENT( IN ) :: COMM, MYID
      INTEGER,    INTENT( IN ) :: LBUFR, LBUFR_BYTES
      INTEGER                  :: BUFR(LBUFR)
      INTEGER(8), INTENT(IN)   :: LRHS_ROOT
      REAL                  :: RHS_ROOT(LRHS_ROOT)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER    :: LPANEL_POS
      INTEGER    :: PANEL_POS(LPANEL_POS)
      LOGICAL, INTENT( in ) :: DO_NBSPARSE
      INTEGER, INTENT( in ) :: LRHS_BOUNDS
      INTEGER, INTENT( in ) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, INTENT( in ) :: PRUN_BELOW_BWD
      INTEGER, INTENT( in ) :: SIZE_TO_PROCESS
      LOGICAL, INTENT( in ) :: TO_PROCESS(SIZE_TO_PROCESS)
      LOGICAL, INTENT( in ) :: FROM_PP
      INTEGER, INTENT( in ) :: LPOOL_B_L0_OMP
      INTEGER, INTENT( in ) :: L_PHYS_L0_OMP
      INTEGER, INTENT( in ) :: PHYS_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT( in ) :: L_VIRT_L0_OMP
      INTEGER, INTENT( in ) :: VIRT_L0_OMP( L_VIRT_L0_OMP )
      INTEGER, INTENT( in ) :: PERM_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT( in ) :: PTR_LEAFS_L0_OMP( L_PHYS_L0_OMP + 1)
      INTEGER, INTENT( in ) :: LL0_OMP_MAPPING
      INTEGER, INTENT( in ) :: L0_OMP_MAPPING( LL0_OMP_MAPPING )
      INTEGER, INTENT( in ) :: LL0_OMP_FACTORS
      TYPE (SMUMPS_L0OMPFAC_T), INTENT(IN) ::
     &                        L0_OMP_FACTORS(LL0_OMP_FACTORS)
      INTEGER :: THREAD_ID, IL0OMPFAC
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOOL_P
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IWCB_P
      REAL, ALLOCATABLE, DIMENSION(:) :: WCB_P
      REAL, ALLOCATABLE, DIMENSION(:) :: W2_P
      INTEGER, ALLOCATABLE, DIMENSION(:) :: PANEL_POS_P
      INTEGER :: LPOOL_P, IIPOOL_P, LIWCB_P, LPANEL_POS_P
      INTEGER :: MYLEAF_LEFT_HUGE_P 
      INTEGER(8) :: LWCB_P
      INTEGER(8) :: POSWCB_P, PLEFTWCB_P
      INTEGER    :: POSIWCB_P
      LOGICAL :: DO_MCAST2_TERMBWD_P
      LOGICAL :: ERROR_WAS_BROADCASTED_P
      INTEGER :: INFO_P(2), allocok
      INTEGER :: VIRTUAL_TASK, PHYSICAL_TASK
      INTEGER :: INODE
      INTEGER :: NEXT_TASK_DYN 
!$    INTEGER :: NOMP_SAVE
      INTEGER :: NBFIN_DUMMY
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: DEJA_SEND_DUMMY
      NBFIN_DUMMY = huge(NBFIN_DUMMY)
      ALLOCATE(DEJA_SEND_DUMMY( 0:SLAVEF-1 ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of DEJA_SEND_DUMMY in '
     &        //'routine SMUMPS_SOL_S '
         INFO(1)=-13
         INFO(2)=SLAVEF
         GOTO 100
      endif
      PTRICB = 0
      NEXT_TASK_DYN = KEEP(400)+1
!$OMP PARALLEL
!$OMP& SHARED ( NEXT_TASK_DYN, LPOOL_B_L0_OMP,
!$OMP&          NBFIN_DUMMY, DEJA_SEND_DUMMY )
!$OMP& PRIVATE ( THREAD_ID, IL0OMPFAC, VIRTUAL_TASK, PHYSICAL_TASK,
!$OMP&           IPOOL_P, LPOOL_P, IIPOOL_P, MYLEAF_LEFT_HUGE_P,
!$OMP&           LIWCB_P, LWCB_P, IWCB_P, WCB_P, W2_P, LPANEL_POS_P,
!$OMP&           PANEL_POS_P,
!$OMP&           PLEFTWCB_P, POSWCB_P, POSIWCB_P,
!$OMP&           INODE,
!$OMP&           INFO_P, DO_MCAST2_TERMBWD_P,
!$OMP&           ERROR_WAS_BROADCASTED_P, NOMP_SAVE, allocok )
!$    NOMP_SAVE = omp_get_max_threads()
      THREAD_ID = 1
!$    THREAD_ID = OMP_GET_THREAD_NUM() + 1
!$OMP BARRIER
!$    CALL omp_set_num_threads(1)
      LPOOL_P = LPOOL_B_L0_OMP
      INFO_P(1:2) = 0
      LWCB_P  = int(KEEP(133),8)*int(NRHS,8)
      LIWCB_P = KEEP(133)
      PLEFTWCB_P = 1_8
      POSWCB_P = LWCB_P
      POSIWCB_P = LIWCB_P
      IF (KEEP(201).EQ.1) THEN
        LPANEL_POS_P = KEEP(228)+1 
        CALL MUMPS_ABORT()
      ELSE
        LPANEL_POS_P = 1
      ENDIF
      ALLOCATE(IPOOL_P(LPOOL_P), IWCB_P(LIWCB_P), WCB_P( LWCB_P),
     &         W2_P(KEEP(133)), PANEL_POS_P(LPANEL_POS_P), stat=allocok)
      IF ( allocok > 0 ) THEN
        INFO_P(1) = -13
        CALL MUMPS_SETI8TOI4(LPOOL_P + LIWCB_P + LWCB_P +
     &                       KEEP(133)+LPANEL_POS_P, INFO(2))
!$OMP CRITICAL(critical_info)
        INFO(1) = -13
        INFO(2) = INFO_P(2)
!$OMP END CRITICAL(critical_info)
      ENDIF
!$OMP BARRIER
      IF (INFO(1) .LT. 0) THEN
        GOTO 50
      ENDIF
      VIRTUAL_TASK = THREAD_ID
 600  CONTINUE
      IF (VIRTUAL_TASK .LT. L_VIRT_L0_OMP) THEN
        DO PHYSICAL_TASK = VIRT_L0_OMP( VIRTUAL_TASK ),
     &       VIRT_L0_OMP ( VIRTUAL_TASK + 1 ) - 1
          INODE = PHYS_L0_OMP( PERM_L0_OMP( PHYSICAL_TASK ) )
          IPOOL_P(1) = INODE
          IIPOOL_P = 2
          MYLEAF_LEFT_HUGE_P = huge(MYLEAF_LEFT_HUGE_P)
          IF ( PRUN_BELOW_BWD ) THEN
            IF ( .NOT. TO_PROCESS(STEP(INODE)) ) THEN
            CYCLE
            ENDIF
          ENDIF
          DO WHILE (IIPOOL_P .NE.1 .AND. INFO_P(1) .GE. 0)
            IIPOOL_P = IIPOOL_P - 1
            INODE = IPOOL_P(IIPOOL_P)
            IL0OMPFAC = L0_OMP_MAPPING(STEP(INODE))
            IF (IL0OMPFAC .NE. THREAD_ID) THEN
            ENDIF
            CALL SMUMPS_SOLVE_NODE_BWD( INODE, N, IPOOL_P, LPOOL_P,
     &      IIPOOL_P,  NBFIN_DUMMY, L0_OMP_FACTORS(IL0OMPFAC)%A(1),
     &      L0_OMP_FACTORS(IL0OMPFAC)%LA, IW, LIW, 
     &      WCB_P, LWCB_P, NRHS, POSWCB_P, PLEFTWCB_P, POSIWCB_P,
     &      RHSINTR, LRHSINTR, POSINRHSINTR_BWD,
     &      PTRICB, PTRACB, IWCB_P, LIWCB_P, W2_P, NE_STEPS, STEP,
     &      FRERE, FILS, PTRIST, PTRFAC, MYLEAF_LEFT_HUGE_P, INFO_P, 
     &      PROCNODE_STEPS,
     &      DEJA_SEND_DUMMY,      
     &      SLAVEF, COMM, MYID, BUFR, LBUFR, LBUFR_BYTES, 
     &      KEEP, KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE,
     &      ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS_P, LPANEL_POS_P,
     &      PRUN_BELOW_BWD, TO_PROCESS, SIZE_TO_PROCESS,
     &      RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &          , ERROR_WAS_BROADCASTED_P
     &          , DO_MCAST2_TERMBWD_P
     &    )
            IF (INFO_P(1) .LT. 0) THEN
!$OMP CRITICAL(critical_info)
              INFO(1) = INFO_P(1)
              INFO(2) = INFO_P(2)
!$OMP END CRITICAL(critical_info)
            ENDIF
            IF ( INFO(1) .LT. 0 ) GOTO 50
            IF (ERROR_WAS_BROADCASTED_P) THEN
              WRITE(*,*) " Internal error 1 in SMUMPS_SOL_L0OMP_R",
     &        ERROR_WAS_BROADCASTED_P
            ENDIF
            IF (DO_MCAST2_TERMBWD_P) THEN
              WRITE(*,*) " Internal error 2 in SMUMPS_SOL_L0OMP_R",
     &        DO_MCAST2_TERMBWD_P
            ENDIF
          ENDDO
        ENDDO
!$OMP ATOMIC CAPTURE
        VIRTUAL_TASK = NEXT_TASK_DYN
        NEXT_TASK_DYN = NEXT_TASK_DYN + 1
!$OMP END ATOMIC
        GOTO 600
      ENDIF
  50  CONTINUE
      IF (allocated(IPOOL_P)) DEALLOCATE(IPOOL_P)
      IF (allocated(IWCB_P)) DEALLOCATE(IWCB_P)
      IF (allocated(WCB_P)) DEALLOCATE(WCB_P)
      IF (allocated(W2_P)) DEALLOCATE(W2_P)
      IF (allocated(PANEL_POS_P)) DEALLOCATE(PANEL_POS_P)
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$    CALL omp_set_num_threads(int(NOMP_SAVE,4))
#else
!$    CALL omp_set_num_threads(NOMP_SAVE)
#endif
!$OMP END PARALLEL
 100  CONTINUE
      IF (allocated(DEJA_SEND_DUMMY)) DEALLOCATE(DEJA_SEND_DUMMY)
      RETURN
      END SUBROUTINE SMUMPS_SOL_L0OMP_S
      END MODULE SMUMPS_SOL_L0OMP_M
