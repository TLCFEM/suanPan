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
      MODULE MUMPS_ANA_OMP_M
      CONTAINS
      SUBROUTINE MUMPS_ANA_L0_OMP( NB_THREADS, N, NSTEPS, SYM, SLAVEF,
     &  DAD, FRERE, FILS, NSTK_STEPS, ND, STEP, PROCNODE_STEPS, KEEP,
     &  KEEP8, MYID_NODES, NA, LNA, ARITH, LPOOL_B_L0_OMP,
     &  IPOOL_B_L0_OMP, LPOOL_A_L0_OMP, IPOOL_A_L0_OMP,
     &  L_VIRT_L0_OMP, VIRT_L0_OMP, VIRT_L0_OMP_MAPPING,
     &  L_PHYS_L0_OMP, PHYS_L0_OMP, PERM_L0_OMP, PTR_LEAFS_L0_OMP,
     & INFO, ICNTL )
      USE MUMPS_IDLL
      USE MUMPS_DDLL
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT ( IN ) :: NB_THREADS, N, NSTEPS, SYM
      INTEGER, INTENT ( IN ) :: SLAVEF, MYID_NODES
      INTEGER, INTENT ( IN ) :: LNA
      INTEGER, INTENT ( IN ) :: DAD (:), FRERE (:) 
      INTEGER, INTENT ( IN ) :: FILS (:)  
      INTEGER, INTENT ( IN ) :: NSTK_STEPS (:) 
      INTEGER, INTENT ( IN ) :: ND (:), STEP (:) 
      INTEGER, INTENT ( IN ) :: PROCNODE_STEPS(:) 
      INTEGER, INTENT ( IN ) :: KEEP ( : ) 
      INTEGER(8), INTENT ( IN ) :: KEEP8(:)  
      INTEGER, INTENT ( IN ) :: NA ( : ) 
      CHARACTER(1), INTENT(IN) :: ARITH 
      INTEGER, INTENT ( OUT ) :: LPOOL_B_L0_OMP
      INTEGER, INTENT ( OUT ) :: LPOOL_A_L0_OMP
      INTEGER, INTENT ( OUT ) :: L_PHYS_L0_OMP
      INTEGER, INTENT ( OUT ) :: L_VIRT_L0_OMP
      INTEGER, DIMENSION(:), POINTER :: IPOOL_B_L0_OMP
      INTEGER, DIMENSION(:), POINTER :: IPOOL_A_L0_OMP
      INTEGER, DIMENSION(:), POINTER :: PHYS_L0_OMP
      INTEGER, DIMENSION(:), POINTER :: VIRT_L0_OMP, VIRT_L0_OMP_MAPPING
      INTEGER, DIMENSION(:), POINTER :: PERM_L0_OMP
      INTEGER, DIMENSION(:), POINTER :: PTR_LEAFS_L0_OMP
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      LOGICAL :: LPOK
      INTEGER :: LP
      INTEGER :: NB_REPEAT_ACCEPTL0, NB_MAX_IN_L0_ACCEPTL0,
     &           NB_IN_L0
      DOUBLE PRECISION :: SMALL_COST
      INTEGER :: THRESH_MEM, SLAVEF_DURING_MAPPING
      REAL    :: THRESH_EQUILIB
      DOUBLE PRECISION, DIMENSION(1,1,1) :: BENCH 
      INTEGER :: INODE
      INTEGER :: NBLEAF_MYID
      DOUBLE PRECISION :: COST_UNDER, COST_ABOVE, COST_TOTAL_BEST
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: THREADS_CHARGE
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COSTS_MONO_THREAD
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COSTS_MULTI_THREAD
      INTEGER         , DIMENSION(:), ALLOCATABLE :: IPOOL_B_INV
      INTEGER(8) :: FACTOR_SIZE_UNDER_L0, FACTOR_SIZE_PER_MPI
      INTEGER, DIMENSION(:), ALLOCATABLE :: CP_NSTK_STEPS
      TYPE ( IDLL_T ), POINTER :: L0_OMP_DLL
      TYPE ( IDLL_T ), POINTER :: LEAFS_ABOVE_L0_OMP_DLL
      INTEGER :: I
      THRESH_EQUILIB = real(KEEP(408))/real(100)
      IF ((THRESH_EQUILIB.GT..99).OR.(THRESH_EQUILIB.LT.0.01)) THEN
         THRESH_EQUILIB = 0.9
      ENDIF
      THRESH_MEM = KEEP(397)
      IF ((THRESH_MEM.LT.-1).OR.(THRESH_MEM.GT.100)) THRESH_MEM=100
      IF (THRESH_MEM.EQ.-1) THEN
        IF (NB_THREADS.EQ.2) THEN
         THRESH_MEM = 50
        ELSEIF (NB_THREADS.LE.4) THEN
         THRESH_MEM = 60
        ELSEIF (NB_THREADS.LE.8) THEN
         THRESH_MEM = 70
        ELSEIF (NB_THREADS.LE.12) THEN
         THRESH_MEM = 80
        ELSEIF (NB_THREADS.LE.20) THEN
         THRESH_MEM = 85
        ELSEIF (NB_THREADS.LE.36) THEN
         THRESH_MEM = 90
        ELSE
         THRESH_MEM = 95
        ENDIF
      ENDIF
        SLAVEF_DURING_MAPPING = SLAVEF
      FACTOR_SIZE_PER_MPI = KEEP8(101) / SLAVEF_DURING_MAPPING
      IF ( KEEP(261) .EQ. 0) THEN
          WRITE(*,*)"KEEP(261) MUST BE SET TO 1 IN ORDER TO USE
     & MULTITHREADED TREE PARALLELISM"
          CALL MUMPS_ABORT()
      END IF
      LP = ICNTL(1)
      LPOK = ( LP .GT. 0 .AND. ICNTL(4) .GE.1 )
      NB_REPEAT_ACCEPTL0 = -1
      NB_MAX_IN_L0_ACCEPTL0 = -1
      CALL MUMPS_ANA_INITIALIZE_L0_OMP ()
      IF (INFO(1) .LT. 0) GOTO 500
      DO WHILE ( .NOT. MUMPS_ANA_ACCEPT_L0_OMP () )
        IF (INFO(1) .LT. 0) GOTO 500
        CALL L0_REMOVE_NODE ( INODE )
          IF (INODE .LT. 0) THEN 
            IPOOL_B_L0_OMP(IPOOL_B_INV(STEP(-INODE))) = INODE
          ELSE
            CALL L0_INSERT_CHILDREN ( INODE )
          END IF
      END DO
      CALL MUMPS_ANA_FINALIZE_L0_OMP ()
 500  CONTINUE
      CALL MUMPS_ANA_FREE_L0_WORKSPACE()
      RETURN
      CONTAINS
      SUBROUTINE MUMPS_ANA_INITIALIZE_L0_OMP ( )
      IMPLICIT NONE
      INTEGER :: INODE, IFATH, IGRANDFATH, SPECIAL_ROOT,
     &           NFRONT, NPIV, LEAF, VARNUM, IERR
      LOGICAL :: INODE_IS_A_LEAF
      INTEGER(8) :: NFRONT8, NPIV8
      DOUBLE PRECISION :: COST_NODE, SMALL_COST_TMP
      LOGICAL :: IN_L0INIT, SKIP_ABOVE
      LOGICAL, EXTERNAL :: MUMPS_ROOTSSARBR, MUMPS_IN_OR_ROOT_SSARBR
      INTEGER, EXTERNAL :: MUMPS_GET_POOL_LENGTH, MUMPS_TYPENODE
      IF (associated(IPOOL_B_L0_OMP)) THEN
        WRITE(*,*) " Internal error 1 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(IPOOL_A_L0_OMP)) THEN
        WRITE(*,*) " Internal error 2 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(VIRT_L0_OMP)) THEN
        WRITE(*,*) " Internal error 3 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(VIRT_L0_OMP_MAPPING)) THEN
        WRITE(*,*) " Internal error 4 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(PERM_L0_OMP)) THEN
        WRITE(*,*) " Internal error 5 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IF (associated(PTR_LEAFS_L0_OMP)) THEN
        WRITE(*,*) " Internal error 6 MUMPS_ANA_INITIALIZE_L0_OMP" 
        CALL MUMPS_ABORT()
      ENDIF
      IERR = IDLL_CREATE ( L0_OMP_DLL )
      NB_IN_L0 = 0
      IF (KEEP(72).eq.1) THEN
        SMALL_COST = 2.0D0
      ELSE
        SMALL_COST = 100.0D0
      ENDIF
      SMALL_COST_TMP = 0.0D0
      IERR = IDLL_CREATE ( LEAFS_ABOVE_L0_OMP_DLL )
      ALLOCATE( THREADS_CHARGE( NB_THREADS ), stat=IERR )
      IF (IERR .GT. 0) THEN
         INFO(1) = -7
         INFO(2) = NB_THREADS
         IF (LPOK) WRITE(LP,150) 'THREADS_CHARGE'
         GOTO 500
      ENDIF
      ALLOCATE( COSTS_MONO_THREAD ( NSTEPS ), stat=IERR )
      IF(IERR.GT.0) THEN
         INFO(1) = -7
         INFO(2) = NSTEPS
         IF (LPOK) WRITE(LP, 150) ' COSTS_MONO_THREAD'
         GOTO 500
      ENDIF
      IF (KEEP(403) .NE. 0) THEN
        ALLOCATE( COSTS_MULTI_THREAD ( NSTEPS ), stat=IERR )
        IF(IERR.GT.0) THEN
          INFO(1) = -7
          INFO(2) = NSTEPS
          IF (LPOK) WRITE(LP, 150) ' COSTS_MULTI_THREAD'
          GOTO 500
        ENDIF
      ENDIF
      ALLOCATE( IPOOL_B_INV ( NSTEPS ), stat=IERR )
      IF(IERR.GT.0) THEN
         INFO(1) = -7
         INFO(2) = NSTEPS
         IF (LPOK) WRITE(LP, 150) ' IPOOL_B_INV'
         GOTO 500
      ENDIF
      ALLOCATE( CP_NSTK_STEPS ( NSTEPS ), stat=IERR )
      IF(IERR.GT.0) THEN
         INFO(1) = -7
         INFO(2) = NSTEPS
         IF (LPOK) WRITE(LP, 150) ' CP_NSTK_STEPS'
         GOTO 500
      ENDIF
      LPOOL_B_L0_OMP=MUMPS_GET_POOL_LENGTH(NA(1),KEEP(1),KEEP8(1))
      ALLOCATE( IPOOL_B_L0_OMP( LPOOL_B_L0_OMP) , stat=IERR )
      IF(IERR.GT.0) THEN
         INFO(1) = -7
         INFO(2) = NSTEPS
         IF (LPOK) WRITE(LP, 150) ' id%IPOOL_B_L0_OMP'
         GOTO 500
      ENDIF
      COSTS_MONO_THREAD = 0.0D0
      IF (KEEP(403) .NE. 0) THEN
        COSTS_MULTI_THREAD = 0.0D0
        COST_UNDER = 0.0D0
        COST_ABOVE = 0.0D0
        COST_TOTAL_BEST = huge(COST_TOTAL_BEST)
      ENDIF
      FACTOR_SIZE_UNDER_L0 = 0_8
      CP_NSTK_STEPS(:) = NSTK_STEPS(:)
      IF (KEEP(403).NE.0) THEN
      CALL READ_BENCH( ARITH, KEEP(50) )
      ENDIF
      CALL MUMPS_INIT_POOL_DIST(N, LEAF,
     &                     MYID_NODES,
     &                     KEEP(199), NA(1), LNA, 
     &                     KEEP(1), KEEP8(1), STEP(1),
     &                     PROCNODE_STEPS(1),
     &                     IPOOL_B_L0_OMP(1), LPOOL_B_L0_OMP)
      DO I = 1, LEAF - 1
        IPOOL_B_INV(STEP(IPOOL_B_L0_OMP(I))) = I
      ENDDO
      LEAF = LEAF - 1
      NBLEAF_MYID = LEAF
      IF (NBLEAF_MYID .EQ. 0) THEN
        RETURN
      ENDIF
 90   CONTINUE
      INODE = IPOOL_B_L0_OMP ( LEAF )
      LEAF = LEAF - 1
      INODE_IS_A_LEAF=.TRUE.
 95   CONTINUE
      NFRONT = ND ( STEP ( INODE ) )
      NFRONT8= int(NFRONT,8)
      NPIV = 0
      VARNUM = INODE
      DO WHILE (VARNUM .GT. 0 )
          NPIV = NPIV + 1
          VARNUM = FILS ( VARNUM )
      END DO
      NPIV8=int(NPIV,8)
      VARNUM = - VARNUM
      IF (KEEP(403) .EQ. 0) THEN
        CALL MUMPS_GET_FLOPS_COST ( NFRONT, NPIV, NPIV,
     &                              SYM, 1, COST_NODE )
        COSTS_MONO_THREAD ( STEP ( INODE ) ) = COST_NODE
      ELSE
        CALL COST_BENCH (NPIV, NFRONT-NPIV, 1, KEEP(50), COST_NODE)
        COSTS_MONO_THREAD ( STEP ( INODE ) ) = COST_NODE
        CALL COST_BENCH (NPIV,NFRONT-NPIV,NB_THREADS,KEEP(50),COST_NODE)
        COSTS_MULTI_THREAD ( STEP ( INODE ) ) = COST_NODE
      END IF
      DO WHILE (VARNUM .GT. 0 )
          COSTS_MONO_THREAD ( STEP ( INODE ) ) =
     &                      COSTS_MONO_THREAD ( STEP ( INODE ) )
     &                      +
     &                      COSTS_MONO_THREAD ( STEP ( VARNUM ) )
          VARNUM = FRERE ( STEP ( VARNUM ) )
      END DO
      IFATH = DAD ( STEP ( INODE ) )
      IF (IFATH .NE. 0) THEN
        IGRANDFATH = DAD( STEP ( IFATH ) )
      ELSE
        IGRANDFATH = 0
      ENDIF
      SPECIAL_ROOT = max(KEEP(38), KEEP(20))
      SKIP_ABOVE = .FALSE.
      IN_L0INIT  = .FALSE.
      IF ( INODE .EQ. SPECIAL_ROOT ) THEN
        IN_L0INIT  = .FALSE.
        IF (INODE_IS_A_LEAF) THEN
          SKIP_ABOVE = .TRUE.
          GOTO 80
        ELSE
          WRITE(*,*) " Internal error 1 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &    INODE, SPECIAL_ROOT
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF ( IFATH .NE. 0 .AND. IFATH .EQ. KEEP(38) ) THEN
        IN_L0INIT  = .FALSE.
        IF (INODE_IS_A_LEAF) THEN
          SKIP_ABOVE = .TRUE.
          GOTO 80
        ELSE
          WRITE(*,*) " Internal error 2 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &    INODE, IFATH, KEEP(38)
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF ( SLAVEF_DURING_MAPPING > 1 ) THEN
        IF (MUMPS_ROOTSSARBR (
     &  PROCNODE_STEPS ( STEP ( INODE ) ), KEEP(199) )
     &  .OR. .NOT. MUMPS_IN_OR_ROOT_SSARBR (
     &  PROCNODE_STEPS ( STEP ( INODE ) ), KEEP(199) )
     &) THEN
          IN_L0INIT = .FALSE.
          IF (INODE_IS_A_LEAF) THEN
            SKIP_ABOVE = .TRUE.
            GOTO 80
          ELSE
            WRITE(*,*)
     &      " Internal error 3 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &      INODE
            CALL MUMPS_ABORT()
          ENDIF
        ENDIF
      ENDIF
      IF (IFATH.NE.0) THEN
        IF ( MUMPS_TYPENODE(STEP(IFATH),KEEP(199)).EQ.2) THEN
          IN_L0INIT = .FALSE.
          IF (INODE_IS_A_LEAF) THEN
            SKIP_ABOVE = .TRUE.
            GOTO 80
          ELSE
            WRITE(*,*)
     &      " Internal error 5 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &      INODE, IFATH
            CALL MUMPS_ABORT()
          ENDIF
        ENDIF
      ENDIF
      IF ( MUMPS_TYPENODE(STEP(INODE),KEEP(199)).EQ.2) THEN
        IN_L0INIT = .FALSE.
        IF (INODE_IS_A_LEAF) THEN
          SKIP_ABOVE = .TRUE.
          GOTO 80
        ELSE
          WRITE(*,*)
     &    " Internal error 6 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &    INODE
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      IF ( IFATH .EQ. 0 ) THEN
        IN_L0INIT = .TRUE.
        GOTO 80
      ELSE
        IF ( IFATH .EQ. KEEP(20) ) THEN
          IN_L0INIT = .TRUE.
          GOTO 80
        ENDIF
        IF ( IGRANDFATH .EQ. KEEP(38) .AND. KEEP(38) .NE. 0 ) THEN
          IN_L0INIT = .TRUE.
          GOTO 80
        ENDIF
        IF ( SLAVEF_DURING_MAPPING > 1 ) THEN
          IF (MUMPS_ROOTSSARBR (
     &    PROCNODE_STEPS ( STEP ( IFATH ) ), KEEP(199) )) THEN
            IN_L0INIT = .TRUE.
            GOTO 80
          ENDIF
        ENDIF
      ENDIF
 80   CONTINUE
      IF (.NOT. SKIP_ABOVE) THEN
        IF (KEEP(50).EQ.0) THEN
          FACTOR_SIZE_UNDER_L0 = FACTOR_SiZE_UNDER_L0 +
     &    NPIV8 * ( NFRONT8 + NFRONT8 - NPIV8 )
        ELSE
          FACTOR_SIZE_UNDER_L0  = FACTOR_SIZE_UNDER_L0  +
     &    NFRONT8 * NPIV8
        ENDIF
      ENDIF
      IF ( IN_L0INIT ) THEN
          SMALL_COST_TMP = max(SMALL_COST_TMP, 
     &                      COSTS_MONO_THREAD ( STEP ( INODE ) ) )
          CALL L0_INSERT_NODE ( L0_OMP_DLL, INODE )
          NB_IN_L0 = NB_IN_L0 + 1
      ELSE IF ( SKIP_ABOVE ) THEN
          IERR = IDLL_PUSH_BACK ( LEAFS_ABOVE_L0_OMP_DLL, INODE )
          IF ( .NOT. INODE_IS_A_LEAF ) THEN
            WRITE(*,*)
     &      " Internal error 7 in MUMPS_ANA_INITIALIZE_L0_OMP",
     &      INODE
            CALL MUMPS_ABORT()
          ENDIF
          IPOOL_B_L0_OMP(LEAF+1) = -INODE
      ELSE
          CP_NSTK_STEPS ( STEP ( IFATH ) ) =
     &                          CP_NSTK_STEPS ( STEP ( IFATH ) ) - 1
          IF ( CP_NSTK_STEPS ( STEP ( IFATH ) ) .EQ. 0 ) THEN
              INODE = IFATH
              INODE_IS_A_LEAF = .FALSE.
              GOTO 95
          ENDIF
      END IF
      IF ( LEAF .GT. 0 ) THEN
        GOTO 90
      END IF
      SMALL_COST = max(SMALL_COST_TMP / 100000d0, SMALL_COST)
      SMALL_COST = min(SMALL_COST, 1D6)
 500  CONTINUE
      RETURN
 150  FORMAT(
     & /' ** ALLOC FAILURE IN MUMPS_ANA_INITIALIZE_L0_OMP FOR ',
     & A30)
      END SUBROUTINE MUMPS_ANA_INITIALIZE_L0_OMP
      SUBROUTINE L0_INSERT_NODE ( DLL, INODE )
      IMPLICIT NONE
      INTEGER, INTENT ( IN ) :: INODE
      TYPE ( IDLL_T ), POINTER :: DLL
      INTEGER :: IERR
      TYPE ( IDLL_NODE_T ), POINTER :: IDLL_NODE
      IF ( COSTS_MONO_THREAD ( STEP ( INODE ) ) .LT. SMALL_COST ) THEN
        IERR = IDLL_PUSH_BACK( DLL, INODE )
        RETURN
      ENDIF
      IERR = IDLL_ITERATOR_BEGIN ( DLL, IDLL_NODE )
      DO WHILE ( associated ( IDLL_NODE ) )
          IF ( COSTS_MONO_THREAD ( STEP ( IDLL_NODE%ELMT ) )
     &         .GT.
     &         COSTS_MONO_THREAD ( STEP ( INODE ) ) ) THEN
              IDLL_NODE => IDLL_NODE%NEXT
          ELSE
              EXIT
          END IF
      END DO
      IF ( .NOT. associated ( IDLL_NODE ) ) THEN
          IERR = IDLL_PUSH_BACK(DLL, INODE)
      ELSE
          IERR = IDLL_INSERT_BEFORE(DLL, IDLL_NODE, INODE)
      ENDIF
      RETURN
      END SUBROUTINE L0_INSERT_NODE
      SUBROUTINE L0_INSERT_CHILDREN ( I_FATHER )
      IMPLICIT NONE
      INTEGER, INTENT ( IN ) :: I_FATHER
      INTEGER :: I_SON, IERR, NB_SONS
      TYPE ( IDLL_T ), POINTER :: SON_DLL
      TYPE ( IDLL_NODE_T ), POINTER :: IDLL_NODE
      I_SON = I_FATHER
      DO WHILE ( I_SON .GT. 0 )
          I_SON = FILS ( I_SON )
      END DO
      I_SON = - I_SON
      IF ( I_SON .EQ. 0 ) THEN
          RETURN
      END IF
      IERR = IDLL_CREATE ( SON_DLL )
      NB_SONS = 0
      DO WHILE ( I_SON .GT. 0 )
          CALL L0_INSERT_NODE ( SON_DLL, I_SON )
          I_SON = FRERE ( STEP ( I_SON ) )
          NB_SONS = NB_SONS + 1
      END DO
      NB_IN_L0 = NB_IN_L0 + NB_SONS
      IERR = IDLL_ITERATOR_BEGIN ( L0_OMP_DLL, IDLL_NODE )
      IERR = IDLL_POP_FRONT ( SON_DLL, I_SON )
      IF ( IERR .NE. 0 ) THEN
          GOTO 190
      END IF
      IF ( .NOT. associated( IDLL_NODE ) ) THEN
              DO
                  IERR = IDLL_PUSH_BACK ( L0_OMP_DLL, I_SON )
                  IERR = IDLL_POP_FRONT ( SON_DLL, I_SON )
                  IF ( IERR .NE. 0 ) THEN
                      GOTO 190
                  END IF
              END DO
      ELSE
          DO
              IF ( COSTS_MONO_THREAD ( STEP ( I_SON ) ) .LT.
     &             SMALL_COST ) THEN
                IERR = IDLL_PUSH_BACK(L0_OMP_DLL, I_SON)
                IF (associated(SON_DLL%FRONT)) THEN
                  L0_OMP_DLL%BACK%NEXT => SON_DLL%FRONT
                  SON_DLL%FRONT%PREV   => L0_OMP_DLL%BACK
                  L0_OMP_DLL%BACK      => SON_DLL%BACK
                  NULLIFY(SON_DLL%FRONT)
                  NULLIFY(SON_DLL%BACK)
                ENDIF
                GOTO 190
              ENDIF
              IF ( COSTS_MONO_THREAD ( STEP ( I_SON )) .LT.
     &             COSTS_MONO_THREAD ( STEP ( IDLL_NODE%ELMT ) ) ) THEN
                  IF ( associated ( IDLL_NODE%NEXT ) ) THEN
                      IDLL_NODE => IDLL_NODE%NEXT
                  ELSE
                      IERR = IDLL_PUSH_BACK(L0_OMP_DLL, I_SON)
                      IERR = IDLL_POP_FRONT ( SON_DLL, I_SON )
                      IF ( IERR .NE. 0 ) THEN
                          GOTO 190
                      END IF
                  END IF
              ELSE
                  IERR = IDLL_INSERT_BEFORE(L0_OMP_DLL, IDLL_NODE,I_SON)
                  IERR = IDLL_POP_FRONT ( SON_DLL, I_SON )
                  IF ( IERR .NE. 0 ) THEN
                      GOTO 190
                  END IF
              END IF
          END DO
      END IF
190   CONTINUE
      IERR = IDLL_DESTROY ( SON_DLL )
      RETURN
      END SUBROUTINE L0_INSERT_CHILDREN
      SUBROUTINE L0_REMOVE_NODE ( INODE )
      IMPLICIT NONE
      INTEGER, INTENT ( OUT ) :: INODE
      INTEGER :: I_SON, IERR, NPIV
      IERR = IDLL_POP_FRONT ( L0_OMP_DLL, INODE )
      NB_IN_L0 = NB_IN_L0 - 1
      I_SON = INODE
      NPIV = 0
      DO WHILE ( I_SON .GT. 0 )
          NPIV = NPIV + 1
          I_SON = FILS ( I_SON )
      END DO
      I_SON = - I_SON
      IF (KEEP(50) .EQ. 0) THEN
        FACTOR_SIZE_UNDER_L0 = FACTOR_SIZE_UNDER_L0 -
     &  int(NPIV, 8) * int(2 *  ND(STEP(INODE)) - NPIV, 8)
      ELSE
        FACTOR_SIZE_UNDER_L0 = FACTOR_SIZE_UNDER_L0 -
     &  int(NPIV, 8) * int(ND(STEP(INODE)), 8)
      ENDIF
      IF ( I_SON .EQ. 0 ) THEN 
          IERR = IDLL_PUSH_BACK ( LEAFS_ABOVE_L0_OMP_DLL, INODE )
          INODE = -INODE
      ELSE IF (INODE .GT. 0 .AND. KEEP(403) .NE. 0) THEN
          COST_ABOVE = COST_ABOVE +
     &                 COSTS_MULTI_THREAD(STEP ( abs(INODE) ))
      END IF
      RETURN
      END SUBROUTINE L0_REMOVE_NODE
      FUNCTION MUMPS_ANA_ACCEPT_L0_OMP ()
      LOGICAL :: MUMPS_ANA_ACCEPT_L0_OMP
      INTEGER :: I, I_LESS_CHARGED, IERR
      DOUBLE PRECISION :: LIGHTEST_CHARGE, HEAVIEST_CHARGE
      TYPE ( IDLL_NODE_T ), POINTER :: IDLL_NODE
      LOGICAL :: DECISION_TAKEN
      NB_MAX_IN_L0_ACCEPTL0 = max(NB_MAX_IN_L0_ACCEPTL0, NB_IN_L0)
      LIGHTEST_CHARGE = -9999.0d0 
      HEAVIEST_CHARGE = -9999.0d0 
      IF ( KEEP(403) .EQ. 0) THEN
        IF (NB_IN_L0 .EQ. 0) THEN
          MUMPS_ANA_ACCEPT_L0_OMP = .TRUE.
          DECISION_TAKEN          = .TRUE.
        ELSE IF ( NB_IN_L0 .LT. NB_MAX_IN_L0_ACCEPTL0 .AND.
     &            NB_IN_L0 .LT. KEEP(400) ) THEN
          MUMPS_ANA_ACCEPT_L0_OMP = .TRUE.
          DECISION_TAKEN = .TRUE.
        ELSE IF ( FACTOR_SIZE_UNDER_L0 .GT.
     &            FACTOR_SIZE_PER_MPI * int(THRESH_MEM,8) / 100_8 ) THEN
          MUMPS_ANA_ACCEPT_L0_OMP= .FALSE.
          DECISION_TAKEN = .TRUE.
        ELSE
          DECISION_TAKEN = .FALSE.
        ENDIF
      ELSE
         DECISION_TAKEN = .FALSE.
      ENDIF
      IF (.NOT. DECISION_TAKEN ) THEN
        THREADS_CHARGE = 0.0D0
        IERR = IDLL_ITERATOR_BEGIN( L0_OMP_DLL, IDLL_NODE )
        DO WHILE ( associated ( IDLL_NODE ) )
          I_LESS_CHARGED = 1
          LIGHTEST_CHARGE = THREADS_CHARGE ( 1 )
          DO I = 2, NB_THREADS
              IF ( THREADS_CHARGE ( I ) .LT. LIGHTEST_CHARGE ) THEN
                  I_LESS_CHARGED = I
                  LIGHTEST_CHARGE = THREADS_CHARGE ( I )
              END IF
          END DO
          THREADS_CHARGE ( I_LESS_CHARGED ) =
     &                   THREADS_CHARGE ( I_LESS_CHARGED )
     &                   +
     &                   COSTS_MONO_THREAD ( STEP ( IDLL_NODE%ELMT ) )
          IDLL_NODE => IDLL_NODE%NEXT
        END DO
        LIGHTEST_CHARGE = THREADS_CHARGE ( 1 )
        HEAVIEST_CHARGE = THREADS_CHARGE ( 1 )
        DO I = 2, NB_THREADS
          IF ( THREADS_CHARGE ( I ) .LT. LIGHTEST_CHARGE ) THEN
              LIGHTEST_CHARGE = THREADS_CHARGE ( I )
          ELSEIF ( THREADS_CHARGE ( I ) .GT. HEAVIEST_CHARGE ) THEN
              HEAVIEST_CHARGE = THREADS_CHARGE ( I )
          END IF
        END DO
        COST_UNDER = HEAVIEST_CHARGE
      ENDIF
      IF (KEEP(403) .EQ. 0) THEN
        IF ( .NOT. DECISION_TAKEN ) THEN
          MUMPS_ANA_ACCEPT_L0_OMP =
     &    (
     &    dble(LIGHTEST_CHARGE)/(dble(HEAVIEST_CHARGE)+1.D-12) 
     &    .GT.THRESH_EQUILIB .AND.
     &
     &    FACTOR_SIZE_UNDER_L0 .LE.
     &    FACTOR_SIZE_PER_MPI * int(THRESH_MEM,8) / 100_8
     &
     &    )
     &    .OR.
     &    ( NB_IN_L0 .LT. NB_MAX_IN_L0_ACCEPTL0 .AND.
     &    LIGHTEST_CHARGE .EQ. 0.0D0 ) 
     &    .OR. ( NB_IN_L0 .EQ. 0 ) 
        ENDIF
        IF (MUMPS_ANA_ACCEPT_L0_OMP) THEN
          IF (associated(PHYS_L0_OMP)) THEN
            DEALLOCATE(PHYS_L0_OMP)
            NULLIFY(PHYS_L0_OMP)
          ENDIF
          IERR = IDLL_2_ARRAY ( L0_OMP_DLL, PHYS_L0_OMP, L_PHYS_L0_OMP )
          IF (IERR .EQ. -2) THEN
            INFO(1) = -7
            INFO(2) = L_PHYS_L0_OMP
            RETURN
          ENDIF
        END IF
      ELSE
        IF (COST_UNDER + COST_ABOVE .LT. COST_TOTAL_BEST) THEN
          IF (associated(PHYS_L0_OMP)) THEN
            DEALLOCATE(PHYS_L0_OMP)
            NULLIFY(PHYS_L0_OMP)
          ENDIF
          IERR = IDLL_2_ARRAY ( L0_OMP_DLL, PHYS_L0_OMP, L_PHYS_L0_OMP )
          COST_TOTAL_BEST = COST_UNDER + COST_ABOVE
          NB_REPEAT_ACCEPTL0 = 100
        END IF
        NB_REPEAT_ACCEPTL0 = NB_REPEAT_ACCEPTL0- 1
        MUMPS_ANA_ACCEPT_L0_OMP = (NB_REPEAT_ACCEPTL0 .EQ. 0)
      END IF
      RETURN
      END FUNCTION MUMPS_ANA_ACCEPT_L0_OMP
      SUBROUTINE MUMPS_ANA_FINALIZE_L0_OMP ()
      IMPLICIT NONE
      INTEGER :: INODE, OLD_INODE, I, J, K, LEAF, IERR
      DOUBLE PRECISION :: LIGHTEST_CHARGE
      INTEGER :: I_LESS_CHARGED
      INTEGER :: MAX_TASK_PER_THREAD
      TYPE ( IDLL_NODE_T ), POINTER :: IDLL_NODE
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: THREADS_TASKS
      INTEGER, DIMENSION(:), ALLOCATABLE :: NB_TASK_PER_THREAD
      INTEGER, DIMENSION(:), ALLOCATABLE :: INV_PERM_L0_OMP
      EXTERNAL :: MUMPS_GET_POOL_LENGTH
      INTEGER :: MUMPS_GET_POOL_LENGTH
      IF (KEEP(402) .EQ. 0) THEN
        L_VIRT_L0_OMP = NB_THREADS + 1
      ELSE
        L_VIRT_L0_OMP = L_PHYS_L0_OMP + 1
      END IF
      LPOOL_A_L0_OMP = MUMPS_GET_POOL_LENGTH(NA(1),KEEP(1),KEEP8(1))
      ALLOCATE ( VIRT_L0_OMP ( max(L_VIRT_L0_OMP,1) ),
     &           VIRT_L0_OMP_MAPPING( max(L_VIRT_L0_OMP,1) ),
     &           STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=2*max(L_VIRT_L0_OMP,1)
         IF (LPOK) WRITE(LP,150) 'id%VIRT_L0_OMP[_MAPPING]'
         GOTO 300
      ENDIF
      ALLOCATE ( PERM_L0_OMP ( max(L_PHYS_L0_OMP,1) ), STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=max(L_PHYS_L0_OMP,1)
         IF (LPOK) WRITE(LP,150) 'id%PERM_L0_OMP'
         GOTO 300
      ENDIF
      ALLOCATE ( PTR_LEAFS_L0_OMP ( L_PHYS_L0_OMP + 1 ), STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=max(L_PHYS_L0_OMP,1)
         IF (LPOK) WRITE(LP,150) 'id%PTR_LEAFS_L0_OMP'
         GOTO 300
      ENDIF
      ALLOCATE ( IPOOL_A_L0_OMP ( LPOOL_A_L0_OMP ), STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=LPOOL_A_L0_OMP
         IF (LPOK) WRITE(LP,150) 'id%IPOOL_A_L0_OMP'
         GOTO 300
      ENDIF
      ALLOCATE ( NB_TASK_PER_THREAD ( NB_THREADS ), STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=NB_THREADS
         IF (LPOK) WRITE(LP,150) 'NB_TASK_PER_THREAD'
         GOTO 300
      ENDIF
      ALLOCATE ( INV_PERM_L0_OMP ( L_PHYS_L0_OMP ), STAT=IERR )
      IF(IERR.GT.0) THEN
         WRITE(*,*) "Allocation Error in MUMPS_ANA_FINALIZE_L0_OMP"
         CALL MUMPS_ABORT()
      ENDIF
      NB_TASK_PER_THREAD = 0
      THREADS_CHARGE = 0.0D0
      DO I = 1, L_PHYS_L0_OMP
          I_LESS_CHARGED = 1
          LIGHTEST_CHARGE = THREADS_CHARGE ( 1 )
          DO J = 2, NB_THREADS
              IF ( THREADS_CHARGE ( J ) .LT. LIGHTEST_CHARGE ) THEN
                  I_LESS_CHARGED = J
                  LIGHTEST_CHARGE = THREADS_CHARGE ( J )
                  IF (THREADS_CHARGE( J ) .EQ. 0) THEN
                    EXIT
                  ENDIF
              END IF
          END DO
          NB_TASK_PER_THREAD ( I_LESS_CHARGED ) =
     &                        NB_TASK_PER_THREAD ( I_LESS_CHARGED ) + 1
          IF (KEEP(402) .NE. 0) THEN
            VIRT_L0_OMP_MAPPING(I) = I_LESS_CHARGED
          ENDIF
          THREADS_CHARGE ( I_LESS_CHARGED ) =
     &                  THREADS_CHARGE ( I_LESS_CHARGED )
     &                  +
     &                  COSTS_MONO_THREAD ( STEP ( PHYS_L0_OMP ( I ) ) )
      END DO
      IF (KEEP(402) .EQ. 0) THEN
        DO I = 1, NB_THREADS
          VIRT_L0_OMP_MAPPING(I) = I
        ENDDO
      ENDIF
      VIRT_L0_OMP_MAPPING(L_VIRT_L0_OMP) = -999999
      MAX_TASK_PER_THREAD = 0
      DO I = 1, NB_THREADS 
          MAX_TASK_PER_THREAD = max (MAX_TASK_PER_THREAD,
     &                               NB_TASK_PER_THREAD ( I ) )
      END DO
      ALLOCATE ( THREADS_TASKS ( NB_THREADS, MAX_TASK_PER_THREAD ),
     &     STAT=IERR )
      IF(IERR.GT.0) THEN
         INFO(1)=-7
         INFO(2)=NB_THREADS*MAX_TASK_PER_THREAD
         IF (LPOK) WRITE(LP,150) 'THREADS_TASK'
         GOTO 300
      ENDIF
      NB_TASK_PER_THREAD = 0
      THREADS_CHARGE = 0.0D0
      THREADS_TASKS = 0 
      DO I = 1, L_PHYS_L0_OMP
          I_LESS_CHARGED = 1
          LIGHTEST_CHARGE = THREADS_CHARGE ( 1 )
          DO J = 2, NB_THREADS
              IF ( THREADS_CHARGE ( J ) .LT. LIGHTEST_CHARGE ) THEN
                  I_LESS_CHARGED = J
                  LIGHTEST_CHARGE = THREADS_CHARGE ( J )
              END IF
          END DO
          NB_TASK_PER_THREAD ( I_LESS_CHARGED ) =
     &                         NB_TASK_PER_THREAD ( I_LESS_CHARGED ) + 1
          THREADS_TASKS ( I_LESS_CHARGED, NB_TASK_PER_THREAD
     &                    ( I_LESS_CHARGED ) ) = PHYS_L0_OMP( I )
          THREADS_CHARGE ( I_LESS_CHARGED ) =
     &                  THREADS_CHARGE ( I_LESS_CHARGED )
     &                  +
     &                  COSTS_MONO_THREAD ( STEP ( PHYS_L0_OMP ( I ) ) )
      END DO
      IF (KEEP(402) .EQ. 0) THEN
        K = 1
        DO I = 1, NB_THREADS
          VIRT_L0_OMP (I) = K
          DO J = 1, NB_TASK_PER_THREAD ( I )
              PHYS_L0_OMP (K) = THREADS_TASKS (I,J)
              K = K + 1
          END DO
        END DO
        VIRT_L0_OMP (NB_THREADS+1) = K
      ELSE
        DO I = 1, L_VIRT_L0_OMP 
            VIRT_L0_OMP (I) = I
        END DO
      END IF
      DO I = 1, L_PHYS_L0_OMP
          INV_PERM_L0_OMP ( I ) = I
      END DO
      IF ( L_PHYS_L0_OMP .GT. 1 ) THEN
        CALL MUMPS_QUICK_SORT_PHYS_L0( N, STEP(1), PHYS_L0_OMP(1),
     &     INV_PERM_L0_OMP, L_PHYS_L0_OMP, 1, L_PHYS_L0_OMP )
      ENDIF
      DO I = 1, L_PHYS_L0_OMP
          PERM_L0_OMP( INV_PERM_L0_OMP ( I ) ) = I
      END DO
      J = NBLEAF_MYID 
      PTR_LEAFS_L0_OMP ( 1 ) = J
      DO I = 1, L_PHYS_L0_OMP 
          OLD_INODE = 0
          INODE = PHYS_L0_OMP ( I )
          DO WHILE ( INODE .NE. 0 )
              OLD_INODE = INODE
              DO WHILE ( INODE .GT. 0 )
                  INODE = FILS ( INODE )
              END DO
              INODE = - INODE
          END DO
          DO WHILE ( IPOOL_B_L0_OMP ( J ) .NE. OLD_INODE )
              J = J - 1
          END DO
          J = J - 1
          PTR_LEAFS_L0_OMP ( I + 1 ) = J
      END DO
      CP_NSTK_STEPS(:) = NSTK_STEPS(:)
      IPOOL_A_L0_OMP = 0
      LEAF = 1
      IERR = IDLL_ITERATOR_BEGIN ( LEAFS_ABOVE_L0_OMP_DLL, IDLL_NODE )
      DO WHILE ( associated( IDLL_NODE ) )
          IPOOL_A_L0_OMP ( LEAF ) = IDLL_NODE%ELMT
          LEAF = LEAF + 1
          IDLL_NODE => IDLL_NODE%NEXT
      END DO
      DO I = 1 , L_PHYS_L0_OMP 
        IF ( DAD ( STEP ( PHYS_L0_OMP (I) ) ) .NE. 0 ) THEN
          CP_NSTK_STEPS ( STEP ( DAD ( STEP ( PHYS_L0_OMP (I) ) ) ) ) =
     &     CP_NSTK_STEPS ( STEP ( DAD ( STEP ( PHYS_L0_OMP (I) ) ) ) )-1
          IF (CP_NSTK_STEPS(STEP(DAD(STEP(PHYS_L0_OMP(I))))) .EQ. 0)THEN
            IPOOL_A_L0_OMP ( LEAF ) = DAD(STEP(PHYS_L0_OMP ( I )))
            LEAF = LEAF + 1
          END IF
        END IF
      END DO
      LEAF = LEAF - 1 
      IPOOL_A_L0_OMP(LPOOL_A_L0_OMP) = LEAF
      IPOOL_A_L0_OMP(LPOOL_A_L0_OMP-1) = 0
      IPOOL_A_L0_OMP(LPOOL_A_L0_OMP-2) = 0
      IF (LEAF .GT. 1) THEN
        CALL MUMPS_QUICK_SORT_IPOOL_PO( N, STEP(1), 
     &            IPOOL_A_L0_OMP(1), LEAF, 1, LEAF )
      ENDIF
 300  CONTINUE
      IF (allocated(NB_TASK_PER_THREAD)) DEALLOCATE (NB_TASK_PER_THREAD)
      IF (allocated(INV_PERM_L0_OMP   )) DEALLOCATE ( INV_PERM_L0_OMP )
      IF (allocated(THREADS_TASKS     )) DEALLOCATE (THREADS_TASKS )
      RETURN
 150  FORMAT(
     & /' ** ALLOC FAILURE IN MUMPS_ANA_FINALIZE_L0_OMP FOR ',
     & A30)
      END SUBROUTINE MUMPS_ANA_FINALIZE_L0_OMP 
      SUBROUTINE MUMPS_ANA_FREE_L0_WORKSPACE()
      INTEGER :: IERR
      IF (allocated(THREADS_CHARGE))     DEALLOCATE(THREADS_CHARGE    )
      IF (allocated(CP_NSTK_STEPS ))     DEALLOCATE(CP_NSTK_STEPS     )
      IF (allocated(COSTS_MONO_THREAD))  DEALLOCATE(COSTS_MONO_THREAD )
      IF (allocated(COSTS_MULTI_THREAD)) DEALLOCATE(COSTS_MULTI_THREAD)
      IF (allocated(IPOOL_B_INV))        DEALLOCATE(IPOOL_B_INV       )
      IERR = IDLL_DESTROY ( LEAFS_ABOVE_L0_OMP_DLL )
      IERR = IDLL_DESTROY ( L0_OMP_DLL )
      RETURN
      END SUBROUTINE MUMPS_ANA_FREE_L0_WORKSPACE
      SUBROUTINE READ_BENCH(ARITH, K50)
      IMPLICIT NONE
      INTEGER,      INTENT(in) :: K50
      CHARACTER(1), INTENT(in) :: ARITH
      INTEGER NLINES, INDEX_NPIV, INDEX_NSCHUR, NB_CORE
      INTEGER V, S, OLD_V, OLD_S, I
      PARAMETER(NLINES=2812)
      DOUBLE PRECISION :: AUX
      CHARACTER(1)             :: K50_STR
      INDEX_NPIV = 0
      INDEX_NSCHUR = 0
      OLD_V = -1
      OLD_S = -1
      WRITE(K50_STR,'(I1)') K50
      OPEN(1,FILE=ARITH//'benchmark_sym_'//K50_STR//'.csv')
      DO I=1,NLINES
        READ(1,*) V, S, NB_CORE, AUX
        IF (V .NE. OLD_V) THEN
            INDEX_NPIV = INDEX_NPIV + 1
            OLD_V = V
        END IF
        IF (S .GT. OLD_S) THEN
            INDEX_NSCHUR = INDEX_NSCHUR + 1
            OLD_S = S
        ELSEIF (S .LT. OLD_S) THEN
            INDEX_NSCHUR = 1
            OLD_S = S
        END IF
        BENCH (INDEX_NPIV, INDEX_NSCHUR, NB_CORE) = AUX
      END DO
      CLOSE(1)
      RETURN
      END SUBROUTINE READ_BENCH
      SUBROUTINE COST_BENCH (NPIV, NSCHUR, NB_CORE, SYM, COST)
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  NPIV, NSCHUR, NB_CORE, SYM
      DOUBLE PRECISION, INTENT(OUT) :: COST
      INTEGER V, VV, S, SS
      INTEGER LOW_INDEX_NPIV, LOW_INDEX_NSCHUR
      INTEGER HIGH_INDEX_NPIV, HIGH_INDEX_NSCHUR
      DOUBLE PRECISION :: APROX_COST_FLOPS, REAL_COST_FLOPS
      IF (NPIV .LE. 10) THEN
        LOW_INDEX_NPIV = NPIV
        V = NPIV
        VV = NPIV + 1
      ELSEIF (NPIV .LE. 100) THEN
        LOW_INDEX_NPIV = 9 + NPIV/10
        V = (NPIV/10)*10
        VV = (NPIV/10+1)*10
      ELSEIF (NPIV .LE. 1000) THEN
        LOW_INDEX_NPIV = 18 + NPIV/100
        V = (NPIV/100)*100
        VV = (NPIV/100+1)*100
      ELSEIF (NPIV .LE. 10000) THEN
        LOW_INDEX_NPIV = 27 + NPIV/1000
        V = (NPIV/1000)*1000
        VV = (NPIV/1000+1)*1000
      ELSE
        LOW_INDEX_NPIV = 37
        V = (NPIV/10000)*10000
        VV = (NPIV/10000+1)*10000
      END IF
      IF (NSCHUR .LE. 10) THEN
        LOW_INDEX_NSCHUR = NSCHUR + 1
        S = NSCHUR
        SS = NSCHUR + 1
      ELSEIF (NSCHUR .LE. 100) THEN
        LOW_INDEX_NSCHUR = 10 + NSCHUR/10
        S = (NSCHUR/10)*10
        SS = (NSCHUR/10+1)*10
      ELSEIF (NSCHUR .LE. 1000) THEN
        LOW_INDEX_NSCHUR = 19 + NSCHUR/100
        S = (NSCHUR/100)*100
        SS = (NSCHUR/100+1)*100
      ELSEIF (NSCHUR .LE. 10000) THEN
        LOW_INDEX_NSCHUR = 28 + NSCHUR/1000
        S = (NSCHUR/1000)*1000
        SS = (NSCHUR/1000+1)*1000
      ELSE
        LOW_INDEX_NSCHUR = 38
        S = (NSCHUR/10000)*10000
        SS = (NSCHUR/10000+1)*10000
      END IF
      IF (V .LT. 10000) THEN
        IF (S .LT. 10000) THEN
          HIGH_INDEX_NPIV = LOW_INDEX_NPIV + 1
          HIGH_INDEX_NSCHUR = LOW_INDEX_NSCHUR + 1
          COST = (BENCH(LOW_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE)
     &             *(VV - NPIV)*(SS - NSCHUR)
     &            +BENCH(LOW_INDEX_NPIV, HIGH_INDEX_NSCHUR, NB_CORE)
     &             *(VV - NPIV)*(NSCHUR - S)
     &            +BENCH(HIGH_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE)
     &             *(NPIV - V)*(SS - NSCHUR)
     &            +BENCH(HIGH_INDEX_NPIV, HIGH_INDEX_NSCHUR, NB_CORE)
     &             *(NPIV - V)*(NSCHUR - S))
     &           /((VV - V)*(SS - S))
        ELSE
          HIGH_INDEX_NPIV = LOW_INDEX_NPIV + 1
          HIGH_INDEX_NSCHUR = LOW_INDEX_NSCHUR
          COST = (BENCH(LOW_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE)
     &             *(VV - NPIV)
     &            +BENCH(HIGH_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE)
     &             *(NPIV - V))
     &           /(VV - V)
          CALL MUMPS_GET_FLOPS_COST ( NPIV+NSCHUR, NPIV, NPIV,
     &                                SYM, 1, REAL_COST_FLOPS )
          CALL MUMPS_GET_FLOPS_COST ( V+S, V, V,
     &                                SYM, 1, APROX_COST_FLOPS )
          COST = COST * (REAL_COST_FLOPS/APROX_COST_FLOPS)
        END IF
      ELSE
        IF (NSCHUR < 10000) THEN
          HIGH_INDEX_NPIV = LOW_INDEX_NPIV
          HIGH_INDEX_NSCHUR = LOW_INDEX_NSCHUR + 1
          COST = (BENCH(LOW_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE)
     &             *(SS - NSCHUR)
     &            +BENCH(LOW_INDEX_NPIV, HIGH_INDEX_NSCHUR, NB_CORE)
     &             *(NSCHUR - S))
     &           /(SS - S)
          CALL MUMPS_GET_FLOPS_COST ( NPIV+NSCHUR, NPIV, NPIV,
     &                                SYM, 1, REAL_COST_FLOPS )
          CALL MUMPS_GET_FLOPS_COST ( V+S, V, V,
     &                                SYM, 1, APROX_COST_FLOPS )
          COST = COST * (REAL_COST_FLOPS/APROX_COST_FLOPS)
        ELSE
          HIGH_INDEX_NPIV = LOW_INDEX_NPIV
          HIGH_INDEX_NSCHUR = LOW_INDEX_NSCHUR
          COST = (BENCH(LOW_INDEX_NPIV, LOW_INDEX_NSCHUR, NB_CORE))
          CALL MUMPS_GET_FLOPS_COST ( NPIV+NSCHUR, NPIV, NPIV,
     &                                SYM, 1, REAL_COST_FLOPS )
          CALL MUMPS_GET_FLOPS_COST ( V+S, V, V,
     &                                SYM, 1, APROX_COST_FLOPS )
          COST = COST * (REAL_COST_FLOPS/APROX_COST_FLOPS)
        END IF
      END IF
      END SUBROUTINE COST_BENCH
      END SUBROUTINE MUMPS_ANA_L0_OMP
      END MODULE MUMPS_ANA_OMP_M
      RECURSIVE SUBROUTINE MUMPS_QUICK_SORT_IPOOL_PO( N, STEP, 
     &            INTLIST, TAILLE, LO, HI )
      IMPLICIT NONE
      INTEGER N, TAILLE
      INTEGER STEP( N )
      INTEGER INTLIST( TAILLE )
      INTEGER LO, HI
      INTEGER I,J
      INTEGER ISWAP, PIVOT
      I = LO
      J = HI
      PIVOT = STEP(INTLIST((I+J)/2))
 10   IF (STEP(INTLIST(I)) > PIVOT) THEN
        I=I+1
        GOTO 10
      ENDIF
 20   IF (STEP(INTLIST(J)) < PIVOT) THEN
        J=J-1
        GOTO 20
      ENDIF
      IF (I < J) THEN
        ISWAP = INTLIST(I)
        INTLIST(I) = INTLIST(J)
        INTLIST(J)=ISWAP
      ENDIF
      IF ( I <= J) THEN
        I = I+1
        J = J-1
      ENDIF
      IF ( I <= J ) GOTO 10
      IF ( LO < J ) CALL MUMPS_QUICK_SORT_IPOOL_PO(N, STEP,
     &              INTLIST, TAILLE, LO, J)
      IF ( I < HI ) CALL MUMPS_QUICK_SORT_IPOOL_PO(N, STEP,
     &              INTLIST, TAILLE, I, HI)
      RETURN
      END SUBROUTINE MUMPS_QUICK_SORT_IPOOL_PO
      RECURSIVE SUBROUTINE MUMPS_QUICK_SORT_PHYS_L0( N, STEP, 
     &            INTLIST, INVPERM, TAILLE, LO, HI )
      IMPLICIT NONE
      INTEGER N, TAILLE
      INTEGER STEP( N )
      INTEGER INTLIST( TAILLE )
      INTEGER INVPERM( TAILLE )
      INTEGER LO, HI
      INTEGER I,J
      INTEGER ISWAP, PIVOT
      INTEGER dswap
      I = LO
      J = HI
      PIVOT = STEP(INTLIST((I+J)/2))
 10   IF (STEP(INTLIST(I)) < PIVOT) THEN
        I=I+1
        GOTO 10
      ENDIF
 20   IF (STEP(INTLIST(J)) > PIVOT) THEN
        J=J-1
        GOTO 20
      ENDIF
      IF (I < J) THEN
        ISWAP = INTLIST(I)
        INTLIST(I) = INTLIST(J)
        INTLIST(J)=ISWAP
        dswap = INVPERM(I)
        INVPERM(I) = INVPERM(J)
        INVPERM(J) = dswap
      ENDIF
      IF ( I <= J) THEN
        I = I+1
        J = J-1
      ENDIF
      IF ( I <= J ) GOTO 10
      IF ( LO < J ) CALL MUMPS_QUICK_SORT_PHYS_L0(N, STEP,
     &              INTLIST, INVPERM, TAILLE, LO, J)
      IF ( I < HI ) CALL MUMPS_QUICK_SORT_PHYS_L0(N, STEP,
     &              INTLIST, INVPERM, TAILLE, I, HI)
      RETURN
      END SUBROUTINE MUMPS_QUICK_SORT_PHYS_L0
      SUBROUTINE MUMPS_ANA_OMP_RETURN()
      RETURN
      END SUBROUTINE MUMPS_ANA_OMP_RETURN
