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
      SUBROUTINE SMUMPS_SOL_C(root, N, A, LA, IW, LIW, W, LWC, 
     & IWCB, LIWW, NRHS, NA, LNA, NE_STEPS, W2, MTYPE, ICNTL, FROM_PP,
     & STEP, FRERE, DAD, FILS, PTRIST, PTRFAC, IW1, LIW1, PTRACB,
     & LIWK_PTRACB, PROCNODE_STEPS, SLAVEF, INFO, KEEP,KEEP8, DKEEP,
     & COMM_NODES, MYID, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES,
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE, IBEG_ROOT_DEF, IEND_ROOT_DEF,
     & IROOT_DEF_RHS_COL1, RHS_ROOT, LRHS_ROOT, SIZE_ROOT, MASTER_ROOT,
     & RHSINTR, LRHSINTR, POSINRHSINTR_FWD, POSINRHSINTR_BWD,
     &           Lnodes_FWD, Lnodes_BWD,
     &           nodes_FWD, nodes_BWD,
     & NZ_RHS, NBCOL_INBLOC, NRHS_ORIG, JBEG_RHS, Step2node, LStep2node,
     & IRHS_SPARSE, IRHS_PTR, SIZE_PERM_RHS, PERM_RHS,
     & SIZE_UNS_PERM_INV, UNS_PERM_INV, NB_FS_IN_RHSINTR_F,
     & NB_FS_IN_RHSINTR_TOT, DO_NBSPARSE , RHS_BOUNDS, LRHS_BOUNDS
#if defined(STAT_ES_SOLVE)
     & , IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING
#endif
     &  ,IPOOL_B_L0_OMP, LPOOL_B_L0_OMP, IPOOL_A_L0_OMP, LPOOL_A_L0_OMP,
     &  L_VIRT_L0_OMP, VIRT_L0_OMP, L_PHYS_L0_OMP, PHYS_L0_OMP,
     &  PERM_L0_OMP, PTR_LEAFS_L0_OMP, L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &  L0_OMP_FACTORS, LL0_OMP_FACTORS
     & )
      USE SMUMPS_OOC
      USE SMUMPS_SOL_ES
      USE SMUMPS_SOL_L0OMP_M, ONLY : SMUMPS_SOL_L0OMP_R,
     &                               SMUMPS_SOL_L0OMP_S
      USE SMUMPS_STRUC_DEF, ONLY : SMUMPS_ROOT_STRUC
     & , SMUMPS_L0OMPFAC_T
      IMPLICIT NONE
#if defined(V_T)
      INCLUDE 'VT.inc'
#endif
      TYPE ( SMUMPS_ROOT_STRUC ) :: root
      INTEGER(8) :: LA
      INTEGER(8) :: LWC
      INTEGER    :: N,LIW,MTYPE,LIW1,LIWW,LNA
      INTEGER ICNTL(60),INFO(80), KEEP(500)
      REAL, intent(inout)     :: DKEEP(230)
      INTEGER(8) KEEP8(150)
      INTEGER IW(LIW),IW1(LIW1),NA(LNA),NE_STEPS(KEEP(28)),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)), FILS(N), PTRIST(KEEP(28)),
     &        DAD(KEEP(28))
      INTEGER(8) ::  PTRFAC(KEEP(28))
      INTEGER    ::  LIWK_PTRACB
      INTEGER(8) ::  PTRACB(LIWK_PTRACB)
      INTEGER NRHS, LRHSINTR, NB_FS_IN_RHSINTR_F, NB_FS_IN_RHSINTR_TOT
      REAL    A(LA), W(LWC), 
     &        W2(KEEP(133))
      REAL ::  RHSINTR(LRHSINTR,NRHS)
      INTEGER SLAVEF, COMM_NODES, MYID, MYID_NODES
      INTEGER PROCNODE_STEPS(KEEP(28)), POSINRHSINTR_FWD(N), 
     &        POSINRHSINTR_BWD(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER IBEG_ROOT_DEF, IEND_ROOT_DEF, IROOT_DEF_RHS_COL1
      INTEGER NRHS_LOC
      INTEGER SIZE_ROOT, MASTER_ROOT
      INTEGER(8) :: LRHS_ROOT
      REAL RHS_ROOT(LRHS_ROOT)
      LOGICAL, intent(in) :: FROM_PP 
      INTEGER, intent(in) :: Lnodes_FWD, Lnodes_BWD
      INTEGER, intent(in) :: nodes_FWD(max(1,Lnodes_FWD)), 
     &                       nodes_BWD(max(1,Lnodes_BWD))
      INTEGER, intent(in)      :: NZ_RHS, NBCOL_INBLOC, NRHS_ORIG
      INTEGER, intent(in)      :: SIZE_UNS_PERM_INV 
      INTEGER, intent(in)      :: SIZE_PERM_RHS 
      INTEGER, intent(in) :: JBEG_RHS
      INTEGER, intent(in) :: IRHS_SPARSE(NZ_RHS)
      INTEGER, intent(in) :: IRHS_PTR(NBCOL_INBLOC+1)
      INTEGER, intent(in) :: PERM_RHS(SIZE_PERM_RHS)
      INTEGER, intent(in) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, intent(in) :: LStep2node
      INTEGER, intent(in) :: Step2node(LStep2node)
      LOGICAL, intent(in)  :: DO_NBSPARSE  
      INTEGER, intent(in)     :: LRHS_BOUNDS
      INTEGER, intent(inout)  :: RHS_BOUNDS (LRHS_BOUNDS) 
#if defined(STAT_ES_SOLVE)
       INTEGER, INTENT(IN) :: SIZE_WORKING, SIZE_IPTR_WORKING,
     &                        IPTR_WORKING(SIZE_IPTR_WORKING),
     &                        WORKING(SIZE_WORKING)
#endif
      INTEGER, INTENT (IN) :: LPOOL_B_L0_OMP
      INTEGER, INTENT (IN) :: IPOOL_B_L0_OMP( LPOOL_B_L0_OMP )
      INTEGER, INTENT (IN) :: LPOOL_A_L0_OMP
      INTEGER, INTENT (IN) :: IPOOL_A_L0_OMP( LPOOL_A_L0_OMP )
      INTEGER, INTENT (IN) :: L_PHYS_L0_OMP
      INTEGER, INTENT (INOUT) :: PHYS_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT (IN) :: L_VIRT_L0_OMP
      INTEGER, INTENT (IN) :: VIRT_L0_OMP( L_VIRT_L0_OMP )
      INTEGER, INTENT (IN) :: PERM_L0_OMP( L_PHYS_L0_OMP )
      INTEGER, INTENT (IN) :: PTR_LEAFS_L0_OMP( L_PHYS_L0_OMP + 1)
      INTEGER, INTENT (IN) :: LL0_OMP_MAPPING
      INTEGER, INTENT (IN) :: L0_OMP_MAPPING( LL0_OMP_MAPPING )
      INTEGER, INTENT (IN) :: LL0_OMP_FACTORS
      TYPE (SMUMPS_L0OMPFAC_T), INTENT(IN) ::
     &                        L0_OMP_FACTORS(LL0_OMP_FACTORS)
      INTEGER MP, LP, LDIAG
      INTEGER K,I,II
      INTEGER allocok
      INTEGER LPOOL,MYLEAF,MYROOT,NBROOT,LPANEL_POS
      INTEGER MYLEAF_NOT_PRUNED
      INTEGER NSTK_S,IPOOL,IPANEL_POS,PTRICB
      INTEGER MTYPE_LOC
      INTEGER MODE_RHS_BOUNDS 
      INTEGER IPT_RHS_ROOT_LOC
      INTEGER IERR
      INTEGER(8) :: IAPOS
      INTEGER       IOLDPS,
     &              LOCAL_M,
     &              LOCAL_N
#if defined(V_T)
      INTEGER soln_c_class, forw_soln, back_soln, root_soln
#endif
      LOGICAL DOFORWARD, DOROOT, DOBACKWARD
      LOGICAL :: DO_L0OMP_FWD, DO_L0OMP_BWD
      LOGICAL I_WORKED_ON_ROOT, SPECIAL_ROOT_REACHED
      INTEGER IROOT
      LOGICAL DOROOT_FWD_OOC, DOROOT_BWD_PANEL
      LOGICAL DUMMY_BOOL
      INTEGER :: IDUMMY
      INTEGER :: NBROOT_UNDER_L0
      REAL, PARAMETER :: ZERO = 0.0E0
      INCLUDE 'mumps_headers.h'
      INTEGER, DIMENSION(:), POINTER             ::  nodes_BWD_PTR
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: Pruned_Roots_FWD
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: Pruned_Roots_NS
      INTEGER  :: Lnodes_BWD_PTR, LPruned_Roots_NS
      INTEGER :: Lnodes_BWD_ROOTS
      INTEGER nb_prun_leaves
      INTEGER, DIMENSION(:), ALLOCATABLE :: Pruned_Leaves
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_List
      INTEGER  nb_prun_nodes
      INTEGER INODE_PRINC, nb_prun_roots
      INTEGER, DIMENSION(:), ALLOCATABLE ::  Pruned_SONS, Pruned_Roots
      INTEGER                            :: SIZE_TO_PROCESS
      LOGICAL, DIMENSION(:), ALLOCATABLE :: TO_PROCESS
      INTEGER ISTEP
      INTEGER :: INODE, ICHILD
      LOGICAL AM1, DO_PRUN_FWD, DO_PRUN_BWD
      LOGICAL Exploit_Sparsity_FWD, Exploit_Sparsity_BWD
      LOGICAL DO_NBSPARSE_BWD, PRUN_BELOW_BWD
      INTEGER :: OOC_FCT_TYPE_TMP
      INTEGER :: MUMPS_OOC_GET_FCT_TYPE
      EXTERNAL :: MUMPS_OOC_GET_FCT_TYPE
      DOUBLE PRECISION TIME_FWD,TIME_BWD,TIME_SpecialRoot
      INTEGER :: nb_sparse
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      LOGICAL, EXTERNAL :: MUMPS_IN_OR_ROOT_SSARBR
      MYLEAF = -1
      LP      = ICNTL(1)
      MP      = ICNTL(2)
      LDIAG   = ICNTL(4)
      CALL MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL, INFO, MYID )
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF(INFO(1).LT.0) GOTO 500
#if defined(V_T)
      CALL VTCLASSDEF( 'Soln_c',soln_c_class,ierr)
      CALL VTFUNCDEF( 'forw_soln',soln_c_class,forw_soln,ierr)
      CALL VTFUNCDEF( 'back_soln',soln_c_class,back_soln,ierr)
      CALL VTFUNCDEF( 'root_soln',soln_c_class,root_soln,ierr)
#endif
      IF (.NOT. FROM_PP) THEN
        CALL MUMPS_SECDEB(TIME_FWD)
      ENDIF
      NSTK_S   = 1
      PTRICB = NSTK_S + KEEP(28)
      IPOOL = PTRICB + KEEP(28)
      LPOOL = NA(1) + 1
      IPANEL_POS = IPOOL + LPOOL
      IF (KEEP(201).EQ.1) THEN
        LPANEL_POS = KEEP(228)+1
      ELSE
        LPANEL_POS = 1
      ENDIF
      IF (IPANEL_POS + LPANEL_POS -1 .ne. LIW1 )  THEN
         WRITE(*,*)  MYID, ": Internal Error 1 in SMUMPS_SOL_C",
     &   IPANEL_POS, LPANEL_POS, LIW1
         CALL MUMPS_ABORT()
      ENDIF
      KEEP(405)=0
      DOFORWARD = .TRUE.
      DOBACKWARD= .TRUE.
      SPECIAL_ROOT_REACHED = .TRUE.
      IF ( KEEP(111).NE.0 .OR. KEEP(252).NE.0 ) THEN
        DOFORWARD = .FALSE.
      ENDIF
      IF (KEEP(221).eq.1) DOBACKWARD = .FALSE.
      IF (KEEP(221).eq.2) DOFORWARD  = .FALSE.
      IF ( KEEP(60).EQ.0 .AND.
     &    ( 
     &      (KEEP(38).NE.0 .AND.  root%yes) 
     &  .OR.
     &      (KEEP(20).NE.0 .AND. MYID_NODES.EQ.MASTER_ROOT)
     &    ) 
     &  .AND. KEEP(252).EQ.0
     &   )
     &THEN
        DOROOT = .TRUE.
      ELSE
        DOROOT = .FALSE.
      ENDIF
      DOROOT_BWD_PANEL = DOROOT .AND. MTYPE.NE.1 .AND. KEEP(50).EQ.0
     &                     .AND. KEEP(201).EQ.1
      DOROOT_FWD_OOC = DOROOT .AND. .NOT.DOROOT_BWD_PANEL
      AM1              = (KEEP(237) .NE. 0)
      Exploit_Sparsity_FWD  = (KEEP(235) .NE. 0) .AND. (.NOT. AM1)
      Exploit_Sparsity_BWD  = (KEEP(212) .NE. 0) .AND. (.NOT. AM1)
      Lnodes_BWD_ROOTS = NA(2)
      DO_PRUN_FWD          = (Exploit_Sparsity_FWD.OR.AM1)
      DO_PRUN_BWD          = (Exploit_Sparsity_BWD.OR.AM1)
      IF (FROM_PP) THEN
        Exploit_Sparsity_FWD = .FALSE.
        DO_PRUN_FWD          = .FALSE.
        Exploit_Sparsity_BWD = .FALSE.
        DO_PRUN_BWD          = .FALSE.
           IF ( AM1 ) THEN
          WRITE(*,*) "Internal error 2 in SMUMPS_SOL_C"
          CALL MUMPS_ABORT()
        ENDIF
      ENDIF
      DO_L0OMP_FWD= ( (KEEP(401).EQ.1).AND.(KEEP(400).GT.0)
     &             .AND.DOFORWARD )
      DO_L0OMP_FWD = DO_L0OMP_FWD .AND. KEEP(201).EQ.0
      DO_L0OMP_BWD = ( (KEEP(401).EQ.1).AND.(KEEP(400).GT.0)
     &              .AND.DOBACKWARD )
      DO_L0OMP_BWD = DO_L0OMP_BWD .AND. KEEP(201).EQ.0
      IF ( DO_PRUN_FWD.OR.DO_PRUN_BWD ) THEN
         ALLOCATE (Pruned_SONS(KEEP(28)), stat=I)
         IF(I.GT.0) THEN
               INFO(1)=-13
               INFO(2)=KEEP(28)
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
      ENDIF
      IF ( DO_PRUN_FWD.OR.DO_PRUN_BWD
     &     .OR. DO_L0OMP_BWD
     &   ) THEN
         SIZE_TO_PROCESS = KEEP(28)
      ELSE
         SIZE_TO_PROCESS = 1
      ENDIF
      ALLOCATE (TO_PROCESS(SIZE_TO_PROCESS), stat=I)
      IF(I.GT.0) THEN
          INFO(1)=-13
          INFO(2)=KEEP(28)
      END IF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF(INFO(1).LT.0) GOTO 500
      IF ( DOFORWARD .AND. DO_PRUN_FWD ) THEN
         CALL SMUMPS_CHAIN_PRUN_NODES( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_FWD, Lnodes_FWD,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots,
     &        nb_prun_leaves )  
         ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
         IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
         END IF
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF(INFO(1).LT.0) GOTO 500
         CALL SMUMPS_CHAIN_PRUN_NODES( 
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_FWD, Lnodes_FWD,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
         CALL SMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
         IF ( KEEP(201) .GT. 0) THEN
           OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                      ('F',MTYPE,KEEP(201),KEEP(50))
         ELSE
           OOC_FCT_TYPE_TMP = -5959 
         ENDIF
         CALL SMUMPS_CHAIN_PRUN_NODES_STATS(
     &        MYID_NODES, N, KEEP(28), KEEP(201), KEEP(485), 
     &        KEEP8(31)+KEEP8(64), 
     &        STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
#if defined(STAT_ES_SOLVE)
     &        , NRHS, COMM_NODES, IW, LIW, PTRIST,KEEP(IXSZ),0,
     &        KEEP(50), KEEP(38)
#endif
     &   )
         IF (DO_NBSPARSE) THEN
           nb_sparse = max(1,KEEP(497))
           MODE_RHS_BOUNDS = 0 
           IF (Exploit_Sparsity_FWD) MODE_RHS_BOUNDS = 2
           CALL SMUMPS_INITIALIZE_RHS_BOUNDS(
     &      STEP, N,
     &      IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &      JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, KEEP(242), KEEP(243),
     &      UNS_PERM_INV, SIZE_UNS_PERM_INV, KEEP(23),
     &      RHS_BOUNDS, KEEP(28),
     &      nb_sparse, MYID_NODES,
     &      MODE_RHS_BOUNDS) 
           CALL SMUMPS_PROPAGATE_RHS_BOUNDS(
     &      Pruned_Leaves, nb_prun_leaves,
     &      STEP, N, Pruned_SONS,
     &      DAD, RHS_BOUNDS, KEEP(28),
     &      MYID_NODES, COMM_NODES, KEEP(485),
#if defined(STAT_ES_SOLVE)
     &      KEEP(46), 
     &      IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING,
#endif
     &      IW, LIW, PTRIST,KEEP(IXSZ),OOC_FCT_TYPE_TMP,0,
     &      KEEP(50), KEEP(38))
         END IF
         SPECIAL_ROOT_REACHED = .FALSE.
         DO I= 1, nb_prun_roots
          IF ( (Pruned_Roots(I).EQ.KEEP(38)).OR.
     &         (Pruned_Roots(I).EQ.KEEP(20)) ) THEN
            SPECIAL_ROOT_REACHED = .TRUE.
            EXIT
          ENDIF
         ENDDO
         DEALLOCATE(Pruned_List)
      ENDIF  
      IF (KEEP(201).GT.0) THEN
        IF (DOFORWARD .OR. DOROOT_FWD_OOC) THEN
           CALL SMUMPS_SOLVE_INIT_OOC_FWD(PTRFAC,KEEP(28),MTYPE,
     &                                A,LA,DOFORWARD,IERR)
           IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            CALL MUMPS_ABORT()
          ENDIF
        ENDIF
      ENDIF
      IF (DOFORWARD) THEN
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 1
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(forw_soln,ierr)
#endif
        IF ( .NOT. DO_PRUN_FWD ) THEN
          CALL MUMPS_INIT_NROOT_DIST(N, NBROOT, MYROOT, MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP, STEP, PROCNODE_STEPS)
          DO ISTEP =1, KEEP(28)
            IW1(NSTK_S+ISTEP-1) = NE_STEPS(ISTEP)
          ENDDO
        ELSE
          CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &         nb_prun_roots, Pruned_Roots,
     &         MYROOT, MYID_NODES, SLAVEF, KEEP, STEP,
     &         PROCNODE_STEPS )
          IF ((Exploit_Sparsity_FWD).AND.(nb_prun_roots.NE.NA(2))) THEN
             Lnodes_BWD_ROOTS = nb_prun_roots
            ALLOCATE(Pruned_Roots_FWD(Lnodes_BWD_ROOTS), STAT=allocok)
            IF(allocok.GT.0) THEN
              WRITE(*,*)'Problem with allocation of nodes_FWD'
              INFO(1) = -13
              INFO(2) = Lnodes_BWD_ROOTS
              CALL MUMPS_ABORT()
            END IF
            Pruned_Roots_FWD(1:Lnodes_BWD_ROOTS)=
     &                        Pruned_Roots(1:Lnodes_BWD_ROOTS)
            DEALLOCATE(Pruned_Roots)
          ELSE
            DEALLOCATE(Pruned_Roots)
          ENDIF
          DO ISTEP = 1, KEEP(28)
            IW1(NSTK_S+ISTEP-1) = Pruned_SONS(ISTEP)
          ENDDO
        ENDIF
        IF ( DO_L0OMP_FWD ) THEN
          KEEP(405)=1
          CALL SMUMPS_SOL_L0OMP_R( N, MTYPE_LOC, NRHS, LIW, IW,
     &    IW1(PTRICB), RHSINTR, LRHSINTR, POSINRHSINTR_FWD,
     &    STEP, FRERE, DAD, FILS, IW1(NSTK_S), 
     &    PTRIST, PTRFAC, INFO,
     &    KEEP, KEEP8, DKEEP, PROCNODE_STEPS, SLAVEF,
     &    COMM_NODES, MYID_NODES,
     &    BUFR, LBUFR, LBUFR_BYTES,
     &    RHS_ROOT, LRHS_ROOT,     
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     &    RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE,
     &    FROM_PP,
     &    NBROOT_UNDER_L0,
     &    LPOOL_B_L0_OMP, IPOOL_B_L0_OMP,
     &    L_VIRT_L0_OMP, VIRT_L0_OMP,
     &    L_PHYS_L0_OMP, PHYS_L0_OMP,
     &    PERM_L0_OMP, PTR_LEAFS_L0_OMP,
     &    L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &    L0_OMP_FACTORS, LL0_OMP_FACTORS,
     &    DO_PRUN_FWD, TO_PROCESS
     &    )
          CALL MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL,
     &                                    INFO, MYID )
          IF (INFO(1).LT.0) THEN
            CALL SMUMPS_BDC_ERROR(MYID_NODES, SLAVEF, COMM_NODES, KEEP)
          ENDIF
          KEEP(405)=0
          MYROOT = MYROOT - NBROOT_UNDER_L0
        ENDIF
        IF ( DO_L0OMP_FWD ) THEN
          IF ( DO_PRUN_FWD ) THEN
            MYLEAF_NOT_PRUNED = IPOOL_A_L0_OMP(LPOOL_A_L0_OMP)
            DO I=1, MYLEAF_NOT_PRUNED
              IF ( TO_PROCESS( STEP( IPOOL_A_L0_OMP(I) ))) THEN
                IW1(IPOOL+MYLEAF-1) = IPOOL_A_L0_OMP(I)
                IW1(NSTK_S+STEP(IPOOL_A_L0_OMP(I))-1) = -99
              ENDIF
            ENDDO
            DO I = 1, nb_prun_leaves
              INODE = Pruned_Leaves(I)
              IF ( MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
     &         .EQ. MYID_NODES ) THEN
                IF (L0_OMP_MAPPING( STEP(INODE) ) .EQ. 0) THEN
                     IW1(NSTK_S+STEP(INODE)-1) = -99
                ENDIF
              ENDIF
            ENDDO
            DO I = 1, L_PHYS_L0_OMP
              INODE = DAD(STEP(PHYS_L0_OMP(I)))
              IF (INODE .NE. 0) THEN
                IF ( TO_PROCESS( STEP( INODE ))) THEN
                  IF ( IW1(NSTK_S+STEP(INODE)-1) .EQ. 0 ) THEN
                       IW1(NSTK_S+STEP(INODE)-1) = -99
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
            MYLEAF = 0
            DO ISTEP = KEEP(28), 1, -1
              INODE=Step2Node(ISTEP)
                IF (IW1(NSTK_S+STEP(INODE)-1).EQ.-99) THEN
                  MYLEAF = MYLEAF + 1
                  IW1(IPOOL+MYLEAF-1) = INODE
                  IW1(NSTK_S+STEP(INODE)-1) = 0
                ENDIF
            ENDDO
            DEALLOCATE(Pruned_Leaves)
          ELSE
            MYLEAF = IPOOL_A_L0_OMP(LPOOL_A_L0_OMP)
            DO I=1, MYLEAF
              IW1(IPOOL+I-1) = IPOOL_A_L0_OMP(I)
            ENDDO
          ENDIF
        ELSE
          IF ( DO_PRUN_FWD ) THEN
            CALL MUMPS_INIT_POOL_DIST_NONA( N, MYLEAF, MYID_NODES,
     &            nb_prun_leaves, Pruned_Leaves, KEEP, KEEP8,
     &            STEP, PROCNODE_STEPS, IW1(IPOOL), LPOOL )
            MYLEAF = MYLEAF - 1
            DEALLOCATE(Pruned_Leaves)
          ELSE
            CALL MUMPS_INIT_POOL_DIST( N, MYLEAF, MYID_NODES,
     &             SLAVEF, NA, LNA, KEEP, KEEP8, STEP,
     &             PROCNODE_STEPS, IW1(IPOOL), LPOOL )
            MYLEAF = MYLEAF - 1 
          ENDIF
        ENDIF
        CALL SMUMPS_SOL_R(N, A(1), LA, IW(1), LIW, W(1),
     &       LWC, NRHS,
     &       IW1(PTRICB), IWCB, LIWW,
     &       RHSINTR,LRHSINTR,POSINRHSINTR_FWD,
     &       STEP, FRERE,DAD,FILS,
     &       IW1(NSTK_S),IW1(IPOOL),LPOOL,PTRIST,PTRFAC,
     &       MYLEAF, MYROOT, INFO,
     &       KEEP, KEEP8, DKEEP,
     &       PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &       BUFR, LBUFR, LBUFR_BYTES,
     &       RHS_ROOT, LRHS_ROOT, MTYPE_LOC, 
     &       
     &       ISTEP_TO_INIV2, TAB_POS_IN_PERE
     &       , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &          , L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &            L0_OMP_FACTORS, LL0_OMP_FACTORS
     &       )
        IF (DO_PRUN_FWD) THEN
          MYLEAF = -1
        ENDIF
#if defined(V_T)
        CALL VTEND(forw_soln,ierr)
#endif
      ENDIF
      CALL MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL, INFO, MYID )
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF(INFO(1).LT.0) GOTO 500
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) THEN
        IF ( LP .GT. 0 ) THEN
          WRITE(LP,*) MYID,
     &    ': ** ERROR RETURN FROM SMUMPS_SOL_R,INFO(1:2)=',
     &    INFO(1:2)
        END IF
        GOTO 500   
      END IF
      CALL MPI_BARRIER( COMM_NODES, IERR )
      IF (.NOT.FROM_PP) THEN
         CALL MUMPS_SECFIN(TIME_FWD)
         DKEEP(117)=real(TIME_FWD) + DKEEP(117)
      ENDIF
      IF ( .NOT.( 
     &       DOBACKWARD.AND.
     &      (DO_PRUN_BWD.OR.(Lnodes_BWD_ROOTS.NE.NA(2)))
     &          )
     &   ) THEN
         IF (.NOT. DO_L0OMP_BWD ) THEN
           IF ( allocated(TO_PROCESS) .AND. SIZE_TO_PROCESS.NE.1 ) THEN
             DEALLOCATE (TO_PROCESS)
             SIZE_TO_PROCESS = 1
             ALLOCATE(TO_PROCESS(SIZE_TO_PROCESS),stat=I)
           ENDIF
         ENDIF
      ENDIF
      IF ( (KEEP(111).NE.0).AND.DOBACKWARD.AND.
     &     (
     &      DO_PRUN_BWD
     &     )
     &   )   THEN
        nb_prun_leaves = 0
        IF ( Lnodes_BWD_ROOTS.NE.NA(2) ) THEN
          nodes_BWD_PTR => Pruned_Roots_FWD
          Lnodes_BWD_PTR = Lnodes_BWD_ROOTS
        ELSE IF ( (Exploit_Sparsity_BWD.AND.KEEP(111).NE.0) 
     &          ) THEN
          LPruned_Roots_NS = 0
          Pruned_SONS(:)   = -1 
          DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)          
            ISTEP = abs(STEP(I))
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               LPruned_Roots_NS   = LPruned_Roots_NS +1
               Pruned_SONS(ISTEP) = 0
            ENDIF
          ENDDO
          ALLOCATE(Pruned_Roots_NS(LPruned_Roots_NS), STAT = allocok) 
          IF(allocok.GT.0) THEN
           WRITE(*,*)'Problem with allocation of nodes_BWD'
           INFO(1) = -13
           INFO(2) = LPruned_Roots_NS
           CALL MUMPS_ABORT()
          END IF
          LPruned_Roots_NS = 0         
          Pruned_SONS(:) = -1  
          DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
            ISTEP = abs(STEP(I))
            INODE_PRINC = Step2node(ISTEP)
            IF ( Pruned_SONS(ISTEP) .eq. -1) THEN
               LPruned_Roots_NS = LPruned_Roots_NS +1
               Pruned_Roots_NS(LPruned_Roots_NS)  = INODE_PRINC
               Pruned_SONS(ISTEP) = 0
            ENDIF
          ENDDO
          nodes_BWD_PTR => Pruned_Roots_NS
          Lnodes_BWD_PTR = LPruned_Roots_NS
        ENDIF
        IF ( 
     &       (Exploit_Sparsity_BWD.AND.KEEP(111).NE.0) .OR.
     &       (Lnodes_BWD_ROOTS.NE.NA(2))
     &     ) THEN
           CALL SMUMPS_TREE_PRUN_NODES( 
     &     .FALSE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_BWD_PTR, Lnodes_BWD_PTR,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves
     &     )
           ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_nodes
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_roots
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
           IF(allocok.GT.0) THEN
              INFO(1)=-13
              INFO(2)=nb_prun_leaves
           END IF
           CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
           IF(INFO(1).LT.0) GOTO 500
           CALL SMUMPS_TREE_PRUN_NODES( 
     &     .TRUE.,
     &     DAD, NE_STEPS, FRERE, KEEP(28),
     &     FILS, STEP, N,
     &     nodes_BWD_PTR, Lnodes_BWD_PTR,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
           IF(allocated(Pruned_Roots_NS)) DEALLOCATE(Pruned_Roots_NS)
           IF(allocated(Pruned_Roots_FWD)) DEALLOCATE(Pruned_Roots_FWD)
           CALL SMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
           IF (KEEP(201).GT.0) THEN
             OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
           ELSE
             OOC_FCT_TYPE_TMP = -5959 
           ENDIF
        ENDIF
      ENDIF
      IF(KEEP(201).EQ.1.AND.DOROOT_BWD_PANEL) THEN
         I_WORKED_ON_ROOT = .FALSE. 
         CALL SMUMPS_SOLVE_INIT_OOC_BWD(PTRFAC,KEEP(28),MTYPE,
     &   I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
         IF (IERR .LT. 0) THEN
           INFO(1) = -90
           INFO(2) = IERR
         ENDIF 
      ENDIF
      IF (KEEP(201).EQ.1) THEN
         CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
         IF ( INFO(1) .LT. 0 ) GOTO 500  
      ENDIF
      IF (KEEP(60).NE.0 .AND. KEEP(221).EQ.0
     &   .AND. MYID_NODES .EQ. MASTER_ROOT) THEN
        RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
      ENDIF
      IF (.NOT. FROM_PP) THEN
        CALL MUMPS_SECDEB(TIME_SpecialRoot)
      ENDIF
      IF ( ( KEEP( 38 ) .NE. 0 ).AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( KEEP(60) .EQ. 0 .AND. KEEP(252) .EQ. 0 ) THEN
          IF ( root%yes ) THEN
            IF (KEEP(201).GT.0) THEN
              IF ( (Exploit_Sparsity_FWD.AND.(KEEP(111).NE.0)) .and.
     &            (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) ) THEN
                  GOTO 1010
              ENDIF
            ENDIF
            IOLDPS = PTRIST(STEP(KEEP(38)))
            LOCAL_M = IW( IOLDPS + 2 + KEEP(IXSZ))
            LOCAL_N = IW( IOLDPS + 1 + KEEP(IXSZ))
            IF (KEEP(201).GT.0) THEN 
              CALL SMUMPS_SOLVE_GET_OOC_NODE(
     &           KEEP(38),PTRFAC,KEEP,A,LA,
     &           STEP,KEEP8,N,DUMMY_BOOL,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                WRITE(*,*) '** ERROR after SMUMPS_SOLVE_GET_OOC_NODE',
     &          INFO(1)
                call MUMPS_ABORT()
              ENDIF
            ENDIF
            IAPOS   = PTRFAC(IW( IOLDPS + 4 + KEEP(IXSZ)))
            IF (LOCAL_M * LOCAL_N .EQ. 0) THEN
              IAPOS = min(IAPOS, LA)
            ENDIF
#if defined(V_T)
            CALL VTBEGIN(root_soln,ierr)
#endif
             CALL SMUMPS_ROOT_SOLVE( NRHS, root%DESCRIPTOR(1), 
     &       root%CNTXT_BLACS, LOCAL_M, LOCAL_N,
     &       root%MBLOCK, root%NBLOCK,
     &       root%IPIV(1), root%LPIV, MASTER_ROOT, MYID_NODES,
     &       COMM_NODES,
     &       RHS_ROOT(1),
     &       root%TOT_ROOT_SIZE, A( IAPOS ),
     &       INFO(1), MTYPE, KEEP(50), FROM_PP)
            IF(KEEP(201).GT.0)THEN 
              CALL SMUMPS_FREE_FACTORS_FOR_SOLVE(KEEP(38),
     &             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
              IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 WRITE(*,*)
     &           '** ERROR after SMUMPS_FREE_FACTORS_FOR_SOLVE ',
     &           INFO(1)
                 call MUMPS_ABORT()
              ENDIF
            ENDIF
          ENDIF  
        ENDIF
      ELSE IF ( ( KEEP(20) .NE. 0) .AND. SPECIAL_ROOT_REACHED ) THEN
        IF ( MYID_NODES .eq.  MASTER_ROOT ) THEN
          IF ( KEEP(60) .eq. 0 ) THEN
            IF (KEEP(201).GT.0) THEN
              CALL SMUMPS_SOLVE_GET_OOC_NODE(
     &              KEEP(20),PTRFAC,KEEP,A,LA,
     &              STEP,KEEP8,N,DUMMY_BOOL,IERR)
              IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                WRITE(*,*) '** ERROR after SMUMPS_SOLVE_GET_OOC_NODE', 
     &          INFO(1)
                call MUMPS_ABORT()
              ENDIF
            END IF
            NRHS_LOC         = NRHS
            IPT_RHS_ROOT_LOC = 1
            IF ( KEEP(111).NE.0 ) THEN
              RHS_ROOT( 1: NRHS*SIZE_ROOT) = ZERO
              NRHS_LOC = IEND_ROOT_DEF - IBEG_ROOT_DEF + 1
              IPT_RHS_ROOT_LOC = IPT_RHS_ROOT_LOC +
     &                           (IROOT_DEF_RHS_COL1-1)*SIZE_ROOT
           ENDIF
           IF (NRHS_LOC .GT. 0) THEN
              CALL SMUMPS_SEQ_SOLVE_ROOT_SVD_QR(NRHS_LOC,
     &             SIZE_ROOT,A( PTRFAC(
     &             IW( PTRIST(STEP(KEEP(20)))+4+KEEP(IXSZ)))),
     &             root, IBEG_ROOT_DEF, IEND_ROOT_DEF,
     &             RHS_ROOT( IPT_RHS_ROOT_LOC ),
     &             KEEP,KEEP8,
     &             MTYPE,INFO,LWC,W(1), LP)
        ENDIF
            IF(KEEP(201).GT.0)THEN
              CALL SMUMPS_FREE_FACTORS_FOR_SOLVE(KEEP(20),
     &             PTRFAC,KEEP(28),A,LA,.FALSE.,IERR)
              IF(IERR.LT.0)THEN
                 INFO(1)=IERR
                 INFO(2)=0
                 WRITE(*,*)
     &           '** ERROR after SMUMPS_FREE_FACTORS_FOR_SOLVE ', 
     &           INFO(1)
                 call MUMPS_ABORT()
              ENDIF
            ENDIF
          ENDIF
        END IF 
      END IF 
      IF (.NOT.FROM_PP) THEN
         CALL MUMPS_SECFIN(TIME_SpecialRoot)
         DKEEP(119)=real(TIME_SpecialRoot) + DKEEP(119)
      ENDIF
#if defined(V_T)
      CALL VTEND(root_soln,ierr)
#endif
 1010 CONTINUE
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF ( INFO(1) .LT. 0 ) RETURN
      IF (DOBACKWARD) THEN
        IF ( KEEP(201).GT.0 .AND.  .NOT. DOROOT_BWD_PANEL )
     &    THEN
          I_WORKED_ON_ROOT = DOROOT
          IF (KEEP(38).gt.0 ) THEN
             IF ( ( Exploit_Sparsity_FWD.AND.(KEEP(111).EQ.0) )
     &            .OR. AM1 ) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   OOC_STATE_NODE(STEP(KEEP(38)))=-4
                ENDIF
             ENDIF
             IF (Exploit_Sparsity_FWD.AND.(KEEP(111).NE.0)) THEN
                IF (OOC_STATE_NODE(STEP(KEEP(38))).eq.-6) THEN
                   I_WORKED_ON_ROOT = .FALSE.
                ENDIF
             ENDIF
          ENDIF
        ENDIF                    
        IF (.NOT.AM1) THEN
          DO_NBSPARSE_BWD = .FALSE.
        ELSE
          DO_NBSPARSE_BWD = DO_NBSPARSE
        ENDIF
        PRUN_BELOW_BWD = DO_PRUN_BWD .AND. KEEP(111).EQ.0
        PRUN_BELOW_BWD = PRUN_BELOW_BWD .OR. DO_L0OMP_BWD
        IF ( DO_PRUN_BWD )  THEN
          CALL SMUMPS_CHAIN_PRUN_NODES( 
     &        .FALSE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_BWD, Lnodes_BWD,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots,
     &        nb_prun_leaves)  
          ALLOCATE(Pruned_List(nb_prun_nodes), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_nodes
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          ALLOCATE(Pruned_Roots(nb_prun_roots), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_roots
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          ALLOCATE(Pruned_Leaves(nb_prun_leaves), STAT=allocok)
          IF(allocok.GT.0) THEN
            INFO(1)=-13
            INFO(2)=nb_prun_leaves
          END IF
          CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
          IF(INFO(1).LT.0) GOTO 500
          CALL SMUMPS_CHAIN_PRUN_NODES(
     &        .TRUE.,
     &        DAD, KEEP(28),
     &        STEP, N,
     &        nodes_BWD, Lnodes_BWD,
     &        Pruned_SONS, TO_PROCESS,
     &        nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &        Pruned_List, Pruned_Roots, Pruned_Leaves )
          CALL SMUMPS_OOC_SET_STATES_ES(N,
     &          KEEP(201), Pruned_List, nb_prun_nodes,
     &          STEP)
          IF (KEEP(201).GT.0) THEN
           OOC_FCT_TYPE_TMP=MUMPS_OOC_GET_FCT_TYPE
     &                    ('B',MTYPE,KEEP(201),KEEP(50))
          ELSE
           OOC_FCT_TYPE_TMP = -5959 
          ENDIF
          CALL SMUMPS_CHAIN_PRUN_NODES_STATS(
     &    MYID_NODES, N, KEEP(28), KEEP(201), KEEP(485), KEEP8(31),
     &    STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_TMP
#if defined(STAT_ES_SOLVE)
     &          , NRHS, COMM_NODES, IW, LIW, PTRIST,KEEP(IXSZ),1,
     &          KEEP(50), KEEP(38)
#endif
     &    )
          IF (DO_NBSPARSE_BWD) THEN
            nb_sparse = max(1,KEEP(497))
            CALL SMUMPS_INITIALIZE_RHS_BOUNDS(
     &       STEP, N,
     &       IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &       JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, KEEP(242), KEEP(243),
     &       UNS_PERM_INV, SIZE_UNS_PERM_INV, KEEP(23),
     &       RHS_BOUNDS, KEEP(28),
     &       nb_sparse, MYID_NODES,
     &       1) 
            CALL SMUMPS_PROPAGATE_RHS_BOUNDS(
     &       Pruned_Leaves, nb_prun_leaves,
     &       STEP, N, Pruned_SONS,
     &       DAD, RHS_BOUNDS, KEEP(28),
     &       MYID_NODES, COMM_NODES, KEEP(485),
#if defined(STAT_ES_SOLVE)
     &       KEEP(46), 
     &       IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING,
#endif     
     &       IW, LIW, PTRIST,KEEP(IXSZ),OOC_FCT_TYPE_TMP,1,
     &       KEEP(50), KEEP(38))
         END IF
        ENDIF
        IF ( KEEP(201).GT.0 ) THEN
          IROOT = max(KEEP(20),KEEP(38)) 
          CALL SMUMPS_SOLVE_INIT_OOC_BWD(PTRFAC,KEEP(28),MTYPE,
     &         I_WORKED_ON_ROOT, IROOT, A, LA, IERR)
        ENDIF
        IF ( KEEP( 50 ) .eq. 0 ) THEN
          MTYPE_LOC = MTYPE
        ELSE
          MTYPE_LOC = 0
        ENDIF
#if defined(V_T)
        CALL VTBEGIN(back_soln,ierr)
#endif
        IF (.NOT.FROM_PP) THEN
          CALL MUMPS_SECDEB(TIME_BWD)
        ENDIF
        IF ( .NOT.SPECIAL_ROOT_REACHED ) THEN
          RHS_ROOT(1:NRHS*SIZE_ROOT) = ZERO
        ENDIF
        IF (AM1.AND.(NB_FS_IN_RHSINTR_F.NE.NB_FS_IN_RHSINTR_TOT)) THEN
         DO I =1, N
           II = POSINRHSINTR_BWD(I)
           IF ((II.GT.0).AND.(II.GT.NB_FS_IN_RHSINTR_F)) THEN
            DO K=1,NRHS
             RHSINTR(II, K) = ZERO
            ENDDO
           ENDIF
         ENDDO
        ENDIF
        IF ( .NOT. DO_PRUN_BWD ) THEN
           IF ( .NOT. DO_L0OMP_BWD ) THEN
             IF (DO_L0OMP_FWD) THEN
               MYLEAF = -1
             ENDIF
           ENDIF
           IF ( DO_L0OMP_BWD ) THEN
             TO_PROCESS(:) = .TRUE.
             DO I=1, L_PHYS_L0_OMP
               TO_PROCESS( STEP(PHYS_L0_OMP( I )))
     &         = .FALSE.
             ENDDO
             IF (MYLEAF .EQ. -1) THEN
               MYLEAF = IPOOL_A_L0_OMP(LPOOL_A_L0_OMP)
             ENDIF
             CALL MUMPS_INIT_POOL_DIST_NA_BWD_L0( N, MYROOT, MYID_NODES,
     &          NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &          IW1(IPOOL), LPOOL, L0_OMP_MAPPING )
           ELSE
             CALL MUMPS_INIT_POOL_DIST_NA_BWD( N, MYROOT, MYID_NODES,
     &          NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &          IW1(IPOOL), LPOOL )
             IF (MYLEAF .EQ. -1) THEN
               CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &         NA(1), 
     &         NA(3), 
     &         MYLEAF, MYID_NODES, SLAVEF, KEEP, STEP,
     &         PROCNODE_STEPS )
             ENDIF
           ENDIF
        ELSE 
          IF ( DO_L0OMP_BWD ) THEN
             DO I=1, L_PHYS_L0_OMP
               IF ( TO_PROCESS( STEP(PHYS_L0_OMP( I ))) ) THEN
                 TO_PROCESS( STEP(PHYS_L0_OMP( I ))) = .FALSE.
                 PHYS_L0_OMP( I ) = -PHYS_L0_OMP( I )
               ENDIF
             ENDDO
             MYLEAF=0
             DO ISTEP = 1, KEEP(28)
               IF ( MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))
     &              .NE. MYID_NODES ) THEN
                  CYCLE
               ENDIF
               IF ( L0_OMP_MAPPING( ISTEP ) .NE. 0 ) THEN
                 CYCLE
               ENDIF
               IF ( .NOT. TO_PROCESS( ISTEP ) ) THEN
                 CYCLE
               ENDIF
               I      = Step2Node( ISTEP )
               ICHILD = FILS     ( I )
               DO WHILE ( ICHILD .GT. 0 )
                 ICHILD = FILS( ICHILD )
               END DO
               IF ( ICHILD .LT. 0 ) THEN
                 ICHILD = -ICHILD
                 DO WHILE ( ICHILD .GT. 0 )
                   IF ( L0_OMP_MAPPING( STEP( ICHILD ) ) .EQ. 0 .AND.
     &                  TO_PROCESS(STEP( ICHILD )) ) THEN
                     GOTO 10
                   ENDIF
                   ICHILD = FRERE( STEP( ICHILD ) )
                 ENDDO
               ENDIF
               MYLEAF = MYLEAF + 1
 10            CONTINUE
             ENDDO
             CALL MUMPS_INIT_POOL_DIST_NA_BWDL0ES( N, MYROOT,
     &          MYID_NODES,
     &          NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &          IW1(IPOOL), LPOOL, L0_OMP_MAPPING, TO_PROCESS )
          ELSE
            CALL MUMPS_INIT_POOL_DIST_BWD(N, nb_prun_roots,
     &      Pruned_Roots,
     &      MYROOT, MYID_NODES, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &      IW1(IPOOL), LPOOL)
            CALL MUMPS_NBLOCAL_ROOTS_OR_LEAVES( N,
     &       nb_prun_leaves, Pruned_Leaves,
     &       MYLEAF, MYID_NODES, SLAVEF, KEEP, STEP,
     &       PROCNODE_STEPS )
          ENDIF
        ENDIF
        IF ( DO_L0OMP_BWD
     &     ) THEN
          KEEP(31) =  1
        ELSE
          KEEP(31) = 0
        ENDIF
        IF (KEEP(31) .EQ. 1) THEN
          DO I = 1, KEEP(28)
            IF (MUMPS_PROCNODE(PROCNODE_STEPS(I),KEEP(199)) .EQ.
     &      MYID_NODES) THEN
              IF ( .NOT. MUMPS_IN_OR_ROOT_SSARBR(PROCNODE_STEPS(I),
     &        KEEP(199)) ) THEN
                IF ( L0_OMP_MAPPING(I) .EQ. 0 ) THEN
                  IF ( DO_PRUN_BWD
     &                 .OR. DO_L0OMP_BWD 
     &               ) THEN
                    IF ( TO_PROCESS(I) ) THEN
                      KEEP(31) = KEEP(31) + 1
                    ENDIF
                  ELSE
                    KEEP(31) = KEEP(31) + 1
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      CALL MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL, INFO, MYID )
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
      IF(INFO(1).LT.0) GOTO 500
        CALL SMUMPS_SOL_S( N, A, LA, IW, LIW, W(1), LWC,
     &          NRHS,
     &          RHSINTR, LRHSINTR, POSINRHSINTR_BWD,
     &          IW1(PTRICB),PTRACB,IWCB,LIWW, W2,
     &          NE_STEPS,
     &          STEP, FRERE,DAD,FILS,
     &          IW1(IPOOL),LPOOL,PTRIST,PTRFAC,MYLEAF,MYROOT,ICNTL,INFO,
     &          PROCNODE_STEPS, SLAVEF, COMM_NODES, MYID_NODES,
     &          BUFR, LBUFR, LBUFR_BYTES, KEEP, KEEP8, DKEEP,
     &          RHS_ROOT, LRHS_ROOT,
     &          MTYPE_LOC, 
     &          ISTEP_TO_INIV2, TAB_POS_IN_PERE, IW1(IPANEL_POS),
     &          LPANEL_POS, PRUN_BELOW_BWD, TO_PROCESS, SIZE_TO_PROCESS
     &        , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE_BWD
     &         , FROM_PP
     &          , L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &            L0_OMP_FACTORS, LL0_OMP_FACTORS
     &          )
        IF ( DO_L0OMP_BWD .AND. DO_PRUN_BWD ) THEN
          DO I = 1, L_PHYS_L0_OMP
            IF ( PHYS_L0_OMP( I ) .LT. 0 ) THEN
              PHYS_L0_OMP( I ) = -PHYS_L0_OMP( I )
              TO_PROCESS(STEP(PHYS_L0_OMP( I ) )) = .TRUE.
            ENDIF
          ENDDO
        ENDIF
        CALL MUMPS_STOP_ON_USER_REQUEST(KEEP, KEEP8, ICNTL, INFO, MYID)
        IF (DO_L0OMP_BWD .AND. INFO(1) .GE. 0) THEN
          KEEP(31) = 0
          PRUN_BELOW_BWD = DO_PRUN_BWD .AND. KEEP(111).EQ.0
          KEEP(405)=1
          CALL SMUMPS_SOL_L0OMP_S(N, MTYPE_LOC, NRHS, LIW, IW,
     &      IW1(PTRICB), PTRACB, RHSINTR, LRHSINTR, POSINRHSINTR_BWD,
     &      STEP, FRERE, FILS, NE_STEPS, PTRIST, PTRFAC, INFO,
     &      KEEP, KEEP8, DKEEP, PROCNODE_STEPS, SLAVEF,
     &      COMM_NODES, MYID_NODES, BUFR, LBUFR, LBUFR_BYTES, 
     &      RHS_ROOT, LRHS_ROOT, ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     &      IW1(IPANEL_POS), LPANEL_POS,
     &      PRUN_BELOW_BWD, TO_PROCESS, SIZE_TO_PROCESS,
     &      RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE_BWD,
     &      FROM_PP,
     &    LPOOL_B_L0_OMP,
     &    L_VIRT_L0_OMP, VIRT_L0_OMP,
     &    L_PHYS_L0_OMP, PHYS_L0_OMP,
     &    PERM_L0_OMP, PTR_LEAFS_L0_OMP,
     &    L0_OMP_MAPPING, LL0_OMP_MAPPING,
     &    L0_OMP_FACTORS, LL0_OMP_FACTORS )
          KEEP(405)=0
        ENDIF
        CALL SMUMPS_CLEAN_PENDING( INFO(1), KEEP,
     &     BUFR, LBUFR,LBUFR_BYTES,
     &     COMM_NODES, IDUMMY,       
     &     SLAVEF, .TRUE., .FALSE. ) 
        CALL MUMPS_STOP_ON_USER_REQUEST(KEEP, KEEP8, ICNTL, INFO, MYID)
        CALL MUMPS_PROPINFO(ICNTL, INFO, COMM_NODES, MYID )
#if defined(V_T)
        CALL VTEND(back_soln,ierr)
#endif
        IF (.NOT.FROM_PP) THEN
          CALL MUMPS_SECFIN(TIME_BWD)
          DKEEP(118)=real(TIME_BWD)+DKEEP(118)
        ENDIF
      ENDIF
      IF (LDIAG.GT.2 .AND. MP.GT.0) THEN
        IF (DOFORWARD) THEN
        K = min(10,size(RHSINTR,1))
        IF (LDIAG.EQ.4) K = size(RHSINTR,1)
        IF ( .NOT. FROM_PP) THEN
          WRITE (MP,99992)
          IF (size(RHSINTR,1).GT.0) 
     &    WRITE (MP,99993) (RHSINTR(I,1),I=1,K)
          IF (size(RHSINTR,1).GT.0.and.NRHS>1) 
     &              WRITE (MP,99994) (RHSINTR(I,2),I=1,K)
          ENDIF
        ENDIF
      ENDIF
500   CONTINUE
      IF ( allocated(TO_PROCESS)) DEALLOCATE (TO_PROCESS)
      IF (DO_PRUN_FWD.OR.DO_PRUN_BWD) THEN
         IF ( allocated(Pruned_Roots_FWD)) 
     &        DEALLOCATE (Pruned_Roots_FWD)
         IF ( allocated(Pruned_Roots_NS)) 
     &        DEALLOCATE (Pruned_Roots_NS)
         IF ( allocated(Pruned_SONS)) DEALLOCATE (Pruned_SONS)
         IF ( allocated(Pruned_Roots)) DEALLOCATE (Pruned_Roots)
         IF ( allocated(Pruned_List)) DEALLOCATE (Pruned_List)
         IF ( allocated(Pruned_Leaves)) DEALLOCATE (Pruned_Leaves)
      ENDIF
      RETURN 
99993 FORMAT (' RHS (internal, first column)'/(1X,1P,5E14.6))
99994 FORMAT (' RHS (internal, 2 nd  column)'/(1X,1P,5E14.6))
99992 FORMAT (//' LEAVING SOLVE (SMUMPS_SOL_C) WITH')
      END SUBROUTINE SMUMPS_SOL_C
      SUBROUTINE SMUMPS_SET_POSTPros (KEEP, ICNTL, NBRHS, MPG, PROKG, 
     &                               ICNTL10, ICNTL11, POSTPros)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: KEEP(500), ICNTL(60), NBRHS, MPG
      LOGICAL, INTENT(IN)  :: PROKG
      INTEGER, INTENT(OUT) :: ICNTL10, ICNTL11
      LOGICAL, INTENT(OUT) :: POSTPros
      POSTPros = .FALSE.
      IF (ICNTL11.NE.0 .OR. ICNTL10.NE.0) THEN
        POSTPros = .TRUE.
        IF (KEEP(111).NE.0) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING: Incompatible features: null space basis',
     &              ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ELSE IF ( KEEP(237) .NE.0 ) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING: Incompatible features: AM1',
     &              ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ELSE IF ( KEEP(252) .NE.0 ) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING: Incompatible features: Fwd in facto ',
     &              ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ELSE IF (KEEP(221).NE.0) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING: Incompatible features: reduced RHS',
     &     ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ELSE IF  (NBRHS.GT. 1 .OR. ICNTL(21) .GT. 0) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING:  Incompatible features: nrhs>1 or distrib sol',
     &     ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ELSE IF ( KEEP(248) .EQ. -1 ) THEN
            IF (PROKG) WRITE(MPG,'(A,A)')
     &     ' WARNING:  Incompatible features: distrib rhs',
     &     ' and Iter. Ref and/or Err. Anal.'
            POSTPros = .FALSE.
        ENDIF  
        IF (.NOT.POSTPros) THEN
          ICNTL11 = 0
          ICNTL10 = 0
        ENDIF
      ENDIF
      IF ((ICNTL(10) .NE. 0) .AND. (ICNTL10 .EQ. 0)) THEN
          IF (PROKG) WRITE(MPG,'(A)')
     &  ' WARNING: ICNTL(10) treated as if set to 0 '
      ENDIF
      IF ((ICNTL(11) .NE. 0)
     &      .AND.(ICNTL11 .EQ. 0)) THEN
          IF (PROKG) WRITE(MPG,'(A)')
     &  ' WARNING: ICNTL(11) treated as if set to 0 '
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SET_POSTPros
      SUBROUTINE SMUMPS_GATHER_SOLUTION( NSLAVES, N, MYID, COMM,
     &           NRHS,
     &           MTYPE, RHS, LRHS, NCOL_RHS, JBEG_RHS, PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, CWORK, LCWORK,
#if defined(USE_OLD_SCALING)
     &           LSCAL, SCALING, LSCALING, 
#else
     &           LSCAL, SCALING_LOC_BWD, LSCALING_LOC_BWD, 
#endif
     &           RHSINTR, LRHSINTR, NCOL_RHSINTR, 
     &           POSINRHSINTR, LPOS_N, PERM_RHS, SIZE_PERM_RHS )
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE, NCOL_RHS
      INTEGER NRHS, LRHS, LCWORK, LPOS_N, NCOL_RHSINTR
      REAL RHS   (LRHS, NCOL_RHS)
      INTEGER, INTENT(in) :: JBEG_RHS
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      REAL :: CWORK(LCWORK)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES
      INTEGER BUFFER(SIZE_BUF)
      INTEGER LRHSINTR, POSINRHSINTR(LPOS_N) 
#if defined(USE_OLD_SCALING)
      REAL, intent(in) :: RHSINTR(LRHSINTR, NCOL_RHSINTR)
      INTEGER, intent(in) :: LSCALING
      REAL, intent(in)    :: SCALING(LSCALING)
#else
      REAL, intent(inout) :: RHSINTR(LRHSINTR, NCOL_RHSINTR)
      INTEGER, intent(in) :: LSCALING_LOC_BWD
      REAL, intent(in) :: SCALING_LOC_BWD(LSCALING_LOC_BWD)
#endif
      LOGICAL, intent(in) :: LSCAL
      INTEGER, INTENT(in) :: SIZE_PERM_RHS
      INTEGER, INTENT(in) :: PERM_RHS(SIZE_PERM_RHS)
      INTEGER I, II, J, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL, N2RECV
      INTEGER LIELL, IPOS, NPIV, MAXNPIV_estim, MAXSurf
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR, allocok
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND, IPOSINRHSINTR
      INTEGER :: JCOL_RHS
      INTEGER :: K242
      LOGICAL :: OMP_FLAG
!$    INTEGER :: CHUNK, NOMP
      INTEGER, PARAMETER :: FIN = -1
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IROWlist
      INCLUDE 'mumps_headers.h'
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN
           IF (LSCAL) THEN 
             OMP_FLAG = .FALSE.
             IF (KEEP(350).EQ.2) THEN
               K242 = KEEP(242)
!$             NOMP = OMP_GET_MAX_THREADS()
!$             CHUNK = max(N/2,1)
!$             IF (int(NRHS,8) * int(N,8) .GE. int(KEEP(363),8)) THEN
!$               OMP_FLAG = .TRUE.
!$           CHUNK=int((int(N,8)*int(NRHS,8)+int(NOMP-1,8))/int(NOMP,8))
!$               CHUNK = min(CHUNK,(N+KEEP(362)-1)/KEEP(362))
!$               CHUNK = max(KEEP(363)/2,CHUNK)
!$             ENDIF
             ENDIF
             IF (OMP_FLAG) THEN
!$OMP PARALLEL FIRSTPRIVATE(JBEG_RHS,N,K242)
!$OMP&  PRIVATE(J,IPOSINRHSINTR,I,JCOL_RHS)
               DO J=1, NRHS
                 IF (K242.EQ.0) THEN
                   JCOL_RHS = J+JBEG_RHS-1
                 ELSE
                   JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
                 ENDIF
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                 DO I=1, N
                   IPOSINRHSINTR = POSINRHSINTR(I)
                   IF (IPOSINRHSINTR.GT.0) THEN 
                     RHS(I,JCOL_RHS) = RHSINTR(IPOSINRHSINTR,J)*
#if defined(USE_OLD_SCALING)
     &               SCALING(I)
#else
     &               SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
                   ELSE
                    RHS(I,JCOL_RHS) = ZERO
                  ENDIF
                 ENDDO
!$OMP END DO NOWAIT
               ENDDO
!$OMP END PARALLEL
             ELSE
               DO J=1, NRHS
                 IF (KEEP(242).EQ.0) THEN
                   JCOL_RHS = J+JBEG_RHS-1
                 ELSE
                   JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
                 ENDIF
                 DO I=1, N
                   IPOSINRHSINTR = POSINRHSINTR(I)
                   IF (IPOSINRHSINTR.GT.0) THEN 
                     RHS(I,JCOL_RHS) = RHSINTR(IPOSINRHSINTR,J)*
#if defined(USE_OLD_SCALING)
     &               SCALING(I)
#else
     &               SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
                   ELSE
                     RHS(I,JCOL_RHS) = ZERO
                   ENDIF
                 ENDDO
               ENDDO
             ENDIF
           ELSE
             OMP_FLAG = .FALSE.
             IF (KEEP(350).EQ.2) THEN
               K242 = KEEP(242)
!$             NOMP = OMP_GET_MAX_THREADS()
!$             OMP_FLAG = .FALSE.
!$             CHUNK = max(N/2,1)
!$             IF (NRHS * N .GE. KEEP(363)) THEN
!$               OMP_FLAG = .TRUE.
!$         CHUNK=int((int(N,8)*int(NRHS,8)+int(NOMP-1,8))/int(NOMP,8))
!$               CHUNK = min(CHUNK,(N+KEEP(362)-1)/KEEP(362))
!$               CHUNK = max(KEEP(363)/2,CHUNK)
!$             ENDIF
             ENDIF
             IF (OMP_FLAG) THEN
!$OMP PARALLEL FIRSTPRIVATE(JBEG_RHS,N,K242)
!$OMP&  PRIVATE(IPOSINRHSINTR,I,JCOL_RHS) IF (OMP_FLAG)
               DO J=1, NRHS
                 IF (K242.EQ.0) THEN
                   JCOL_RHS = J+JBEG_RHS-1
                 ELSE
                   JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
                 ENDIF
!$OMP DO SCHEDULE(DYNAMIC,CHUNK)
                 DO I=1, N
                   IPOSINRHSINTR = POSINRHSINTR(I)
                   IF (IPOSINRHSINTR.GT.0) THEN
                     RHS(I,JCOL_RHS) = RHSINTR(IPOSINRHSINTR,J)
                   ELSE
                     RHS(I,JCOL_RHS) = ZERO
                   ENDIF
                 ENDDO
!$OMP END DO NOWAIT
               ENDDO
!$OMP END PARALLEL
             ELSE
               DO J=1, NRHS
                 IF (KEEP(242).EQ.0) THEN
                   JCOL_RHS = J+JBEG_RHS-1
                 ELSE
                   JCOL_RHS = PERM_RHS(J+JBEG_RHS-1)
                 ENDIF
                 DO I=1, N
                   IPOSINRHSINTR = POSINRHSINTR(I)
                   IF (IPOSINRHSINTR.GT.0) THEN
                     RHS(I,JCOL_RHS) = RHSINTR(IPOSINRHSINTR,J)
                   ELSE
                     RHS(I,JCOL_RHS) = ZERO
                   ENDIF
                 ENDDO
               ENDDO
             ENDIF
           ENDIF
        RETURN
      ENDIF
      MAXNPIV_estim = max(KEEP(246), KEEP(247))
      MAXSurf       = MAXNPIV_estim*NRHS
      IF (LCWORK .LT. MAXNPIV_estim) THEN
        WRITE(*,*) MYID, 
     &  ": Internal error 2 in SMUMPS_GATHER_SOLUTION:",
     &  TYPE_PARAL, LCWORK, KEEP(247), NRHS
        CALL MUMPS_ABORT()
      ENDIF
      IF (MYID.EQ.MASTER) THEN
         ALLOCATE(IROWlist(KEEP(247)),stat=allocok)
         IF(allocok.GT.0) THEN
            WRITE(*,*)'Problem with allocation of IROWlist'
            CALL MUMPS_ABORT()
         ENDIF
      ENDIF
      IF (NSLAVES .EQ. 1 .AND. TYPE_PARAL .EQ. 1) THEN
         CALL MUMPS_ABORT()
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(MAXNPIV_estim+2,MPI_INTEGER, COMM, 
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(MAXSurf,MPI_REAL, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  SMUMPS_GATHER_SOLUTION '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =N
      POS_BUF  =0
      IF (I_AM_SLAVE) THEN
        POS_BUF = 0
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &          KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N)
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              IF (MYID .EQ. MASTER) THEN
                   N2RECV=N2RECV-NPIV
                   IF (NPIV.GT.0) 
     &             CALL SMUMPS_NPIV_BLOCK_ADD ( .TRUE. )
              ELSE
                   IF (NPIV.GT.0) 
     &             CALL SMUMPS_NPIV_BLOCK_ADD ( .FALSE.)
              ENDIF
          ENDIF
        ENDDO
        CALL SMUMPS_NPIV_BLOCK_SEND()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (NPIV.NE.FIN)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             IROWlist, NPIV, MPI_INTEGER, COMM, IERR)
          DO J=1, NRHS
            IF (KEEP(242).EQ.0) THEN
              JCOL_RHS=J+JBEG_RHS-1
            ELSE
              JCOL_RHS=PERM_RHS(J+JBEG_RHS-1)
            ENDIF
            CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                 CWORK, NPIV, MPI_REAL,
     &                 COMM, IERR)
#if defined(USE_OLD_SCALING)
            IF (LSCAL) THEN
              DO I=1,NPIV
                RHS(IROWlist(I),JCOL_RHS)=CWORK(I)*SCALING(IROWlist(I))
              ENDDO
            ELSE
#else
#endif
              DO I=1,NPIV
                RHS(IROWlist(I),JCOL_RHS)=CWORK(I)
              ENDDO
#if defined(USE_OLD_SCALING)
            ENDIF
#endif
          ENDDO
          N2RECV=N2RECV-NPIV
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       DEALLOCATE(IROWlist)
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE SMUMPS_NPIV_BLOCK_ADD ( ON_MASTER )
        LOGICAL, intent(in) ::  ON_MASTER     
        INTEGER :: JPOS, K242
        LOGICAL :: PRIV_LSCAL
        IF (ON_MASTER) THEN
        IF (KEEP(350).EQ.2
     &  .AND. (NRHS.EQ.1.OR.((NPIV*NRHS*2*KEEP(16)).GE.KEEP(364)))) THEN
          PRIV_LSCAL = LSCAL
          K242 = KEEP(242)
          DO J=1, NRHS
            IF (K242.EQ.0) THEN
              JPOS = J+JBEG_RHS-1
            ELSE
              JPOS = PERM_RHS(J+JBEG_RHS-1)
            ENDIF
            DO II=1,NPIV
              I=IW(J1+II-1)
              IPOSINRHSINTR= POSINRHSINTR(I) 
              IF (PRIV_LSCAL) THEN
                RHS(I,JPOS) = RHSINTR(IPOSINRHSINTR,J)*
#if defined(USE_OLD_SCALING)
     &          SCALING(I)
#else
     &          SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
              ELSE
                RHS(I,JPOS) = RHSINTR(IPOSINRHSINTR,J)
              ENDIF
            ENDDO
          ENDDO
        ELSE
         IF (KEEP(242).EQ.0) THEN
           IF (LSCAL) THEN
             DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSINTR= POSINRHSINTR(I) 
                DO J=1, NRHS
                  RHS(I,J+JBEG_RHS-1) = RHSINTR(IPOSINRHSINTR,J)*
#if defined(USE_OLD_SCALING)
     &            SCALING(I)
#else
     &            SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
                ENDDO
             ENDDO
           ELSE
             DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSINTR= POSINRHSINTR(I)
                DO J=1, NRHS
                  RHS(I,J+JBEG_RHS-1) = RHSINTR(IPOSINRHSINTR,J)
                ENDDO
             ENDDO
           ENDIF
         ELSE
           IF (LSCAL) THEN
              DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSINTR= POSINRHSINTR(I)
!DIR$ NOVECTOR
                DO J=1, NRHS
                  RHS(I,PERM_RHS(J+JBEG_RHS-1)) =
     &               RHSINTR(IPOSINRHSINTR,J)*
#if defined(USE_OLD_SCALING)
     &               SCALING(I)
#else
     &               SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
                ENDDO
              ENDDO
           ELSE
              DO II=1,NPIV
                I=IW(J1+II-1)
                IPOSINRHSINTR= POSINRHSINTR(I)
!DIR$ NOVECTOR
                DO J=1, NRHS
                  RHS(I,PERM_RHS(J+JBEG_RHS-1)) =
     &               RHSINTR(IPOSINRHSINTR,J)
                ENDDO
              ENDDO
           ENDIF
         ENDIF
        ENDIF
         RETURN
        ENDIF
        CALL MPI_PACK(NPIV, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(IW(J1), NPIV, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        IPOSINRHSINTR= POSINRHSINTR(IW(J1)) 
        DO J=1,NRHS
#if ! defined(USE_OLD_SCALING)
          IF (LSCAL) THEN
            DO II=IPOSINRHSINTR, IPOSINRHSINTR+NPIV-1
              RHSINTR(II,J)=
     &        RHSINTR(II,J)*SCALING_LOC_BWD(II)
            ENDDO
          ENDIF
#endif
          CALL MPI_PACK(RHSINTR(IPOSINRHSINTR,J), NPIV,
     &              MPI_REAL,
     &              BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &              IERR)
        ENDDO
        N2SEND=N2SEND+NPIV  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL SMUMPS_NPIV_BLOCK_SEND()
        END IF
        RETURN
        END SUBROUTINE SMUMPS_NPIV_BLOCK_ADD
        SUBROUTINE SMUMPS_NPIV_BLOCK_SEND()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE SMUMPS_NPIV_BLOCK_SEND
      END SUBROUTINE SMUMPS_GATHER_SOLUTION
      SUBROUTINE SMUMPS_GATHER_SOLUTION_AM1(NSLAVES, N, MYID, COMM,
     &           NRHS, RHSINTR,  LRHSINTR, NRHSINTR_COL,
     &           KEEP, BUFFER,
     &           SIZE_BUF, SIZE_BUF_BYTES, 
#if defined(USE_OLD_SCALING)
     &           LSCAL, SCALING, LSCALING,
#else
     &           LSCAL, SCALING_LOC_BWD, LSCALING_LOC_BWD, 
#endif
     &          IRHS_PTR_COPY, LIRHS_PTR_COPY, 
     &          IRHS_SPARSE_COPY, LIRHS_SPARSE_COPY,
     &          RHS_SPARSE_COPY, LRHS_SPARSE_COPY,
     &          UNS_PERM_INV, LUNS_PERM_INV,
     &          POSINRHSINTR, LPOS_ROW, NB_FS_IN_RHSINTR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM
      INTEGER NRHS, LRHSINTR, NRHSINTR_COL 
      REAL, intent(in) :: RHSINTR (LRHSINTR, NRHSINTR_COL)
      INTEGER KEEP(500)
      INTEGER SIZE_BUF, SIZE_BUF_BYTES, LPOS_ROW
      INTEGER BUFFER(SIZE_BUF)
      INTEGER, intent(in) :: LIRHS_PTR_COPY, LIRHS_SPARSE_COPY, 
     &                       LRHS_SPARSE_COPY, LUNS_PERM_INV, 
     &                       NB_FS_IN_RHSINTR
      INTEGER :: IRHS_SPARSE_COPY(LIRHS_SPARSE_COPY), 
     &           IRHS_PTR_COPY(LIRHS_PTR_COPY), 
     &           UNS_PERM_INV(LUNS_PERM_INV), 
     &           POSINRHSINTR(LPOS_ROW) 
      REAL :: RHS_SPARSE_COPY(LRHS_SPARSE_COPY)
#if defined(USE_OLD_SCALING)
      INTEGER, intent(in) :: LSCALING
      REAL, intent(in)    :: SCALING(LSCALING)
#else
      INTEGER, intent(in) :: LSCALING_LOC_BWD
      REAL, intent(in) :: SCALING_LOC_BWD(LSCALING_LOC_BWD)
#endif
      LOGICAL, intent(in) :: LSCAL
      INTEGER COLSIZE, K, IZ, IPREV, NBCOL_INBLOC
      INTEGER I, II, J, MASTER,
     &         TYPE_PARAL, N2RECV, IPOSINRHSINTR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
      INTEGER RECORD_SIZE_P_1, SIZE1, SIZE2
      INTEGER POS_BUF, N2SEND
      INTEGER, PARAMETER :: FIN = -1
      INCLUDE 'mumps_headers.h'
      TYPE_PARAL = KEEP(46)  
      I_AM_SLAVE = MYID .ne. MASTER .OR. TYPE_PARAL .eq. 1
      NBCOL_INBLOC = size(IRHS_PTR_COPY)-1
      IF (NSLAVES.EQ.1 .AND. TYPE_PARAL.EQ.1) THEN 
        K=1              
        DO J = 1,  NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IPOSINRHSINTR = POSINRHSINTR(I)
             IF (IPOSINRHSINTR.GT.0) THEN 
                IF (LSCAL) THEN
                  RHS_SPARSE_COPY(IZ)=
     &                      RHSINTR(IPOSINRHSINTR,K)
#if defined(USE_OLD_SCALING)
     &                      * SCALING(I)
#else
     &                      * SCALING_LOC_BWD(IPOSINRHSINTR)
#endif
                ELSE
                  RHS_SPARSE_COPY(IZ)=RHSINTR(IPOSINRHSINTR,K)
                ENDIF
             ENDIF
           ENDDO
           K = K + 1
        ENDDO
        RETURN
      ENDIF
      IF (I_AM_SLAVE) THEN
        K=1              
        DO J = 1, NBCOL_INBLOC
           COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
           IF (COLSIZE.EQ.0) CYCLE
           DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
             IPOSINRHSINTR = POSINRHSINTR(I)
             IF (IPOSINRHSINTR.GT.0) THEN 
#if ! defined(USE_OLD_SCALING)
               IF (LSCAL) THEN
                 RHS_SPARSE_COPY(IZ)=RHSINTR(IPOSINRHSINTR,K)
     &           * SCALING_LOC_BWD(IPOSINRHSINTR)
               ELSE
#endif
                 RHS_SPARSE_COPY(IZ)=RHSINTR(IPOSINRHSINTR,K)
#if ! defined(USE_OLD_SCALING)
               ENDIF
#endif
             ENDIF          
           ENDDO
           K = K + 1
        ENDDO
      ENDIF
      SIZE1=0
      CALL MPI_PACK_SIZE(3,MPI_INTEGER, COMM,  
     &          SIZE1, IERR)
      SIZE2=0
      CALL MPI_PACK_SIZE(1,MPI_REAL, COMM,
     &                   SIZE2, IERR)
      RECORD_SIZE_P_1= SIZE1+SIZE2
      IF (RECORD_SIZE_P_1.GT.SIZE_BUF_BYTES) THEN
         write(6,*) MYID, 
     &    ' Internal error 3 in  SMUMPS_GATHER_SOLUTION_AM1 '
         write(6,*) MYID, ' RECORD_SIZE_P_1, SIZE_BUF_BYTES=', 
     &                 RECORD_SIZE_P_1, SIZE_BUF_BYTES
         CALL MUMPS_ABORT()
      ENDIF
      N2SEND   =0
      N2RECV   =size(IRHS_SPARSE_COPY)
      POS_BUF  =0
      IF (I_AM_SLAVE) THEN
        DO J = 1,  NBCOL_INBLOC
            COLSIZE = IRHS_PTR_COPY(J+1) - IRHS_PTR_COPY(J)
            IF (COLSIZE.LE.0) CYCLE
            K = 0 
            DO IZ=IRHS_PTR_COPY(J), IRHS_PTR_COPY(J+1)-1
             I = IRHS_SPARSE_COPY(IZ)
             II = I
             IF  (KEEP(23).NE.0) II = UNS_PERM_INV(I)
             IPOSINRHSINTR = POSINRHSINTR(II)
             IF (IPOSINRHSINTR.GT.0) THEN
               IF (MYID .EQ. MASTER) THEN
                  N2RECV=N2RECV-1
#if defined(USE_OLD_SCALING)
                  IF (LSCAL) 
     &            CALL SMUMPS_AM1_BLOCK_ADD ( .TRUE. )
#endif
                  IRHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &               I
                  RHS_SPARSE_COPY( IRHS_PTR_COPY(J) + K) =
     &                RHS_SPARSE_COPY(IZ)
                  K = K+1 
               ELSE
#if defined(USE_OLD_SCALING)
                  CALL SMUMPS_AM1_BLOCK_ADD ( .FALSE. )
#else
                  CALL SMUMPS_AM1_BLOCK_ADD ()
#endif
               ENDIF
              ENDIF          
            ENDDO
            IF (MYID.EQ.MASTER) 
     &             IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + K
        ENDDO
        CALL SMUMPS_AM1_BLOCK_SEND()   
      ENDIF
      IF ( MYID .EQ. MASTER ) THEN
       DO WHILE (N2RECV .NE. 0)
        CALL MPI_RECV( BUFFER, SIZE_BUF_BYTES, MPI_PACKED,
     &                 MPI_ANY_SOURCE,
     &                 GatherSol, COMM, STATUS, IERR )
        POS_BUF = 0
        CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        DO WHILE (J.NE.FIN)
          IZ = IRHS_PTR_COPY(J)
          CALL MPI_UNPACK( BUFFER,SIZE_BUF_BYTES, POS_BUF,
     &             I, 1, MPI_INTEGER, COMM, IERR)
          IRHS_SPARSE_COPY(IZ) = I
          CALL MPI_UNPACK(BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &             RHS_SPARSE_COPY(IZ), 1, MPI_REAL,
     &             COMM, IERR)
#if defined(USE_OLD_SCALING)
          IF (LSCAL) THEN
              IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
              RHS_SPARSE_COPY(IZ) = RHS_SPARSE_COPY(IZ)*SCALING(I)    
          ENDIF
#endif
          N2RECV=N2RECV-1
          IRHS_PTR_COPY(J) = IRHS_PTR_COPY(J) + 1
          CALL MPI_UNPACK( BUFFER, SIZE_BUF_BYTES, POS_BUF,
     &                   J, 1, MPI_INTEGER, COMM, IERR)
        ENDDO
       ENDDO
       IPREV = 1
       DO J=1, size(IRHS_PTR_COPY)-1
         I= IRHS_PTR_COPY(J) 
         IRHS_PTR_COPY(J) = IPREV
         IPREV = I
       ENDDO
      ENDIF
      RETURN
      CONTAINS
        SUBROUTINE SMUMPS_AM1_BLOCK_ADD (
#if defined(USE_OLD_SCALING)
     &                                    SCALE_ONLY 
#endif
     &                                               )
#if defined(USE_OLD_SCALING)
        LOGICAL, intent(in) ::  SCALE_ONLY    
#endif
#if defined(USE_OLD_SCALING)
        INTEGER III
#endif
#if defined(USE_OLD_SCALING)
        IF (SCALE_ONLY) THEN
          WRITE(*,*) "SMUMPS_AM1_BLOCK_ADD(true) should not be called"
          CALL MUMPS_ABORT()
         III = I
         IF (KEEP(23).NE.0) III = UNS_PERM_INV(I)
         IF (LSCAL) THEN
            RHS_SPARSE_COPY(IZ)=RHS_SPARSE_COPY(IZ)*SCALING(III)
         ENDIF
         RETURN
        ENDIF
#endif
        CALL MPI_PACK(J, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(I, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
        CALL MPI_PACK(RHS_SPARSE_COPY(IZ), 1, MPI_REAL,
     &                BUFFER, SIZE_BUF_BYTES, POS_BUF, COMM,
     &                IERR)
        N2SEND=N2SEND+1  
        IF ( POS_BUF + RECORD_SIZE_P_1 > SIZE_BUF_BYTES ) THEN
          CALL SMUMPS_AM1_BLOCK_SEND()
        END IF
        RETURN
        END SUBROUTINE SMUMPS_AM1_BLOCK_ADD
        SUBROUTINE SMUMPS_AM1_BLOCK_SEND()
        IF (N2SEND .NE. 0) THEN
         CALL MPI_PACK(FIN, 1, MPI_INTEGER, BUFFER,
     &                SIZE_BUF_BYTES, POS_BUF, COMM, IERR )
         CALL MPI_SEND(BUFFER, POS_BUF, MPI_PACKED, MASTER, 
     &                 GatherSol, COMM, IERR)
        ENDIF
        POS_BUF=0
        N2SEND=0
        RETURN
        END SUBROUTINE SMUMPS_AM1_BLOCK_SEND
      END SUBROUTINE SMUMPS_GATHER_SOLUTION_AM1
      SUBROUTINE SMUMPS_DISTSOL_INDICES(MTYPE, ISOL_LOC,
     &             PTRIST, KEEP,KEEP8,
     &             IW, LIW_PASSED, MYID_NODES, N, STEP,
     &             PROCNODE, NSLAVES,
#if defined(USE_OLD_SCALING)
     &             scaling_data, LSCAL,
#endif
     &             IRHS_loc_MEANINGFUL, IRHS_loc, Nloc_RHS
     &             )
      IMPLICIT NONE
      INTEGER MTYPE, MYID_NODES, N, NSLAVES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE(KEEP(28))
      INTEGER ISOL_LOC(KEEP(89))
      INTEGER LIW_PASSED
      INTEGER IW(LIW_PASSED)
      INTEGER STEP(N)
#if defined(USE_OLD_SCALING)
      LOGICAL LSCAL
#endif
      LOGICAL :: IRHS_loc_MEANINGFUL
      INTEGER :: Nloc_RHS
      INTEGER :: IRHS_loc(Nloc_RHS)
#if defined(USE_OLD_SCALING)
      type scaling_data_t
        SEQUENCE
        REAL, dimension(:), pointer :: SCALING
        REAL, dimension(:), pointer :: SCALING_LOC
        INTEGER, dimension(:), pointer :: SCALING_IND 
      end type scaling_data_t
      type (scaling_data_t) :: scaling_data
#endif
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER ISTEP, K
      INTEGER J1, IPOS, LIELL, NPIV, JJ
      LOGICAL :: CHECK_IRHS_loc
      INTEGER(8) :: DIFF_ADDR
      INCLUDE 'mumps_headers.h'
      CHECK_IRHS_loc=.FALSE.
      IF ( IRHS_loc_MEANINGFUL ) THEN
        IF (Nloc_RHS .GT. 0) THEN
          CALL MUMPS_SIZE_C( IRHS_loc(1), ISOL_loc(1),
     &                       DIFF_ADDR )
          IF (DIFF_ADDR .EQ. 0_8) THEN
            CHECK_IRHS_loc=.TRUE.
          ENDIF
        ENDIF
      ENDIF
      K=0
      DO ISTEP=1, KEEP(28)
          IF ( MYID_NODES == MUMPS_PROCNODE( PROCNODE(ISTEP),
     &                   KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW_PASSED, PTRIST, STEP, N)
              IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
                   J1=IPOS+1+LIELL
              ELSE
                   J1=IPOS+1
              END IF
              DO JJ=J1,J1+NPIV-1
                  K=K+1
                  IF (CHECK_IRHS_loc) THEN
                    IF (K.LE.Nloc_RHS) THEN
                      IF ( IW(JJ) .NE.IRHS_LOC(K) ) THEN
                      ENDIF
                    ENDIF
                  ENDIF
                  ISOL_LOC(K)=IW(JJ)
#if defined(USE_OLD_SCALING)
                  IF (LSCAL) THEN
                    scaling_data%SCALING_LOC(K)=
     &              scaling_data%SCALING(IW(JJ))
                  ENDIF
#endif
              ENDDO
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_DISTSOL_INDICES
#if ! defined(USE_OLD_SCALING)
      SUBROUTINE SMUMPS_SCALINGRHSINTR(LSCAL, N,
     &           SCALING_LOC, SCALING_RHSINTR,
     &           L, POSINRHSINTR, KEEP, ROWORCOL, PTRIST,
     &           IW, LIW_PASSED, MYID_NODES, STEP,
     &           PROCNODE, NSLAVES)
      IMPLICIT NONE
      INTEGER              :: KEEP(500)
      LOGICAL, INTENT(IN)  :: LSCAL
      INTEGER, INTENT(IN)  :: N, L
      INTEGER, INTENT(IN)  :: POSINRHSINTR(N)
      REAL   , INTENT(IN)  :: SCALING_LOC(max(KEEP(89),1))
      REAL   , INTENT(OUT) :: SCALING_RHSINTR(L)
      INTEGER, INTENT(IN)  :: ROWORCOL, NSLAVES, LIW_PASSED, MYID_NODES
      INTEGER, INTENT(IN)  :: STEP(KEEP(28)),
     &                        PROCNODE(KEEP(28)),
     &                        PTRIST(KEEP(28)),
     &                        IW(LIW_PASSED)
      INTEGER :: IPOSINRHSINTR
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: ISTEP
      INTEGER :: KLOC, J1, JJ, LIELL, IPOS, NPIV
      IF (.NOT. LSCAL) THEN
        WRITE(*,*) "Internal error 1 in SMUMPS_DS_SCALINGRHSINTR"
        CALL MUMPS_ABORT()
      ENDIF
      IF (ROWORCOL .NE. 1 .AND. ROWORCOL.NE.2) THEN
        WRITE(*,*) "Internal error 2 in SMUMPS_DS_SCALINGRHSINTR",
     &  ROWORCOL
      ENDIF
      IF (KEEP(89).EQ.0) RETURN
      KLOC = 1
      DO ISTEP=1, KEEP(28)
        IF ( MYID_NODES == MUMPS_PROCNODE( PROCNODE(ISTEP),
     &               KEEP(199))) THEN
          CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &    NPIV, LIELL, IPOS, IW, LIW_PASSED, PTRIST, STEP, N)
          IF (ROWORCOL .EQ. 1) THEN
               J1=IPOS+1
          ELSE
               J1=IPOS+1+LIELL
          END IF
          IPOSINRHSINTR = POSINRHSINTR(IW(J1)) 
          IF ( IPOSINRHSINTR .GT. 0 ) THEN
            DO JJ=1, NPIV
             SCALING_RHSINTR(IPOSINRHSINTR+JJ-1) =
     &       SCALING_LOC(KLOC+JJ-1)
            ENDDO
          ENDIF
          KLOC = KLOC + NPIV
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_SCALINGRHSINTR
#endif
      SUBROUTINE SMUMPS_DISTRIBUTED_SOLUTION(
     &           SLAVEF, N, MYID_NODES,
     &           MTYPE, RHSINTR, LRHSINTR, NBRHS_EFF,
     &           POSINRHSINTR,
     &           ISOL_LOC, 
     &           SOL_LOC, NRHS, BEG_RHS, LSOL_LOC,
     &           PTRIST,
     &           PROCNODE_STEPS, KEEP,KEEP8, IW, LIW, STEP,
     &           SCALING_LOC_BWD, LSCALING_LOC_BWD,
     &           LSCAL, NB_RHSSKIPPED,
     &           PERM_RHS, SIZE_PERM_RHS)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      LOGICAL, intent(in) :: LSCAL
      INTEGER, intent(in) :: LSCALING_LOC_BWD
      REAL, intent(in) :: SCALING_LOC_BWD(LSCALING_LOC_BWD)
      INTEGER SLAVEF, N, MYID_NODES, LIW, MTYPE, NBRHS_EFF, LRHSINTR
      INTEGER POSINRHSINTR(N), NB_RHSSKIPPED
      INTEGER LSOL_LOC, BEG_RHS
      INTEGER ISOL_LOC(LSOL_LOC)
      INTEGER, INTENT(in) :: NRHS 
      REAL SOL_LOC( LSOL_LOC, NRHS )
      REAL RHSINTR( LRHSINTR, NBRHS_EFF )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N)
      INTEGER, INTENT(in) :: SIZE_PERM_RHS
      INTEGER, INTENT(in) :: PERM_RHS( SIZE_PERM_RHS )
      INTEGER :: JJ, J1, ISTEP, K, KLOC, IPOSINRHSINTR, JEMPTY
      INTEGER :: JCOL, JCOL_PERM
      INTEGER :: IPOS, LIELL, NPIV, JEND
      LOGICAL :: ROOT
!$    LOGICAL :: OMP_FLAG
      REAL, PARAMETER :: ZERO = 0.0E0
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      K=0
      JEMPTY = BEG_RHS+NB_RHSSKIPPED-1
      JEND   = BEG_RHS+NB_RHSSKIPPED+NBRHS_EFF-1
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &      KEEP(199))) THEN
          ROOT=.false.
          IF (KEEP(38).ne.0) ROOT = STEP(KEEP(38))==ISTEP
          IF (KEEP(20).ne.0) ROOT = STEP(KEEP(20))==ISTEP
          IF ( ROOT ) THEN
                IPOS = PTRIST(ISTEP) + KEEP(IXSZ)
                LIELL = IW(IPOS+3)
                NPIV = LIELL
                IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
          ELSE
              IPOS = PTRIST(ISTEP) + 2 +KEEP(IXSZ)
              LIELL = IW(IPOS-2)+IW(IPOS+1)
              IPOS= IPOS+1
              NPIV = IW(IPOS)
              IPOS= IPOS+1
              IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
          END IF
          IF (MTYPE.eq.1 .AND. KEEP(50).EQ.0) THEN
               J1=IPOS+1+LIELL
          ELSE
               J1=IPOS+1
          END IF
            IF (NB_RHSSKIPPED.GT.0) THEN
              DO JCOL = BEG_RHS, JEMPTY
                IF (KEEP(242) .NE. 0) THEN
                JCOL_PERM = PERM_RHS(JCOL)
                ELSE
                 JCOL_PERM = JCOL
                ENDIF
                KLOC=K
                DO JJ=J1,J1+NPIV-1
                  KLOC=KLOC+1
                  SOL_LOC(KLOC, JCOL_PERM) = ZERO
                ENDDO
              ENDDO
            ENDIF
!$            OMP_FLAG = ( JEND-JEMPTY.GE.KEEP(362) .AND.
!$   &                     (NPIV*(JEND-JEMPTY) .GE. KEEP(363)/2 ) )
!$OMP PARALLEL DO PRIVATE(JCOL,JCOL_PERM,KLOC,JJ,IPOSINRHSINTR) 
!$OMP&         IF(OMP_FLAG)
            DO JCOL = JEMPTY+1, JEND
              IF (KEEP(242) .NE. 0) THEN
                JCOL_PERM = PERM_RHS(JCOL)
              ELSE
                 JCOL_PERM = JCOL
              ENDIF
              DO JJ=J1,J1+NPIV-1
                 KLOC=K + JJ-J1 + 1
                IF (LSCAL) THEN
                  SOL_LOC(KLOC,JCOL_PERM) =
     &            SCALING_LOC_BWD(KLOC)*
     &            RHSINTR(KLOC,JCOL-JEMPTY)
                ELSE
                  SOL_LOC(KLOC,JCOL_PERM) =
     &            RHSINTR(KLOC,JCOL-JEMPTY)
                ENDIF
              ENDDO
            ENDDO
!$OMP END PARALLEL DO
          K=K+NPIV
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_DISTRIBUTED_SOLUTION
      SUBROUTINE SMUMPS_SCATTER_RHS
     &           (NSLAVES, N, MYID, COMM,
     &           LSCAL, SCALING_LOC_FWD,
     &           MTYPE, RHS, LRHS, NCOL_RHS, NRHS,
     &           RHSINTR, LRHSINTR, NCOL_RHSINTR,
     &           POSINRHSINTR_FWD, NB_FS_IN_RHSINTR_F,
     &           PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           ICNTL, INFO)
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER NSLAVES, N, MYID, COMM, LIW, MTYPE
      INTEGER NRHS, LRHS, NCOL_RHS,  LRHSINTR, NCOL_RHSINTR
      INTEGER ICNTL(60), INFO(80)
      REAL, intent(in)  :: RHS (LRHS, NCOL_RHS)
      REAL, intent(out) :: RHSINTR(LRHSINTR, NCOL_RHSINTR)
      INTEGER, intent(in)  :: POSINRHSINTR_FWD(N), NB_FS_IN_RHSINTR_F
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER IW(LIW), STEP(N) 
      LOGICAL, intent(in) :: LSCAL
      REAL, intent(in) :: SCALING_LOC_FWD(max(1,KEEP(89)))
      INTEGER BUF_MAXSIZE, BUF_MAXREF
      PARAMETER (BUF_MAXREF=200000)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUF_INDX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BUF_RHS
      REAL, ALLOCATABLE, DIMENSION(:) :: BUF_RHS_2
      INTEGER ENTRIES_2_PROCESS, PROC_WHO_ASKS, BUF_EFFSIZE
      INTEGER INDX 
      INTEGER allocok
      REAL ZERO
      PARAMETER( ZERO = 0.0E0 )
      INTEGER I, J, K, JJ, J1, ISTEP, MASTER,
     &        MYID_NODES, TYPE_PARAL
      INTEGER LIELL, IPOS, NPIV
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      PARAMETER(MASTER=0)
      LOGICAL I_AM_SLAVE
!$    INTEGER :: CHUNK, NOMP
!$    LOGICAL :: OMP_FLAG
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      TYPE_PARAL = KEEP(46)
      I_AM_SLAVE = MYID .ne. 0 .OR. TYPE_PARAL .eq. 1
      IF ( TYPE_PARAL == 1 ) THEN
        MYID_NODES = MYID
      ELSE
        MYID_NODES = MYID-1
      ENDIF
      BUF_EFFSIZE = 0
      BUF_MAXSIZE = max(min(BUF_MAXREF,int(2000000/NRHS)),2000)
      IF ( KEEP(350).EQ.2 ) THEN
!$       NOMP = OMP_GET_MAX_THREADS()
         ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     &        BUF_RHS_2(BUF_MAXSIZE*NRHS),
     &        stat=allocok)
      ELSE
         ALLOCATE (BUF_INDX(BUF_MAXSIZE),
     &        BUF_RHS(NRHS,BUF_MAXSIZE),
     &        stat=allocok)
      END IF
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)=BUF_MAXSIZE*(NRHS+1)
      ENDIF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM, MYID )
      IF (INFO(1).LT.0) RETURN
      IF (MYID.EQ.MASTER) THEN
        ENTRIES_2_PROCESS = N - KEEP(89)
        IF (TYPE_PARAL.EQ.1.AND.ENTRIES_2_PROCESS.NE.0) THEN
         IF (NB_FS_IN_RHSINTR_F.LT.LRHSINTR) THEN
           DO K=1, NCOL_RHSINTR
            DO I = NB_FS_IN_RHSINTR_F +1, LRHSINTR
             RHSINTR (I, K) = ZERO
            ENDDO
           ENDDO
         ENDIF
        ENDIF
        IF ( KEEP(350).EQ.2 ) THEN
           DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
              CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     &             MPI_ANY_SOURCE,
     &             ScatterRhsI, COMM, STATUS, IERR )
              CALL MPI_GET_COUNT(STATUS,MPI_INTEGER,BUF_EFFSIZE,IERR)
              PROC_WHO_ASKS = STATUS(MPI_SOURCE)
!$            OMP_FLAG = .FALSE.
!$            CHUNK = NRHS
!$            IF (BUF_EFFSIZE*NRHS .GE. KEEP(363)) THEN
!$              OMP_FLAG = .TRUE.
!$              CHUNK = max((BUF_EFFSIZE*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$            ENDIF
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(INDX)
!$OMP&  IF (OMP_FLAG)
              DO K = 1, NRHS
                 DO I = 1, BUF_EFFSIZE
                    INDX = BUF_INDX( I )
                    BUF_RHS_2( I+(K-1)*BUF_EFFSIZE) = RHS( INDX, K )
                 ENDDO
              ENDDO
!$OMP END PARALLEL DO
              CALL MPI_SEND( BUF_RHS_2,
     &             NRHS*BUF_EFFSIZE,
     &             MPI_REAL, PROC_WHO_ASKS,
     &             ScatterRhsR, COMM, IERR)
              ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
           ENDDO
           BUF_EFFSIZE= 0
        ELSE
            DO WHILE ( ENTRIES_2_PROCESS .NE. 0)
              CALL MPI_RECV( BUF_INDX, BUF_MAXSIZE, MPI_INTEGER,
     &             MPI_ANY_SOURCE,
     &             ScatterRhsI, COMM, STATUS, IERR )
              CALL MPI_GET_COUNT( STATUS, MPI_INTEGER,BUF_EFFSIZE,IERR)
              PROC_WHO_ASKS = STATUS(MPI_SOURCE)
              DO I = 1, BUF_EFFSIZE
                 INDX = BUF_INDX( I )
                 DO K = 1, NRHS
                    BUF_RHS( K, I ) = RHS( INDX, K )
                 ENDDO
              ENDDO
              CALL MPI_SEND( BUF_RHS, NRHS*BUF_EFFSIZE,
     &             MPI_REAL, PROC_WHO_ASKS,
     &             ScatterRhsR, COMM, IERR)
              ENTRIES_2_PROCESS = ENTRIES_2_PROCESS - BUF_EFFSIZE
           ENDDO
           BUF_EFFSIZE= 0
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (MYID.NE.MASTER) THEN
         IF (NB_FS_IN_RHSINTR_F.LT.LRHSINTR) THEN
           DO K=1, NCOL_RHSINTR
            DO I = NB_FS_IN_RHSINTR_F +1, LRHSINTR
             RHSINTR (I, K) = ZERO
            ENDDO
           ENDDO
         ENDIF
        ENDIF
        DO ISTEP = 1, KEEP(28)
          IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &          KEEP(199))) THEN
             CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &       NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
              ELSE
                   J1=IPOS+1+LIELL
              END IF
              IF (MYID.EQ.MASTER) THEN
                INDX = POSINRHSINTR_FWD(IW(J1))
                IF (KEEP(350).EQ.2 .AND.
     &   (NRHS.EQ.1.OR.((NPIV*NRHS*2*KEEP(16)).GE.KEEP(364)))) THEN
!$                 OMP_FLAG = .FALSE.
!$                 CHUNK = NRHS
!$                 IF (NPIV*NRHS .GE. KEEP(363)) THEN
!$                   OMP_FLAG = .TRUE.
!$                   CHUNK = max((NPIV*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$                 ENDIF
#if ! defined(USE_OLD_SCALING)
                  IF (LSCAL) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(J,JJ)
!$OMP&  FIRSTPRIVATE(INDX) IF (OMP_FLAG)
                     DO K = 1, NRHS
                       DO JJ=J1,J1+NPIV-1
                         J=IW(JJ)
                         RHSINTR( INDX+JJ-J1, K ) = RHS( J, K ) *
     &                   SCALING_LOC_FWD( INDX+JJ-J1 )
                       ENDDO
                     ENDDO
!$OMP END PARALLEL DO
                  ELSE
#endif
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(J,JJ)
!$OMP&  FIRSTPRIVATE(INDX) IF (OMP_FLAG)
                    DO K = 1, NRHS
                      DO JJ=J1,J1+NPIV-1
                         J=IW(JJ)
                         RHSINTR( INDX+JJ-J1, K ) = RHS( J, K )
                      ENDDO
                    ENDDO
!$OMP END PARALLEL DO
#if ! defined(USE_OLD_SCALING)
                  ENDIF
#endif
                ELSE 
#if ! defined(USE_OLD_SCALING)
                  IF (LSCAL) THEN
                   DO JJ=J1,J1+NPIV-1
                      J=IW(JJ)
                      DO K = 1, NRHS
                         RHSINTR( INDX+JJ-J1, K ) = RHS( J, K )
     &                 * SCALING_LOC_FWD( INDX + JJ - J1 )
                      ENDDO
                   ENDDO
                  ELSE
#endif
                   DO JJ=J1,J1+NPIV-1
                      J=IW(JJ)
                      DO K = 1, NRHS
                         RHSINTR( INDX+JJ-J1, K ) = RHS( J, K )
                      ENDDO
                   ENDDO
#if ! defined(USE_OLD_SCALING)
                  ENDIF
#endif
                END IF 
              ELSE
                DO JJ=J1,J1+NPIV-1
                  BUF_EFFSIZE = BUF_EFFSIZE + 1
                  BUF_INDX(BUF_EFFSIZE) = IW(JJ)
                  IF (BUF_EFFSIZE + 1 .GT. BUF_MAXSIZE) THEN
                   CALL SMUMPS_GET_BUF_INDX_RHS()
                  ENDIF
                ENDDO
              ENDIF
          ENDIF
        ENDDO
        IF ( BUF_EFFSIZE .NE. 0 .AND. MYID.NE.MASTER ) 
     &              CALL SMUMPS_GET_BUF_INDX_RHS()
      ENDIF
      IF (KEEP(350).EQ.2) THEN
        DEALLOCATE (BUF_INDX, BUF_RHS_2)
      ELSE
        DEALLOCATE (BUF_INDX, BUF_RHS)
      ENDIF
      RETURN
      CONTAINS
                  SUBROUTINE SMUMPS_GET_BUF_INDX_RHS()
                  CALL MPI_SEND(BUF_INDX, BUF_EFFSIZE, MPI_INTEGER,
     &            MASTER, ScatterRhsI, COMM, IERR )
                  IF (KEEP(350).EQ.2) THEN
                    CALL MPI_RECV(BUF_RHS_2, BUF_EFFSIZE*NRHS,
     &                 MPI_REAL,
     &                 MASTER,
     &                 ScatterRhsR, COMM, STATUS, IERR )
!$                  OMP_FLAG = .FALSE.
!$                  CHUNK = NRHS
!$                  IF (BUF_EFFSIZE*NRHS .GE. KEEP(363)) THEN
!$                    OMP_FLAG = .TRUE.
!$              CHUNK = max((BUF_EFFSIZE*NRHS+NOMP-1)/NOMP,KEEP(363)/2)
!$                  ENDIF
#if ! defined(USE_OLD_SCALING)
                    IF (LSCAL) THEN
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(INDX)
!$OMP&  IF (OMP_FLAG)
                      DO K = 1, NRHS
                        DO I = 1, BUF_EFFSIZE
                           INDX = POSINRHSINTR_FWD(BUF_INDX(I))
                          RHSINTR( INDX, K ) =
     &                           BUF_RHS_2( I+(K-1)*BUF_EFFSIZE) *
     &                           SCALING_LOC_FWD( INDX )
                        ENDDO
                      ENDDO
!$OMP END PARALLEL DO
                    ELSE
#endif
!$OMP PARALLEL DO COLLAPSE(2) SCHEDULE(STATIC,CHUNK) PRIVATE(INDX)
!$OMP&  IF (OMP_FLAG)
                      DO K = 1, NRHS
                        DO I = 1, BUF_EFFSIZE
                           INDX = POSINRHSINTR_FWD(BUF_INDX(I))
                          RHSINTR( INDX, K ) =
     &                           BUF_RHS_2( I+(K-1)*BUF_EFFSIZE)
                        ENDDO
                      ENDDO
!$OMP END PARALLEL DO
#if ! defined(USE_OLD_SCALING)
                    ENDIF
#endif
                   ELSE
                     CALL MPI_RECV(BUF_RHS, BUF_EFFSIZE*NRHS,
     &                 MPI_REAL,
     &                 MASTER,
     &                 ScatterRhsR, COMM, STATUS, IERR )
#if ! defined(USE_OLD_SCALING)
                     IF (LSCAL) THEN
                     DO I = 1, BUF_EFFSIZE
                        INDX = POSINRHSINTR_FWD(BUF_INDX(I))
                        DO K = 1, NRHS
                           RHSINTR( INDX, K ) = BUF_RHS( K, I )
     &                     * SCALING_LOC_FWD( INDX )
                        ENDDO
                     ENDDO
                     ELSE
#endif
                     DO I = 1, BUF_EFFSIZE
                        INDX = POSINRHSINTR_FWD(BUF_INDX(I))
                        DO K = 1, NRHS
                           RHSINTR( INDX, K ) = BUF_RHS( K, I )
                        ENDDO
                     ENDDO
#if ! defined(USE_OLD_SCALING)
                     ENDIF
#endif
                  END IF 
                  BUF_EFFSIZE = 0
                  RETURN
                  END SUBROUTINE SMUMPS_GET_BUF_INDX_RHS
      END SUBROUTINE SMUMPS_SCATTER_RHS
      SUBROUTINE SMUMPS_BUILD_GLOB2LOC
     &           (NSLAVES, N, MYID_NODES,
     &           PTRIST,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP, 
     &           GLOB2LOC_RHS, GLOB2LOC_SOL,
     &           GLOB2LOC_SOL_ALLOC,
     &           MTYPE,
     &           NBENT_RHSINTR, NB_FS_IN_RHSINTR )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, intent(in) :: NSLAVES, N, MYID_NODES, LIW
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER, intent(in) :: IW(LIW), STEP(N)
      INTEGER, intent(in) :: MTYPE 
      LOGICAL, intent(in) :: GLOB2LOC_SOL_ALLOC
      INTEGER, intent(out):: GLOB2LOC_RHS(N), GLOB2LOC_SOL(N) 
      INTEGER, intent(out):: NBENT_RHSINTR, NB_FS_IN_RHSINTR
      INTEGER ISTEP
      INTEGER NPIV
      INTEGER IPOS, LIELL
      INTEGER JJ, J1, JCOL
      INTEGER IPOSINRHSINTR, IPOSINRHSINTR_SOL
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      GLOB2LOC_RHS = 0
      IF (GLOB2LOC_SOL_ALLOC) GLOB2LOC_SOL = 0
      IPOSINRHSINTR   = 1     
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &     KEEP(199))) THEN
           CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP, NPIV, LIELL,
     &     IPOS, IW, LIW, PTRIST, STEP, N )
           IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
           ELSE
                   J1=IPOS+1+LIELL
           END IF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                   JCOL = IPOS+1+LIELL
           ELSE
                   JCOL = IPOS+1
           ENDIF
           DO JJ = J1, J1+NPIV-1
                GLOB2LOC_RHS(IW(JJ)) = IPOSINRHSINTR+JJ-J1
           ENDDO
           IF (GLOB2LOC_SOL_ALLOC) THEN
             DO JJ = JCOL, JCOL+NPIV-1
               GLOB2LOC_SOL(IW(JJ)) = IPOSINRHSINTR+JJ-JCOL
             ENDDO
           ENDIF
           IPOSINRHSINTR       = IPOSINRHSINTR + NPIV
        ENDIF
      ENDDO
      NB_FS_IN_RHSINTR = IPOSINRHSINTR -1
      IF (GLOB2LOC_SOL_ALLOC) IPOSINRHSINTR_SOL=IPOSINRHSINTR
      IF (IPOSINRHSINTR.GT.N) GOTO 500 
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &     KEEP(199))) THEN
           CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &     NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
           IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                   J1=IPOS+1
           ELSE
                   J1=IPOS+1+LIELL
           END IF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                   JCOL = IPOS+1+LIELL
           ELSE
                   JCOL = IPOS+1
           ENDIF
           IF (GLOB2LOC_SOL_ALLOC) THEN
            DO JJ = NPIV, LIELL-1-KEEP(253)
              IF (GLOB2LOC_RHS(IW(J1+JJ)).EQ.0) THEN
               GLOB2LOC_RHS(IW(J1+JJ)) = - IPOSINRHSINTR
               IPOSINRHSINTR = IPOSINRHSINTR + 1
              ENDIF
              IF (GLOB2LOC_SOL(IW(JCOL+JJ)).EQ.0) THEN
               GLOB2LOC_SOL(IW(JCOL+JJ)) = - IPOSINRHSINTR_SOL
               IPOSINRHSINTR_SOL = IPOSINRHSINTR_SOL + 1
              ENDIF
             ENDDO
           ELSE
             DO JJ = J1+NPIV, J1+LIELL-1-KEEP(253)
              IF (GLOB2LOC_RHS(IW(JJ)).EQ.0) THEN
               GLOB2LOC_RHS(IW(JJ)) = - IPOSINRHSINTR
               IPOSINRHSINTR = IPOSINRHSINTR + 1
              ENDIF
             ENDDO
           ENDIF
        ENDIF
      ENDDO
 500  NBENT_RHSINTR = IPOSINRHSINTR - 1
      IF (GLOB2LOC_SOL_ALLOC) 
     &     NBENT_RHSINTR = max(NBENT_RHSINTR, IPOSINRHSINTR_SOL-1)
      RETURN
      END SUBROUTINE SMUMPS_BUILD_GLOB2LOC
      SUBROUTINE SMUMPS_NODES_FWD_BWD_SIZE_FILL (
     &           fill, ICNTL,
     &           N, NSTEPS, KEEP, STEP, Step2node,
     &           IRHS_loc, Nloc_RHS, 
     &           IRHS_PTR, NBCOL_INBLOC, IRHS_SPARSE, NZ_RHS,
     &           PERM_RHS, SIZE_PERM_RHS, JBEG_RHS,
     &           UNS_PERM_INV, SIZE_UNS_PERM_INV,
     &           ICNTL21, 
     &           MYID, COMM,
     &           INFO, 
     &           Pruned_Sons_FWD, Pruned_Sons_BWD, 
     &           Lnodes_FWD, Lnodes_BWD
     &           , nodes_FWD, nodes_BWD
     &           , Lnodes_FWD_in, Lnodes_BWD_in
     &           )
      USE SMUMPS_SOL_ES, ONLY : SMUMPS_ES_NODES_SIZE_AND_FILL
      IMPLICIT NONE
      LOGICAL, INTENT(IN)  :: fill
      INTEGER, INTENT(IN)  :: ICNTL(60),N, NSTEPS, KEEP(500) 
      INTEGER, INTENT(IN)  :: STEP(N), Step2node(NSTEPS)
      INTEGER, INTENT(IN)  :: Nloc_RHS, 
     &                        IRHS_loc(max(1,Nloc_RHS))
      INTEGER, intent(in) :: NBCOL_INBLOC, IRHS_PTR(NBCOL_INBLOC+1)
      INTEGER, intent(in) :: NZ_RHS, IRHS_SPARSE(NZ_RHS)
      INTEGER, intent(in) :: SIZE_PERM_RHS, PERM_RHS(SIZE_PERM_RHS)
      INTEGER, intent(in) :: JBEG_RHS, SIZE_UNS_PERM_INV
      INTEGER, intent(in) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, intent(IN)  :: ICNTL21
      INTEGER, intent(in)  :: MYID, COMM
      INTEGER, intent(inout) :: INFO(80)
      INTEGER, intent(inout) :: Pruned_Sons_FWD(NSTEPS), 
     &                          Pruned_Sons_BWD(NSTEPS)
      INTEGER, intent(inout) :: Lnodes_FWD, Lnodes_BWD
      INTEGER, intent(in) :: Lnodes_FWD_in, Lnodes_BWD_in
      INTEGER, intent(out) :: nodes_FWD(Lnodes_FWD_in), 
     &                        nodes_BWD(Lnodes_BWD_in)
      INCLUDE 'mpif.h'
      LOGICAL :: DO_PRUN_FWD, AM1, Exploit_Sparsity_FWD, 
     &           Exploit_Sparsity_BWD
      INTEGER :: Lnodes_FWD_loc, Lnodes_BWD_loc, ISTEP,
     &           INODE_PRINC, I, II, JAM1
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER, DIMENSION(:), ALLOCATABLE ::  TMP_INT_ARRAY
      INTEGER :: allocok
#endif
#if defined(AVOID_MPI_IN_PLACE)
      ALLOCATE(TMP_INT_ARRAY(KEEP(28)), STAT = allocok)
      IF(allocok.GT.0) THEN
        INFO(1)=-13
        INFO(2)=KEEP(28)
      END IF
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM, MYID )
      IF(INFO(1).LT.0) GOTO 500
#endif
      AM1                  = (KEEP(237) .NE. 0)
      Exploit_Sparsity_FWD = (KEEP(235) .NE. 0) .AND. (.NOT. AM1)
      DO_PRUN_FWD          = (Exploit_Sparsity_FWD.OR.AM1)
      Exploit_Sparsity_BWD = (KEEP(212) .NE. 0) .AND. (.NOT. AM1)
      IF (.NOT.fill) Lnodes_FWD=-1
      IF (.NOT.fill) Lnodes_BWD=-1
      IF (.NOT.fill.AND.KEEP(252).NE.0) THEN
        Lnodes_FWD = 0
      ENDIF
      IF ( KEEP(252).NE.0 ) DO_PRUN_FWD = .FALSE.  
      IF ( DO_PRUN_FWD ) THEN
        IF ( Exploit_Sparsity_FWD.AND.KEEP(248).EQ.-1 ) THEN
         IF (.NOT.fill) THEN
           CALL SMUMPS_ES_NODES_SIZE_AND_FILL ( fill,
     &           N, KEEP(28), KEEP, STEP, Step2node,
     &           IRHS_loc, Nloc_RHS, MYID, COMM,
     &           Pruned_Sons_FWD, Lnodes_FWD 
#if defined(AVOID_MPI_IN_PLACE)
     &           , TMP_INT_ARRAY
#endif
     &           )
         ELSE IF (Lnodes_FWD.GT.0) THEN
            CALL SMUMPS_ES_NODES_SIZE_AND_FILL ( fill,
     &           N, KEEP(28), KEEP, STEP, Step2node,
     &           IRHS_loc, Nloc_RHS, MYID, COMM,
     &           Pruned_Sons_FWD, Lnodes_FWD, 
#if defined(AVOID_MPI_IN_PLACE)
     &            TMP_INT_ARRAY, 
#endif
     &           nodes_FWD
     &           )
         ENDIF
        ELSE IF ( Exploit_Sparsity_FWD.AND.KEEP(248).NE.-1 ) THEN
           IF (.NOT.fill) THEN
            Lnodes_FWD = 0
            Pruned_Sons_FWD(:) = -1
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               IF ( Pruned_Sons_FWD(ISTEP) .eq. -1) THEN
                  Lnodes_FWD = Lnodes_FWD +1
                  Pruned_Sons_FWD(ISTEP) = 0 
               ENDIF
            ENDDO
           ELSE IF (Lnodes_FWD.GT.0) THEN
            Lnodes_FWD_loc = 0
            Pruned_Sons_FWD(:) = -1
            DO I = 1, NZ_RHS
               ISTEP       = abs( STEP(IRHS_SPARSE(I)) )
               INODE_PRINC = Step2node( ISTEP )
               IF ( Pruned_Sons_FWD(ISTEP) .eq. -1) THEN
                  Lnodes_FWD_loc = Lnodes_FWD_loc +1
                  nodes_FWD(Lnodes_FWD_loc)  = INODE_PRINC
                  Pruned_Sons_FWD(ISTEP) = 0 
               ENDIF
            ENDDO
           ENDIF
        ELSE IF ( AM1 ) THEN  
          IF (.NOT.fill) THEN
            Lnodes_FWD = 0
            Pruned_Sons_FWD(:) = -1
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF       
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_Sons_FWD(ISTEP) .eq. -1) THEN
                 Lnodes_FWD = Lnodes_FWD +1
                 Pruned_Sons_FWD(ISTEP) = 0                 
              ENDIF
            ENDDO
          ELSE IF (Lnodes_FWD.GT.0) THEN
            Lnodes_FWD_loc = 0
            Pruned_Sons_FWD = -1
            DO I = 1, NBCOL_INBLOC
              IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
              IF ( (KEEP(242) .NE. 0 ).OR. (KEEP(243).NE.0) ) THEN
                   JAM1 = PERM_RHS(JBEG_RHS+I-1)
              ELSE
                   JAM1 = JBEG_RHS+I-1
              ENDIF
              ISTEP = abs(STEP(JAM1))
              INODE_PRINC = Step2node(ISTEP)
              IF ( Pruned_Sons_FWD(ISTEP) .eq. -1) THEN
                 Lnodes_FWD_loc = Lnodes_FWD_loc +1
                 nodes_FWD(Lnodes_FWD_loc)  = INODE_PRINC
                 Pruned_Sons_FWD(ISTEP) = 0
              ENDIF
            ENDDO            
          ENDIF
         ENDIF                  
      ENDIF
      IF (AM1) THEN
        IF (.NOT.fill) THEN
         Lnodes_BWD = 0
         Pruned_Sons_BWD(:) = -1 
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)          
            ISTEP = abs(STEP(I))
            IF ( Pruned_Sons_BWD(ISTEP) .eq. -1) THEN
               Lnodes_BWD = Lnodes_BWD +1
               Pruned_Sons_BWD(ISTEP) = 0
            ENDIF
         ENDDO
        ELSE IF (Lnodes_BWD.GT.0) THEN
         Lnodes_BWD_loc = 0         
         Pruned_Sons_BWD(:) = -1  
         DO II = 1, NZ_RHS
            I = IRHS_SPARSE(II)
            IF (KEEP(23).NE.0) I = UNS_PERM_INV(I)
            ISTEP = abs(STEP(I))
            INODE_PRINC = Step2node(ISTEP)
            IF ( Pruned_Sons_BWD(ISTEP) .eq. -1) THEN
               Lnodes_BWD_loc = Lnodes_BWD_loc +1
               nodes_BWD(Lnodes_BWD_loc)  = INODE_PRINC
               Pruned_Sons_BWD(ISTEP) = 0
            ENDIF
         ENDDO
        ENDIF
      ENDIF
#if defined(AVOID_MPI_IN_PLACE)
      GOTO 600  
  500 CONTINUE
      Lnodes_FWD = -1
      Lnodes_BWD = -1
  600 CONTINUE
#endif
#if defined(AVOID_MPI_IN_PLACE)
      IF ( allocated(TMP_INT_ARRAY)) DEALLOCATE(TMP_INT_ARRAY)
#endif
      RETURN
      END SUBROUTINE SMUMPS_NODES_FWD_BWD_SIZE_FILL
      SUBROUTINE SMUMPS_BUILD_GLOB2LOC_NODES_ES (
     &           NSLAVES, N, MYID_NODES,
     &           PTRIST, DAD,
     &           KEEP,KEEP8, PROCNODE_STEPS, IW, LIW, STEP,
     &           Lnodes_FWD, Lnodes_BWD,
     &           nodes_FWD, nodes_BWD,
     &           GLOB2LOC_RHS, GLOB2LOC_SOL,
     &           GLOB2LOC_SOL_ALLOC,
     &           MTYPE,
     &           NBENT_RHSINTR,
     &           NB_FS_IN_RHSINTR_FWD, NB_FS_IN_RHSINTR_TOT
     &            )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER, intent(in)    :: NSLAVES, N, MYID_NODES, LIW
      INTEGER, intent(in)    :: KEEP(500)
      INTEGER(8), intent(in) :: KEEP8(150)
      INTEGER, intent(in) :: PTRIST(KEEP(28)), PROCNODE_STEPS(KEEP(28))
      INTEGER, intent(in) :: Lnodes_FWD, Lnodes_BWD
      INTEGER, intent(in) :: nodes_FWD(max(1,Lnodes_FWD)), 
     &                       nodes_BWD(max(1,Lnodes_BWD))
      INTEGER, intent(inout) :: DAD(KEEP(28))
      INTEGER, intent(in) :: IW(LIW), STEP(N)
      INTEGER, intent(in) :: MTYPE
      LOGICAL, intent(in) :: GLOB2LOC_SOL_ALLOC
      INTEGER, intent(out):: GLOB2LOC_RHS(N), GLOB2LOC_SOL(N)
      INTEGER, intent(out):: NBENT_RHSINTR
      INTEGER, intent(out):: NB_FS_IN_RHSINTR_FWD, NB_FS_IN_RHSINTR_TOT
      INTEGER I
      INTEGER ISTEP, OLDISTEP
      INTEGER NPIV
      INTEGER IPOS, LIELL
      INTEGER JJ, J1, JCOL, ABSJCOL
      INTEGER IPOSINRHSINTR_RHS, IPOSINRHSINTR_SOL
      INTEGER NBENT_RHSINTR_ROW, NBENT_RHSINTR_COL
      LOGICAL GO_UP
      INCLUDE 'mumps_headers.h'
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      GLOB2LOC_RHS = 0
      IF (GLOB2LOC_SOL_ALLOC) GLOB2LOC_SOL = 0
      IPOSINRHSINTR_RHS = 0
      IPOSINRHSINTR_SOL = 0
      DO I = 1, Lnodes_FWD
        ISTEP = STEP(nodes_FWD(I))
        GO_UP = .TRUE.
        DO WHILE(GO_UP) 
          IF(MYID_NODES.EQ.
     &      MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
            CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &      NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
            IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
              J1=IPOS+1
            ELSE
              J1=IPOS+1+LIELL
            END IF
            IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
              JCOL = IPOS+1+LIELL
            ELSE
              JCOL = IPOS+1
            ENDIF
            IF(NPIV.GT.0) THEN 
              IF(GLOB2LOC_RHS(IW(J1)).EQ.0) THEN 
                DO JJ = J1, J1+NPIV-1
                  GLOB2LOC_RHS(IW(JJ))
     &            = IPOSINRHSINTR_RHS + JJ - J1 + 1
                ENDDO
                IPOSINRHSINTR_RHS = IPOSINRHSINTR_RHS + NPIV
                IF (GLOB2LOC_SOL_ALLOC) THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    GLOB2LOC_SOL(IW(JJ))
     &              = - (IPOSINRHSINTR_SOL + JJ - JCOL + 1)
                  ENDDO
                  IPOSINRHSINTR_SOL = IPOSINRHSINTR_SOL + NPIV
                ENDIF
              ELSE
                GO_UP = .FALSE. 
              END IF
            END IF 
          END IF 
          IF(DAD(ISTEP).NE.0) THEN
            ISTEP = STEP(DAD(ISTEP))
          ELSE 
            GO_UP = .FALSE.
          END IF
        END DO 
      END DO 
      NB_FS_IN_RHSINTR_FWD = IPOSINRHSINTR_RHS
      IF(GLOB2LOC_SOL_ALLOC) THEN 
        DO I=1, Lnodes_BWD
          ISTEP = STEP(nodes_BWD(I))
          GO_UP = .TRUE.
          DO WHILE(GO_UP) 
            IF(MYID_NODES.EQ.
     &        MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                J1=IPOS+1
              ELSE
                J1=IPOS+1+LIELL
              END IF
              IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                JCOL = IPOS+1+LIELL
              ELSE
                JCOL = IPOS+1
              ENDIF
              ABSJCOL = abs(IW(JCOL))
              IF(NPIV.GT.0) THEN 
                IF(GLOB2LOC_SOL(ABSJCOL).EQ.0)  THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    GLOB2LOC_SOL(abs(IW(JJ))) = 
     &                     IPOSINRHSINTR_SOL+JJ-JCOL+1
                  END DO
                  IPOSINRHSINTR_SOL = IPOSINRHSINTR_SOL + NPIV
                ELSE IF (GLOB2LOC_SOL(ABSJCOL).LT.0) THEN
                  DO JJ = JCOL, JCOL+NPIV-1
                    GLOB2LOC_SOL(abs(IW(JJ)))=
     &                -(GLOB2LOC_SOL(abs(IW(JJ))))
                  END DO
                ELSE
                  GO_UP = .FALSE.
                END IF
              END IF 
            END IF 
            IF(DAD(ISTEP).NE.0) THEN
              ISTEP = STEP(DAD(ISTEP))
            ELSE 
              GO_UP = .FALSE.
            END IF
          END DO 
        END DO 
      END IF 
      NB_FS_IN_RHSINTR_TOT = IPOSINRHSINTR_SOL
      IF (NSLAVES.NE.1) THEN
      DO I = 1, Lnodes_FWD
        ISTEP = STEP(nodes_FWD(I))
        GO_UP = .TRUE.
        DO WHILE(GO_UP) 
          IF(MYID_NODES.EQ.
     &      MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
            CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &      NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
            IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
              J1=IPOS+1
            ELSE
              J1=IPOS+1+LIELL
            END IF
            IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
              JCOL = IPOS+1+LIELL
            ELSE
              JCOL = IPOS+1
            ENDIF
            DO JJ = NPIV, LIELL-1-KEEP(253)
              IF(GLOB2LOC_RHS(IW(J1+JJ)).EQ.0) THEN 
                IPOSINRHSINTR_RHS = IPOSINRHSINTR_RHS + 1
                GLOB2LOC_RHS(IW(JJ+J1)) = -IPOSINRHSINTR_RHS
              END IF
            END DO
          END IF
          IF(DAD(ISTEP).GT.0) THEN
            OLDISTEP=ISTEP
            ISTEP = STEP(DAD(ISTEP))
            DAD(OLDISTEP)=-DAD(OLDISTEP)
          ELSE 
            GO_UP = .FALSE.
          END IF
        END DO 
      END DO 
      DAD=ABS(DAD)
      IF(GLOB2LOC_SOL_ALLOC) THEN
        DO I=1, Lnodes_BWD
        ISTEP = STEP(nodes_BWD(I))
          GO_UP = .TRUE.
          DO WHILE(GO_UP) 
            IF(MYID_NODES.EQ.
     &        MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),KEEP(199))) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, KEEP,
     &        NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
              IF (MTYPE.eq.1 .OR. KEEP(50).NE.0) THEN
                J1=IPOS+1
              ELSE
                J1=IPOS+1+LIELL
              END IF
              IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
                JCOL = IPOS+1+LIELL
              ELSE
                JCOL = IPOS+1
              ENDIF
              DO JJ = NPIV, LIELL-1-KEEP(253)
                IF(GLOB2LOC_SOL(IW(JCOL+JJ)).EQ.0) THEN 
                  IPOSINRHSINTR_SOL = IPOSINRHSINTR_SOL + 1
                  GLOB2LOC_SOL(IW(JCOL+JJ)) = -IPOSINRHSINTR_SOL
                END IF
              END DO
            END IF
            IF(DAD(ISTEP).GT.0) THEN
              OLDISTEP=ISTEP
              ISTEP = STEP(DAD(ISTEP))
              DAD(OLDISTEP)=-DAD(OLDISTEP)
            ELSE 
              GO_UP = .FALSE.
            END IF
          END DO 
        END DO 
       DAD=ABS(DAD)
      END IF
      ENDIF
      NBENT_RHSINTR_ROW = IPOSINRHSINTR_RHS
      NBENT_RHSINTR_COL = IPOSINRHSINTR_SOL
      NBENT_RHSINTR = max(NBENT_RHSINTR_ROW,NBENT_RHSINTR_COL)
      RETURN
      END SUBROUTINE SMUMPS_BUILD_GLOB2LOC_NODES_ES
