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
      MODULE SMUMPS_SOL_ES
      PRIVATE
      PUBLIC:: SMUMPS_CHAIN_PRUN_NODES
      PUBLIC:: SMUMPS_CHAIN_PRUN_NODES_STATS
      PUBLIC:: SMUMPS_INITIALIZE_RHS_BOUNDS
      PUBLIC:: SMUMPS_PROPAGATE_RHS_BOUNDS
      PUBLIC:: SMUMPS_TREE_PRUN_NODES
      PUBLIC:: SMUMPS_SOL_ES_INIT
#     if defined(STAT_ES_SOLVE)
      PUBLIC:: SMUMPS_SOL_ES_PRINT_STATS
#     endif
      PUBLIC:: SMUMPS_ES_GET_SUM_Nloc
      PUBLIC:: SMUMPS_ES_NODES_SIZE_AND_FILL
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK
#     if defined(STAT_ES_SOLVE)
      DOUBLE PRECISION :: nb_flops,
     &                    nb_sparse_flops,
     &                    total_efficiency
      INTEGER :: total_procs, total_blocks
#endif
      INCLUDE 'mumps_headers.h'
      CONTAINS
      SUBROUTINE SMUMPS_SOL_ES_INIT(SIZE_OF_BLOCK_ARG, KEEP201)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KEEP201
      INTEGER(8), POINTER, DIMENSION(:,:) :: SIZE_OF_BLOCK_ARG
      IF (KEEP201 > 0) THEN
        SIZE_OF_BLOCK => SIZE_OF_BLOCK_ARG
      ELSE
        NULLIFY(SIZE_OF_BLOCK)
      ENDIF
#if defined(STAT_ES_SOLVE)
      nb_flops=0.0d0
      nb_sparse_flops=0.0d0
      total_efficiency=0.0d0
      total_procs=0
      total_blocks=0
#endif
      RETURN
      END SUBROUTINE SMUMPS_SOL_ES_INIT
      SUBROUTINE SMUMPS_TREE_PRUN_NODES( 
     &     fill,
     &     DAD, NE_STEPS, FRERE, KEEP28,
     &     FILS, STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     TO_PROCESS,
     &     nb_prun_nodes, nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N, KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28),NE_STEPS(KEEP28),FRERE(KEEP28)
      INTEGER, INTENT(IN) :: FILS(N), STEP(N)
      INTEGER, INTENT(IN) :: nodes_RHS(:),  nb_nodes_RHS
      INTEGER :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes)
      INTEGER :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots)
      INTEGER :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves)
      LOGICAL :: TO_PROCESS(KEEP28) 
      INTEGER :: IN, I, ISTEP, TMP, TMPsave
      LOGICAL :: FILS_VISITED
      nb_prun_nodes = 0
      nb_prun_leaves = 0
      TO_PROCESS(:) = .FALSE.
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         TMPsave = TMP
         ISTEP = STEP(TMP)
         DO WHILE(.NOT.TO_PROCESS(ISTEP))
            TO_PROCESS(ISTEP) = .TRUE.
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = TMP
            END IF
            IN = FILS(TMP) 
            DO WHILE(IN.GT.0) 
               IN = FILS(IN)
            END DO
            FILS_VISITED = .FALSE.
            IF (IN.LT.0) THEN 
             FILS_VISITED = TO_PROCESS(STEP(-IN))
            ENDIF
            IF ( IN.LT.0.and..NOT.FILS_VISITED)
     &            THEN 
               TMP = -IN
               ISTEP = STEP(TMP)
            ELSE 
               IF (IN.EQ.0) THEN
                 nb_prun_leaves = nb_prun_leaves + 1
                 IF (fill) THEN
                    Pruned_Leaves(nb_prun_leaves) = TMP
                 END IF
               ELSE 
                 TMP = -IN
                 ISTEP = STEP(TMP)
               ENDIF
               DO WHILE (TMP.NE.TMPsave) 
                  TMP = abs(FRERE(ISTEP))
                  IF(TMP.NE.0) THEN 
                     ISTEP = STEP(TMP) 
                  ELSE 
                     exit
                  END IF
                  IF (.NOT.TO_PROCESS(ISTEP)) exit
               END DO
            END IF
         END DO
      END DO
      nb_prun_roots = 0
      DO I=1,nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF(DAD(ISTEP).NE.0) THEN 
            IF(.NOT.TO_PROCESS(STEP(DAD(ISTEP)))) THEN
               nb_prun_roots = nb_prun_roots + 1
               IF(fill) THEN
                  Pruned_Roots(nb_prun_roots) = TMP
               END IF
            END IF
         ELSE 
            nb_prun_roots = nb_prun_roots + 1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF          
         END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_TREE_PRUN_NODES
      SUBROUTINE SMUMPS_CHAIN_PRUN_NODES(
     &     fill,
     &     DAD, KEEP28,
     &     STEP, N,
     &     nodes_RHS, nb_nodes_RHS,
     &     Pruned_SONS, TO_PROCESS,
     &     nb_prun_nodes,nb_prun_roots, nb_prun_leaves,
     &     Pruned_List, Pruned_Roots, Pruned_Leaves
     &     )
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: fill
      INTEGER, INTENT(IN) :: N
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: KEEP28
      INTEGER, INTENT(IN) :: DAD(KEEP28)
      INTEGER, INTENT(IN) :: nb_nodes_RHS
      INTEGER, INTENT(IN) :: nodes_RHS(max(nb_nodes_RHS,1))
      INTEGER, INTENT(INOUT) :: nb_prun_nodes
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_List(nb_prun_nodes)
      INTEGER, INTENT(INOUT) :: nb_prun_roots
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Roots(nb_prun_roots)
      INTEGER, INTENT(INOUT) :: nb_prun_leaves
      INTEGER, OPTIONAL, INTENT(INOUT):: Pruned_Leaves(nb_prun_leaves)
      INTEGER, INTENT(OUT) :: Pruned_SONS(KEEP28)
      LOGICAL, INTENT(OUT) :: TO_PROCESS(KEEP28)
      INTEGER :: IN, I, ISTEP, TMP
      nb_prun_nodes = 0
      nb_prun_roots = 0
      TO_PROCESS(:) = .FALSE.
      Pruned_SONS(:) = -1
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         TO_PROCESS(ISTEP) = .TRUE.
         IF (Pruned_SONS(ISTEP) .eq. -1) THEN
            Pruned_SONS(ISTEP) = 0
            nb_prun_nodes = nb_prun_nodes + 1
            IF(fill) THEN
               Pruned_List(nb_prun_nodes) = nodes_RHS(I)
            END IF
            IN = nodes_RHS(I)
            IN = DAD(STEP(IN))
            DO WHILE (IN.NE.0)
               TO_PROCESS(STEP(IN)) = .TRUE.
               IF (Pruned_SONS(STEP(IN)).eq.-1) THEN 
                  nb_prun_nodes = nb_prun_nodes + 1
                  IF(fill) THEN
                     Pruned_List(nb_prun_nodes) = IN
                  END IF
                  Pruned_SONS(STEP(IN)) = 1
                  TMP = IN
                  IN = DAD(STEP(IN))
               ELSE 
                  Pruned_SONS(STEP(IN)) = Pruned_SONS(STEP(IN)) + 1
                  GOTO 201
               ENDIF
            ENDDO
            nb_prun_roots = nb_prun_roots +1
            IF(fill) THEN
               Pruned_Roots(nb_prun_roots) = TMP
            END IF
         ENDIF
  201    CONTINUE
      ENDDO
      nb_prun_leaves = 0
      DO I = 1, nb_nodes_RHS
         TMP = nodes_RHS(I)
         ISTEP = STEP(TMP)
         IF (Pruned_SONS(ISTEP).EQ.0) THEN
            nb_prun_leaves = nb_prun_leaves +1
            IF(fill) THEN
              Pruned_Leaves(nb_prun_leaves) = TMP
            END IF
         END IF
      ENDDO
      RETURN
      END SUBROUTINE SMUMPS_CHAIN_PRUN_NODES
      SUBROUTINE SMUMPS_INITIALIZE_RHS_BOUNDS(
     & STEP, N,
     & IRHS_PTR, NBCOL, IRHS_SPARSE, NZ_RHS,
     & JBEG_RHS, PERM_RHS, SIZE_PERM_RHS, K242, K243,
     & UNS_PERM_INV, SIZE_UNS_PERM_INV, K23,
     & RHS_BOUNDS, NSTEPS,
     & nb_sparse, MYID,
     & mode)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: MYID, N, NSTEPS, K242, K243, K23
      INTEGER, INTENT(IN) :: JBEG_RHS, SIZE_PERM_RHS, nb_sparse
      INTEGER, INTENT(IN) :: NBCOL, NZ_RHS, SIZE_UNS_PERM_INV
      INTEGER, INTENT(IN) :: STEP(N), PERM_RHS(SIZE_PERM_RHS)
      INTEGER, INTENT(IN) :: IRHS_PTR(NBCOL+1),IRHS_SPARSE(NZ_RHS)
      INTEGER, INTENT(IN) :: UNS_PERM_INV(SIZE_UNS_PERM_INV)
      INTEGER, INTENT(INOUT):: RHS_BOUNDS(2*NSTEPS)
      INTEGER, INTENT(IN) :: mode 
      INTEGER :: I, ICOL, JPTR, J, JAM1, node, bound
      RHS_BOUNDS = 0
      ICOL = 0 
      DO I = 1, NBCOL
        IF ( (IRHS_PTR(I+1)-IRHS_PTR(I)).EQ.0) CYCLE
        ICOL = ICOL + 1
        bound = ICOL - mod(ICOL, nb_sparse) + 1
        IF(mod(ICOL, nb_sparse).EQ.0) bound = bound - nb_sparse
        IF(mode.EQ.0) THEN 
          IF ((K242.NE.0).OR.(K243.NE.0)) THEN
            JAM1 = PERM_RHS(JBEG_RHS+I-1)
          ELSE
            JAM1 = JBEG_RHS+I-1
          ENDIF
          node = abs(STEP(JAM1))
          IF(RHS_BOUNDS(2*node - 1).EQ.0) THEN 
            RHS_BOUNDS(2*node - 1) = bound                 
            RHS_BOUNDS(2*node)     = bound + nb_sparse - 1 
          ELSE
            RHS_BOUNDS(2*node) = bound + nb_sparse - 1
          END IF
        ELSE  
          DO JPTR = IRHS_PTR(I), IRHS_PTR(I+1)-1
            J = IRHS_SPARSE(JPTR)
            IF ( mode .EQ. 1 ) THEN
              IF (K23.NE.0) J = UNS_PERM_INV(J)
            ENDIF
            node = abs(STEP(J))
            IF(RHS_BOUNDS(2*node - 1).EQ.0) THEN
              RHS_BOUNDS(2*node - 1) = bound
              RHS_BOUNDS(2*node)     = bound + nb_sparse - 1
            ELSE
              RHS_BOUNDS(2*node) = bound + nb_sparse - 1
            END IF
          END DO
        END IF
      END DO
      RETURN
      END SUBROUTINE SMUMPS_INITIALIZE_RHS_BOUNDS
      SUBROUTINE SMUMPS_PROPAGATE_RHS_BOUNDS(
     & pruned_leaves, nb_pruned_leaves,
     & STEP, N, Pruned_SONS,
     & DAD, RHS_BOUNDS, NSTEPS,
     & MYID, COMM, KEEP485,
#if defined(STAT_ES_SOLVE)
     & KEEP46,
     & IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING,
#endif
     & IW, LIW, PTRIST, KIXSZ,OOC_FCT_LOC, PHASE, LDLT, K38)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_headers.h'
      INTEGER, INTENT(IN) :: nb_pruned_leaves, N, NSTEPS
      INTEGER, INTENT(IN) :: STEP(N), DAD(NSTEPS), Pruned_SONS(NSTEPS)
      INTEGER, INTENT(IN) :: MYID, COMM, KEEP485
      INTEGER, INTENT(IN) :: pruned_leaves(nb_pruned_leaves)
      INTEGER, INTENT(IN) :: LIW, IW(LIW), PTRIST(NSTEPS)
      INTEGER, INTENT(IN) :: KIXSZ, OOC_FCT_LOC, PHASE, LDLT, K38
#     if defined(STAT_ES_SOLVE)
      INTEGER, INTENT(IN) :: KEEP46 
      INTEGER, INTENT(IN) :: SIZE_IPTR_WORKING, SIZE_WORKING, 
     &                       IPTR_WORKING(SIZE_IPTR_WORKING),
     &                       WORKING(SIZE_WORKING)
#endif
      INTEGER, INTENT(INOUT):: RHS_BOUNDS(2*NSTEPS)
      INTEGER :: I, node, father, size_pool, next_size_pool
      INTEGER :: IERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: POOL, NBSONS
#if defined(STAT_ES_SOLVE)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: isleaf
      INTEGER :: J, NPROCS, proc, allocok
      LOGICAL :: found
      DOUBLE PRECISION :: avg_load, efficiency, max_load, effmax
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: proc_flops_buf
      DOUBLE PRECISION :: proc_block_flops, block_flops
      INTEGER :: SK38  
      INTEGER, PARAMETER :: MASTER  = 0
#endif
      ALLOCATE(POOL(nb_pruned_leaves),
     &         NBSONS(NSTEPS),
     &         STAT=IERR)
      IF (IERR.NE.0) THEN
         WRITE(6,*)'Allocation problem in SMUMPS_PROPAGATE_RHS_BOUNDS'
         CALL MUMPS_ABORT()
      END IF
      size_pool = nb_pruned_leaves
      POOL = pruned_leaves
      NBSONS = Pruned_SONS
#     if defined(STAT_ES_SOLVE)
      NPROCS = SIZE_IPTR_WORKING-1
      IF((MYID.EQ.MASTER).AND.(KEEP46.EQ.1)) THEN
        ALLOCATE(isleaf(NSTEPS), STAT=IERR)
        IF (IERR.NE.0) THEN
          WRITE(6,*)'Allocation problem in SMUMPS_PROPAGATE_RHS_BOUNDS'
          CALL MUMPS_ABORT()
        END IF
        isleaf = .FALSE.
        DO I=1,nb_pruned_leaves
          isleaf(abs(STEP(pruned_leaves(I)))) = .true.
        END DO
        proc = 0
        DO I=1,NPROCS
          found = .FALSE.
          J = IPTR_WORKING(I)
          DO WHILE((J.LE.IPTR_WORKING(I+1)-1).AND.(.NOT.found))
            IF (isleaf(WORKING(J)))THEN
              found = .TRUE.            
            END IF
            J = J + 1
          END DO
          IF(found) THEN
            proc = proc + 1
          END IF
        END DO
        total_procs = total_procs + proc
        total_blocks = total_blocks + 1
        DEALLOCATE(isleaf)
      END IF
#     endif
      DO WHILE (size_pool.ne.0)
        next_size_pool =0 
        DO I=1, size_pool 
          node = STEP(POOL(I)) 
          IF (DAD(node).NE.0) THEN
            father = STEP(DAD(node))
            NBSONS(father) = NBSONS(father)-1
            IF (RHS_BOUNDS(2*father-1).EQ.0) THEN
              RHS_BOUNDS(2*father-1) = RHS_BOUNDS(2*node-1)
              RHS_BOUNDS(2*father)   = RHS_BOUNDS(2*node)
            ELSE
              RHS_BOUNDS(2*father-1) = min(RHS_BOUNDS(2*father-1),
     &                                     RHS_BOUNDS(2*node-1))
              RHS_BOUNDS(2*father) = max(RHS_BOUNDS(2*father),
     &                                     RHS_BOUNDS(2*node))
            END IF
            IF(NBSONS(father).EQ.0) THEN 
              next_size_pool = next_size_pool+1
              POOL(next_size_pool) = DAD(node)
            END IF
          END IF
        END DO
        size_pool = next_size_pool 
      END DO
      DEALLOCATE(POOL, NBSONS)
#     if defined(STAT_ES_SOLVE)
      IF (KEEP46.EQ.1) THEN
       IF(MYID.EQ.MASTER) THEN
        block_flops = 0D0
       END IF
       proc_block_flops = 0D0
       IF (K38 .GT. 0) THEN
        SK38 = STEP(K38)
       ELSE
        SK38 = 0
       END IF
       DO I=1,NSTEPS
        IF (RHS_BOUNDS(2*I).NE.0) THEN 
          IF(PTRIST(I).GT.0) THEN
            proc_block_flops = proc_block_flops
     &           + dble(2*(RHS_BOUNDS(2*I) - RHS_BOUNDS(2*I-1) +1))
     &           * dble(SMUMPS_LOCAL_FACTOR_SIZE(IW, LIW,
     &             PTRIST(I)+KIXSZ, 
     &             PHASE,LDLT,I.EQ.SK38))
          END IF
        END IF
       END DO
       IF(MYID.EQ.MASTER) THEN
          ALLOCATE(proc_flops_buf(SIZE_IPTR_WORKING-1),stat=allocok)
          IF(allocok.GT.0) THEN
             WRITE(6,*)'Allocation problem of proc_flops_buf'
     &            ,' in SMUMPS_PROPAGATE_RHS_BOUNDS'
             CALL MUMPS_ABORT()
          ENDIF
          proc_flops_buf=0.0d0
       ELSE 
          ALLOCATE(proc_flops_buf(1),stat=allocok)
          IF(allocok.GT.0) THEN
             WRITE(6,*)'Allocation problem of proc_flops_buf'
     &            ,' in SMUMPS_PROPAGATE_RHS_BOUNDS'
             CALL MUMPS_ABORT()
          ENDIF
          proc_flops_buf=0.0d0
       END IF
       CALL MPI_GATHER(proc_block_flops, 1, MPI_DOUBLE_PRECISION,
     &       proc_flops_buf, 1, MPI_DOUBLE_PRECISION,
     &       0, COMM, IERR)
       CALL MPI_REDUCE(proc_block_flops, block_flops, 1,
     &  MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERR)  
       IF(MYID.EQ.MASTER) THEN                         
        nb_sparse_flops = nb_sparse_flops+block_flops
        avg_load = sum(proc_flops_buf)/dble(NPROCS)
        max_load = maxval(proc_flops_buf)
        efficiency = 0D0
        effmax = 0D0
        DO I=1,NPROCS
          efficiency= efficiency + (proc_flops_buf(I)-avg_load)**2
          IF (proc_flops_buf(I)-avg_load.GT.0.0D0) THEN
            effmax = effmax + (max_load-avg_load)**2
          ELSE IF (proc_flops_buf(I)-avg_load.LT.0.0D0) THEN
            effmax = effmax + avg_load**2
          END IF
        END DO
        efficiency = sqrt(efficiency/dble(NPROCS))
        effmax = sqrt(effmax/dble(NPROCS))
        IF(effmax.ne.0.0d0) efficiency = efficiency / effmax
        efficiency = 1.0d0 - efficiency
        efficiency = efficiency * block_flops 
        total_efficiency = total_efficiency + efficiency
        DEALLOCATE(proc_flops_buf)
       ELSE
        DEALLOCATE(proc_flops_buf)
       END IF
      END IF
#endif
      RETURN
      END SUBROUTINE SMUMPS_PROPAGATE_RHS_BOUNDS
      INTEGER(8) FUNCTION SMUMPS_LOCAL_FACTOR_SIZE(IW,LIW,PTR,
     &                                 PHASE, LDLT, IS_ROOT)
        INTEGER, INTENT(IN) :: LIW, PTR, PHASE, LDLT
        INTEGER, INTENT(IN) :: IW(LIW)
        LOGICAL, INTENT(IN) :: IS_ROOT
        INTEGER(8) :: NCB, NELIM, LIELL, NPIV, NROW
        NCB   = int(IW(PTR),8)     
        NELIM = int(IW(PTR+1),8)   
        NROW  = int(IW(PTR+2),8)
        NPIV  = int(IW(PTR+3),8)
        LIELL = NPIV + NCB
        IF (IS_ROOT) THEN
          SMUMPS_LOCAL_FACTOR_SIZE = int(IW(PTR+1),8) *     
     &                               int(IW(PTR+2),8) / 2_8 
          RETURN
        ENDIF
        IF (NCB.GE.0_8) THEN 
          IF (PHASE.EQ.0   
     &      .OR. (PHASE.EQ.1.AND.LDLT.NE.0) 
     &  ) THEN
            SMUMPS_LOCAL_FACTOR_SIZE =
     &            NPIV*(NPIV-1_8)/2_8 + (NROW-NPIV)*NPIV
          ELSE
            SMUMPS_LOCAL_FACTOR_SIZE =
     &      NPIV*(NPIV-1_8)/2_8 + (LIELL-NPIV)*NPIV
          ENDIF
        ELSE
          SMUMPS_LOCAL_FACTOR_SIZE =
     &      -NCB*NELIM
        END IF
      RETURN
      END FUNCTION SMUMPS_LOCAL_FACTOR_SIZE
      SUBROUTINE SMUMPS_CHAIN_PRUN_NODES_STATS
     &                (MYID, N, KEEP28, KEEP201, KEEP485, FR_FACT,
     &                STEP, Pruned_List, nb_prun_nodes, OOC_FCT_TYPE_LOC
#     if defined(STAT_ES_SOLVE)
     &                , NRHS, COMM, IW, LIW, PTRIST, KIXSZ, PHASE,
     &                  LDLT, K38
#endif     
     & )
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP28, KEEP201, OOC_FCT_TYPE_LOC, N,
     &                       KEEP485
      INTEGER(8), intent(in) :: FR_FACT
      INTEGER, intent(in) :: nb_prun_nodes, MYID
      INTEGER, intent(in) :: Pruned_List(nb_prun_nodes)
      INTEGER, intent(in) :: STEP(N)
#if defined(STAT_ES_SOLVE)
      INTEGER, INTENT(IN):: LIW, COMM, NRHS, LDLT, K38
      INTEGER, INTENT(IN):: IW(LIW), PTRIST(KEEP28), KIXSZ, PHASE  
      DOUBLE PRECISION :: proc_block_flops, block_flops
      INTEGER(8) :: Pruned_Size_ic 
      INTEGER :: IERR
      INTEGER :: SK38 
#endif
      INCLUDE 'mpif.h'
      INTEGER I, ISTEP
      INTEGER(8) :: Pruned_Size
#if defined(STAT_ES_SOLVE)
      Pruned_Size_ic = 0_8
#endif
      Pruned_Size = 0_8
#if defined(STAT_ES_SOLVE)
      IF (K38 .GT. 0) THEN
        SK38 = STEP(K38)
      ELSE
        SK38 = 0
      END IF
#endif
      DO I = 1, nb_prun_nodes
        ISTEP = STEP(Pruned_List(I))
        IF (KEEP201 .GT. 0) THEN
            Pruned_Size = Pruned_Size + SIZE_OF_BLOCK
     &                    (ISTEP, OOC_FCT_TYPE_LOC)
        ENDIF
#if       defined(STAT_ES_SOLVE)
          IF (PTRIST(ISTEP) .GT. 0) THEN
            Pruned_Size_ic = Pruned_Size_ic + 
     &               SMUMPS_LOCAL_FACTOR_SIZE(IW, LIW,
     &               PTRIST(ISTEP)+KIXSZ, 
     &               PHASE, LDLT,
     &               ISTEP.EQ.SK38)
          ENDIF
#         endif
      ENDDO
#if defined(STAT_ES_SOLVE)
      proc_block_flops = dble(2_8*Pruned_Size_ic)*dble(NRHS)
      CALL MPI_REDUCE(proc_block_flops, block_flops, 1,
     & MPI_DOUBLE_PRECISION, MPI_SUM, 0, COMM, IERR)
      IF(MYID.EQ.0) nb_flops = nb_flops + block_flops
#endif
      RETURN
      END SUBROUTINE SMUMPS_CHAIN_PRUN_NODES_STATS
#if   defined(STAT_ES_SOLVE)
      SUBROUTINE SMUMPS_SOL_ES_PRINT_STATS(
     &                  K212, K235, K237, K485, K497, NZLU8,
     &                  NRHS, ICNTL27, N, K50, DKEEPS, RINFOGS, MPG)
      IMPLICIT NONE
      INTEGER, intent(in)    ::  K212, K235, K237, K485, K497, 
     &                           NRHS, MPG, ICNTL27, N, K50
      INTEGER(8), intent(in) :: NZLU8
      REAL, intent(out)    :: DKEEPS(5), RINFOGS(5)
      LOGICAL              :: AM1, ES_FWD, ES_BWD, DO_NBSPARSE
      IF (MPG.LE.0) RETURN
      AM1         = (K237 .NE. 0)
      ES_FWD      = (K235 .NE. 0) .AND. (.NOT. AM1)
      ES_BWD      = (K212 .NE. 0) .AND. (.NOT. AM1)
      DO_NBSPARSE = (K497.NE.0).AND.(NRHS.GT.1).AND.(ICNTL27.GT.1)
      IF (AM1) 
     & WRITE(MPG,'(/A)') ' ** FLOPS SUMMARY during SOLVE AM1 ** '
      IF ((ES_FWD).AND. (.NOT.ES_BWD)) 
     & WRITE(MPG,'(/A,A)') ' ** FLOPS SUMMARY during fwd step',
     &                  ' (exploit RHS sparsity) ** '
      IF ((.NOT.ES_FWD).AND. (ES_BWD))
     & WRITE(MPG,'(/A,A)') ' ** FLOPS SUMMARY during bwd step',
     &                  ' (selected entries in solution) ** '
      IF ((ES_FWD).AND. (ES_BWD)) 
     & WRITE(MPG,'(/A,/A)') 
     & ' ** FLOPS SUMMARY during SOLVE (fwd+bwd steps)',
     & '    (sparse RHS and selected entries in solution) **'
      IF (
     &     (ES_FWD) .AND. (.NOT.ES_BWD)
     &    .OR.
     &     (.NOT.ES_FWD) .AND. (ES_BWD)
     &   ) THEN
       IF (K50.NE.0) THEN
         DKEEPS(1)=(real(NZLU8)-real(N))*real(2*NRHS)
       ELSE
         DKEEPS(1)=(real(NZLU8)-real(N))*real(NRHS)
       ENDIF
      ELSE IF ((ES_FWD).AND.(ES_BWD)) THEN
       IF (K50.NE.0) THEN
         DKEEPS(1) = (real(NZLU8)-real(N))*real(4*NRHS)
       ELSE
         DKEEPS(1)=(real(NZLU8)-real(N))*real(2*NRHS)
       ENDIF
      ENDIF
      RINFOGS(1) = DKEEPS(1)
      IF (.NOT.AM1) THEN
       WRITE(MPG,'(A,F25.1)') 
     & ' RINFOG(24) FLOPS with dense full rank format  =', DKEEPS(1)
      ENDIF
      DKEEPS(2)=real(nb_flops) 
      IF (DO_NBSPARSE) DKEEPS(4)=real(nb_sparse_flops)
      IF (DO_NBSPARSE) THEN
       RINFOGS(2)= DKEEPS(4)
      ELSE
       RINFOGS(2)= DKEEPS(2)
      ENDIF
      WRITE(MPG,'(A,F25.1)') 
     & ' RINFOG(25) FLOPS with exploit sparsity (ES)   =', RINFOGS(2)
      RETURN
      END SUBROUTINE SMUMPS_SOL_ES_PRINT_STATS
#endif
      SUBROUTINE SMUMPS_ES_GET_SUM_Nloc (
     &           N, Nloc_ITAB, ITAB_loc, COMM,
     &           SUM_idNloc_8 )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N
#if defined(MUMPS_NOF2003)
        INTEGER, POINTER :: ITAB_loc (:)
#else
        INTEGER, INTENT( IN ), POINTER :: ITAB_loc (:)
#endif
      INTEGER, INTENT(IN)  :: Nloc_ITAB
      INTEGER, INTENT(IN)  :: COMM
      INTEGER(8)   :: SUM_idNloc_8
      INCLUDE 'mpif.h'
      INTEGER I, II, IERR_MPI 
      INTEGER(8) :: idNloc_8
      idNloc_8 = 0_8
      DO I= 1, Nloc_ITAB
          II = ITAB_loc(I)
          IF (II.GE.1 .and. II.LE.N) 
     &    idNloc_8 = idNloc_8 + 1_8
      ENDDO
      CALL MPI_ALLREDUCE (idNloc_8, SUM_idNloc_8, 1, 
     &                  MPI_INTEGER8,
     &                  MPI_SUM, COMM, IERR_MPI )
      RETURN
      END SUBROUTINE SMUMPS_ES_GET_SUM_Nloc
      SUBROUTINE SMUMPS_ES_NODES_SIZE_AND_FILL (
     &           fill, 
     &           N, NSTEPS, KEEP, STEP, Step2node,
     &           ITAB_loc, Nloc_ITAB, 
     &           MYID, COMM,
     &           Pruned_Sons, Lnodes_ITAB
#if defined(AVOID_MPI_IN_PLACE)
     &           , TMP_INT_ARRAY
#endif
     &           , nodes_ITAB
     &           )
      IMPLICIT NONE
      LOGICAL, INTENT(IN)  :: fill
      INTEGER, INTENT(IN)  :: N, NSTEPS, KEEP(500) 
      INTEGER, INTENT(IN)  :: Nloc_ITAB
      INTEGER, INTENT(IN)  :: STEP(N), Step2node(NSTEPS), 
     &                        ITAB_loc(max(1,Nloc_ITAB))
      INTEGER, INTENT(IN)  :: MYID, COMM
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER              :: TMP_INT_ARRAY(NSTEPS) 
#endif
      INTEGER, INTENT(INOUT) :: Pruned_Sons(NSTEPS), Lnodes_ITAB
      INTEGER, OPTIONAL, INTENT(OUT) :: nodes_ITAB(max(1,Lnodes_ITAB))
      INCLUDE 'mpif.h'
      INTEGER I, II, ISTEP, IERR_MPI, Lnodes_ITAB_loc, INODE_PRINC
      IF (.NOT.fill) THEN
        Pruned_SONS = 0
        DO I= 1, Nloc_ITAB
          II = ITAB_loc(I)
          IF (II.GE.1 .and. II.LE.N) THEN
            ISTEP = abs(STEP(II))
            IF ( Pruned_SONS(ISTEP) .eq. 0 ) THEN
              Pruned_SONS(ISTEP) = 1
            ENDIF
           ENDIF
        ENDDO
#if defined(AVOID_MPI_IN_PLACE)
        TMP_INT_ARRAY = Pruned_Sons
#endif
        CALL MPI_ALLREDUCE(
#if defined(AVOID_MPI_IN_PLACE)
     &         TMP_INT_ARRAY, 
#else
     &         MPI_IN_PLACE,
#endif
     &         Pruned_Sons, NSTEPS, 
     &         MPI_INTEGER, MPI_SUM, COMM, IERR_MPI)
        Lnodes_ITAB = 0
        DO ISTEP=1,NSTEPS
          if (Pruned_SONS(ISTEP) .NE.0) Lnodes_ITAB=Lnodes_ITAB+1
        ENDDO
      ELSE IF (Lnodes_ITAB.GT.0) THEN
        Lnodes_ITAB_loc = 0
        DO ISTEP=1,NSTEPS
          if (Pruned_SONS(ISTEP) .GT. 0) then
            Lnodes_ITAB_loc=Lnodes_ITAB_loc+1
            INODE_PRINC = Step2node( ISTEP )
            nodes_ITAB(Lnodes_ITAB_loc)  = INODE_PRINC
          endif
        ENDDO
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_ES_NODES_SIZE_AND_FILL
      END MODULE SMUMPS_SOL_ES
      SUBROUTINE SMUMPS_PERMUTE_RHS_GS
     &          (LP, LPOK, PROKG, MPG, PERM_STRAT, 
     &           SYM_PERM, N, NRHS,
     &           IRHS_PTR, SIZE_IRHS_PTR, 
     &           IRHS_SPARSE, NZRHS, 
     &           PERM_RHS, IERR
     &         )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LP, MPG, PERM_STRAT, N, NRHS, 
     &                       SIZE_IRHS_PTR,
     &                       NZRHS
      LOGICAL, INTENT(IN) :: LPOK, PROKG
      INTEGER, INTENT(IN) :: SYM_PERM(N)
      INTEGER, INTENT(IN) :: IRHS_PTR(SIZE_IRHS_PTR)
      INTEGER, INTENT(IN) :: IRHS_SPARSE(NZRHS)
      INTEGER, INTENT(OUT) :: PERM_RHS(NRHS)
      INTEGER, INTENT(OUT) :: IERR
      INTEGER :: I,J,K, POSINPERMRHS, JJ,
     &           KPOS
      INTEGER, ALLOCATABLE :: ROW_REFINDEX(:)
      IERR = 0
      IF ((PERM_STRAT.NE.-1).AND.(PERM_STRAT.NE.1)) THEN
       IERR=-1
       IF (LPOK)
     & WRITE(LP,*) " INTERNAL ERROR -1 in ",
     &       " SMUMPS_PERMUTE_RHS_GS, PERM_STRAT =", PERM_STRAT, 
     &       " is out of range "
       RETURN
      ENDIF
      IF (PERM_STRAT.EQ.-1) THEN
       DO I=1,NRHS
        PERM_RHS(I) = I
       END DO
       GOTO 490
      ENDIF
      ALLOCATE(ROW_REFINDEX(NRHS), STAT=IERR)
      IF (IERR.GT.0) THEN
       IERR=-1
       IF (LPOK) THEN
          WRITE(LP,*) " ERROR -2 : ", 
     &         " ALLOCATE IN SMUMPS_PERMUTE_RHS_GS OF SIZE :",
     &         NRHS
       ENDIF
       RETURN
      ENDIF
      DO I=1,NRHS
        IF (IRHS_PTR(I+1)-IRHS_PTR(I).LE.0) THEN
          IERR =  1
          IF (I.EQ.1) THEN
            ROW_REFINDEX(I) = IRHS_SPARSE(IRHS_PTR(I))
          ELSE
            ROW_REFINDEX(I) = ROW_REFINDEX(I-1)
          ENDIF
        ELSE
          ROW_REFINDEX(I) = IRHS_SPARSE(IRHS_PTR(I))
        ENDIF
      END DO
      POSINPERMRHS = 0
      DO I=1,NRHS
       KPOS = N+1 
       JJ   = 0   
       DO J=1,NRHS
        K = ROW_REFINDEX(J)
        IF (K.LE.0) CYCLE 
        IF (SYM_PERM(K).LT.KPOS) THEN
         KPOS = SYM_PERM(K)
         JJ   = J
        ENDIF
       END DO
       IF (JJ.EQ.0) THEN
         IERR = -3 
         IF (LPOK)
     &   WRITE(LP,*) " INTERNAL ERROR -3 in ",
     &       " SMUMPS_PERMUTE_RHS_GS "
         GOTO 500
       ENDIF
       POSINPERMRHS           = POSINPERMRHS + 1
       PERM_RHS(POSINPERMRHS) = JJ
       ROW_REFINDEX(JJ)       = -ROW_REFINDEX(JJ)
      END DO
      IF (POSINPERMRHS.NE.NRHS) THEN
         IF (LPOK)
     &   WRITE(LP,*) " INTERNAL ERROR -4 in ",
     &       " SMUMPS_PERMUTE_RHS_GS ", maxval(ROW_REFINDEX)
         IERR = -4
         GOTO 500
      ENDIF
  490 CONTINUE
 500  CONTINUE
      IF (allocated(ROW_REFINDEX)) DEALLOCATE(ROW_REFINDEX)
      END SUBROUTINE SMUMPS_PERMUTE_RHS_GS
      SUBROUTINE SMUMPS_PERMUTE_RHS_AM1
     &          (PERM_STRAT, SYM_PERM,
     &           IRHS_PTR, NHRS,
     &           PERM_RHS, SIZEPERM, IERR
     &         )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: PERM_STRAT, NHRS, SIZEPERM
      INTEGER, INTENT(IN) :: SYM_PERM(SIZEPERM)
      INTEGER, INTENT(IN) :: IRHS_PTR(NHRS)
      INTEGER, INTENT(OUT):: IERR   
      INTEGER, INTENT(OUT):: PERM_RHS(SIZEPERM)
      DOUBLE PRECISION :: RAND_NUM
      INTEGER  I, J, STRAT
      IERR = 0
      STRAT = PERM_STRAT 
      IF( (STRAT.NE.-3).AND.
     &    (STRAT.NE.-2).AND.
     &    (STRAT.NE.-1).AND.
     &    (STRAT.NE. 1).AND.
     &    (STRAT.NE. 2).AND.
     &    (STRAT.NE. 6) ) THEN
        WRITE(*,*)"Warning: incorrect value for the RHS permutation; ",
     &            "defaulting to post-order"
        STRAT = 1
      END IF
      IF (STRAT .EQ. -3) THEN
         PERM_RHS(1:SIZEPERM)=0
         DO I=1, SIZEPERM 
           CALL random_number(RAND_NUM) 
           RAND_NUM = RAND_NUM*dble(SIZEPERM) 
           J = ceiling(RAND_NUM) 
           DO WHILE (PERM_RHS(J).NE.0) 
             CALL random_number(RAND_NUM)
             RAND_NUM = RAND_NUM*dble(SIZEPERM)
             J = ceiling(RAND_NUM)
           ENDDO
           PERM_RHS(J)=I
         ENDDO
      ELSEIF (STRAT .EQ. -2) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM -I +1) = I
         ENDDO
      ELSEIF (STRAT .EQ. -1) THEN
         DO I=1, SIZEPERM
            PERM_RHS(I) = I
         ENDDO
      ELSEIF (STRAT .EQ.  1) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SYM_PERM(I)) = I
         ENDDO
      ELSEIF (STRAT .EQ.  2) THEN
         DO I=1, SIZEPERM
            PERM_RHS(SIZEPERM-SYM_PERM(I)+1) = I
         ENDDO
      ENDIF
      END SUBROUTINE SMUMPS_PERMUTE_RHS_AM1
      SUBROUTINE SMUMPS_INTERLEAVE_RHS_AM1(
     &  PERM_RHS, SIZE_PERM,
     &  IPTR_WORKING, SIZE_IPTR_WORKING, WORKING, SIZE_WORKING,
     &  IRHS_PTR,
     &  STEP, SYM_PERM, N, NBRHS,
     &  PROCNODE, NSTEPS, SLAVEF, KEEP199,
     &  behaviour_L0, reorder, n_select, PROKG, MPG 
     &  )
      IMPLICIT NONE
      INTEGER, INTENT(IN) ::  SIZE_PERM,
     &                        SIZE_IPTR_WORKING,
     &                        IPTR_WORKING(SIZE_IPTR_WORKING),
     &                        SIZE_WORKING,
     &                        WORKING(SIZE_WORKING),
     &                        N,
     &                        IRHS_PTR(N+1),
     &                        STEP(N),
     &                        SYM_PERM(N),
     &                        NBRHS,
     &                        NSTEPS,
     &                        PROCNODE(NSTEPS),
     &                        SLAVEF, KEEP199,
     &                        n_select, MPG
      LOGICAL, INTENT(IN) :: behaviour_L0,
     &                        reorder, PROKG
      INTEGER, INTENT(INOUT) :: PERM_RHS(SIZE_PERM)
      INTEGER :: I, J, K, IVAR, IBLOCK,
     &           entry,            
     &           node,             
     &           SIZE_PERM_WORKING,
     &           NB_NON_EMPTY,     
     &           to_be_found,      
     &           posintmprhs,      
     &           selected,         
     &           local_selected,   
     &           current_proc,     
     &           NPROCS,           
     &           n_pass,           
     &           pass,             
     &           nblocks,          
     &           n_select_loc,     
     &           IERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: TMP_RHS,    
     &                                      PTR_PROCS,  
     &                                      IPTR_PERM_WORKING,
     &                                      PERM_WORKING,
     &                                      MYTYPENODE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: 
     &         PERM_PO,     
     &         ISTEP2BLOCK, 
     &         NEXTINBLOCK  
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: USED
      LOGICAL :: allow_above_L0
      INTEGER, EXTERNAL :: MUMPS_TYPENODE_ROUGH
      NPROCS = SIZE_IPTR_WORKING - 1
      ALLOCATE(TMP_RHS(SIZE_PERM),
     &         PTR_PROCS(NPROCS),
     &         USED(SIZE_PERM),
     &         IPTR_PERM_WORKING(NPROCS+1),
     &         MYTYPENODE(NSTEPS),
     &         STAT=IERR)
      IF(IERR.GT.0) THEN
        WRITE(*,*)'Allocation error in SMUMPS_INTERLEAVE_RHS_AM1'
        CALL MUMPS_ABORT()
      END IF
      DO I=1, NSTEPS
        MYTYPENODE(I) = MUMPS_TYPENODE_ROUGH( PROCNODE(I), KEEP199 )
      ENDDO
      NB_NON_EMPTY = 0
      DO I=1,SIZE_PERM
        IF(IRHS_PTR(I+1)-IRHS_PTR(I).NE.0) THEN 
          NB_NON_EMPTY = NB_NON_EMPTY + 1
        END IF
      END DO
      K = 0
      IPTR_PERM_WORKING(1)=1
      DO I=1,NPROCS
        USED = .FALSE.
        DO J=IPTR_WORKING(I),IPTR_WORKING(I+1)-1
          USED(WORKING(J)) = .TRUE.        
        END DO
        DO J=1,N
          IF (USED(abs(STEP(PERM_RHS(J)))).AND.
     &      ((IRHS_PTR(PERM_RHS(J)+1)-IRHS_PTR(PERM_RHS(J))).NE.0))
     &    THEN
            K = K + 1
          END IF
        END DO
        IPTR_PERM_WORKING(I+1) = K+1
      END DO
      SIZE_PERM_WORKING = K
      ALLOCATE(PERM_WORKING(SIZE_PERM_WORKING),
     &         STAT=IERR)
      IF(IERR.GT.0) THEN
        WRITE(*,*)'Allocation error in SMUMPS_INTERLEAVE_RHS_AM1'
        CALL MUMPS_ABORT()
      END IF
      K = 0
      DO I=1,NPROCS
        USED = .FALSE.
        DO J=IPTR_WORKING(I),IPTR_WORKING(I+1)-1
          USED(WORKING(J)) = .TRUE.        
        END DO
        DO J=1,N
          IF (USED(abs(STEP(PERM_RHS(J)))).AND.
     &      ((IRHS_PTR(PERM_RHS(J)+1)-IRHS_PTR(PERM_RHS(J))).NE.0)) 
     &    THEN
            K = K + 1
            PERM_WORKING(K) = PERM_RHS(J)
          END IF
        END DO
      END DO      
      IF(behaviour_L0) THEN
        n_pass = 2
        allow_above_L0 = .false.
        to_be_found = 0
        DO I=1,SIZE_PERM
          IF((MYTYPENODE(abs(STEP(I))).LE.1).AND. 
     &    (IRHS_PTR(I+1)-IRHS_PTR(I).NE.0))       
     &    THEN
            to_be_found = to_be_found + 1
          END IF
        END DO
      ELSE
        n_pass = 1
        allow_above_L0 = .true.
        to_be_found = NB_NON_EMPTY
      END IF
      PTR_PROCS(1:NPROCS) = IPTR_PERM_WORKING(1:NPROCS)
      USED = .FALSE.
      current_proc = 1
      n_select_loc = n_select
      IF (n_select_loc.LE.0) THEN
       n_select_loc = 1
      ENDIF
      posintmprhs = 0
      DO pass=1,n_pass
        selected = 0
        DO WHILE(selected.LT.to_be_found)
          local_selected = 0
          DO WHILE(local_selected.LT.n_select_loc)
            IF(PTR_PROCS(current_proc).EQ.
     &        IPTR_PERM_WORKING(current_proc+1))
     &      THEN 
              EXIT
            ELSE 
              entry = PERM_WORKING(PTR_PROCS(current_proc))
              node  = abs(STEP(entry))
              IF(.NOT.USED(entry)) THEN
                IF(allow_above_L0.OR.(MYTYPENODE(node).LE.1)) THEN
                  USED(entry) = .TRUE.
                  selected = selected + 1
                  local_selected = local_selected + 1
                  posintmprhs = posintmprhs + 1
                  TMP_RHS(posintmprhs) = entry
                  IF(selected.EQ.to_be_found) EXIT
                END IF
              END IF
              PTR_PROCS(current_proc) = PTR_PROCS(current_proc) + 1
            END IF
          END DO
          current_proc = mod(current_proc,NPROCS)+1
        END DO
        to_be_found = NB_NON_EMPTY - to_be_found
        allow_above_L0 = .true.
        PTR_PROCS(1:NPROCS) = IPTR_PERM_WORKING(1:NPROCS)
      END DO
      DO I=1,SIZE_PERM
        IF(IRHS_PTR(PERM_RHS(I)+1)-IRHS_PTR(PERM_RHS(I)).EQ.0) THEN
          posintmprhs = posintmprhs+1
          TMP_RHS(posintmprhs) = PERM_RHS(I)
          IF(posintmprhs.EQ.SIZE_PERM) EXIT
        END IF
      END DO
      DEALLOCATE(PTR_PROCS, USED,
     &           IPTR_PERM_WORKING,
     &           PERM_WORKING, MYTYPENODE)
      IF(reorder) THEN
        nblocks = (N+NBRHS-1)/NBRHS 
        ALLOCATE(PERM_PO(N), ISTEP2BLOCK(N), NEXTINBLOCK(nblocks),
     &           stat=IERR)
        IF(IERR.GT.0) THEN
          IF (PROKG ) WRITE(MPG,*)
     &    'Warning: reorder not done in SMUMPS_INTERLEAVE_RHS_AM1'
          PERM_RHS = TMP_RHS
          GOTO 500
        ENDIF
        DO IVAR = 1, N
          K = SYM_PERM( IVAR )
          PERM_PO( K ) = IVAR
        END DO
        DO I = 1, N
          IBLOCK = 1 + ( I - 1 ) / NBRHS
          IVAR = TMP_RHS( I )
          K = SYM_PERM( IVAR ) 
          ISTEP2BLOCK( K ) = IBLOCK
        END DO
        DO IBLOCK = 1, NBLOCKS
          NEXTINBLOCK(IBLOCK) = 1 + (IBLOCK-1)*NBRHS
        ENDDO
        DO K = 1, N
          IBLOCK = ISTEP2BLOCK(K)
          IVAR = PERM_PO(K)
          PERM_RHS(NEXTINBLOCK(IBLOCK)) = IVAR
          NEXTINBLOCK(IBLOCK) = NEXTINBLOCK(IBLOCK) + 1
        ENDDO
      ELSE
        PERM_RHS = TMP_RHS
      END IF
 500  CONTINUE
      DEALLOCATE(TMP_RHS)
      IF (allocated(PERM_PO    )) DEALLOCATE(PERM_PO    )
      IF (allocated(ISTEP2BLOCK)) DEALLOCATE(ISTEP2BLOCK)
      IF (allocated(NEXTINBLOCK)) DEALLOCATE(NEXTINBLOCK)
      RETURN
      END SUBROUTINE SMUMPS_INTERLEAVE_RHS_AM1
