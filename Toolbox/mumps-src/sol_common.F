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
      SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS ( ISTEP, KEEP,
     &           NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ISTEP, LIW, KEEP(500), N
      INTEGER, INTENT(IN)  :: IW( LIW )
      INTEGER, INTENT(IN)  :: STEP( N ), PTRIST( KEEP(28) )
      INTEGER, INTENT(OUT) :: NPIV, LIELL, IPOS
      INCLUDE 'mumps_headers.h'
      INTEGER :: SROOT
      IF (KEEP(38) .NE. 0) THEN
        SROOT = STEP(KEEP(38))
      ELSE IF (KEEP(20) .NE. 0) THEN
        SROOT = STEP(KEEP(20))
      ELSE
        SROOT = 0
      ENDIF
      IPOS = PTRIST(ISTEP)
      IF (IPOS .LE. 0) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_SOL_GET_NPIV_LIELL_IPOS",
     &  ISTEP 
        CALL MUMPS_ABORT()
      ENDIF
      NPIV = IW(IPOS+3+KEEP(IXSZ))
      IF ( ISTEP.EQ.SROOT ) THEN
             IPOS = PTRIST(ISTEP)
             LIELL = IW(IPOS+3+KEEP(IXSZ))
             NPIV = LIELL
             IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
      ELSE
             IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
             LIELL = IW(IPOS-2)+IW(IPOS+1)
             IPOS= IPOS+1
             NPIV = IW(IPOS)
             IPOS= IPOS+1
             IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS
      SUBROUTINE MUMPS_GET_INDICES(MYID_NODES, NSLAVES, N,
     &           PTRIST, KEEP,KEEP8, IW, LIW, STEP, PROCNODE_STEPS,
     &           INDICES, ROW_OR_COL_INDICES)
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: KEEP(500)
      INTEGER(8), INTENT(IN) :: KEEP8(150)
      INTEGER, INTENT(IN)    :: MYID_NODES, NSLAVES, N, LIW
      INTEGER, INTENT(IN)    :: PTRIST(KEEP(28))
      INTEGER, INTENT(IN)    :: IW(LIW), STEP(N)
      INTEGER, INTENT(IN)    :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT(OUT)   :: INDICES(KEEP(89))
      INTEGER, INTENT(IN)    :: ROW_OR_COL_INDICES
      INTEGER :: ISTEP
      INTEGER :: NPIV, LIELL, IPOS
      INTEGER :: IINDICES 
      INTEGER :: J1 
      INCLUDE 'mumps_headers.h'
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      IINDICES = 0
      DO ISTEP = 1, KEEP(28)
        IF (MYID_NODES == MUMPS_PROCNODE(PROCNODE_STEPS(ISTEP),
     &                                   KEEP(199))) THEN
          CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS ( ISTEP, KEEP,
     &           NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
          IF ( ROW_OR_COL_INDICES .EQ. 0 ) THEN
            J1 = IPOS + 1
          ELSE IF (ROW_OR_COL_INDICES .EQ. 1 ) THEN
            J1 = IPOS + LIELL + 1
          ELSE
            WRITE(*,*) "Internal error 1 in MUMPS_GET_INDICES",
     &      ROW_OR_COL_INDICES
            CALL MUMPS_ABORT()
          ENDIF
          IF (IINDICES+NPIV .GT. KEEP(89)) THEN
            WRITE(*,*) "Internal error 2 in MUMPS_GET_INDICES",
     &      IINDICES, KEEP(89)
            CALL MUMPS_ABORT()
          ENDIF
          INDICES(IINDICES+1:IINDICES+NPIV)=IW(J1:J1+NPIV-1)
          IINDICES=IINDICES+NPIV
        ENDIF
      ENDDO
      IF (IINDICES .NE. KEEP(89)) THEN
        WRITE(*,*) "Internal error 3 in MUMPS_GET_INDICES",
     &  IINDICES, KEEP(89)
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_GET_INDICES
      SUBROUTINE MUMPS_SOL_RHSMAPINFO( N, Nloc_RHS, INFO23,
     &                           IRHS_loc, MAP_RHS_loc,
     &                           POSINRHSINTR_FWD,
     &                           NSLAVES, MYID_NODES, COMM_NODES,
     &                           ICNTL, INFO )
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, Nloc_RHS
      INTEGER, INTENT(IN)    :: INFO23
      INTEGER, INTENT(IN)    :: IRHS_loc   (max(1,Nloc_RHS))     
      INTEGER, INTENT(OUT)   :: MAP_RHS_loc(max(1,Nloc_RHS)) 
      INTEGER, INTENT(IN)    :: POSINRHSINTR_FWD (N)
      INTEGER, INTENT(IN)    :: NSLAVES, MYID_NODES, COMM_NODES
      INTEGER, INTENT(INOUT) :: INFO(80)
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INCLUDE 'mpif.h'
      INTEGER :: I, NFS_LOC, NFS_TOT, IERR_MPI, allocok
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER :: allocoktmp
#endif
      INTEGER, ALLOCATABLE, DIMENSION(:) :: GLOBAL_MAPPING 
#if defined(AVOID_MPI_IN_PLACE)
      ALLOCATE(GLOBAL_MAPPING(int(N,8)+int(N,8)), stat=allocok)
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        CALL MUMPS_SET_IERROR(int(N,8)+int(N,8), INFO(2))
      ENDIF
#else
      ALLOCATE(GLOBAL_MAPPING(N), stat=allocok)
      IF (allocok .GT. 0) THEN
        INFO(1)=-13
        INFO(2)= N
      ENDIF
#endif
#if defined(AVOID_MPI_IN_PLACE)
      allocoktmp = allocok
      CALL MPI_ALLREDUCE(allocoktmp, allocok, 1,
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, allocok, 1,
#endif
     &                   MPI_INTEGER, MPI_SUM,
     &                   COMM_NODES, IERR_MPI)
      IF (allocok .NE. 0) RETURN
      NFS_LOC = 0
      NFS_TOT = 0
      DO I = 1, N
        IF (POSINRHSINTR_FWD(I) .LE. 0) THEN
          GLOBAL_MAPPING(I) = 0
        ELSE
          GLOBAL_MAPPING(I) = MYID_NODES
          NFS_LOC = NFS_LOC + 1
        ENDIF
      ENDDO
      CALL MPI_ALLREDUCE(NFS_LOC, NFS_TOT, 1, MPI_INTEGER,
     &                   MPI_SUM, COMM_NODES, IERR_MPI)
#if defined(AVOID_MPI_IN_PLACE)
      DO I = 1, N
        GLOBAL_MAPPING(int(N,8)+int(I,8)) = GLOBAL_MAPPING(I)
      ENDDO
      CALL MUMPS_BIGALLREDUCE( .FALSE.,
     &                    GLOBAL_MAPPING(N+1), GLOBAL_MAPPING,
     &                    N, MPI_INTEGER,
     &                    MPI_SUM, COMM_NODES, IERR_MPI )
#else
      CALL MUMPS_BIGALLREDUCE( .TRUE.,
     &                    MPI_IN_PLACE, GLOBAL_MAPPING, N, MPI_INTEGER,
     &                    MPI_SUM, COMM_NODES, IERR_MPI )
#endif
      DO I = 1, Nloc_RHS
        IF (IRHS_loc(I) .GE.1 .AND. IRHS_loc(I) .LE. N) THEN
          MAP_RHS_loc(I) = GLOBAL_MAPPING(IRHS_loc(I))
        ELSE
          MAP_RHS_loc(I) = -87878787
        ENDIF
      ENDDO
      DEALLOCATE(GLOBAL_MAPPING)
      RETURN
      END SUBROUTINE MUMPS_SOL_RHSMAPINFO
      SUBROUTINE MUMPS_COMPUTE_LASTFS_DYN( INODE, LASTFSSBTR_DYN,
     &MTYPE, KEEP, IW, LIW, N, STEP, PTRIST, FILS, FRERE )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: INODE
      INTEGER, INTENT(OUT) :: LASTFSSBTR_DYN
      INTEGER, INTENT(IN)  :: N, MTYPE, LIW, KEEP(500)
      INTEGER, INTENT(IN)  :: IW(LIW), STEP( N ), PTRIST( KEEP(28) )
      INTEGER, INTENT(IN)  :: FILS(N), FRERE( KEEP(28) )
      INTEGER :: NPIV, LIELL, IPOS, INODE_CUR, IN
      INODE_CUR = INODE
 10   CONTINUE
      CALL  MUMPS_SOL_GET_NPIV_LIELL_IPOS ( STEP(INODE_CUR), KEEP,
     &     NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
      IF ( NPIV .NE. 0 ) THEN
        IF (MTYPE .EQ. 1 .OR. KEEP(50) .NE. 0) THEN
          LASTFSSBTR_DYN = IW( IPOS + NPIV )
        ELSE
          LASTFSSBTR_DYN = IW( IPOS+NPIV+LIELL )
        ENDIF
      ELSE
        IN = INODE_CUR
        DO WHILE (IN.GT. 0)
          IN = FILS(IN)
        ENDDO
        IF (IN .LT. 0) THEN
          INODE_CUR = -IN
          GOTO 10
        ELSE
          DO
            IF (INODE_CUR .EQ. INODE) THEN
              LASTFSSBTR_DYN = 0
              EXIT
            ENDIF
            INODE_CUR = FRERE(STEP(INODE_CUR))
            IF (INODE_CUR .GT. 0) THEN
              GOTO 10
            ELSE
              INODE_CUR = -INODE_CUR
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_COMPUTE_LASTFS_DYN
