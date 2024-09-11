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
      SUBROUTINE MUMPS_ADJUST_SIZE_LRGROUPS ( STEP, FILS, N, 
     &           ND_STEPS, NSTEPS, KEEP, 
     &           LRGROUPS, IFLAG, IERROR)
!$    USE OMP_LIB
      USE MUMPS_LR_COMMON,  ONLY : COMPUTE_BLR_VCS
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, NSTEPS
      INTEGER, INTENT(IN)    :: FILS(N), STEP(N), ND_STEPS(NSTEPS)
      INTEGER, INTENT(INOUT) :: KEEP(500)
      INTEGER, INTENT(INOUT) :: LRGROUPS(KEEP(280))
      INTEGER, INTENT(INOUT) :: IFLAG, IERROR
      INTEGER :: I, NASS, NFRONT, IN, 
     &           GROUP_SIZE2, NBGROUPS_local, LRGROUPS_SIGN, IG, 
     &           NBGROUPS, MAXSIZE_GROUPS,
     &           CurrentGroup, CurrentSize 
      INTEGER :: OriginalSize, NBSPLIT, SZ_FINAL
      INTEGER, DIMENSION(:), ALLOCATABLE :: GROUP_SIZE
      INTEGER  :: allocok
      IF (KEEP(494).EQ.0) RETURN
      IF (KEEP(280).NE.N) THEN
        WRITE(6,*) " Internal error in MUMPS_ADJUST_SIZE_LRGROUPS ",
     &             "N, KEEP(280) =", N, KEEP(280)
        RETURN
      ENDIF
      NBGROUPS = 0
      DO I=1, N
        IN = abs(LRGROUPS(I))
        IF (IN.GT.NBGROUPS) NBGROUPS=IN
      END DO
      ALLOCATE( GROUP_SIZE(NBGROUPS), stat=allocok )
      IF (allocok > 0) THEN
        IFLAG  = -7
        IERROR = NBGROUPS
        RETURN
      ENDIF
      GROUP_SIZE = 0
      MAXSIZE_GROUPS   = 0
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1) 
!$OMP&         PRIVATE( I, NBGROUPS_local, NASS, IN, NFRONT,
!$OMP&                  GROUP_SIZE2, CurrentGroup, LRGROUPS_SIGN,
!$OMP&                  IG, CurrentSize, OriginalSize, NBSPLIT, 
!$OMP&                  SZ_FINAL
!$OMP&                )
!$OMP& REDUCTION(max:MAXSIZE_GROUPS
!$OMP&                )
      DO I= 1, N
        IF (STEP(I).LE.0) CYCLE 
        NASS = 0
        IN   = I
        DO WHILE (IN.GT.0)
          CurrentGroup = abs(LRGROUPS(IN))
          GROUP_SIZE(CurrentGroup) = GROUP_SIZE(CurrentGroup) + 1
          NASS = NASS + 1
          IN = FILS(IN)
        ENDDO
        NFRONT = ND_STEPS(STEP(I))
        CALL COMPUTE_BLR_VCS(KEEP(472), GROUP_SIZE2, KEEP(488) , 
     &                       NASS, NFRONT, KEEP(35))
        IN = I
        DO WHILE (IN.GT.0)
          CurrentGroup  = LRGROUPS(IN)
          OriginalSize  = GROUP_SIZE(abs(CurrentGroup))
          NBSPLIT = (OriginalSize+GROUP_SIZE2-1)/ GROUP_SIZE2
          SZ_FINAL = (OriginalSize+NBSPLIT-1)/NBSPLIT
          MAXSIZE_GROUPS = max(SZ_FINAL, MAXSIZE_GROUPS)
          IF (CurrentGroup.LT.0) THEN 
             LRGROUPS_SIGN=-1
          ELSE
             LRGROUPS_SIGN=1
          ENDIF
          IG             = CurrentGroup 
          NBGROUPS_local = IG
          CurrentSize    = 0
          DO WHILE (IG .eq. CurrentGroup)
            LRGROUPS(IN) = NBGROUPS_local
            CurrentSize = CurrentSize + 1
            IF ( CurrentSize .GT. SZ_FINAL ) THEN
!$OMP ATOMIC CAPTURE
              NBGROUPS = NBGROUPS + 1
              NBGROUPS_local = NBGROUPS
!$OMP END ATOMIC 
              NBGROUPS_local = NBGROUPS_local*LRGROUPS_SIGN
              CurrentSize = 0
            ENDIF
            IN = FILS(IN)
            IF (IN.LE.0) EXIT
            IG =  LRGROUPS(IN)
           END DO
        END DO
      END DO
      KEEP(142) = MAXSIZE_GROUPS
      DEALLOCATE(GROUP_SIZE)
      RETURN
      END  SUBROUTINE MUMPS_ADJUST_SIZE_LRGROUPS
      SUBROUTINE MUMPS_STOP_ON_USER_REQUEST( KEEP, KEEP8, ICNTL, INFO,
     &                                       MYID )
      INTEGER    :: INFO(80), KEEP(500), ICNTL(60), MYID
      INTEGER(8) :: KEEP8(150)
      INTEGER :: ICNTL50
        CALL MUMPS_GETVAL_ADDR_C(ICNTL50, KEEP8(83))
      IF (ICNTL50 .eq. 1) THEN
        INFO(1)=-80
        INFO(2)=MYID
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_STOP_ON_USER_REQUEST
      SUBROUTINE MUMPS_BUILD_COMM_PARA_ANA (
     &      OPTION, N, COMM, MYID, COMM_NODES, MYID_NODES,
     &      NPROCS, NSLAVES,
     &      KEEP, 
     &      COMM_PARAORD, NPROCS_PARAORD,  
     &      COMM_PARAORD_ALLOCATED,
     &      COMM_PARASYMB, NPROCS_PARASYMB, 
     &      COMM_PARASYMB_ALLOCATED, 
     &      ICNTL, INFO
     &           )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT(IN) :: OPTION, N, COMM, COMM_NODES, 
     &                       MYID, MYID_NODES, NPROCS, NSLAVES
      INTEGER, INTENT(IN) :: ICNTL(60)
      INTEGER, INTENT(INOUT) :: KEEP(500)  
      INTEGER, INTENT(INOUT) :: INFO(80)  
      INTEGER, INTENT(OUT) :: COMM_PARAORD, NPROCS_PARAORD, 
     &                        COMM_PARASYMB, NPROCS_PARASYMB
      LOGICAL, INTENT(OUT) :: COMM_PARAORD_ALLOCATED, 
     &                        COMM_PARASYMB_ALLOCATED
      INTEGER  :: IERR, BASE, I, COLOR
      INTEGER  :: WORKERS, WTMP
      LOGICAL  :: IDO
      INTEGER :: MINPROCSPERNODE, MAXPROCSPERNODE
      INTEGER :: MYNODEID, NNODES
      INTEGER :: WORKERSPERNODE, WORKERSLEFT, WORKERSONMYNODE
      INTEGER, PARAMETER :: MAXNODE_SCOTCH=16, MAXNODE_METIS=64,  
     &                      ROOT=0
      IF (KEEP(339).NE.0) THEN
        IF (NPROCS.EQ.1) GOTO 100
        IF(N.LE.100) THEN
           WORKERS = 2
        ELSE
           WORKERS = min(NPROCS,N/16)
        END IF
        I    = 1
        WTMP = 1
        DO
           IF (I .GT. WORKERS) EXIT
           WTMP = I
           I = I*2
        END DO
        WORKERS = WTMP
        IF (WORKERS.LE.1) GOTO 100
        CALL MPI_ALLREDUCE( KEEP(412), MINPROCSPERNODE, 1, MPI_INTEGER,
     &                      MPI_MIN, COMM, IERR )
        CALL MPI_ALLREDUCE( KEEP(412), MAXPROCSPERNODE, 1, MPI_INTEGER,
     &                      MPI_MAX, COMM, IERR )
        IF (MINPROCSPERNODE .NE. MAXPROCSPERNODE ) THEN
         GOTO  100
        ENDIF
        IF ( KEEP(410) .NE. MPI_COMM_NULL ) THEN
          CALL MPI_COMM_RANK(KEEP(410), MYNODEID, IERR)
          CALL MPI_COMM_SIZE(KEEP(410), NNODES, IERR)
        ENDIF
        CALL MPI_BCAST( MYNODEID, 1, MPI_INTEGER,
     &                  ROOT,  KEEP(411), IERR )
        CALL MPI_BCAST( NNODES, 1, MPI_INTEGER,
     &                  ROOT,  KEEP(411), IERR )
        IF (WORKERS .LT. NNODES ) THEN
          IF (MYNODEID .LT. WORKERS ) THEN
            WORKERSONMYNODE = 1
          ELSE
            WORKERSONMYNODE = 0
          ENDIF
        ELSE
          WORKERSPERNODE = WORKERS / NNODES
          WORKERSLEFT    = WORKERS - WORKERSPERNODE*NNODES
          WORKERSONMYNODE = WORKERSPERNODE
          IF (NNODES - MYNODEID .LE. WORKERSLEFT ) THEN
              WORKERSONMYNODE = WORKERSONMYNODE+1
          ENDIF
        ENDIF
        NPROCS_PARAORD = WORKERS
        IF ( KEEP(413) .LE. WORKERSONMYNODE - 1 ) THEN
          IDO = .TRUE.
          COLOR = 1
          COMM_PARAORD_ALLOCATED = .TRUE.
        ELSE
          IDO = .FALSE.
          COLOR = MPI_UNDEFINED
          COMM_PARAORD_ALLOCATED = .FALSE.
        ENDIF
        CALL MPI_COMM_SPLIT( COMM, COLOR, 0, COMM_PARAORD, IERR)
        COMM_PARASYMB           = COMM_PARAORD
        COMM_PARASYMB_ALLOCATED = .FALSE.
        NPROCS_PARASYMB         = NPROCS_PARAORD
        GOTO 500
      ENDIF
  100 CONTINUE
      BASE = NPROCS-NSLAVES
      COMM_PARAORD    = MPI_COMM_NULL
      NPROCS_PARAORD  = 0
      COMM_PARAORD_ALLOCATED = .FALSE.
      NPROCS_PARASYMB = NPROCS
      IF (OPTION.EQ.0) THEN
        IF (KEEP(245).EQ.1) THEN
#if defined(ptscotch)
          COMM_PARAORD    = COMM_NODES
          NPROCS_PARAORD  = NSLAVES
          COMM_PARAORD_ALLOCATED = .FALSE.
#else
         INFO(1)= -999
         GOTO 600
#endif
        ELSE IF (KEEP(245) .EQ. 2) THEN
#if defined(parmetis) || defined(parmetis3)
          IF(N.LE.100) THEN
             WORKERS = 2
          ELSE
             WORKERS = min(NSLAVES,N/16)
          END IF
          I=1
          DO
             IF (I .GT. WORKERS) EXIT
             NPROCS_PARAORD = I
             I = I*2
          END DO
          IDO = (MYID .GE. BASE) .AND.
     &         (MYID .LE. BASE+NPROCS_PARAORD-1)
          IF ( IDO ) THEN
             COLOR   = 1
             COMM_PARAORD_ALLOCATED = .TRUE.
          ELSE
             COLOR = MPI_UNDEFINED
             COMM_PARAORD_ALLOCATED = .FALSE.
          END IF
          CALL MPI_COMM_SPLIT( COMM, COLOR, 0, COMM_PARAORD, IERR )
          COMM_PARASYMB           = COMM_PARAORD
          COMM_PARASYMB_ALLOCATED = .FALSE.
          NPROCS_PARASYMB         = NPROCS_PARAORD
#else
         INFO(1)= -999
         GOTO 600
#endif
        ENDIF
      ELSE 
        call MUMPS_ABORT()
        GOTO 600
      ENDIF
      NPROCS_PARASYMB         = NPROCS_PARAORD+BASE
      IF (BASE.EQ.0) THEN
        COMM_PARASYMB           = COMM_PARAORD
        COMM_PARASYMB_ALLOCATED = .FALSE.
      ELSE
         IF ((MYID.EQ.0).OR.COMM_PARAORD.NE.MPI_COMM_NULL) THEN
          COLOR = 1
          COMM_PARASYMB_ALLOCATED = .TRUE.
         ELSE
          COLOR = MPI_UNDEFINED
          COMM_PARASYMB_ALLOCATED = .FALSE.
         ENDIF
        CALL MPI_COMM_SPLIT( COMM, COLOR, 0, COMM_PARASYMB, IERR )
      ENDIF
 500  CONTINUE
  600 CONTINUE
      RETURN
      END SUBROUTINE MUMPS_BUILD_COMM_PARA_ANA
      SUBROUTINE MUMPS_BUILD_PARAORD_to_idCOMM (
     &      COMM, MYID, KEEP, 
     &      COMM_PARASYMB, NPROCS_PARASYMB,
     &      COMM_PARAORD,
     &      NPROCS_PARAORD,  
     &      PARAORD_to_idCOMM, 
#if defined(AVOID_MPI_IN_PLACE)
     &      TMP,
#endif
     &      RKinSYMB_PROC0ORD, 
     &      RKinidCOMM_PROC0SYMB, NPROCS)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT(IN)  :: COMM, MYID 
      INTEGER, INTENT(IN)  :: COMM_PARAORD, NPROCS_PARAORD,
     &                        COMM_PARASYMB, NPROCS_PARASYMB 
      INTEGER, INTENT(IN)  :: KEEP(500) 
      INTEGER, INTENT(OUT) :: PARAORD_to_idCOMM(NPROCS_PARAORD)
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER              :: TMP(NPROCS_PARAORD) 
#endif
      INTEGER, INTENT(OUT) :: RKinSYMB_PROC0ORD, RKinidCOMM_PROC0SYMB
      INTEGER  :: idPARAORD, idPARASYMB
      INTEGER  :: IERR
      INTEGER  :: NPROCS
#if defined(AVOID_MPI_IN_PLACE)
      INTEGER :: TMP_INT
#endif
#if defined(AVOID_MPI_IN_PLACE)
      TMP(1:NPROCS_PARAORD)= -1
#else
      PARAORD_to_idCOMM(1:NPROCS_PARAORD) = -1
#endif
      IF (COMM_PARAORD.NE.MPI_COMM_NULL) THEN
        CALL MPI_COMM_RANK (COMM_PARAORD, idPARAORD, IERR)
#if defined(AVOID_MPI_IN_PLACE)
        TMP(idPARAORD+1) = MYID
#else
        PARAORD_to_idCOMM(idPARAORD+1) = MYID
#endif
      ENDIF
      CALL MPI_ALLREDUCE(
#if defined(AVOID_MPI_IN_PLACE)
     &     TMP,
#else
     &     MPI_IN_PLACE,
#endif
     &     PARAORD_to_idCOMM(1), 
     &     NPROCS_PARAORD, MPI_INTEGER, MPI_MAX, COMM, IERR)
      RKinSYMB_PROC0ORD = -1
      IF (COMM_PARAORD.NE.MPI_COMM_NULL) THEN
        CALL MPI_COMM_RANK (COMM_PARAORD, idPARAORD, IERR)
        CALL MPI_COMM_RANK (COMM_PARASYMB, idPARASYMB, IERR)
        IF (idPARAORD.EQ.0) RKinSYMB_PROC0ORD = idPARASYMB
      ENDIF
#if defined(AVOID_MPI_IN_PLACE)
      TMP_INT = RKinSYMB_PROC0ORD
      CALL MPI_ALLREDUCE(TMP_INT, RKinSYMB_PROC0ORD, 1,
     &     MPI_INTEGER, MPI_MAX, COMM, IERR)
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, RKinSYMB_PROC0ORD, 1,
     &     MPI_INTEGER, MPI_MAX, COMM, IERR)
#endif
      RKinidCOMM_PROC0SYMB=-1
      IF (COMM_PARASYMB.NE.MPI_COMM_NULL) THEN
        CALL MPI_COMM_RANK (COMM_PARASYMB, idPARASYMB, IERR)
        IF (idPARASYMB.EQ.0) RKinidCOMM_PROC0SYMB = MYID
      ENDIF
#if defined(AVOID_MPI_IN_PLACE)
      TMP_INT = RKinidCOMM_PROC0SYMB
      CALL MPI_ALLREDUCE(TMP_INT, RKinidCOMM_PROC0SYMB, 1,
     &     MPI_INTEGER, MPI_MAX, COMM, IERR)
#else
      CALL MPI_ALLREDUCE(MPI_IN_PLACE, RKinidCOMM_PROC0SYMB, 1,
     &     MPI_INTEGER, MPI_MAX, COMM, IERR)
#endif
      RETURN
      END SUBROUTINE MUMPS_BUILD_PARAORD_to_idCOMM 
      SUBROUTINE MUMPS_FIND_UNIT(IUNIT)
      IMPLICIT NONE
      INTEGER :: IUNIT
      INTEGER, PARAMETER :: UNIT_MIN = 10 
      INTEGER, PARAMETER :: UNIT_MAX = 500
      INTEGER :: I
      LOGICAL :: BUSY
      IUNIT = -1
      DO I = UNIT_MIN, UNIT_MAX
        INQUIRE(UNIT=I, OPENED=BUSY)
        IF ( .NOT. BUSY ) THEN
          IUNIT = I
          RETURN
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_FIND_UNIT
      SUBROUTINE MUMPS_MAKE1ROOT( N, FRERE, FILS, NFSIZ, THEROOT )
      IMPLICIT NONE
      INTEGER, intent( in    )  :: N
      INTEGER, intent( in    )  :: NFSIZ( N )
      INTEGER, intent( inout )  :: FRERE( N ), FILS( N )
      INTEGER, intent( out   )  :: THEROOT
      INTEGER INODE, IROOT, IFILS, IN, IROOTLAST, SIZE
      IROOT = -9999
      SIZE  = 0
      DO INODE = 1, N
        IF ( FRERE( INODE ) .EQ. 0 )  THEN
          IF ( NFSIZ( INODE ) .GT. SIZE ) THEN
            SIZE  = NFSIZ( INODE )
            IROOT = INODE
          END IF
        ENDIF
      END DO
      IN = IROOT
      DO WHILE ( FILS( IN ) .GT. 0 )
        IN = FILS( IN )
      END DO
      IROOTLAST = IN
      IFILS     = - FILS ( IN )
      DO INODE = 1, N
        IF ( FRERE( INODE ) .eq. 0 .and. INODE .ne. IROOT ) THEN
          IF ( IFILS .eq. 0 ) THEN
            FILS( IROOTLAST ) = - INODE
            FRERE( INODE )    = -IROOT
            IFILS             = INODE
          ELSE
            FRERE( INODE ) = -FILS( IROOTLAST )
            FILS( IROOTLAST ) = - INODE
          END IF
        END IF
      END DO
      THEROOT = IROOT
      RETURN
      END SUBROUTINE MUMPS_MAKE1ROOT
      INTEGER FUNCTION MUMPS_ENCODE_TPN_IPROC(TPN,IPROC,K199)
      INTEGER, INTENT(IN) :: TPN, IPROC, K199
      IF (K199 < 0) THEN
        MUMPS_ENCODE_TPN_IPROC = IPROC + ISHFT(TPN+1, 24)
      ELSE
        MUMPS_ENCODE_TPN_IPROC = (TPN-1)*K199+IPROC+1
      ENDIF
      RETURN
      END FUNCTION MUMPS_ENCODE_TPN_IPROC
      INTEGER FUNCTION MUMPS_TYPENODE_ROUGH(PROCINFO_INODE, K199)
      IMPLICIT NONE
      INTEGER K199 
      INTEGER PROCINFO_INODE
      IF (K199 < 0) THEN
        MUMPS_TYPENODE_ROUGH = ISHFT(PROCINFO_INODE,-24) - 1
      ELSE
        MUMPS_TYPENODE_ROUGH = (PROCINFO_INODE-1+2*K199)/K199 - 1
      ENDIF
      RETURN 
      END FUNCTION MUMPS_TYPENODE_ROUGH
      INTEGER FUNCTION MUMPS_TYPENODE(PROCINFO_INODE, K199)
      IMPLICIT NONE
      INTEGER K199 
      INTEGER PROCINFO_INODE, TPN
      IF (K199 < 0) THEN
        TPN = ISHFT(PROCINFO_INODE,-24) - 1
        IF (TPN < 1 ) THEN
          TPN = 1
        ELSE IF (TPN.GE.4) THEN
          TPN = 2
        ENDIF
      ELSE
        IF (PROCINFO_INODE <= K199 ) THEN
          TPN = 1
        ELSE
          TPN = (PROCINFO_INODE-1+2*K199)/K199 - 1
          IF ( TPN .LT. 1 ) TPN = 1
          IF (TPN.EQ.4.OR.TPN.EQ.5.OR.TPN.EQ.6) TPN = 2
        END IF
      END IF
      MUMPS_TYPENODE = TPN
      RETURN 
      END FUNCTION MUMPS_TYPENODE
      SUBROUTINE MUMPS_TYPEANDPROCNODE( TPN,
     &  MUMPS_PROCNODE, PROCINFO_INODE, K199 )
      INTEGER, INTENT(IN) :: K199, PROCINFO_INODE
      INTEGER, intent(out) :: TPN, MUMPS_PROCNODE
      IF (K199 < 0 ) THEN
        MUMPS_PROCNODE=iand(PROCINFO_INODE,
#if defined(MUMPS_NOF2003)
     &         16777215
#else
     &         int(B"111111111111111111111111")
#endif
     &  )
        TPN = ISHFT(PROCINFO_INODE,-24) - 1
        IF (TPN < 1 ) THEN
          TPN = 1
        ELSE IF (TPN.GE.4) THEN
          TPN = 2
        ENDIF
      ELSE
        IF (K199 == 1) THEN
          MUMPS_PROCNODE = 0
          IF (PROCINFO_INODE <= K199) THEN
            TPN = 1
          ELSE
            TPN = 3
          ENDIF
        ELSE
          TPN = (PROCINFO_INODE-1+2*K199)/K199-1
          MUMPS_PROCNODE = (PROCINFO_INODE-1+2*K199)-
     &                      (TPN+1)*K199
          IF (TPN .LT. 1) THEN
            TPN = 1
          ELSE IF (TPN .ge. 4) THEN
            TPN = 2
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_TYPEANDPROCNODE
      INTEGER FUNCTION MUMPS_PROCNODE(PROCINFO_INODE, K199)
      IMPLICIT NONE
      INTEGER K199 
      INTEGER PROCINFO_INODE
      IF ( K199 < 0 ) THEN
        MUMPS_PROCNODE=iand(PROCINFO_INODE,
#if defined(MUMPS_NOF2003)
     &          16777215 
#else
     &         int(B"111111111111111111111111")
#endif
     &  )
      ELSE
        IF (K199 == 1) THEN
          MUMPS_PROCNODE = 0
        ELSE
          MUMPS_PROCNODE=mod(2*K199+PROCINFO_INODE-1,K199)
        END IF
      ENDIF
      RETURN
      END FUNCTION MUMPS_PROCNODE
      INTEGER FUNCTION MUMPS_TYPESPLIT (PROCINFO_INODE, K199)
      IMPLICIT NONE
      INTEGER, intent(in) ::  K199 
      INTEGER PROCINFO_INODE, TPN
      IF (K199 < 0) THEN
        TPN = ishft(PROCINFO_INODE,-24) - 1
        IF (TPN < 1 ) TPN = 1
      ELSE
        IF (PROCINFO_INODE <= K199 ) THEN
           TPN = 1
        ELSE
          TPN = (PROCINFO_INODE-1+2*K199)/K199 - 1
          IF ( TPN .LT. 1 ) TPN = 1
        ENDIF
      ENDIF
      MUMPS_TYPESPLIT = TPN
      RETURN
      END FUNCTION MUMPS_TYPESPLIT
      LOGICAL FUNCTION MUMPS_ROOTSSARBR( PROCINFO_INODE, K199 )
      IMPLICIT NONE
      INTEGER K199
      INTEGER TPN, PROCINFO_INODE
      IF (K199 < 0) THEN
        TPN = ishft(PROCINFO_INODE,-24) - 1
      ELSE
        TPN = (PROCINFO_INODE-1+2*K199)/K199 - 1
      ENDIF
      MUMPS_ROOTSSARBR = ( TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_ROOTSSARBR
      LOGICAL FUNCTION MUMPS_INSSARBR( PROCINFO_INODE, K199 )
      IMPLICIT NONE
      INTEGER K199
      INTEGER TPN, PROCINFO_INODE
      IF (K199 < 0) THEN
        TPN = ishft(PROCINFO_INODE,-24) - 1
      ELSE
        TPN = (PROCINFO_INODE-1+K199+K199)/K199 - 1
      ENDIF
      MUMPS_INSSARBR = ( TPN .eq. -1 )
      RETURN 
      END FUNCTION MUMPS_INSSARBR
      LOGICAL FUNCTION MUMPS_IN_OR_ROOT_SSARBR
     &        ( PROCINFO_INODE, K199 )
      IMPLICIT NONE
      INTEGER K199
      INTEGER TPN, PROCINFO_INODE
      IF (K199 < 0) THEN
        TPN = ishft(PROCINFO_INODE,-24) - 1
      ELSE
        TPN = (PROCINFO_INODE-1+K199+K199)/K199 - 1
      ENDIF
      MUMPS_IN_OR_ROOT_SSARBR =
     &           ( TPN .eq. -1 .OR. TPN .eq. 0 )
      RETURN
      END FUNCTION MUMPS_IN_OR_ROOT_SSARBR
            SUBROUTINE MUMPS_SET_SSARBR_DAD(
     &           SSARBR, INODE, DAD, N,
     &           KEEP28,
     &           STEP, PROCNODE_STEPS, K199)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: N, KEEP28, K199, INODE
            INTEGER, INTENT(IN) :: DAD(KEEP28), PROCNODE_STEPS(KEEP28)
            INTEGER, INTENT(IN) :: STEP(N)
            LOGICAL, INTENT(OUT) :: SSARBR
            INTEGER :: DADINODE, TYPEDAD
            LOGICAL, EXTERNAL :: MUMPS_INSSARBR
            INTEGER, EXTERNAL :: MUMPS_TYPENODE
            SSARBR   = .FALSE.
            DADINODE = DAD(STEP(INODE))
            IF (DADINODE .NE. 0) THEN
              TYPEDAD  = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(DADINODE)),
     &                                  K199)
              IF (TYPEDAD.EQ.1) THEN
                SSARBR=MUMPS_INSSARBR(PROCNODE_STEPS(STEP(DADINODE)),
     &                                K199)
              ENDIF
            ENDIF
            RETURN
            END SUBROUTINE MUMPS_SET_SSARBR_DAD
      LOGICAL FUNCTION MUMPS_I_AM_CANDIDATE( MYID, SLAVEF, INODE,
     &                 NMB_PAR2, ISTEP_TO_INIV2 , K71, STEP, N, 
     &                 CANDIDATES, KEEP24 )
      IMPLICIT NONE
      INTEGER MYID, SLAVEF, INODE, NMB_PAR2, KEEP24, I
      INTEGER K71, N
      INTEGER ISTEP_TO_INIV2 ( K71 ), STEP ( N )
      INTEGER CANDIDATES(SLAVEF+1, max(NMB_PAR2,1))
      INTEGER NCAND, POSINODE
      MUMPS_I_AM_CANDIDATE = .FALSE.
      IF (KEEP24 .eq. 0) RETURN
      POSINODE = ISTEP_TO_INIV2 ( STEP (INODE) )
      NCAND = CANDIDATES( SLAVEF+1, POSINODE )
      DO I = 1, NCAND
        IF (MYID .EQ. CANDIDATES( I, POSINODE ))
     &     MUMPS_I_AM_CANDIDATE = .TRUE.
      END DO
      RETURN
      END FUNCTION MUMPS_I_AM_CANDIDATE
      SUBROUTINE MUMPS_SECDEB(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()
      RETURN
      END SUBROUTINE MUMPS_SECDEB
      SUBROUTINE MUMPS_SECFIN(T)
      DOUBLE PRECISION T
      DOUBLE PRECISION MPI_WTIME
      EXTERNAL MPI_WTIME
      T=MPI_WTIME()-T
      RETURN
      END SUBROUTINE MUMPS_SECFIN
      SUBROUTINE MUMPS_SORT_DOUBLES( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_DOUBLES
      SUBROUTINE MUMPS_SORT_DOUBLES_DEC( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      DOUBLE PRECISION VAL( N )
      INTEGER I, ISWAP
      DOUBLE PRECISION SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
          IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
            DONE = .FALSE.
            ISWAP = ID( I )
            ID ( I ) = ID ( I + 1 )
            ID ( I + 1 ) = ISWAP
            SWAP = VAL( I )
            VAL( I ) = VAL( I + 1 )
            VAL( I + 1 ) = SWAP
          END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_DOUBLES_DEC
      SUBROUTINE MUMPS_MEM_CENTRALIZE(MYID, COMM, INFO, INFOG, IRANK)
      IMPLICIT NONE
      INTEGER MYID, COMM, IRANK, INFO, INFOG(2)
      INCLUDE 'mpif.h'
      INTEGER IERR_MPI, MASTER
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      INTEGER(4) :: TEMP1(2),TEMP2(2)
#else
      INTEGER :: TEMP1(2),TEMP2(2)
#endif
      PARAMETER( MASTER = 0 )
      CALL MPI_REDUCE( INFO, INFOG(1), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, COMM, IERR_MPI )
      CALL MPI_REDUCE( INFO, INFOG(2), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, COMM, IERR_MPI )
      TEMP1(1) = INFO
      TEMP1(2) = MYID
      CALL MPI_REDUCE( TEMP1, TEMP2, 1, MPI_2INTEGER,
     &                 MPI_MAXLOC, MASTER, COMM, IERR_MPI )
      IF ( MYID.eq. MASTER ) THEN
         IF ( INFOG(1) .ne. TEMP2(1) ) THEN
          write(*,*) 'Error in MUMPS_MEM_CENTRALIZE'
          CALL MUMPS_ABORT()
        END IF
        IRANK    = TEMP2(2)
      ELSE
        IRANK    = -1
      END IF
      RETURN
      END SUBROUTINE MUMPS_MEM_CENTRALIZE
      INTEGER FUNCTION MUMPS_GET_POOL_LENGTH
     &        (MAX_ACTIVE_NODES,KEEP,KEEP8)
      IMPLICIT NONE
      INTEGER MAX_ACTIVE_NODES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      MUMPS_GET_POOL_LENGTH = MAX_ACTIVE_NODES + 1 + 3
      RETURN
      END FUNCTION MUMPS_GET_POOL_LENGTH
      SUBROUTINE MUMPS_INIT_POOL_DIST_BWD(N,
     &           nb_prun_roots, Pruned_Roots,
     &           MYROOT, MYID_NODES,
     &           KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &           IPOOL, LPOOL )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N, MYID_NODES, LPOOL, nb_prun_roots
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT(IN)  :: STEP(N)
      INTEGER, INTENT(IN)  :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT(IN)  :: Pruned_Roots(nb_prun_roots)
      INTEGER, INTENT(OUT) :: MYROOT
      INTEGER, INTENT(OUT) :: IPOOL(LPOOL)
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: I, INODE
      MYROOT = 0
      DO I = nb_prun_roots, 1, -1
        INODE = Pruned_Roots(I)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      KEEP(199)) .EQ. MYID_NODES) THEN
          MYROOT = MYROOT + 1
          IPOOL(MYROOT) = INODE
        ENDIF
      END DO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_BWD
      SUBROUTINE MUMPS_INIT_POOL_DIST_BWD_L0(N,
     &           nb_prun_roots, Pruned_Roots,
     &           MYROOT, MYID_NODES,
     &           KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &           IPOOL, LPOOL, TO_PROCESS )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: N, MYID_NODES, LPOOL, nb_prun_roots
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT(IN)  :: STEP(N)
      INTEGER, INTENT(IN)  :: PROCNODE_STEPS(KEEP(28))
      LOGICAL, INTENT(IN)  :: TO_PROCESS(KEEP(28))
      INTEGER, INTENT(IN)  :: Pruned_Roots(nb_prun_roots)
      INTEGER, INTENT(OUT) :: MYROOT
      INTEGER, INTENT(OUT) :: IPOOL(LPOOL)
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: I, INODE
      MYROOT = 0
      DO I = nb_prun_roots, 1, -1
        INODE = Pruned_Roots(I)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      KEEP(199)) .EQ. MYID_NODES) THEN
          IF ( TO_PROCESS(STEP(INODE)) ) THEN
            MYROOT = MYROOT + 1
            IPOOL(MYROOT) = INODE
          ENDIF
        ENDIF
      END DO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_BWD_L0
      SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWD(N, MYROOT, MYID_NODES,
     &           NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &           IPOOL, LPOOL )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, MYID_NODES, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(KEEP(28)), NA(LNA)
      INTEGER, INTENT(OUT) :: IPOOL(LPOOL)
      INTEGER, INTENT(OUT) :: MYROOT
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: NBLEAF, NBROOT, I, INODE
      NBLEAF = NA(1)
      NBROOT = NA(2)
      MYROOT = 0
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      KEEP(199)) .EQ. MYID_NODES) THEN
          MYROOT = MYROOT + 1
          IPOOL(MYROOT) = INODE
        ENDIF
      END DO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWD
      SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWD_L0(N, MYROOT, MYID_NODES,
     &           NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &           IPOOL, LPOOL, L0_OMP_MAPPING )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, MYID_NODES, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT(IN) :: STEP(N)
      INTEGER, INTENT(IN) :: PROCNODE_STEPS(KEEP(28)), NA(LNA)
      INTEGER, INTENT(IN) :: L0_OMP_MAPPING(KEEP(28))
      INTEGER, INTENT(OUT) :: IPOOL(LPOOL)
      INTEGER, INTENT(OUT) :: MYROOT
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: NBLEAF, NBROOT, I, INODE
      NBLEAF = NA(1)
      NBROOT = NA(2)
      MYROOT = 0
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      KEEP(199)) .EQ. MYID_NODES) THEN
          IF ( L0_OMP_MAPPING(STEP(INODE)).EQ.0 ) THEN
            MYROOT = MYROOT + 1
            IPOOL(MYROOT) = INODE
          ENDIF
        ENDIF
      END DO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWD_L0
      SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWDL0ES(N, MYROOT,
     &           MYID_NODES,
     &           NA, LNA, KEEP, KEEP8, STEP, PROCNODE_STEPS,
     &           IPOOL, LPOOL, L0_OMP_MAPPING, TO_PROCESS )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, MYID_NODES, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT(IN)  :: STEP(N)
      INTEGER, INTENT(IN)  :: PROCNODE_STEPS(KEEP(28)), NA(LNA)
      INTEGER, INTENT(IN)  :: L0_OMP_MAPPING(KEEP(28))
      INTEGER, INTENT(OUT) :: IPOOL(LPOOL)
      INTEGER, INTENT(OUT) :: MYROOT
      LOGICAL, INTENT(IN)  :: TO_PROCESS( KEEP(28) )
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: NBLEAF, NBROOT, I, INODE
      NBLEAF = NA(1)
      NBROOT = NA(2)
      MYROOT = 0
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      KEEP(199)) .EQ. MYID_NODES) THEN
          IF ( L0_OMP_MAPPING(STEP(INODE)).EQ.0 ) THEN
            IF ( TO_PROCESS( STEP(INODE) ) ) THEN
              MYROOT = MYROOT + 1
              IPOOL(MYROOT) = INODE
            ENDIF
          ENDIF
        ENDIF
      END DO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_NA_BWDL0ES
      SUBROUTINE MUMPS_INIT_POOL_DIST(N, LEAF,
     &           MYID_NODES,
     &           K199, NA, LNA, KEEP,KEEP8, STEP,
     &           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, MYID_NODES,
     &        K199, LPOOL, LNA
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), NA(LNA),
     &        IPOOL(LPOOL)
      INTEGER NBLEAF, INODE, I
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      NBLEAF = NA(1)
      LEAF = 1
      DO I = 1, NBLEAF
        INODE = NA(I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
     &   .EQ.MYID_NODES) THEN
           IPOOL(LEAF) = INODE
           LEAF        = LEAF + 1
          ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST
      SUBROUTINE MUMPS_INIT_POOL_DIST_NONA
     &           (N, LEAF, MYID_NODES,
     &           LLEAVES, LEAVES, KEEP,KEEP8, STEP,
     &           PROCNODE_STEPS, IPOOL, LPOOL)
      IMPLICIT NONE
      INTEGER N, LEAF, MYID_NODES,
     &        LPOOL, LLEAVES
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER STEP(N)
      INTEGER PROCNODE_STEPS(KEEP(28)), LEAVES(LLEAVES),
     &        IPOOL(LPOOL)
      INTEGER I, INODE
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      LEAF = 1
      DO I = 1, LLEAVES
        INODE = LEAVES(I)
        IF ( MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),KEEP(199))
     &   .EQ.MYID_NODES ) THEN
          IPOOL( LEAF ) = INODE
          LEAF = LEAF + 1
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_INIT_POOL_DIST_NONA
      SUBROUTINE MUMPS_INIT_NROOT_DIST(N, NBROOT,
     &           NROOT_LOC, MYID_NODES,
     &           SLAVEF, NA, LNA, KEEP, STEP,
     &           PROCNODE_STEPS)
      IMPLICIT NONE
      INTEGER, INTENT( OUT ) :: NROOT_LOC 
      INTEGER, INTENT( OUT ) :: NBROOT 
      INTEGER, INTENT( IN ) :: KEEP( 500 )
      INTEGER, INTENT( IN ) :: SLAVEF
      INTEGER, INTENT( IN ) :: N
      INTEGER, INTENT( IN ) :: STEP(N)
      INTEGER, INTENT( IN ) :: LNA
      INTEGER, INTENT( IN ) :: NA(LNA)
      INTEGER, INTENT( IN ) :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT( IN ) :: MYID_NODES
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      INTEGER :: INODE, I, NBLEAF
      NBLEAF = NA(1)
      NBROOT = NA(2)
      NROOT_LOC = 0
      DO I = 1, NBROOT
        INODE = NA(I+2+NBLEAF)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &    KEEP(199)).EQ.MYID_NODES) THEN
            NROOT_LOC = NROOT_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_INIT_NROOT_DIST
      SUBROUTINE MUMPS_NBLOCAL_ROOTS_OR_LEAVES
     &           (N, NBRORL, RORL_LIST,
     &           NRORL_LOC, MYID_NODES,
     &           SLAVEF, KEEP, STEP,
     &           PROCNODE_STEPS)
      IMPLICIT NONE
      INTEGER, INTENT( OUT ) :: NRORL_LOC 
      INTEGER, INTENT( IN ) :: NBRORL 
      INTEGER, INTENT( IN ) :: RORL_LIST(NBRORL)
      INTEGER, INTENT( IN ) :: KEEP( 500 )
      INTEGER, INTENT( IN ) :: SLAVEF
      INTEGER, INTENT( IN ) :: N
      INTEGER, INTENT( IN ) :: STEP(N)
      INTEGER, INTENT( IN ) :: PROCNODE_STEPS(KEEP(28))
      INTEGER, INTENT( IN ) :: MYID_NODES
      INTEGER I, INODE
      INTEGER, EXTERNAL :: MUMPS_PROCNODE
      NRORL_LOC = 0
      DO I = 1, NBRORL
        INODE = RORL_LIST(I)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &    KEEP(199)).EQ.MYID_NODES) THEN
            NRORL_LOC = NRORL_LOC + 1
        END IF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_NBLOCAL_ROOTS_OR_LEAVES
      LOGICAL FUNCTION MUMPS_COMPARE_TAB(TAB1,TAB2,LEN1,LEN2)
      IMPLICIT NONE
      INTEGER LEN1 , LEN2 ,I
      INTEGER TAB1(LEN1)
      INTEGER TAB2(LEN2)
      MUMPS_COMPARE_TAB=.FALSE.
      IF(LEN1 .NE. LEN2) THEN
         RETURN
      ENDIF
      DO I=1 , LEN1
         IF(TAB1(I) .NE. TAB2(I)) THEN
            RETURN
         ENDIF
      ENDDO
      MUMPS_COMPARE_TAB=.TRUE.
      RETURN
      END FUNCTION MUMPS_COMPARE_TAB
      SUBROUTINE MUMPS_SORT_INT( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_INT
      SUBROUTINE MUMPS_SORT_INT_DEC( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER VAL( N )
      INTEGER I, ISWAP
      INTEGER SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .LT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_INT_DEC
      SUBROUTINE MUMPS_SORT_INT8( N, VAL, ID )
      INTEGER N
      INTEGER ID( N )
      INTEGER(8) :: VAL( N )
      INTEGER I, ISWAP
      INTEGER(8) SWAP
      LOGICAL DONE
      DONE = .FALSE.
      DO WHILE ( .NOT. DONE )
        DONE = .TRUE.
        DO I = 1, N - 1
           IF ( VAL( I ) .GT. VAL( I + 1 ) ) THEN
              DONE = .FALSE.
              ISWAP = ID( I )
              ID ( I ) = ID ( I + 1 )
              ID ( I + 1 ) = ISWAP
              SWAP = VAL( I )
              VAL( I ) = VAL( I + 1 )
              VAL( I + 1 ) = SWAP
           END IF
        END DO
      END DO
      RETURN
      END SUBROUTINE MUMPS_SORT_INT8
      SUBROUTINE MUMPS_ABORT()
#if defined(PRINT_BACKTRACE_ON_ABORT) 
#if defined(__INTEL_COMPILER)
      USE IFCORE
#endif
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
#if defined(__GFORTRAN__)
      CALL BACKTRACE()
#endif
#if defined(__INTEL_COMPILER)
!$OMP CRITICAL(MUMPS_TRACEBACKQQ)
      CALL TRACEBACKQQ("MUMPS_ABORT calls TRACEBACKQQ:",
     &                 user_exit_code=-1)
!$OMP END CRITICAL(MUMPS_TRACEBACKQQ)
#endif
#else
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER IERR, IERRCODE
#endif
      IERRCODE = -99
      CALL MPI_ABORT(MPI_COMM_WORLD, IERRCODE, IERR)
      RETURN
      END SUBROUTINE MUMPS_ABORT
      SUBROUTINE MUMPS_GET_PERLU(KEEP12,ICNTL14,
     &     KEEP50,KEEP54,ICNTL6,ICNTL8)
      IMPLICIT NONE
      INTEGER, intent(out)::KEEP12
      INTEGER, intent(in)::ICNTL14,KEEP50,KEEP54,ICNTL6,ICNTL8
      KEEP12 = ICNTL14 
      RETURN
      END SUBROUTINE MUMPS_GET_PERLU
      SUBROUTINE MUMPS_REDUCEI8( IN, OUT, MPI_OP, ROOT, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER ROOT, COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_REDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, ROOT, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_REDUCEI8
      SUBROUTINE MUMPS_ALLREDUCEI8( IN, OUT, MPI_OP, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER COMM, MPI_OP
      INTEGER(8) IN, OUT
      INTEGER IERR
      DOUBLE PRECISION DIN, DOUT
      DIN =dble(IN)
      DOUT=0.0D0
      CALL MPI_ALLREDUCE(DIN, DOUT, 1, MPI_DOUBLE_PRECISION,
     &                   MPI_OP, COMM, IERR)
      OUT=int(DOUT,kind=8)
      RETURN
      END SUBROUTINE MUMPS_ALLREDUCEI8
      SUBROUTINE MUMPS_SETI8TOI4(I8, I)
      IMPLICIT NONE
      INTEGER   , INTENT(OUT) :: I
      INTEGER(8), INTENT(IN)  :: I8
      IF ( I8 .GT. int(huge(I),8) ) THEN
        I = -int(I8/1000000_8,kind(I))
      ELSE
        I = int(I8,kind(I))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SETI8TOI4
      SUBROUTINE MUMPS_ABORT_ON_OVERFLOW(I8, STRING)
      IMPLICIT NONE
      INTEGER(8), INTENT(IN) :: I8
      CHARACTER(*), INTENT(IN) :: STRING
      INTEGER I
      IF ( I8 .GT. int(huge(I),8)) THEN
        WRITE(*,*) STRING
        CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_ABORT_ON_OVERFLOW
      SUBROUTINE MUMPS_SET_IERROR( SIZE8, IERROR  )
      INTEGER(8), INTENT(IN) :: SIZE8
      INTEGER, INTENT(OUT) :: IERROR
      CALL MUMPS_SETI8TOI4(SIZE8, IERROR)
      RETURN
      END SUBROUTINE MUMPS_SET_IERROR
      SUBROUTINE MUMPS_STOREI8(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(in)  :: I8
      INTEGER,    intent(out) :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( I8 .LT. IDIV ) THEN
        INT_ARRAY(1) = 0
        INT_ARRAY(2) = int(I8)
      ELSE
        INT_ARRAY(1) = int(I8 / IDIV)
        INT_ARRAY(2) = int(mod(I8,IDIV))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_STOREI8
      SUBROUTINE MUMPS_GETI8(I8, INT_ARRAY)
      IMPLICIT NONE
      INTEGER(8), intent(out)  :: I8
      INTEGER,    intent(in)  :: INT_ARRAY(2)
      INTEGER(kind(0_4)) :: I32
      INTEGER(8) :: IDIV, IPAR
      PARAMETER (IPAR=int(huge(I32),8))
      PARAMETER (IDIV=IPAR+1_8)
      IF ( INT_ARRAY(1) .EQ. 0 ) THEN
        I8=int(INT_ARRAY(2),8)
      ELSE
        I8=int(INT_ARRAY(1),8)*IDIV+int(INT_ARRAY(2),8)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_GETI8
      SUBROUTINE MUMPS_ADDI8TOARRAY( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_GETI8(I8TMP, INT_ARRAY)
      I8TMP = I8TMP + I8
      CALL MUMPS_STOREI8(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_ADDI8TOARRAY
      SUBROUTINE MUMPS_SUBTRI8TOARRAY( INT_ARRAY, I8 )
      IMPLICIT NONE
      INTEGER(8), intent(in) :: I8
      INTEGER, intent(inout) :: INT_ARRAY(2)
      INTEGER(8) :: I8TMP
      CALL MUMPS_GETI8(I8TMP, INT_ARRAY)
      I8TMP = I8TMP - I8
      CALL MUMPS_STOREI8(I8TMP, INT_ARRAY)
      RETURN
      END SUBROUTINE MUMPS_SUBTRI8TOARRAY
      FUNCTION MUMPS_SEQANA_AVAIL(ICNTL7)
      LOGICAL :: MUMPS_SEQANA_AVAIL
      INTEGER, INTENT(IN) :: ICNTL7
      LOGICAL :: SCOTCH=.FALSE.
      LOGICAL :: METIS =.FALSE.
#if defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
      METIS = .TRUE.
#endif
#if defined(scotch) || defined(ptscotch)
      SCOTCH = .TRUE.
#endif
      IF ( ICNTL7 .LT. 0 .OR. ICNTL7 .GT. 7 ) THEN
        MUMPS_SEQANA_AVAIL = .FALSE.
      ELSE
        MUMPS_SEQANA_AVAIL = .TRUE.
      ENDIF
      IF ( ICNTL7 .EQ. 5 ) MUMPS_SEQANA_AVAIL = METIS
      IF ( ICNTL7 .EQ. 3 ) MUMPS_SEQANA_AVAIL = SCOTCH
      RETURN
      END FUNCTION MUMPS_SEQANA_AVAIL
      FUNCTION MUMPS_PARANA_AVAIL(WHICH)
      LOGICAL :: MUMPS_PARANA_AVAIL
      CHARACTER :: WHICH*(*)
      LOGICAL :: PTSCOTCH=.FALSE., PARMETIS=.FALSE.
#if defined(ptscotch)
      PTSCOTCH = .TRUE.
#endif
#if defined(parmetis) || defined(parmetis3)
      PARMETIS = .TRUE.
#endif
      SELECT CASE(WHICH)
      CASE('ptscotch','PTSCOTCH')
         MUMPS_PARANA_AVAIL = PTSCOTCH
      CASE('parmetis','PARMETIS')
         MUMPS_PARANA_AVAIL = PARMETIS
      CASE('both','BOTH')
         MUMPS_PARANA_AVAIL = PTSCOTCH .AND. PARMETIS
      CASE('any','ANY')
         MUMPS_PARANA_AVAIL = PTSCOTCH .OR. PARMETIS
      CASE default
         write(*,'("Invalid input in MUMPS_PARANA_AVAIL")')
      END SELECT
      RETURN
      END FUNCTION MUMPS_PARANA_AVAIL
      SUBROUTINE MUMPS_SORT_STEP(N,FRERE,STEP,FILS,
     &     NA,LNA,NE,ND,DAD,LDAD,USE_DAD,
     &     NSTEPS,INFO,LP,
     &     PROCNODE,SLAVEF
     &     )
      IMPLICIT NONE
      INTEGER N, NSTEPS, LNA, LP,LDAD
      INTEGER FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER NA(LNA), NE(NSTEPS), ND(NSTEPS)
      INTEGER DAD(LDAD)
      LOGICAL USE_DAD
      INTEGER INFO(80)
      INTEGER SLAVEF,PROCNODE(NSTEPS)
      INTEGER  POSTORDER,TMP_SWAP
      INTEGER, DIMENSION (:), ALLOCATABLE :: STEP_TO_NODE
      INTEGER, DIMENSION (:), ALLOCATABLE :: IPOOL,TNSTK
      INTEGER I,II,allocok
      INTEGER NBLEAF,NBROOT,LEAF,IN,INODE,IFATH
      POSTORDER=1
      NBLEAF = NA(1)
      NBROOT = NA(2)
      ALLOCATE( IPOOL(NBLEAF), TNSTK(NSTEPS), stat=allocok )
      IF (allocok > 0) THEN
        IF ( LP .GT. 0 )
     &    WRITE(LP,*)'Memory allocation error in MUMPS_SORT_STEP'
        INFO(1)=-7
        INFO(2)=NSTEPS
        RETURN
      ENDIF
      DO I=1,NSTEPS
         TNSTK(I) = NE(I)
      ENDDO
      ALLOCATE(STEP_TO_NODE(NSTEPS),stat=allocok)
      IF (allocok > 0) THEN
         IF ( LP .GT. 0 )
     &        WRITE(LP,*)'Memory allocation error in
     &MUMPS_SORT_STEP'
         INFO(1)=-7
         INFO(2)=NSTEPS
         RETURN
      ENDIF
      DO I=1,N
         IF(STEP(I).GT.0)THEN
            STEP_TO_NODE(STEP(I))=I
         ENDIF
      ENDDO
      IPOOL(1:NBLEAF)=NA(3:2+NBLEAF)
      LEAF = NBLEAF + 1
 91   CONTINUE
      IF (LEAF.NE.1) THEN
         LEAF = LEAF -1
         INODE = IPOOL(LEAF)
      ENDIF
 96   CONTINUE
      IF (USE_DAD) THEN
         IFATH = DAD( STEP(INODE) )
      ELSE
         IN = INODE
 113     IN = FRERE(IN)
         IF (IN.GT.0) GO TO 113
         IFATH = -IN
      ENDIF
      TMP_SWAP=FRERE(STEP(INODE))
      FRERE(STEP(INODE))=FRERE(POSTORDER)
      FRERE(POSTORDER)=TMP_SWAP
      TMP_SWAP=ND(STEP(INODE))
      ND(STEP(INODE))=ND(POSTORDER)
      ND(POSTORDER)=TMP_SWAP
      TMP_SWAP=NE(STEP(INODE))
      NE(STEP(INODE))=NE(POSTORDER)
      NE(POSTORDER)=TMP_SWAP
      TMP_SWAP=PROCNODE(STEP(INODE))
      PROCNODE(STEP(INODE))=PROCNODE(POSTORDER)
      PROCNODE(POSTORDER)=TMP_SWAP
      IF(USE_DAD)THEN
         TMP_SWAP=DAD(STEP(INODE))
         DAD(STEP(INODE))=DAD(POSTORDER)
         DAD(POSTORDER)=TMP_SWAP
      ENDIF
      TMP_SWAP=TNSTK(STEP(INODE))
      TNSTK(STEP(INODE))=TNSTK(POSTORDER)
      TNSTK(POSTORDER)=TMP_SWAP
      II=STEP_TO_NODE(POSTORDER)
      TMP_SWAP=STEP(INODE)
      STEP(STEP_TO_NODE(POSTORDER))=TMP_SWAP
      STEP(INODE)=POSTORDER
      STEP_TO_NODE(POSTORDER)=INODE
      STEP_TO_NODE(TMP_SWAP)=II
      IN=II
 101  IN = FILS(IN)
      IF (IN .GT. 0 ) THEN
         STEP(IN)=-STEP(II)
         GOTO 101
      ENDIF
      IN=INODE
 102  IN = FILS(IN)
      IF (IN .GT. 0 ) THEN
         STEP(IN)=-STEP(INODE)
         GOTO 102
      ENDIF
      POSTORDER = POSTORDER + 1
      IF (IFATH.EQ.0) THEN
         NBROOT = NBROOT - 1
         IF (NBROOT.EQ.0) GOTO 116
         GOTO 91
      ENDIF
      TNSTK(STEP(IFATH)) = TNSTK(STEP(IFATH)) - 1
      IF ( TNSTK(STEP(IFATH)) .EQ. 0 ) THEN      
         INODE = IFATH
         GOTO 96
      ELSE
         GOTO 91
      ENDIF
 116  CONTINUE
      DEALLOCATE(STEP_TO_NODE)
      DEALLOCATE(IPOOL,TNSTK)
      RETURN
      END SUBROUTINE MUMPS_SORT_STEP
      SUBROUTINE MUMPS_CHECK_COMM_NODES(COMM_NODES, EXIT_FLAG)
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: COMM_NODES
      LOGICAL, INTENT(OUT) :: EXIT_FLAG
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE), IERR
      CALL MPI_IPROBE( MPI_ANY_SOURCE, TERREUR, COMM_NODES,
     &            EXIT_FLAG, STATUS, IERR)
      RETURN
      END SUBROUTINE MUMPS_CHECK_COMM_NODES
      SUBROUTINE MUMPS_GET_PROC_PER_NODE(K414, MyID, NbProcs, COMM)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER :: K414, MyID, NbProcs, COMM, ALLOCOK
      INTEGER :: ierr,MyNAME_length,MyNAME_length_RCV,i,j
      CHARACTER(len=MPI_MAX_PROCESSOR_NAME) :: MyNAME
      CHARACTER, dimension(:), allocatable :: MyNAME_TAB,MyNAME_TAB_RCV
      logical :: SAME_NAME
      call MPI_GET_PROCESSOR_NAME(MyNAME, MyNAME_length, ierr)
      allocate(MyNAME_TAB(MyNAME_length), STAT=ALLOCOK)
      IF(ALLOCOK.LT.0) THEN
         write(*,*) "Allocation error in MUMPS_GET_PROC_PER_NODE"
         call MUMPS_ABORT()
      ENDIF
      DO i=1, MyNAME_length
         MyNAME_TAB(i) = MyNAME(i:i)
      ENDDO
      K414=0
      do i=0, NbProcs-1
         if(MyID .eq. i) then
            MyNAME_length_RCV  = MyNAME_length
         else
            MyNAME_length_RCV = 0
         endif
         call MPI_BCAST(MyNAME_length_RCV,1,MPI_INTEGER,
     &        i,COMM,ierr)
         allocate(MyNAME_TAB_RCV(MyNAME_length_RCV), STAT=ALLOCOK)
         IF(ALLOCOK.LT.0) THEN
            write(*,*) "Allocation error in MUMPS_GET_PROC_PER_NODE"
            call MUMPS_ABORT()
         ENDIF
         if(MyID .eq. i) then
            MyNAME_TAB_RCV = MyNAME_TAB
         endif
         call MPI_BCAST(MyNAME_TAB_RCV,MyNAME_length_RCV,MPI_CHARACTER,
     &        i,COMM,ierr)
         SAME_NAME=.FALSE.
         IF(MyNAME_length .EQ. MyNAME_length_RCV) THEN
            DO j=1, MyNAME_length
               IF(MyNAME_TAB(j) .NE. MyNAME_TAB_RCV(j)) THEN
                  goto 100
               ENDIF
            ENDDO
            SAME_NAME=.TRUE.
         ENDIF
 100     continue
         IF(SAME_NAME) K414=K414+1
         deallocate(MyNAME_TAB_RCV)
      enddo
      deallocate(MyNAME_TAB)
      END SUBROUTINE MUMPS_GET_PROC_PER_NODE
      SUBROUTINE MUMPS_ICOPY_32TO64 (INTAB, SIZETAB, OUTTAB8)
      IMPLICIT NONE
      INTEGER, intent(in)     ::  SIZETAB
      INTEGER, intent(in)     ::  INTAB(SIZETAB)
      INTEGER(8), intent(out) ::  OUTTAB8(SIZETAB)
      INTEGER :: I
      DO I=1,SIZETAB
       OUTTAB8(I) = int(INTAB(I),8)
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_ICOPY_32TO64
      SUBROUTINE MUMPS_ICOPY_32TO64_64C(INTAB, SIZETAB8, OUTTAB8)
      IMPLICIT NONE
      INTEGER(8), intent(in)  ::  SIZETAB8
      INTEGER, intent(in)     ::  INTAB(SIZETAB8)
      INTEGER(8), intent(out) ::  OUTTAB8(SIZETAB8)
      INTEGER(8) :: I8
      LOGICAL    :: OMP_FLAG
      OMP_FLAG = (SIZETAB8 .GE.500000_8 )
!$OMP PARALLEL DO PRIVATE(I8)
!$OMP&         IF(OMP_FLAG)
      DO I8=1_8, SIZETAB8
        OUTTAB8(I8) = int(INTAB(I8),8)
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE MUMPS_ICOPY_32TO64_64C
      SUBROUTINE MUMPS_ICOPY_32TO64_64C_IP(IN_OUT_TAB48, SIZETAB,
     &    SIZETABX2 ) 
      IMPLICIT NONE
      INTEGER(8), intent(in) :: SIZETAB, SIZETABX2
      INTEGER, intent(inout) :: IN_OUT_TAB48(SIZETABX2)
      CALL MUMPS_ICOPY_32TO64_64C_IP_REC(IN_OUT_TAB48, SIZETAB,
     &    SIZETABX2)
      RETURN
      END SUBROUTINE MUMPS_ICOPY_32TO64_64C_IP
      RECURSIVE SUBROUTINE MUMPS_ICOPY_32TO64_64C_IP_REC(
     &                     IN_OUT_TAB48, SIZETAB, SIZETABX2)
      IMPLICIT NONE
      INTEGER(8), intent(in) :: SIZETAB, SIZETABX2
      INTEGER :: IN_OUT_TAB48(SIZETABX2)
      INTEGER(8) :: IBEG24, IBEG28, SIZE1, SIZE2
      IF (SIZETAB.LE. 1000_8) THEN
        CALL MUMPS_ICOPY_32TO64_64C_IP_C(IN_OUT_TAB48,
     &       SIZETAB)
      ELSE
        SIZE2  = SIZETAB / 2
        SIZE1  = SIZETAB - SIZE2
        IBEG24 = SIZE1+1                 
        IBEG28 = 2*SIZE1+1_8             
        CALL MUMPS_ICOPY_32TO64_64C(IN_OUT_TAB48(IBEG24),
     &                           SIZE2, IN_OUT_TAB48(IBEG28))
        CALL MUMPS_ICOPY_32TO64_64C_IP_REC(IN_OUT_TAB48,
     &                           SIZE1, 2_8*SIZE1)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_ICOPY_32TO64_64C_IP_REC
      SUBROUTINE MUMPS_ICOPY_64TO32(INTAB8, SIZETAB, OUTTAB)
      IMPLICIT NONE
      INTEGER, intent(in) ::  SIZETAB
      INTEGER(8), intent(in) ::  INTAB8(SIZETAB)
      INTEGER, intent(out)   ::  OUTTAB(SIZETAB)
      INTEGER :: I
      DO I=1,SIZETAB
       OUTTAB(I) = int(INTAB8(I))
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_ICOPY_64TO32
      SUBROUTINE MUMPS_ICOPY_64TO32_64C (INTAB8, SIZETAB, OUTTAB)
      IMPLICIT NONE
      INTEGER(8), intent(in)    ::  SIZETAB
      INTEGER(8), intent(in) ::  INTAB8(SIZETAB)
      INTEGER, intent(out)   ::  OUTTAB(SIZETAB)
      INTEGER(8) :: I8
      DO I8=1_8,SIZETAB
       OUTTAB(I8) = int(INTAB8(I8))
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_ICOPY_64TO32_64C
      SUBROUTINE MUMPS_ICOPY_64TO32_64C_IP(IN_OUT_TAB48, SIZETAB,
     &                                     SIZETABX2)
      IMPLICIT NONE
      INTEGER(8), intent(in) :: SIZETAB, SIZETABX2
      INTEGER, intent(inout) :: IN_OUT_TAB48(SIZETABX2)
      CALL MUMPS_ICOPY_64TO32_64C_IP_REC(IN_OUT_TAB48, SIZETAB,
     &                                   SIZETABX2)
      RETURN
      END SUBROUTINE MUMPS_ICOPY_64TO32_64C_IP
      RECURSIVE SUBROUTINE MUMPS_ICOPY_64TO32_64C_IP_REC(
     &                     IN_OUT_TAB48, SIZETAB, SIZETABX2)
      IMPLICIT NONE
      INTEGER(8), intent(in) :: SIZETAB, SIZETABX2
      INTEGER :: IN_OUT_TAB48(SIZETABX2)
      INTEGER(8) :: IBEG24, IBEG28, SIZE1, SIZE2
      IF (SIZETAB.LE. 1000_8) THEN
        CALL MUMPS_ICOPY_64TO32_64C_IP_C(IN_OUT_TAB48,
     &       SIZETAB)
      ELSE
        SIZE2  = SIZETAB / 2
        SIZE1  = SIZETAB - SIZE2
        IBEG24 = SIZE1 + 1
        IBEG28 = SIZE1 + SIZE1 + 1_8
        CALL MUMPS_ICOPY_64TO32_64C_IP_REC(IN_OUT_TAB48,
     &                                     SIZE1, 2_8*SIZE1)
        CALL MUMPS_ICOPY_64TO32_64C(IN_OUT_TAB48(IBEG28),
     &                           SIZE2, IN_OUT_TAB48(IBEG24))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_ICOPY_64TO32_64C_IP_REC
      SUBROUTINE MUMPS_GET_NNZ_INTERNAL( NNZ, NZ, NNZ_i )
      INTEGER   , INTENT(IN)  :: NZ
      INTEGER(8), INTENT(IN)  :: NNZ
      INTEGER(8), INTENT(OUT) :: NNZ_i
      IF (NNZ > 0_8) THEN
        NNZ_i = NNZ
      ELSE
        NNZ_i = int(NZ, 8)
      ENDIF
      END SUBROUTINE MUMPS_GET_NNZ_INTERNAL
      SUBROUTINE MUMPS_NPIV_CRITICAL_PATH(
     &     N, NSTEPS, STEP, FRERE, FILS,
     &     NA, LNA, NE, MAXNPIVTREE )
      IMPLICIT NONE
      INTEGER, intent(in) :: N, NSTEPS, LNA
      INTEGER, intent(in) :: FRERE(NSTEPS), FILS(N), STEP(N)
      INTEGER, intent(in) :: NA(LNA), NE(NSTEPS)
      INTEGER, intent(out) :: MAXNPIVTREE
      INTEGER :: IFATH,INODE,ISON
      INTEGER :: NPIV,ILEAF,NBLEAF
      INTEGER, DIMENSION(:) , ALLOCATABLE :: MAXNPIV
      INTEGER :: I, allocok
      MAXNPIVTREE = -9999
      ALLOCATE ( MAXNPIV(NSTEPS), stat=allocok)
      IF (allocok .gt.0) THEN
         WRITE(*, *) 'Allocation error in MUMPS_NPIV_CRITICAL_PATH' 
     &           ,NSTEPS
         CALL MUMPS_ABORT()
      ENDIF
      MAXNPIV = 0
      NBLEAF = NA(1)
      DO ILEAF = 1, NBLEAF
        INODE = NA(2+ILEAF)
 95     CONTINUE
        NPIV = 0
        ISON = INODE
 100    NPIV = NPIV + 1
        ISON = FILS(ISON)
        IF (ISON .GT. 0 ) GOTO 100
        ISON = -ISON
        MAXNPIV( STEP(INODE) ) = NPIV
        DO I = 1, NE(STEP(INODE))
          MAXNPIV(STEP(INODE)) = max( MAXNPIV(STEP(INODE)),
     &                                NPIV + MAXNPIV(STEP(ISON)) )
          ISON     = FRERE(STEP(ISON))
        ENDDO
        IFATH = INODE
        DO WHILE (IFATH .GT. 0)
          IFATH = FRERE(STEP(IFATH))
        ENDDO
        IFATH = -IFATH
        IF (IFATH.EQ.0) THEN 
          MAXNPIVTREE = max(MAXNPIVTREE, MAXNPIV(STEP(INODE)))
        ELSE 
          IF (FRERE(STEP(INODE)) .LT. 0) THEN
            INODE = IFATH
            GOTO 95
          ENDIF
        ENDIF
      ENDDO
      DEALLOCATE( MAXNPIV )
      RETURN
      END SUBROUTINE MUMPS_NPIV_CRITICAL_PATH
      SUBROUTINE MUMPS_LDLTPANEL_NBTARGET( NPIV, NB_TARGET, KEEP )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPIV
      INTEGER, INTENT(IN) :: KEEP(500)
      INTEGER, INTENT(OUT) :: NB_TARGET
      INTEGER :: NBPANELS, NBCOLMIN, NBPANELSMAX
      IF (NPIV .EQ. 0) THEN
        NB_TARGET = 0
      ELSE
        NBCOLMIN    = KEEP(460)
        NBPANELSMAX = KEEP(459)
        NBPANELS = min( (NPIV+NBCOLMIN-1) / NBCOLMIN, NBPANELSMAX  )
        NB_TARGET = ( NPIV+NBPANELS-1 ) / NBPANELS
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_LDLTPANEL_NBTARGET
      SUBROUTINE MUMPS_LDLTPANEL_STORAGE
     &           ( NPIV, KEEP, IW, NB_ENTRIES )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPIV 
      INTEGER, INTENT(IN) :: KEEP(500), IW(*) 
      INTEGER(8), INTENT(OUT) :: NB_ENTRIES 
      INTEGER :: NB_TARGET, NBCOLS_PANEL, NBROWS_PANEL
      INTEGER :: ICOL_BEG, ICOL_END, NBPANELS
      CALL MUMPS_LDLTPANEL_NBTARGET( NPIV, NB_TARGET, KEEP )
      NB_ENTRIES = 0_8
      NBROWS_PANEL = NPIV
      ICOL_BEG = 1
      NBPANELS = 0
      DO WHILE ( ICOL_BEG .LE. NPIV )
        NBPANELS = NBPANELS + 1
        ICOL_END = min(NB_TARGET * NBPANELS, NPIV)
        IF (IW(1) .NE. 0) THEN
          IF ( IW( ICOL_END ) < 0 ) THEN 
            ICOL_END = ICOL_END + 1
          ENDIF
        ENDIF
        NBCOLS_PANEL = ICOL_END - ICOL_BEG + 1
        NB_ENTRIES = NB_ENTRIES + int(NBCOLS_PANEL,8) *
     &                            int(NBROWS_PANEL,8)
        NBROWS_PANEL = NBROWS_PANEL - NBCOLS_PANEL
        ICOL_BEG = ICOL_END + 1
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_LDLTPANEL_STORAGE
      SUBROUTINE MUMPS_LDLTPANEL_PANELINFOS( NPIV, KEEP, IW,
     &       NB_TARGET, NBPANELS, PANEL_COL, PANEL_POS, PANEL_TABSIZE,
     &       IGNORE_K459 )
      IMPLICIT NONE
      INTEGER,    INTENT(IN)  :: NPIV
      INTEGER,    INTENT(IN)  :: IW( NPIV )
      INTEGER,    INTENT(IN)  :: KEEP(500)
      INTEGER,    INTENT(IN)  :: PANEL_TABSIZE
      INTEGER,    INTENT(OUT) :: NB_TARGET, NBPANELS
      INTEGER,    INTENT(OUT) :: PANEL_COL( PANEL_TABSIZE )
      INTEGER(8), INTENT(OUT) :: PANEL_POS( PANEL_TABSIZE )
      LOGICAL,    INTENT(IN)  :: IGNORE_K459
      INTEGER :: IPANEL, ICOL_END, NBROWS_PANEL, NBCOLS_PANEL
      IF ( IGNORE_K459 ) THEN
        NB_TARGET = NPIV
      ELSE
        CALL MUMPS_LDLTPANEL_NBTARGET( NPIV, NB_TARGET, KEEP )
      ENDIF
      PANEL_POS(1) = 1_8
      PANEL_COL(1) = 1
      NBROWS_PANEL = NPIV
      NBPANELS     = 1
      IF ( KEEP(459) .GT. 1 .AND. KEEP(50) .NE. 0 .AND.
     &   NB_TARGET.NE.NPIV ) THEN
        NBPANELS = ( NPIV + NB_TARGET -1 ) / NB_TARGET
        IF ( PANEL_TABSIZE .LT. NBPANELS + 1 ) THEN
          WRITE(*,*) " Internal error in MUMPS_LDLTPANEL_PANELINFOS",
     &    PANEL_TABSIZE, NBPANELS
          CALL MUMPS_ABORT()
        ENDIF
        DO IPANEL=1, NBPANELS
          ICOL_END = min(IPANEL*NB_TARGET, NPIV)
          IF ( IW(ICOL_END) .LT. 0 ) THEN
            ICOL_END = ICOL_END + 1
          ENDIF
          NBCOLS_PANEL = ICOL_END - PANEL_COL(IPANEL) + 1
          PANEL_POS(IPANEL+1) = PANEL_POS(IPANEL) + 
     &    int(NBROWS_PANEL,8)*int(NBCOLS_PANEL,8)
          PANEL_COL(IPANEL+1) = PANEL_COL(IPANEL) + NBCOLS_PANEL
          NBROWS_PANEL = NBROWS_PANEL - NBCOLS_PANEL
        ENDDO
      ELSE
          PANEL_POS(2)  = int(NPIV,8)*int(NPIV,8)+1_8
          PANEL_COL(2) = NPIV+1
      ENDIF
      END SUBROUTINE MUMPS_LDLTPANEL_PANELINFOS
      SUBROUTINE MUMPS_LDLTPANEL_SIZES
     &           ( NPIV, KEEP, IW, PANEL_SIZES, NBPANELS )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPIV 
      INTEGER, INTENT(IN) :: KEEP(500), IW(NPIV) 
      INTEGER(8), INTENT(OUT) :: PANEL_SIZES( KEEP(459) ) 
      INTEGER, INTENT(OUT) :: NBPANELS
      INTEGER :: NB_TARGET
      INTEGER :: ICOL_BEG, ICOL_END
      NBPANELS = 0
      CALL MUMPS_LDLTPANEL_NBTARGET( NPIV, NB_TARGET, KEEP )
      ICOL_BEG = 1
      NBPANELS = 0
      DO WHILE ( ICOL_BEG .LE. NPIV )
        NBPANELS = NBPANELS + 1
        ICOL_END = min(NB_TARGET * NBPANELS, NPIV)
        IF ( IW( ICOL_END ) < 0 ) THEN 
          ICOL_END = ICOL_END + 1
        ENDIF
        PANEL_SIZES(NBPANELS) = ICOL_END-ICOL_BEG+1
        ICOL_BEG = ICOL_END + 1
      ENDDO
      PANEL_SIZES(NBPANELS+1:KEEP(459))=0
      RETURN
      END SUBROUTINE MUMPS_LDLTPANEL_SIZES
      SUBROUTINE MUMPS_BUILD_ARCH_NODE_COMM
     &           ( COMM, NEWCOMM, NEWSIZE, NEWRANK, COMM0 )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER, INTENT(IN)  :: COMM
      INTEGER, INTENT(OUT) :: NEWCOMM, NEWSIZE, NEWRANK, COMM0
      INTEGER :: SMALLEST_ID_ON_SAME_NODE, IPROC, MYID, IERR, NPROCS
      INTEGER :: TMPNAME_LENGTH, MYNAME_LENGTH
      CHARACTER(len=MPI_MAX_PROCESSOR_NAME) :: MYNAME, TMPNAME
      INTEGER :: COLOR
      SMALLEST_ID_ON_SAME_NODE = -1
      CALL MPI_COMM_RANK( COMM, MYID, IERR )
      CALL MPI_COMM_SIZE( COMM, NPROCS, IERR )
      CALL MPI_GET_PROCESSOR_NAME(MYNAME, MYNAME_LENGTH, IERR )
      DO IPROC = 0, NPROCS - 1
        IF (MYID .EQ. IPROC) THEN
          TMPNAME = MYNAME
          TMPNAME_LENGTH = MYNAME_LENGTH
        ENDIF
        CALL MPI_BCAST( TMPNAME_LENGTH, 1, MPI_INTEGER,
     &                  IPROC, COMM, IERR )
        CALL MPI_BCAST( TMPNAME, TMPNAME_LENGTH, MPI_CHARACTER,
     &                  IPROC, COMM, IERR)
        IF (SMALLEST_ID_ON_SAME_NODE .LT. 0) THEN
         IF ( TMPNAME_LENGTH .EQ. MYNAME_LENGTH ) THEN
          IF ( TMPNAME(1:TMPNAME_LENGTH) .EQ. MYNAME(1:MYNAME_LENGTH) )
     &    THEN
            SMALLEST_ID_ON_SAME_NODE = IPROC
          ENDIF
         ENDIF
        ENDIF
      ENDDO
      CALL MPI_COMM_SPLIT( COMM, SMALLEST_ID_ON_SAME_NODE, 0,
     &                     NEWCOMM, IERR )
      CALL MPI_COMM_RANK( NEWCOMM, NEWRANK, IERR )
      CALL MPI_COMM_SIZE( NEWCOMM, NEWSIZE, IERR )
      IF (NEWRANK .EQ.0) THEN
        COLOR = 0
      ELSE
        COLOR = MPI_UNDEFINED
      ENDIF
      CALL MPI_COMM_SPLIT( COMM, COLOR, 0, COMM0, IERR )
      RETURN
      END SUBROUTINE MUMPS_BUILD_ARCH_NODE_COMM
      SUBROUTINE MUMPS_DESTROY_ARCH_NODE_COMM( ARCH_NODE_COMM,
     &                                         COMM0, RK )
      IMPLICIT NONE
      INTEGER :: ARCH_NODE_COMM, COMM0, RK, IERR
      INCLUDE 'mpif.h'
      CALL MPI_COMM_FREE( ARCH_NODE_COMM, IERR )
      IF ( RK .EQ. 0 ) CALL MPI_COMM_FREE( COMM0, IERR )
      RETURN
      END SUBROUTINE MUMPS_DESTROY_ARCH_NODE_COMM
      SUBROUTINE  MUMPS_DM_FAC_UPD_DYN_MEMCNTS
     &                    ( MEM_COUNT_ALLOCATED, ATOMIC_UPDATES, KEEP8,
     &                      IFLAG, IERROR, K69UPD, K71UPD )
      IMPLICIT NONE
      INTEGER(8), INTENT(IN)        :: MEM_COUNT_ALLOCATED
      INTEGER(8), INTENT(INOUT)     :: KEEP8(150)
      LOGICAL, INTENT(IN)           :: ATOMIC_UPDATES
      INTEGER, INTENT(INOUT)        :: IFLAG, IERROR
      LOGICAL, INTENT(IN) :: K69UPD
      LOGICAL, INTENT(IN) :: K71UPD
      INTEGER(8) :: KEEP8TMPCOPY
      IF (MEM_COUNT_ALLOCATED.GT.0) THEN
           IF (ATOMIC_UPDATES ) THEN
!$OMP        ATOMIC CAPTURE
             KEEP8(73) = KEEP8(73) + MEM_COUNT_ALLOCATED
             KEEP8TMPCOPY = KEEP8(73)
!$OMP        END ATOMIC
!$OMP        ATOMIC UPDATE
             KEEP8(74) = max(KEEP8(74), KEEP8TMPCOPY)
!$OMP        END ATOMIC
           ELSE
             KEEP8(73) = KEEP8(73) + MEM_COUNT_ALLOCATED
             KEEP8TMPCOPY = KEEP8(73)
             KEEP8(74) = max(KEEP8(74), KEEP8(73))
           ENDIF
           IF ( KEEP8TMPCOPY .GT. KEEP8(75) ) THEN
             IFLAG = -19
             CALL MUMPS_SET_IERROR(
     &            (KEEP8TMPCOPY-KEEP8(75)), IERROR)
           ENDIF
           IF ( K69UPD ) THEN
             IF ( ATOMIC_UPDATES ) THEN
!$OMP          ATOMIC CAPTURE
               KEEP8(69) = KEEP8(69) + MEM_COUNT_ALLOCATED
               KEEP8TMPCOPY = KEEP8(69)
!$OMP          END ATOMIC
!$OMP          ATOMIC UPDATE
               KEEP8(68) = max(KEEP8(68), KEEP8TMPCOPY)
!$OMP          END ATOMIC
             ELSE
               KEEP8(69) = KEEP8(69) + MEM_COUNT_ALLOCATED
               KEEP8(68) = max(KEEP8(69), KEEP8(68))
             ENDIF
           ENDIF
           IF ( K71UPD ) THEN
             IF ( ATOMIC_UPDATES ) THEN
!$OMP          ATOMIC CAPTURE
               KEEP8(71) = KEEP8(71) + MEM_COUNT_ALLOCATED
               KEEP8TMPCOPY = KEEP8(71)
!$OMP          END ATOMIC
!$OMP          ATOMIC UPDATE
               KEEP8(70) = max(KEEP8(70), KEEP8TMPCOPY)
!$OMP          END ATOMIC
             ELSE
               KEEP8(71) = KEEP8(71) + MEM_COUNT_ALLOCATED
               KEEP8(70) = max(KEEP8(71), KEEP8(70))
             ENDIF
           ENDIF
      ELSE
           IF (ATOMIC_UPDATES) THEN
!$OMP ATOMIC UPDATE
             KEEP8(73) = KEEP8(73) + MEM_COUNT_ALLOCATED
!$OMP END ATOMIC
             IF ( K69UPD ) THEN
!$OMP ATOMIC UPDATE
               KEEP8(69) = KEEP8(69) + MEM_COUNT_ALLOCATED
!$OMP END ATOMIC
             ENDIF
             IF ( K71UPD ) THEN
!$OMP ATOMIC UPDATE
               KEEP8(71) = KEEP8(71) + MEM_COUNT_ALLOCATED
!$OMP END ATOMIC
             ENDIF
           ELSE
             KEEP8(73) = KEEP8(73) + MEM_COUNT_ALLOCATED
             IF ( K69UPD ) THEN
               KEEP8(69) = KEEP8(69) + MEM_COUNT_ALLOCATED
             ENDIF
             IF ( K71UPD ) THEN
               KEEP8(71) = KEEP8(71) + MEM_COUNT_ALLOCATED
             ENDIF
           ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_DM_FAC_UPD_DYN_MEMCNTS
      SUBROUTINE MUMPS_SET_PARTI_REGULAR(
     &     SLAVEF,
     &     KEEP,KEEP8,
     &     PROCS,
     &     MEM_DISTRIB, NCB, NFRONT, NSLAVES_NODE,
     &     TAB_POS, SLAVES_LIST, SIZE_SLAVES_LIST,MYID,INODE,
     &     TAB_MAXS_ARG,SUP_PROC_ARG,MAX_SURF,NB_ROW_MAX
     &     )
      IMPLICIT NONE
      INTEGER, intent(in) :: KEEP(500),SIZE_SLAVES_LIST
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in) :: SLAVEF, NFRONT, NCB,MYID
      INTEGER, intent(in) :: PROCS(SLAVEF+1)
      INTEGER(8), intent(in) :: TAB_MAXS_ARG(0:SLAVEF-1)
      INTEGER, intent(in) :: SUP_PROC_ARG(2)
      INTEGER, intent(in) :: MEM_DISTRIB(0:SLAVEF-1),INODE
      INTEGER, intent(out):: SLAVES_LIST(SIZE_SLAVES_LIST)
      INTEGER, intent(out):: TAB_POS(SLAVEF+2)
      INTEGER, intent(out):: NSLAVES_NODE,NB_ROW_MAX
      INTEGER(8), intent(out):: MAX_SURF
      LOGICAL :: FORCE_LDLTRegular_NIV2
      INTEGER NSLAVES,ACC
      INTEGER i,J,NELIM,NB_SUP,K50,NB_ROWS(PROCS(SLAVEF+1))
      INTEGER TMP_NROW,X,K
      LOGICAL SUP,MEM_CSTR
      DOUBLE PRECISION MAX_LOAD,TOTAL_LOAD,VAR,TMP,A,B,C,DELTA,
     &     LOAD_CORR
      INTEGER IDWLOAD(SLAVEF)
      INTEGER(8) MEM_CONSTRAINT(2)
      K50=KEEP(50)
      FORCE_LDLTRegular_NIV2 = .FALSE.
      MAX_SURF=0
      NB_ROW_MAX=0
      NELIM=NFRONT-NCB
      NB_SUP=0
      TOTAL_LOAD=0.0D0
      SUP=.FALSE.
      IF(SUP_PROC_ARG(1).NE.
     &     0)THEN
         MEM_CONSTRAINT(1)=TAB_MAXS_ARG(PROCS(1))
         TOTAL_LOAD=TOTAL_LOAD+dble(SUP_PROC_ARG(1))/100.0D0
         NB_SUP=NB_SUP+1
      ENDIF
      IF(SUP_PROC_ARG(2).NE.
     &     0)THEN
         MEM_CONSTRAINT(2)=TAB_MAXS_ARG(PROCS(PROCS(SLAVEF+1)))
         TOTAL_LOAD=TOTAL_LOAD+dble(SUP_PROC_ARG(2))/100.0D0
         NB_SUP=NB_SUP+1
      ENDIF
      TOTAL_LOAD=TOTAL_LOAD+(PROCS(SLAVEF+1)-NB_SUP)
      IF(K50.EQ.0)THEN
         MAX_LOAD=dble( NELIM ) * dble ( NCB ) +
     *        dble(NCB) * dble(NELIM)*dble(2*NFRONT-NELIM-1)
      ELSE
         MAX_LOAD=dble(NELIM) * dble ( NCB ) *
     *        dble(NFRONT+1)
      ENDIF
      TMP=min(MAX_LOAD,MAX_LOAD/TOTAL_LOAD)
      J=1
      DO i=1,PROCS(SLAVEF+1)
         IF((NB_SUP.GT.0).AND.(i.EQ.1))THEN
            CYCLE
         ELSEIF((NB_SUP.EQ.2).AND.(i.EQ.PROCS(SLAVEF+1)))THEN
            CYCLE
         ENDIF
         IDWLOAD(J)=PROCS(i)
         J=J+1
      ENDDO
      DO i=1,NB_SUP
         IF(i.EQ.1)THEN
            IDWLOAD(J)=PROCS(1)
         ELSE
            IDWLOAD(J)=PROCS(PROCS(SLAVEF+1))
         ENDIF
         J=J+1
      ENDDO
      IF ((K50.EQ.0).OR.FORCE_LDLTRegular_NIV2) THEN
         ACC=0
         J=PROCS(SLAVEF+1)-NB_SUP+1
         DO i=1,NB_SUP
            VAR=dble(SUP_PROC_ARG(i))/100.0D0
            TMP_NROW=int(dble(MEM_CONSTRAINT(i))/dble(NFRONT))
            NB_ROWS(J)=int(max((VAR*dble(TMP))/
     &           (dble(NELIM)*dble(2*NFRONT-NELIM)),
     &           dble(1)))
            IF(NB_ROWS(J).GT.TMP_NROW)THEN
               NB_ROWS(J)=TMP_NROW
            ENDIF
            IF(NCB-ACC.LT.NB_ROWS(J)) THEN 
               NB_ROWS(J)=NCB-ACC
               ACC=NCB
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(J)
            J=J+1
         ENDDO
         IF(ACC.EQ.NCB)THEN
            GOTO 777
         ENDIF
         DO i=1,PROCS(SLAVEF+1)-NB_SUP
            VAR=1.0D0
            TMP_NROW=int((dble(TAB_MAXS_ARG(IDWLOAD(i))))/dble(NFRONT))
            NB_ROWS(i)=int((dble(VAR)*dble(TMP))/
     &           (dble(NELIM)*dble(2*NFRONT-NELIM)))
            IF(NB_ROWS(i).GT.TMP_NROW)THEN
               NB_ROWS(i)=TMP_NROW
            ENDIF
            IF(NCB-ACC.LT.NB_ROWS(i)) THEN 
               NB_ROWS(i)=NCB-ACC
               ACC=NCB
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(i)
         ENDDO
         IF(ACC.NE.NCB)THEN
            IF(PROCS(SLAVEF+1).EQ.NB_SUP)THEN
               TMP_NROW=(NCB-ACC)/PROCS(SLAVEF+1)+1
               DO i=1,PROCS(SLAVEF+1)
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  IF(ACC+TMP_NROW.GT.NCB)THEN
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+NCB-ACC
                     ACC=NCB
                     EXIT
                  ENDIF
                  ACC=ACC+TMP_NROW
               ENDDO
            ELSE
               TMP_NROW=(NCB-ACC)/(PROCS(SLAVEF+1)-NB_SUP)+1
               DO i=1,PROCS(SLAVEF+1)-NB_SUP
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  ACC=ACC+TMP_NROW
                  IF(ACC.GT.NCB) THEN 
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+
     &                    (NCB-(ACC-TMP_NROW))
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ELSE
         ACC=0
         i=PROCS(SLAVEF+1)-NB_SUP+1
         X=NCB
         LOAD_CORR=0.0D0
         MEM_CSTR=.FALSE.
         DO J=1,NB_SUP
            VAR=DBLE(SUP_PROC_ARG(J))/DBLE(100)
            A=1.0D0
            B=dble(X+NELIM)
            C=-dble(max(MEM_CONSTRAINT(J),0_8))
            DELTA=((B*B)-(4*A*C))
            TMP_NROW=int((-B+sqrt(DELTA))/(2*A))
            A=dble(-NELIM)
            B=dble(NELIM)*(dble(-NELIM)+dble(2*(X+NELIM)+1))
            C=-(VAR*TMP)
            DELTA=(B*B-(4*A*C))
            NB_ROWS(i)=int((-B+sqrt(DELTA))/(2*A))
            IF(NB_ROWS(i).GT.TMP_NROW)THEN
               NB_ROWS(i)=TMP_NROW
               MEM_CSTR=.TRUE.               
            ENDIF
            IF(ACC+NB_ROWS(i).GT.NCB)THEN
               NB_ROWS(i)=NCB-ACC
               ACC=NCB
               X=0
               EXIT
            ENDIF
            X=X-NB_ROWS(i)
            ACC=ACC+NB_ROWS(i) 
            LOAD_CORR=LOAD_CORR+(dble(NELIM) * dble (NB_ROWS(i)) *
     *           dble(2*(X+NELIM) - NELIM - NB_ROWS(i) + 1))
            i=i+1
         ENDDO
         IF(ACC.EQ.NCB)THEN
            GOTO 777
         ENDIF
            IF((PROCS(SLAVEF+1).NE.NB_SUP).AND.MEM_CSTR)THEN
               TMP=(MAX_LOAD-LOAD_CORR)/(PROCS(SLAVEF+1)-NB_SUP)
            ENDIF
         X=ACC
         ACC=0
         DO i=1,PROCS(SLAVEF+1)-NB_SUP
            IF (KEEP(375) .EQ. 1) THEN 
              VAR=1.0D0
              A=dble(NELIM)
              B=dble(NELIM)*(dble(NELIM)+dble(2*ACC+1))
              C=-(VAR*TMP)
            ELSE    
              A=1.0D0
              B=dble(ACC+NELIM)
              C=-TMP
            ENDIF
            DELTA=((B*B)-(4*A*C))
            NB_ROWS(i)=int((-B+sqrt(DELTA))/(2*A))
            IF(NCB-ACC-X.LT.NB_ROWS(i))THEN
               NB_ROWS(i)=NCB-ACC-X
               ACC=NCB-X
               EXIT
            ENDIF
            ACC=ACC+NB_ROWS(i)
         ENDDO
         ACC=ACC+X
         IF(ACC.NE.NCB)THEN
            IF(PROCS(SLAVEF+1).EQ.NB_SUP)THEN
               TMP_NROW=(NCB-ACC)/PROCS(SLAVEF+1)+1
               DO i=1,PROCS(SLAVEF+1)
                  NB_ROWS(i)=NB_ROWS(i)+TMP_NROW
                  IF(ACC+TMP_NROW.GT.NCB)THEN
                     NB_ROWS(i)=NB_ROWS(i)-TMP_NROW+NCB-ACC
                     ACC=NCB
                     EXIT
                  ENDIF
                  ACC=ACC+TMP_NROW
               ENDDO
            ELSE
               NB_ROWS(PROCS(SLAVEF+1)-NB_SUP)=
     &              NB_ROWS(PROCS(SLAVEF+1)
     &              -NB_SUP)+NCB-ACC
            ENDIF
         ENDIF
      ENDIF
 777  CONTINUE
      NSLAVES=0
      ACC=1
      J=1
      K=1
      DO i=1,PROCS(SLAVEF+1)
         IF(NB_ROWS(i).NE.0)THEN
            SLAVES_LIST(J)=IDWLOAD(i)
            TAB_POS(J)=ACC
            ACC=ACC+NB_ROWS(i)
            NB_ROW_MAX=max(NB_ROW_MAX,NB_ROWS(i))
            IF(K50.EQ.0)THEN
               MAX_SURF=max(int(NB_ROWS(i),8)*int(NCB,8),int(0,8))
            ELSE
               MAX_SURF=max(int(NB_ROWS(i),8)*int(ACC,8),int(0,8))
            ENDIF
            NSLAVES=NSLAVES+1
            J=J+1
         ELSE
            SLAVES_LIST(PROCS(SLAVEF+1)-K+1)=IDWLOAD(i)
            K=K+1
         ENDIF
      ENDDO
      TAB_POS(SLAVEF+2) = NSLAVES
      TAB_POS(NSLAVES+1)= NCB+1
      NSLAVES_NODE=NSLAVES
      END SUBROUTINE MUMPS_SET_PARTI_REGULAR
      SUBROUTINE MUMPS_BIGALLREDUCE(IS_IN_PLACE,
     &           B1, B2, N,
     &           MPI_TYPE, MPI_OP, COMM, IERR_MPI)
      IMPLICIT NONE
      LOGICAL :: IS_IN_PLACE
      INTEGER :: N, MPI_TYPE, MPI_OP, COMM, IERR_MPI
      INTEGER :: B1(*), B2(*) 
      INTEGER, PARAMETER :: CHUNK = 250 000 000
      INCLUDE 'mpif.h'
      INTEGER :: I, THISCHUNK
      INTEGER(8) :: IPOS
      IF ( MPI_TYPE .NE. MPI_INTEGER .AND. MPI_TYPE .NE. MPI_2INTEGER )
     &  THEN
        WRITE(*,*) "Internal error MUMPS_BIGALLREDUCE",MPI_TYPE
      ENDIF
      DO I = 1, N, CHUNK
        THISCHUNK=min(CHUNK, N-I+1)
        IF ( MPI_TYPE.EQ.MPI_INTEGER) THEN
                 IPOS = int(I,8)
        ELSE
                 IPOS =int(I,8)+int(I,8)-1_8
        ENDIF
        IF ( IS_IN_PLACE ) THEN
          CALL MPI_ALLREDUCE( MPI_IN_PLACE, B2(IPOS), THISCHUNK,
     &         MPI_TYPE, MPI_OP, COMM, IERR_MPI )
        ELSE
          CALL MPI_ALLREDUCE( B1(IPOS), B2(IPOS), THISCHUNK,
     &         MPI_TYPE, MPI_OP, COMM, IERR_MPI )
        ENDIF
      ENDDO
      RETURN
      END SUBROUTINE MUMPS_BIGALLREDUCE
      INTEGER FUNCTION MUMPS_NUMROC(ROOT_SIZE, BLOCK_SIZE,
     &           MYPROC_ROW, IPROC_ROW1, NPROCS_ROW)
      IMPLICIT NONE
      INTEGER :: ROOT_SIZE, BLOCK_SIZE,
     &           MYPROC_ROW, IPROC_ROW1, NPROCS_ROW
      INTEGER :: NBLOCKS, EXTRA_FULL_BLOCKS, DIST_IPROC_ROW1, MYNBROWS
      NBLOCKS = ROOT_SIZE / BLOCK_SIZE
      MYNBROWS = (NBLOCKS/NPROCS_ROW)*BLOCK_SIZE
      EXTRA_FULL_BLOCKS = mod(NBLOCKS, NPROCS_ROW)
      IF (MYPROC_ROW .GE. IPROC_ROW1) THEN
        DIST_IPROC_ROW1 = MYPROC_ROW - IPROC_ROW1
      ELSE
        DIST_IPROC_ROW1 = NPROCS_ROW + MYPROC_ROW - IPROC_ROW1
      ENDIF
      IF ( DIST_IPROC_ROW1 .LT. EXTRA_FULL_BLOCKS ) THEN
        MYNBROWS = MYNBROWS + BLOCK_SIZE
      ELSE IF ( DIST_IPROC_ROW1 .EQ. EXTRA_FULL_BLOCKS ) THEN
        MYNBROWS = MYNBROWS + mod(ROOT_SIZE,BLOCK_SIZE)
      ENDIF
      MUMPS_NUMROC = MYNBROWS
      RETURN
      END FUNCTION MUMPS_NUMROC
