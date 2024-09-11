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
#if ! defined(PCPRET)
      SUBROUTINE DMUMPS_FAC_DIST_ARROWHEADS_OMP  (      
     & N, NZ_loc8, 
     &    A_loc, IRN_loc, JCN_loc, 
     &    SIZESCAL, LSCAL, COLSCA, ROWSCA, 
     & DBLARR, LDBLARR, INTARR, LINTARR,
     & PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     & FILS, KEEP,KEEP8, MYID, COMM, NBRECORDS,
     &
     & S, LA, root, PROCNODE_STEPS, NPROCS, SLAVEF, PERM, STEP,
     & ICNTL, INFO, NSEND8, NLOCAL8,
     & ISTEP_TO_INIV2, CANDIDATES
     & )
!$    USE OMP_LIB
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, SIZESCAL
      LOGICAL, INTENT(IN)    :: LSCAL
      INTEGER(8), INTENT(IN) :: NZ_loc8
      INTEGER, INTENT(IN)    :: IRN_LOC(max(1_8,NZ_loc8)), 
     &                          JCN_LOC(max(1_8,NZ_loc8))
      DOUBLE PRECISION, INTENT(IN)    :: A_loc(max(1_8,NZ_loc8))
      DOUBLE PRECISION, INTENT(IN)       :: ROWSCA(SIZESCAL),
     &                          COLSCA(SIZESCAL)
      INTEGER(8), INTENT(IN) :: LDBLARR, LINTARR
      DOUBLE PRECISION, INTENT(OUT)   :: DBLARR( LDBLARR )
      INTEGER, INTENT(OUT)   :: INTARR( LINTARR )
      INTEGER, INTENT(INOUT) :: KEEP(500)   
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER(8), INTENT(INOUT):: KEEP8(150) 
      INTEGER, INTENT(IN)    :: FILS( N )
      INTEGER, INTENT(IN)    :: MYID, COMM, NBRECORDS
      INTEGER(8), INTENT(IN) :: LA
      INTEGER, INTENT(IN)    :: NPROCS, SLAVEF
      INTEGER(8), INTENT(OUT):: NSEND8, NLOCAL8
      INTEGER, INTENT(IN)    :: ISTEP_TO_INIV2(KEEP(71))
      INTEGER, INTENT(IN)    :: CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      DOUBLE PRECISION, INTENT(INOUT) :: S( LA )
      TYPE (DMUMPS_ROOT_STRUC), INTENT(INOUT) :: root
      INTEGER, INTENT(IN)    :: PROCNODE_STEPS(KEEP(28)), 
     &                         PERM( N ), STEP( N )
      INTEGER, INTENT(INOUT) :: INFO( 80 ) 
      INTEGER, INTENT(IN)    :: ICNTL(60)
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, ALLOCATABLE, DIMENSION(:,:)  :: IW4
      INTEGER(8), ALLOCATABLE, DIMENSION(:) :: PTRAW
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFSENDI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: BUFSENDR
      INTEGER, ALLOCATABLE, DIMENSION(:,:)   :: BUFSEND_POSRESERVED
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BUFRECVI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: BUFRECVR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IACT, ISENDREQI, ISENDREQR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEND_ACTIVE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IRECVREQI, IRECVREQR
      INTEGER, ALLOCATABLE, DIMENSION(:):: RECV_BUF_STATUS
         INTEGER, PARAMETER :: Processed_IrecNeeded          = 1
         INTEGER, PARAMETER :: IrecPosted                    = 2
         INTEGER, PARAMETER :: Received_NotProcessed         = 3
         INTEGER, PARAMETER :: Processed_IrecNotneeded       = 4
         INTEGER, PARAMETER :: Treating                      = 5
         INTEGER, PARAMETER :: BeingTreatednotbyme           = 6 
      INTEGER(8) :: ILOC8
      INTEGER :: EndNZloc, NB_END_MSG_2_RECV
      LOGICAL :: MPI_End_Send, End_TreatRecvBuf, MPI_InvolvedinSend, 
     &           MPI_InvolvedinRecv, TH_InvolvedinComm, 
     &           NO_ATOMIC_Wsendbuf, NO_ATOMIC_Warrow, FINISHED,
     &           TH_InvolvedinArrange, TH_InvolvedinTreatRecv
      INTEGER(8) :: PTR_ROOT
      INTEGER    :: LOCAL_M, LOCAL_N, ARROW_ROOT
      LOGICAL :: EARLYT3ROOTINS
      LOGICAL :: I_AM_SLAVE, OneMPI
      INTEGER :: IARR1, IORG, NOMP, NOMP_MAX
      INTEGER :: ISTEP, ISLAVE_MAIN, IMAIN, JMAIN
      INTEGER :: allocok 
      LOGICAL :: OMP_FLAG
      INTEGER(8) :: IS8MAIN
      INTEGER :: TYPE_NODE_P, MASTER_NODE_P, NBJ_P
      INTEGER(8) :: IS8_P
      INTEGER :: LP, MP
      LOGICAL :: LPOK, PROK
      INTEGER(8) :: NB_RANGE_8
      INTEGER :: SHIFT_PID
      INTEGER :: NOMP_SHARED
      LOGICAL    :: NOTHINGTOARRANGE_P
      INTEGER    :: IOMP, NB_RANGE_P, EndNZloc_P
      LOGICAL    :: ThWorking
      INTEGER(8) :: ILOC8_P
      INTEGER, PARAMETER :: MPI_MASTER = 0  
      NSEND8  = 0_8
      NLOCAL8 = 0_8
      LP = ICNTL(1)
      MP = ICNTL(2)
      LPOK    = ((LP.GT.0).AND.(ICNTL(4).GE.1))
      PROK    = ((MP.GT.0).AND.(ICNTL(4).GE.2))
      NB_RANGE_8        = int(max(NBRECORDS/10, 1), 8)
      IF (KEEP(46).EQ.0) THEN
       SHIFT_PID = 1
      ELSE
       SHIFT_PID = 0
      ENDIF
      I_AM_SLAVE = (MYID.NE.0.OR.KEEP(46).EQ.1)
      OneMPI     = NPROCS.EQ.1 
      IF ( OneMPI.OR. 
     &     (KEEP(54).EQ.0.AND.(MYID.NE.MPI_MASTER)) 
     &   ) THEN
        MPI_InvolvedinSend = .FALSE.  
        MPI_End_Send       = .TRUE.
      ELSE
        MPI_InvolvedinSend = .TRUE.   
        MPI_End_Send       = .FALSE.
      ENDIF
      IF (.NOT.OneMPI) THEN
        ALLOCATE( 
     &    BUFSENDI(NBRECORDS * 2 + 1, 2, NPROCS),  
     &    BUFSENDR(NBRECORDS, 2, NPROCS),
     &    IACT(NPROCS), SEND_ACTIVE(NPROCS), 
     &    ISENDREQI(NPROCS), ISENDREQR(NPROCS),
     &    BUFSEND_POSRESERVED(2, NPROCS),
     &    stat=allocok)
        IF ( allocok .GT. 0 ) THEN
          IF ( LPOK ) THEN
            WRITE(LP,*)
     &     '** Error allocating SEND buffers for matrix distribution'
          END IF
          INFO(1) = -13
          INFO(2) =  ( NBRECORDS * 2 + 1 ) * NPROCS * 2 +
     &             NBRECORDS * NPROCS * 2 + 
     &             NPROCS*6
          GOTO 20
        END IF
        DO ISLAVE_MAIN=1, NPROCS
          IACT (ISLAVE_MAIN)         = 1
          ISENDREQI(ISLAVE_MAIN) = MPI_REQUEST_NULL
          ISENDREQR(ISLAVE_MAIN) = MPI_REQUEST_NULL
          BUFSENDI(1, 1, ISLAVE_MAIN)       = 0
          BUFSEND_POSRESERVED(1,ISLAVE_MAIN)= 0
          BUFSENDI(1, 2, ISLAVE_MAIN)       = NBRECORDS  
          BUFSEND_POSRESERVED(2,ISLAVE_MAIN)= NBRECORDS
          SEND_ACTIVE(ISLAVE_MAIN)   =  .FALSE.
        ENDDO
      ENDIF
      IF (OneMPI.OR. 
     &    (KEEP(54).EQ.0.AND.(MYID.EQ.MPI_MASTER))
     &   ) THEN
        NB_END_MSG_2_RECV  = 0
        MPI_InvolvedinRecv = .FALSE. 
        End_TreatRecvBuf   = .TRUE.
      ELSE IF (KEEP(54).EQ.0.AND.MYID.NE.MPI_MASTER) THEN
        NB_END_MSG_2_RECV  =  1
        MPI_InvolvedinRecv = .TRUE.  
        End_TreatRecvBuf   = .FALSE.
      ELSE
        NB_END_MSG_2_RECV  = NPROCS-1
        MPI_InvolvedinRecv = .TRUE.  
        End_TreatRecvBuf   = .FALSE.
      ENDIF
      IF (.NOT.OneMPI) THEN
        ALLOCATE( 
     &  BUFRECVI(NBRECORDS * 2 + 1, NPROCS),
     &  BUFRECVR(NBRECORDS, NPROCS),
     &  IRECVREQI(NPROCS), IRECVREQR(NPROCS), 
     &  RECV_BUF_STATUS(NPROCS),
     &            stat=allocok)
        IF ( allocok .GT. 0 ) THEN
           IF ( LPOK ) THEN
             WRITE(LP,*)
     &      '** Error allocating RECV buffers for matrix distribution'
           END IF
           INFO(1) = -13
           INFO(2) =  ( NBRECORDS * 2 + 1 ) * NPROCS +
     &              NBRECORDS * NPROCS  +      
     &              NPROCS*3
           GOTO 20
        ENDIF
        BUFRECVI(1, 1:NPROCS) = 0
        IRECVREQI(1:NPROCS)   = MPI_REQUEST_NULL
        IRECVREQR(1:NPROCS)   = MPI_REQUEST_NULL
        RECV_BUF_STATUS (1:NPROCS)= Processed_IrecNeeded
        RECV_BUF_STATUS (MYID+1)  = Processed_IrecNotneeded
        IF (KEEP(54).EQ.0) THEN
         DO ISLAVE_MAIN=1, NPROCS
           RECV_BUF_STATUS (ISLAVE_MAIN)= Processed_IrecNotneeded
         ENDDO
         IF (MYID.NE.MPI_MASTER) THEN
           RECV_BUF_STATUS(MPI_MASTER+1)  = Processed_IrecNeeded
         ENDIF
        ENDIF
      ENDIF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        CALL MUMPS_SET_IERROR( int(N,8) * 2_8, INFO(2) )
        GOTO 20
      END IF
      ALLOCATE( PTRAW( N ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        CALL MUMPS_SET_IERROR( int(N,8) * KEEP(10), INFO(2) )
        GOTO 20
      END IF
 20   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 100
      ARROW_ROOT = 0
      EARLYT3ROOTINS = KEEP(200) .EQ.0
     &        .OR.  ( KEEP(200) .LT. 0 .AND. KEEP(400) .EQ. 0 )
      IF (I_AM_SLAVE) THEN
       DO JMAIN = 1, N
        ISTEP=STEP(JMAIN)
        IF (ISTEP .GT. 0) THEN
          IARR1 = PTRDEBARR( ISTEP )
          IF ( IARR1 .GT. 0 ) THEN
            IMAIN = JMAIN
            IORG = 0
            DO WHILE ( IMAIN .GT. 0 )
              IORG = IORG + 1
              IW4(IMAIN, 1)  = NINCOLARR( IARR1 + IORG - 1 )
              IW4(IMAIN, 2)  = NINROWARR( IARR1 + IORG - 1 ) +
     &                     NINCOLARR( IARR1 + IORG - 1 )
              IS8MAIN        = PTR8ARR( IARR1 + IORG - 1 )
              PTRAW( IMAIN ) = IS8MAIN
              INTARR( IS8MAIN ) = IMAIN
              DBLARR( IS8MAIN ) = ZERO
              IMAIN = FILS(IMAIN)
            ENDDO
          ENDIF
        ENDIF
       ENDDO
       IF ( KEEP(38) .NE. 0 .AND. EARLYT3ROOTINS ) THEN
         CALL DMUMPS_GET_ROOT_INFO(root,LOCAL_M, LOCAL_N, PTR_ROOT, LA)
         CALL DMUMPS_SET_ROOT_TO_ZERO(root, KEEP, S, LA)
       ELSE
         LOCAL_M = -19999; LOCAL_N = -29999; PTR_ROOT = -99999_8
       END IF
      ELSE 
         LOCAL_M = -19999; LOCAL_N = -29999; PTR_ROOT = -99999_8
      ENDIF
      NOMP=1
!$    NOMP=omp_get_max_threads()
      NOMP_MAX = NOMP
      IF (NOMP_MAX.GT.2.AND.KEEP(399).EQ.2) THEN
       IF (.NOT.OneMPI) THEN
        NOMP_MAX = 2
       ENDIF
      ENDIF
      IF (NOMP_MAX.GT.3.AND.KEEP(399).EQ.3) THEN
       IF (.NOT.OneMPI) THEN
        NOMP_MAX = 3
       ENDIF
      ENDIF
      ILOC8      = 1_8  
      OMP_FLAG   = ((NOMP .GE.2).AND.(KEEP(399).NE.99))
      FINISHED = .FALSE.
      NOMP_SHARED = 1
!$OMP PARALLEL 
!$OMP&  NUM_THREADS(NOMP_MAX)
!$OMP&    PRIVATE(
!$OMP&         IOMP, ThWorking, ILOC8_P, NB_RANGE_P,
!$OMP&         NOTHINGTOARRANGE_P, EndNZloc_P, TH_InvolvedinComm,
!$OMP&         TH_InvolvedinTreatRecv,  TH_InvolvedinArrange )
!$OMP& REDUCTION(+:NSEND8, NLOCAL8, ARROW_ROOT) 
!$OMP& IF (OMP_FLAG)
      IOMP        = 0
!$    IOMP=omp_get_thread_num()
!$OMP SINGLE
!$      NOMP_SHARED= omp_get_num_threads()  
        IF (OneMPI) THEN
         EndNZloc = NOMP_SHARED
        ELSE
         EndNZloc = max(1,NOMP_SHARED -1)
         IF ( KEEP(399).EQ.2.OR.KEEP(399).EQ.3 ) THEN
          EndNZloc = min(EndNZloc,1)
         ENDIF
        ENDIF
        IF (NZ_loc8.EQ.0_8) EndNZloc = 0  
        IF (.NOT.MPI_InvolvedinSend.AND.(EndNZloc.EQ.0)) EndNZloc=-1
        NO_ATOMIC_Wsendbuf = ( NOMP_SHARED.EQ.1 )
        NO_ATOMIC_Warrow   = ( NOMP_SHARED.EQ.1 )
        IF (NPROCS.GT.1) THEN
          NO_ATOMIC_Warrow = (NOMP_SHARED.LE.2)
          IF ( KEEP(399).EQ.2 .OR. KEEP(399).EQ.3) THEN
           NO_ATOMIC_Wsendbuf =  .TRUE.
           IF (.NOT.MPI_InvolvedinSend) NO_ATOMIC_Warrow=.TRUE.
           IF (.NOT.MPI_InvolvedinRecv) NO_ATOMIC_Warrow=.TRUE.
          ENDIF
        ENDIF
!$OMP END SINGLE
      ThWorking          =  OneMPI.OR.
     &                     (NOMP_SHARED.EQ.1) .OR. (IOMP.NE.0)
      TH_InvolvedinTreatRecv = (MPI_InvolvedinRecv.AND.ThWorking)
      IF ( TH_InvolvedinTreatRecv.AND.
     &     (NOMP_SHARED.EQ.3).AND.(KEEP(399).EQ.3) ) THEN
          IF (IOMP.NE.2) TH_InvolvedinTreatRecv = .FALSE.
      ENDIF
      TH_InvolvedinArrange = ThWorking
      IF (.NOT.OneMPI.AND.ThWorking) THEN
        IF (KEEP(399).EQ.2.OR.KEEP(399).EQ.3) 
     &    THEN
          IF ((NOMP_SHARED.NE.1).AND.(IOMP.NE.1)) 
     &       TH_InvolvedinArrange = .FALSE.
        ENDIF
      ENDIF
      TH_InvolvedinComm  = ((.NOT.OneMPI).AND.(IOMP.EQ.0))
      NOTHINGTOARRANGE_P = (NZ_loc8.EQ.0_8)
      ILOC8_P            = 0_8
       DO WHILE ( .NOT.FINISHED )
          IF (TH_InvolvedinComm) THEN
            CALL DMUMPS_ARROW_TRY_PROGRESS_COMM ( IOMP, MYID,
     &       NPROCS, NBRECORDS, MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &       COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &       IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &       IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &       RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &       EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &        )
          ENDIF
          IF (.NOT.ThWorking) THEN
            CALL MUMPS_USLEEP(20)
            GOTO  50
          ENDIF
          IF (TH_InvolvedinTreatRecv) THEN
            CALL DMUMPS_ARROW_TRY_TREAT_RECV_BUF ( IOMP,
     &             BUFRECVI, BUFRECVR, RECV_BUF_STATUS, 
     &             NPROCS, NBRECORDS, N, IW4,
     &             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT, S, LA,
     &             NB_END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF,
     &             PTRAW, 
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR, 
     &             NOMP_SHARED, EARLYT3ROOTINS, NO_ATOMIC_Warrow )
            IF (NOMP_SHARED.EQ.1) THEN
              CALL DMUMPS_ARROW_TRY_PROGRESS_COMM ( IOMP, MYID,
     &        NPROCS, NBRECORDS, 
     &        MPI_InvolvedinSend, MPI_InvolvedinRecv,
     &        COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &        IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &        IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &        RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &        EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &        )
            ENDIF
          ENDIF
          IF (.NOT. NOTHINGTOARRANGE_P.AND.TH_InvolvedinArrange) THEN
!$OMP ATOMIC CAPTURE
            ILOC8_P = ILOC8
            ILOC8   = ILOC8 + NB_RANGE_8
!$OMP END ATOMIC 
            IF (ILOC8_P.LE.NZ_loc8) THEN
              NB_RANGE_P = int(min(NB_RANGE_8, NZ_loc8-ILOC8_P+1))
              CALL DMUMPS_FAC_ARROW_ARRANGE ( MYID, IOMP, N, SHIFT_PID,
     &   SLAVEF, LSCAL, NSEND8, NLOCAL8, ILOC8_P, NB_RANGE_P, 
     &   NZ_loc8, IRN_LOC, JCN_LOC, A_loc, SIZESCAL, ROWSCA, COLSCA,
     &   ISTEP_TO_INIV2, CANDIDATES, KEEP, MPI_MASTER, NO_ATOMIC_Warrow,
     &   NO_ATOMIC_Wsendbuf, TH_InvolvedinTreatRecv,
     &   NPROCS, NBRECORDS, MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &   COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, IACT,
     &   SEND_ACTIVE, BUFSEND_POSRESERVED, IRECVREQI, IRECVREQR,
     &   BUFRECVI, BUFRECVR, RECV_BUF_STATUS, EndNZloc, MPI_End_Send,
     &   End_TreatRecvBuf,
     &   root, PTR_ROOT, ARROW_ROOT, EARLYT3ROOTINS, LOCAL_M, LOCAL_N,
     &   S, LA, NB_END_MSG_2_RECV, PROCNODE_STEPS, IW4, PTRAW, PTR8ARR,
     &   NINCOLARR, NINROWARR, PTRDEBARR, PERM, STEP, INTARR, LINTARR,
     &   DBLARR, LDBLARR, NOMP_SHARED )
            ENDIF
            IF (ILOC8_P+NB_RANGE_8.GT.NZ_loc8) THEN
              IF (.NOT. NOTHINGTOARRANGE_P) THEN
               NOTHINGTOARRANGE_P=.TRUE. 
!$OMP          ATOMIC CAPTURE
               EndNZloc   = EndNZloc-1
               EndNZloc_P = EndNZloc
!$OMP          END ATOMIC
               IF (MPI_End_Send.AND.EndNZloc_P.EQ.0) THEN
!$OMP           ATOMIC WRITE
                EndNZloc=-1
!$OMP           END ATOMIC
               ENDIF
              ENDIF
            ENDIF
          ENDIF
  50      CONTINUE  
!$OMP     MASTER
!$OMP     ATOMIC WRITE
          FINISHED = ( (EndNZloc.EQ.-1) 
     &            .AND.(MPI_End_Send.OR.(.not.MPI_InvolvedinSend))
     &            .AND. End_TreatRecvBuf
     &               )
!$OMP     END ATOMIC
!$OMP     END MASTER
       ENDDO
!$OMP END PARALLEL
!$OMP PARALLEL DO SCHEDULE(DYNAMIC,1)
!$OMP&         PRIVATE(ISTEP, TYPE_NODE_P, MASTER_NODE_P, NBJ_P, 
!$OMP&                 IARR1, IS8_P  )
!$OMP&        IF (OMP_FLAG)
      DO ISTEP=1, KEEP(28)
        CALL MUMPS_TYPEANDPROCNODE( TYPE_NODE_P, MASTER_NODE_P,
     &    PROCNODE_STEPS(ISTEP), KEEP(199) )
        MASTER_NODE_P = MASTER_NODE_P + SHIFT_PID
        IF ( MASTER_NODE_P.NE.MYID.OR.
     &      (  (TYPE_NODE_P.NE.1) .AND. (TYPE_NODE_P.NE.2) )
     &     ) CYCLE
        IARR1  = PTRDEBARR( ISTEP )
        NBJ_P  = NINCOLARR( IARR1)
        IF (NBJ_P.LE.0) CYCLE
        IS8_P  = PTR8ARR( IARR1) + 1_8
        CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &              INTARR( IS8_P ),
     &              DBLARR( IS8_P ),
     &              NBJ_P, 1, NBJ_P )
      ENDDO
!$OMP END PARALLEL DO
 100  CONTINUE
      IF (ALLOCATED(IW4))         DEALLOCATE( IW4 )
      IF (ALLOCATED(PTRAW))       DEALLOCATE( PTRAW )
      IF (ALLOCATED(BUFSENDI))    DEALLOCATE( BUFSENDI )
      IF (ALLOCATED(BUFSENDR))    DEALLOCATE( BUFSENDR )
      IF (ALLOCATED(BUFRECVI))    DEALLOCATE( BUFRECVI )
      IF (ALLOCATED(BUFRECVR))    DEALLOCATE( BUFRECVR )
      IF (ALLOCATED(IACT))        DEALLOCATE( IACT )
      IF (ALLOCATED(ISENDREQI))   DEALLOCATE( ISENDREQI )
      IF (ALLOCATED(ISENDREQR))   DEALLOCATE( ISENDREQR )
      IF (ALLOCATED(IRECVREQI))   DEALLOCATE( IRECVREQI )
      IF (ALLOCATED(IRECVREQR))   DEALLOCATE( IRECVREQR )
      IF (ALLOCATED(SEND_ACTIVE)) DEALLOCATE( SEND_ACTIVE )
      IF (ALLOCATED(BUFSEND_POSRESERVED))   
     &                      DEALLOCATE( BUFSEND_POSRESERVED )
      IF (ALLOCATED(RECV_BUF_STATUS)) DEALLOCATE( RECV_BUF_STATUS )
      RETURN
      END SUBROUTINE DMUMPS_FAC_DIST_ARROWHEADS_OMP 
      SUBROUTINE DMUMPS_ARROW_TRY_PROGRESS_COMM ( IOMP, MYID,
     &     NPROCS, NBRECORDS, MPI_InvolvedinSend, MPI_InvolvedinRecv,
     &     COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &     IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &     IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &     RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &     EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &        )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IOMP, MYID, NPROCS, NBRECORDS, COMM
      LOGICAL, INTENT(IN) :: MPI_InvolvedinSend, MPI_InvolvedinRecv
      INTEGER, INTENT(IN) :: NB_END_MSG_2_RECV
      INTEGER, INTENT(INOUT) :: EndNZloc
      LOGICAL, INTENT(INOUT) :: MPI_End_Send, End_TreatRecvBuf
      INTEGER, INTENT(INOUT) :: 
     &     ISENDREQI(NPROCS), ISENDREQR(NPROCS), 
     &     BUFSENDI(NBRECORDS * 2 + 1, 2, NPROCS), 
     &     IACT(NPROCS), BUFSEND_POSRESERVED(2,NPROCS)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFSENDR(NBRECORDS, 2, NPROCS)
      LOGICAL, INTENT(INOUT) :: SEND_ACTIVE(NPROCS)
      INTEGER, INTENT(INOUT) :: 
     &     IRECVREQI(NPROCS), IRECVREQR(NPROCS), 
     &     BUFRECVI(NBRECORDS * 2 + 1, NPROCS),
     &     RECV_BUF_STATUS(NPROCS)
      DOUBLE PRECISION, INTENT(INOUT) ::  BUFRECVR(NBRECORDS, NPROCS)
         INTEGER, PARAMETER :: Processed_IrecNeeded          = 1
         INTEGER, PARAMETER :: IrecPosted                    = 2
         INTEGER, PARAMETER :: Received_NotProcessed         = 3
         INTEGER, PARAMETER :: Processed_IrecNotneeded       = 4
         INTEGER, PARAMETER :: Treating                      = 5
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER :: EndNZloc_copy, ISLAVE, NBREC, STATE, 
     &           NB_END_MSG_2_RECV_COPY, ISLAVE_RECV 
      INTEGER :: IERR, IACT_P, NEXT_IACT
      INTEGER :: TAILLE_SEND_I, TAILLE_SEND_R
      LOGICAL :: FLAG, FLAGRECV, ALL_LAST_MESS_SENT
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      IF (MPI_InvolvedinSend.and.(.NOT.MPI_End_Send)) THEN
        DO ISLAVE = 1, NPROCS
          IF (ISLAVE-1.EQ.MYID) CYCLE
          IF (SEND_ACTIVE( ISLAVE )) THEN
              CALL MPI_TEST( ISENDREQR( ISLAVE ), FLAG, STATUS, IERR )
              IF (FLAG) THEN
                CALL MPI_WAIT( ISENDREQI( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
              ENDIF
          ENDIF
        ENDDO
!$OMP   ATOMIC READ
        EndNZloc_copy = EndNZloc
!$OMP   END ATOMIC 
        ALL_LAST_MESS_SENT = (EndNZloc_copy.EQ.0)
        IF (EndNZloc_copy.NE.-1) THEN
          DO ISLAVE=1, NPROCS
            IF (ISLAVE-1.EQ.MYID) CYCLE
            IF (EndNZloc_copy .EQ. 0) THEN
             NBREC = 
     &         min(BUFSEND_POSRESERVED(IACT(ISLAVE),ISLAVE),NBRECORDS)
             IF (NBREC.EQ.-99) CYCLE
             BUFSENDI(1,IACT(ISLAVE),ISLAVE) = - NBREC
            ELSE
!$OMP         ATOMIC READ
              NBREC = BUFSENDI(1,IACT(ISLAVE),ISLAVE)
!$OMP         END ATOMIC
            ENDIF
            IF ((EndNZloc_copy.EQ.0).OR.(NBREC.EQ.NBRECORDS)) THEN
              IF (.NOT.SEND_ACTIVE(ISLAVE)) THEN
                TAILLE_SEND_I = NBREC * 2 + 1
                TAILLE_SEND_R = NBREC
                IACT_P = IACT(ISLAVE)
                CALL MPI_ISEND( BUFSENDI(1, IACT_P, ISLAVE ),
     &              TAILLE_SEND_I,
     &              MPI_INTEGER, ISLAVE - 1, 
     &              ARR_INT, COMM,
     &              ISENDREQI( ISLAVE ), IERR )
                CALL MPI_ISEND( BUFSENDR(1, IACT_P, ISLAVE ),
     &              TAILLE_SEND_R,
     &              MPI_DOUBLE_PRECISION, ISLAVE - 1, 
     &              ARR_REAL, COMM,
     &              ISENDREQR( ISLAVE ), IERR )
                SEND_ACTIVE( ISLAVE ) = .TRUE.
                IF (EndNZloc_copy.NE.0) THEN
                  NEXT_IACT = 3-IACT_P
!$OMP             ATOMIC WRITE
                  BUFSEND_POSRESERVED(NEXT_IACT,ISLAVE) = 0
!$OMP             END ATOMIC
!$OMP             ATOMIC WRITE
                  BUFSENDI(1,NEXT_IACT,ISLAVE) = 0
!$OMP             END ATOMIC
!$OMP             ATOMIC WRITE
                  IACT( ISLAVE ) = NEXT_IACT
!$OMP             END ATOMIC
                ELSE
                  BUFSEND_POSRESERVED(IACT_P,ISLAVE) = -99
                ENDIF
              ELSE
                ALL_LAST_MESS_SENT=.FALSE.
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        IF (EndNZloc_copy.EQ.0.AND.ALL_LAST_MESS_SENT) THEN 
            EndNZloc      = -1  
            EndNZloc_copy = -1  
        ENDIF
        IF (.NOT.MPI_End_Send) THEN
          IF ( (EndNZloc_copy.EQ.-1) )  THEN
            MPI_End_Send = .TRUE.
            DO ISLAVE = 1, NPROCS
               IF (ISLAVE-1.EQ.MYID) CYCLE
               IF (SEND_ACTIVE( ISLAVE )) THEN 
                 MPI_End_Send=.FALSE.
                 EXIT
               ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
      IF (MPI_InvolvedinRecv.AND.(.NOT.End_TreatRecvBuf)) THEN
        CALL MPI_TESTANY(NPROCS, IRECVREQR, ISLAVE_RECV, 
     &                   FLAGRECV, STATUS,IERR)
        IF (FLAGRECV.AND.(ISLAVE_RECV.NE.MPI_UNDEFINED)) 
     &            THEN
            CALL MPI_WAIT(IRECVREQI(ISLAVE_RECV),STATUS,IERR)
!$OMP       ATOMIC WRITE
            RECV_BUF_STATUS(ISLAVE_RECV)=Received_NotProcessed
!$OMP       END ATOMIC
        ENDIF
        DO  ISLAVE = 1, NPROCS
          IF (ISLAVE - 1 .EQ. MYID) CYCLE
!$OMP     ATOMIC READ
          STATE = RECV_BUF_STATUS(ISLAVE)
!$OMP     END ATOMIC
          IF (STATE.EQ.Processed_IrecNeeded) THEN
            CALL MPI_IRECV ( BUFRECVI(1,ISLAVE), NBRECORDS * 2 + 1,  
     &              MPI_INTEGER, ISLAVE-1, ARR_INT, COMM, 
     &              IRECVREQI(ISLAVE), IERR)
            CALL MPI_IRECV ( BUFRECVR(1,ISLAVE), NBRECORDS,
     &              MPI_DOUBLE_PRECISION, ISLAVE-1, 
     &              ARR_REAL, COMM,
     &              IRECVREQR(ISLAVE), IERR)
!$OMP       ATOMIC WRITE
            RECV_BUF_STATUS(ISLAVE) = IrecPosted
!$OMP       END ATOMIC
          ENDIF
        ENDDO
!$OMP   ATOMIC READ
        NB_END_MSG_2_RECV_COPY = NB_END_MSG_2_RECV
!$OMP   END ATOMIC
        IF (NB_END_MSG_2_RECV_COPY.EQ.0) THEN
          End_TreatRecvBuf = .TRUE.
          DO ISLAVE = 1, NPROCS
            IF (ISLAVE - 1 .EQ. MYID) CYCLE
            IF (RECV_BUF_STATUS(ISLAVE).NE.Processed_IrecNotneeded) THEN
               End_TreatRecvBuf = .FALSE.
               EXIT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_ARROW_TRY_PROGRESS_COMM 
      SUBROUTINE DMUMPS_ARROW_TRY_TREAT_RECV_BUF (IOMP,
     &             BUFRECVI, BUFRECVR, RECV_BUF_STATUS, 
     &             NPROCS, NBRECORDS, N, IW4,
     &             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT, S, LA,
     &             NB_END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF,
     &             PTRAW, 
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR, 
     &             NOMP_SHARED, EARLYT3ROOTINS, NO_ATOMIC_Warrow )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER, INTENT(IN):: IOMP, NPROCS, NBRECORDS, N, MYID, SLAVEF,
     &                      NOMP_SHARED
      LOGICAL, INTENT(IN):: EARLYT3ROOTINS
      INTEGER, INTENT(IN):: BUFRECVI( NBRECORDS * 2 + 1, NPROCS )
      DOUBLE PRECISION, INTENT(IN):: BUFRECVR( NBRECORDS, NPROCS )
      INTEGER, INTENT(INOUT) :: RECV_BUF_STATUS(NPROCS)
      INTEGER, INTENT(INOUT):: IW4( N, 2 )
      INTEGER, INTENT(IN):: KEEP(500)
      INTEGER, INTENT(INOUT):: NB_END_MSG_2_RECV
      INTEGER(8), INTENT(IN) :: PTRAW( N )
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER, INTENT(IN) :: PERM( N ), STEP( N )
      INTEGER, INTENT(IN) :: PROCNODE_STEPS( KEEP(28) )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER, INTENT(INOUT) ::  INTARR( LINTARR )
      INTEGER, INTENT(IN)::  LOCAL_M, LOCAL_N
      INTEGER(8), INTENT(IN) :: PTR_ROOT, LA
      DOUBLE PRECISION, INTENT(INOUT) :: S( LA ), DBLARR( LDBLARR )
      LOGICAL, INTENT(IN)    :: NO_ATOMIC_Warrow
         INTEGER, PARAMETER :: Processed_IrecNeeded          = 1
         INTEGER, PARAMETER :: IrecPosted                    = 2
         INTEGER, PARAMETER :: Received_NotProcessed         = 3
         INTEGER, PARAMETER :: Processed_IrecNotneeded       = 4
         INTEGER, PARAMETER :: Treating                      = 5
         INTEGER, PARAMETER :: BeingTreatednotbyme           = 6
      INTEGER MUMPS_TYPENODE
      EXTERNAL MUMPS_TYPENODE
      INTEGER STATE, ISLAVE
      DO ISLAVE =1, NPROCS
         IF (MYID.EQ.ISLAVE-1) CYCLE
!$OMP    ATOMIC READ
         STATE = RECV_BUF_STATUS(ISLAVE)
!$OMP    END ATOMIC
         IF (STATE.EQ.Treating) CYCLE
         IF (STATE.EQ.Received_NotProcessed) THEN
          IF (NOMP_SHARED.EQ.1) THEN
            RECV_BUF_STATUS(ISLAVE) = Treating
            STATE                   = Treating
          ELSE IF (KEEP(399).LE.3) THEN
!$OMP       ATOMIC WRITE
            RECV_BUF_STATUS(ISLAVE) = Treating
!$OMP       END ATOMIC
            STATE                   = Treating
          ELSE
!$OMP      CRITICAL(ARROW_RECV_BUF_STATUS)
!$OMP        ATOMIC READ
             STATE = RECV_BUF_STATUS(ISLAVE)
!$OMP        END ATOMIC
             IF (STATE.EQ.Received_NotProcessed) THEN
!$OMP          ATOMIC WRITE
               RECV_BUF_STATUS(ISLAVE) = Treating
!$OMP          END ATOMIC 
               STATE                   = Treating
             ELSE
               STATE                   = BeingTreatednotbyme
             ENDIF
!$OMP      END CRITICAL(ARROW_ RECV_BUF_STATUS)
          ENDIF
         ENDIF
        IF (STATE.NE.Treating) CYCLE
        IF (NO_ATOMIC_Warrow) THEN
          CALL DMUMPS_ARROW_TREAT_RECV_BUF_1TH()
        ELSE
          CALL DMUMPS_ARROW_TREAT_RECV_BUF()
        ENDIF
      ENDDO 
      RETURN
      CONTAINS
      SUBROUTINE  DMUMPS_ARROW_TREAT_RECV_BUF()
      INTEGER :: IREC, NB_REC, TYPE_NODE
      INTEGER :: IPOSROOT, JPOSROOT, ILOCROOT, JLOCROOT
      INTEGER(8) :: IS8
      INTEGER    :: IARR, JARR, IW4_CAPTURED
      INTEGER    :: NB_END_MSG_2_RECV_COPY
      DOUBLE PRECISION    :: VAL
      LOGICAL    :: LAST_MESSAGE
      LAST_MESSAGE = .FALSE.
      NB_REC = BUFRECVI( 1, ISLAVE )
      TYPE_NODE = -998
      IF ( NB_REC .LE. 0 ) THEN
        LAST_MESSAGE = .TRUE.
!$OMP ATOMIC CAPTURE
        NB_END_MSG_2_RECV_COPY = NB_END_MSG_2_RECV
        NB_END_MSG_2_RECV      = NB_END_MSG_2_RECV - 1   
!$OMP END ATOMIC
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFRECVI( IREC * 2, ISLAVE )
        JARR = BUFRECVI( IREC * 2 + 1, ISLAVE )
        VAL  = BUFRECVR( IREC, ISLAVE )
        IF (EARLYT3ROOTINS) THEN
          TYPE_NODE = MUMPS_TYPENODE( 
     &              PROCNODE_STEPS(abs(STEP(abs( IARR )))),
     &              KEEP(199) )
        ENDIF
        IF ( TYPE_NODE .eq. 3 .AND. EARLYT3ROOTINS ) THEN
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L( IARR )
            JPOSROOT = root%RG2L( JARR )
          ELSE
            IPOSROOT = root%RG2L( JARR )
            JPOSROOT = root%RG2L( -IARR )
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
!$OMP ATOMIC UPDATE
            S( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT-1,8)) =  S( PTR_ROOT
     &        + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT - 1,8) )
     &      + VAL
!$OMP END ATOMIC
          ELSE
!$OMP ATOMIC UPDATE
            root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8) )
     &      = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8))
     &      + VAL
!$OMP END ATOMIC
          ENDIF
        ELSE IF (IARR.GE.0) THEN
          IF (IARR.EQ.JARR) THEN
            IS8 = PTRAW(IARR)
!$OMP ATOMIC UPDATE
            DBLARR(IS8) = DBLARR(IS8) + VAL
!$OMP END ATOMIC
          ELSE
!$OMP ATOMIC CAPTURE
            IW4_CAPTURED= IW4(IARR,2)
            IW4(IARR,2) = IW4(IARR,2) - 1
!$OMP END ATOMIC
            IS8         = PTRAW(IARR)+int(IW4_CAPTURED,8)
            INTARR(IS8) = JARR
            DBLARR(IS8) = VAL
          ENDIF
        ELSE
          IARR         = -IARR
!$OMP ATOMIC CAPTURE
          IW4_CAPTURED   = IW4(IARR,1)
          IW4(IARR,1)    = IW4(IARR,1) - 1
!$OMP END ATOMIC
          IS8          = PTRAW(IARR) + int(IW4_CAPTURED,8)
          INTARR(IS8)  = JARR
          DBLARR(IS8)  = VAL
        ENDIF
      ENDDO
 100  CONTINUE
      IF (LAST_MESSAGE) THEN
!$OMP ATOMIC WRITE
         RECV_BUF_STATUS(ISLAVE) = Processed_IrecNotneeded
!$OMP END ATOMIC
      ELSE
!$OMP ATOMIC WRITE
         RECV_BUF_STATUS(ISLAVE) = Processed_IrecNeeded
!$OMP END ATOMIC
      ENDIF
      RETURN
      END SUBROUTINE  DMUMPS_ARROW_TREAT_RECV_BUF
      SUBROUTINE  DMUMPS_ARROW_TREAT_RECV_BUF_1TH()
      INTEGER :: IREC, NB_REC, TYPE_NODE
      INTEGER :: IPOSROOT, JPOSROOT, ILOCROOT, JLOCROOT
      INTEGER(8) :: IS8
      INTEGER    :: IARR, JARR
      INTEGER    :: NB_END_MSG_2_RECV_COPY
      DOUBLE PRECISION    :: VAL
      LOGICAL    :: LAST_MESSAGE
      LAST_MESSAGE = .FALSE.
      NB_REC = BUFRECVI( 1, ISLAVE )
      IF ( NB_REC .LE. 0 ) THEN
        LAST_MESSAGE = .TRUE.
!$OMP ATOMIC CAPTURE
        NB_END_MSG_2_RECV_COPY = NB_END_MSG_2_RECV
        NB_END_MSG_2_RECV      = NB_END_MSG_2_RECV - 1   
!$OMP END ATOMIC
        NB_REC = - NB_REC
      END IF 
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFRECVI( IREC * 2, ISLAVE )
        JARR = BUFRECVI( IREC * 2 + 1, ISLAVE )
        VAL  = BUFRECVR( IREC, ISLAVE )
        IF (EARLYT3ROOTINS) THEN
          TYPE_NODE = MUMPS_TYPENODE( 
     &              PROCNODE_STEPS(abs(STEP(abs( IARR )))),
     &              KEEP(199) )
        ENDIF
        IF ( TYPE_NODE .eq. 3 .AND. EARLYT3ROOTINS ) THEN
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L( IARR )
            JPOSROOT = root%RG2L( JARR )
          ELSE
            IPOSROOT = root%RG2L( JARR )
            JPOSROOT = root%RG2L( -IARR )
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
            S( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT-1,8)) =  S( PTR_ROOT
     &        + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT - 1,8) )
     &      + VAL
          ELSE
            root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8) )
     &      = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                      * int(root%SCHUR_LLD,8)
     &                      + int(ILOCROOT,8))
     &      + VAL
          ENDIF
        ELSE IF (IARR.GE.0) THEN
          IF (IARR.EQ.JARR) THEN
            IS8 = PTRAW(IARR)
            DBLARR(IS8) = DBLARR(IS8) + VAL
          ELSE
            IS8         = PTRAW(IARR)+IW4(IARR,2)
            IW4(IARR,2) = IW4(IARR,2) - 1
            INTARR(IS8) = JARR
            DBLARR(IS8) = VAL
          ENDIF
        ELSE
          IARR         = -IARR
          IS8          = PTRAW(IARR)+IW4(IARR,1)
          IW4(IARR,1)  = IW4(IARR,1) - 1
          INTARR(IS8)  = JARR
          DBLARR(IS8)  = VAL
        ENDIF
      ENDDO
 100  CONTINUE
      IF (LAST_MESSAGE) THEN
!$OMP ATOMIC WRITE
         RECV_BUF_STATUS(ISLAVE) = Processed_IrecNotneeded
!$OMP END ATOMIC
      ELSE
!$OMP ATOMIC WRITE
         RECV_BUF_STATUS(ISLAVE) = Processed_IrecNeeded
!$OMP END ATOMIC
      ENDIF
      RETURN
      END SUBROUTINE  DMUMPS_ARROW_TREAT_RECV_BUF_1TH
      END SUBROUTINE DMUMPS_ARROW_TRY_TREAT_RECV_BUF 
      SUBROUTINE DMUMPS_FAC_ARROW_ARRANGE ( 
     &  MYID, IOMP, N, SHIFT_PID, SLAVEF, LSCAL, NSEND8, NLOCAL8,
     &  ILOC8_P, NB_RANGE_P, NZ_loc8, IRN_LOC, JCN_LOC, A_loc, SIZESCAL,
     &  ROWSCA, COLSCA, ISTEP_TO_INIV2, CANDIDATES, KEEP, MPI_MASTER, 
     &  NO_ATOMIC_Warrow, NO_ATOMIC_Wsendbuf, TH_InvolvedinTreatRecv,
     &  NPROCS, NBRECORDS, MPI_InvolvedinSend, MPI_InvolvedinRecv,
     &  COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, IACT,
     &  SEND_ACTIVE, BUFSEND_POSRESERVED, IRECVREQI, IRECVREQR,
     &  BUFRECVI, BUFRECVR, RECV_BUF_STATUS, EndNZloc, MPI_End_Send,
     &  End_TreatRecvBuf,
     &  root, PTR_ROOT, ARROW_ROOT, EARLYT3ROOTINS, LOCAL_M, LOCAL_N,
     &  S, LA, NB_END_MSG_2_RECV, PROCNODE_STEPS, IW4, PTRAW, 
     &  PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &  PERM, STEP, INTARR, LINTARR, DBLARR, LDBLARR, NOMP_SHARED )
!$    USE OMP_LIB
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: MYID, IOMP, N, SHIFT_PID, 
     &                            SLAVEF, SIZESCAL, MPI_MASTER
      INTEGER, INTENT(IN)      :: NB_RANGE_P, KEEP(500)
      INTEGER(8), INTENT(IN)   :: NZ_loc8, ILOC8_P
      INTEGER(8), INTENT(INOUT):: NSEND8, NLOCAL8
      INTEGER, INTENT(IN)      :: IRN_LOC(max(1_8,NZ_loc8)), 
     &                            JCN_LOC(max(1_8,NZ_loc8))
      INTEGER, INTENT(IN):: ISTEP_TO_INIV2(KEEP(71))
      INTEGER, INTENT(IN):: CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      DOUBLE PRECISION, INTENT(IN):: A_loc(max(1_8,NZ_loc8))
      DOUBLE PRECISION, INTENT(IN)   :: ROWSCA(SIZESCAL),
     &                      COLSCA(SIZESCAL)
      LOGICAL, INTENT(IN):: NO_ATOMIC_Warrow, NO_ATOMIC_Wsendbuf, 
     &                      TH_InvolvedinTreatRecv
      INTEGER, INTENT(IN) :: NPROCS, NBRECORDS, COMM, NOMP_SHARED
      LOGICAL, INTENT(IN) :: MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &                       LSCAL
      INTEGER, INTENT(INOUT) :: EndNZloc
      LOGICAL, INTENT(INOUT) :: MPI_End_Send, End_TreatRecvBuf
      INTEGER, INTENT(INOUT) :: ISENDREQI(NPROCS), ISENDREQR(NPROCS), 
     &     BUFSENDI(NBRECORDS * 2 + 1, 2, NPROCS), 
     &     IACT(NPROCS), BUFSEND_POSRESERVED(2,NPROCS)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFSENDR(NBRECORDS, 2, NPROCS)
      LOGICAL, INTENT(INOUT) :: SEND_ACTIVE(NPROCS)
      INTEGER, INTENT(INOUT) :: IRECVREQI(NPROCS), IRECVREQR(NPROCS), 
     &     BUFRECVI(NBRECORDS * 2 + 1, NPROCS),
     &     RECV_BUF_STATUS(NPROCS)
      DOUBLE PRECISION, INTENT(INOUT) :: BUFRECVR(NBRECORDS, NPROCS)
         INTEGER, PARAMETER :: Processed_IrecNeeded          = 1
         INTEGER, PARAMETER :: IrecPosted                    = 2
         INTEGER, PARAMETER :: Received_NotProcessed         = 3
         INTEGER, PARAMETER :: Processed_IrecNotneeded       = 4
         INTEGER, PARAMETER :: Treating                      = 5
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER, INTENT(IN)      :: LOCAL_M, LOCAL_N
      LOGICAL, INTENT(IN)      :: EARLYT3ROOTINS
      INTEGER, INTENT(INOUT)   :: ARROW_ROOT
      INTEGER, INTENT(INOUT):: IW4( N, 2 )
      INTEGER, INTENT(INOUT):: NB_END_MSG_2_RECV
      INTEGER(8), INTENT(IN) :: PTRAW( N )
      INTEGER(8), INTENT(IN) :: PTR8ARR(KEEP(193))
      INTEGER, INTENT(IN) :: NINCOLARR(KEEP(194))
      INTEGER, INTENT(IN) :: NINROWARR(KEEP(195))
      INTEGER, INTENT(IN) :: PTRDEBARR(KEEP(196))
      INTEGER, INTENT(IN) :: PERM( N ), STEP( N )
      INTEGER, INTENT(IN) :: PROCNODE_STEPS( KEEP(28) )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER, INTENT(INOUT) ::  INTARR( LINTARR )
      INTEGER(8), INTENT(IN) :: PTR_ROOT, LA
      DOUBLE PRECISION, INTENT(INOUT) ::  S( LA ), DBLARR( LDBLARR )
      INTEGER MUMPS_TYPESPLIT
      EXTERNAL MUMPS_TYPESPLIT
      INTEGER     :: DEST, JSEND, ISEND , DEST_SAVE
      INTEGER     :: I, INIV2, NCAND, T4MASTER
      INTEGER     :: IOLD, JOLD, IARR, TYPESPLIT
      INTEGER(8)  :: IS8, IZ8, LAST8
      LOGICAL     :: T4_MASTER_CONCERNED
      INTEGER     :: MASTER_NODE, TYPE_NODE, ISTEP_P
      INTEGER     :: IPOSROOT, JPOSROOT, ILOCROOT, JLOCROOT
      INTEGER     :: IROW_GRID, JCOL_GRID, IW4_CAPTURED
      LOGICAL     :: LOCAL_ASSEMBLY, LOCAL
      DOUBLE PRECISION     :: VAL
      INTEGER     :: ISTEP_T3_1PROC
      LAST8 = ILOC8_P + int(NB_RANGE_P-1,8)
      LOCAL_ASSEMBLY  =  (NPROCS.EQ.1)
      IF (NPROCS.EQ.1 .AND. KEEP(38).EQ.0) THEN
        TYPE_NODE = 1
        ISTEP_T3_1PROC  = -9999
      ELSE IF (NPROCS.EQ.1 .AND. KEEP(38).NE.0) THEN
        ISTEP_T3_1PROC = STEP(KEEP(38))
      ELSE
        ISTEP_T3_1PROC = -99999
      ENDIF
      DO IZ8=ILOC8_P, LAST8
          IOLD = IRN_loc(IZ8)
          JOLD = JCN_loc(IZ8)
          IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) THEN
            CYCLE
          ENDIF
          IF (IOLD.EQ.JOLD) THEN
            ISEND = IOLD
            JSEND = IOLD
            IARR  = IOLD
          ELSE IF (PERM(IOLD).LT.PERM(JOLD)) THEN
            IARR  = IOLD
            IF ( KEEP(50) .NE. 0 ) THEN
              ISEND = -IOLD
            ELSE
              ISEND = IOLD
            ENDIF
            JSEND = JOLD
          ELSE
            IARR  = JOLD
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
          VAL = A_loc(IZ8)
          IF (LSCAL) THEN
            VAL = VAL * ROWSCA(IOLD)*COLSCA(JOLD)
          ENDIF
          IF (NPROCS.GT.1 .OR. KEEP(38).NE.0) THEN
            ISTEP_P = abs(STEP(IARR))
            IF ( NPROCS.EQ.1 .AND. ISTEP_P.NE.ISTEP_T3_1PROC ) THEN
              TYPE_NODE=1
            ELSE IF (NPROCS.EQ.1) THEN
              TYPE_NODE=3
              IF (EARLYT3ROOTINS) THEN
                IF ( ISEND < 0 ) THEN
                  IPOSROOT = root%RG2L(JSEND)
                  JPOSROOT = root%RG2L(IARR )
                ELSE
                  IPOSROOT = root%RG2L(IARR )
                  JPOSROOT = root%RG2L(JSEND)
                END IF
              END IF
            ELSE
            ISTEP_P = abs(STEP(IARR))
            CALL MUMPS_TYPEANDPROCNODE( TYPE_NODE, MASTER_NODE,
     &      PROCNODE_STEPS(ISTEP_P), KEEP(199) )
            MASTER_NODE = MASTER_NODE + SHIFT_PID
            T4_MASTER_CONCERNED = .FALSE.
            T4MASTER            = -9999
            IF ( TYPE_NODE .eq. 1 ) THEN
              DEST = MASTER_NODE
            ELSE IF ( TYPE_NODE .eq. 2 ) THEN
              IF ( ISEND .LT. 0 ) THEN
                DEST = -1
              ELSE
                DEST = MASTER_NODE
              END IF
              INIV2         = ISTEP_TO_INIV2(ISTEP_P)
              IF ( KEEP(79) .GT. 0) THEN
                TYPESPLIT  = MUMPS_TYPESPLIT( PROCNODE_STEPS(ISTEP_P),
     &                                        KEEP(199) )
                IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
                  T4_MASTER_CONCERNED = .TRUE.
                  T4MASTER=
     &                   CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
     &                 + SHIFT_PID
                ENDIF
              ENDIF
            ELSE 
              ARROW_ROOT = ARROW_ROOT + 1
              IF (EARLYT3ROOTINS) THEN
                IF ( ISEND < 0 ) THEN
                  IPOSROOT = root%RG2L(JSEND)
                  JPOSROOT = root%RG2L(IARR )
                ELSE
                  IPOSROOT = root%RG2L(IARR )
                  JPOSROOT = root%RG2L(JSEND)
                END IF
                IROW_GRID = mod( ( IPOSROOT-1 )/root%MBLOCK, root%NPROW)
                JCOL_GRID = mod( ( JPOSROOT-1 )/root%NBLOCK, root%NPCOL)
                DEST = IROW_GRID * root%NPCOL + JCOL_GRID + SHIFT_PID
              ELSE
                DEST = -2
              ENDIF
            ENDIF
            IF (DEST .eq. -1) THEN
              NLOCAL8 = NLOCAL8 + 1_8
              NSEND8  = NSEND8 + int(SLAVEF -1,8)
            ELSE IF (DEST .EQ. -2) THEN
              NLOCAL8 = NLOCAL8 + 1_8
              NSEND8  = NSEND8 + int(SLAVEF -1,8)
            ELSE
              IF (DEST .eq.MYID ) THEN
                NLOCAL8 = NLOCAL8 + 1_8
              ELSE
                NSEND8 = NSEND8 + 1_8
              ENDIF
            ENDIF
            LOCAL_ASSEMBLY = .FALSE.
            IF ( DEST.EQ.-1) THEN
              INIV2 = ISTEP_TO_INIV2(ISTEP_P)
              NCAND = CANDIDATES(SLAVEF+1,INIV2)
              IF (KEEP(79) .GT. 0) THEN
                DO I=1, SLAVEF
                  DEST=CANDIDATES(I,INIV2) + SHIFT_PID
                  IF (DEST.LT.0) EXIT 
                  LOCAL = (DEST.EQ.MYID)
                  IF (LOCAL) LOCAL_ASSEMBLY = .TRUE.
                  IF (LOCAL) CYCLE
                  IF (I.EQ.NCAND+1) CYCLE
                  CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                ENDDO
              ELSE
                DO I=1, NCAND
                  DEST=CANDIDATES(I,INIV2) + SHIFT_PID
                  LOCAL = (DEST.EQ.MYID)
                  IF (LOCAL) LOCAL_ASSEMBLY = .TRUE.
                  IF (LOCAL) CYCLE
                  CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                ENDDO
              ENDIF
              IF ( LOCAL_ASSEMBLY ) THEN
                DEST_SAVE = DEST
                DEST      = MASTER_NODE
                IF (DEST.NE.MYID)
     &           CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                IF (T4_MASTER_CONCERNED) THEN
                  DEST   = T4MASTER
                  IF (DEST.NE.MYID)
     &             CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                ENDIF
                DEST = DEST_SAVE
              ELSE
                DEST=MASTER_NODE
                LOCAL_ASSEMBLY = (DEST.EQ.MYID)
                IF (.NOT.LOCAL_ASSEMBLY) 
     &           CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                IF (T4_MASTER_CONCERNED) THEN
                  DEST = T4MASTER
                  LOCAL_ASSEMBLY = (DEST.EQ.MYID)
                  IF (.NOT.LOCAL_ASSEMBLY) 
     &             CALL DMUMPS_DIST_FILL_SEND_BUFFER()
                ENDIF
              ENDIF
            ELSE IF (DEST .GE. 0) THEN
              LOCAL_ASSEMBLY = (DEST.EQ.MYID)
              IF (.NOT.LOCAL_ASSEMBLY) 
     &             CALL DMUMPS_DIST_FILL_SEND_BUFFER()
              IF (T4_MASTER_CONCERNED) THEN
                DEST = T4MASTER
                LOCAL_ASSEMBLY = (DEST.EQ.MYID)
                IF (.NOT.LOCAL_ASSEMBLY) 
     &             CALL DMUMPS_DIST_FILL_SEND_BUFFER()
              ENDIF
            ELSE IF (DEST .EQ. -2) THEN
              DO I = 0, SLAVEF-1
                DEST=I + SHIFT_PID
                IF (DEST.NE.MYID) 
     &             CALL DMUMPS_DIST_FILL_SEND_BUFFER()
              ENDDO
              IF (SHIFT_PID.EQ.1.AND.MYID.EQ.MPI_MASTER) THEN
                LOCAL_ASSEMBLY=.FALSE.
              ELSE
                LOCAL_ASSEMBLY=.TRUE. 
              ENDIF
            ENDIF
          ENDIF
          ENDIF
          IF (LOCAL_ASSEMBLY) THEN
            IF ( TYPE_NODE .eq. 3 .AND. EARLYT3ROOTINS ) THEN
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                     ( root%MBLOCK * root%NPROW ) )
     &                   + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                     ( root%NBLOCK * root%NPCOL ) )
     &                   + mod( JPOSROOT - 1, root%NBLOCK ) + 1
              IF (KEEP(60)==0) THEN
!$OMP ATOMIC UPDATE
                S( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &            + int(ILOCROOT-1,8)) =  S( PTR_ROOT
     &            + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &            + int(ILOCROOT - 1,8) )
     &          + VAL
!$OMP END ATOMIC
              ELSE
!$OMP ATOMIC UPDATE
                root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                          * int(root%SCHUR_LLD,8)
     &                          + int(ILOCROOT,8) )
     &          = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                          * int(root%SCHUR_LLD,8)
     &                          + int(ILOCROOT,8))
     &          + VAL
!$OMP END ATOMIC
              ENDIF
            ELSE
      IF (NO_ATOMIC_Warrow) THEN
             IF (ISEND.EQ.JSEND) THEN
               IS8         = PTRAW(ISEND)
               DBLARR(IS8) = DBLARR(IS8) + VAL
             ELSE IF (ISEND.GE.0) THEN 
               IS8         = PTRAW(IARR)+IW4(IARR,2)
               IW4(IARR,2) = IW4(IARR,2) - 1
               INTARR(IS8) = JSEND
               DBLARR(IS8) = VAL
             ELSE 
               IS8         = PTRAW(IARR)+IW4(IARR,1)
               IW4(IARR,1) = IW4(IARR,1) - 1
               INTARR(IS8) = JSEND
               DBLARR(IS8) = VAL
             ENDIF
      ELSE
             IF (ISEND.EQ.JSEND) THEN
               IS8         = PTRAW(ISEND)
!$OMP ATOMIC UPDATE
               DBLARR(IS8) = DBLARR(IS8) + VAL
!$OMP END ATOMIC
             ELSE IF (ISEND.GE.0) THEN 
!$OMP ATOMIC CAPTURE
               IW4_CAPTURED = IW4(IARR,2)
               IW4(IARR,2)  = IW4(IARR,2) - 1
!$OMP END ATOMIC
               IS8         = PTRAW(IARR) + int(IW4_CAPTURED,8)
               INTARR(IS8) = JSEND
               DBLARR(IS8) = VAL
             ELSE 
!$OMP ATOMIC CAPTURE
               IW4_CAPTURED = IW4(IARR,1)
               IW4(IARR,1)  = IW4(IARR,1) - 1
!$OMP END ATOMIC
               IS8         = PTRAW(IARR) + int(IW4_CAPTURED,8)
               INTARR(IS8) = JSEND
               DBLARR(IS8) = VAL
             ENDIF
      ENDIF
            ENDIF
          ENDIF
      ENDDO
      RETURN
      CONTAINS
      SUBROUTINE DMUMPS_DIST_FILL_SEND_BUFFER( )
      INTEGER IREQ, IACT_P, ISLAVE
      ISLAVE = DEST+1
 100  CONTINUE
!$OMP ATOMIC READ
      IACT_P = IACT(ISLAVE)
!$OMP END ATOMIC
      IF (NO_ATOMIC_Wsendbuf) THEN
        BUFSEND_POSRESERVED(IACT_P,ISLAVE) = 
     &         BUFSEND_POSRESERVED(IACT_P, ISLAVE) + 1
        IREQ = BUFSEND_POSRESERVED(IACT_P, ISLAVE)
       IF (IREQ > NBRECORDS ) THEN
         IF (IREQ > huge(NBRECORDS)-1000 - NOMP_SHARED-2) THEN
          BUFSEND_POSRESERVED(IACT_P, ISLAVE) =  min(NBRECORDS+1, 
     &                     BUFSEND_POSRESERVED(IACT_P, ISLAVE) )
         ENDIF
           IF (NOMP_SHARED.EQ.1)
     &      CALL DMUMPS_ARROW_TRY_PROGRESS_COMM (IOMP,
     &       MYID, NPROCS, NBRECORDS, 
     &       MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &       COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &       IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &       IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &       RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &       EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &       )
          IF (TH_InvolvedinTreatRecv) 
     &    CALL DMUMPS_ARROW_TRY_TREAT_RECV_BUF (IOMP,
     &             BUFRECVI, BUFRECVR, RECV_BUF_STATUS, 
     &             NPROCS, NBRECORDS, N, IW4,
     &             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT, S, LA,
     &             NB_END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF,
     &             PTRAW, 
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR, 
     &             NOMP_SHARED, EARLYT3ROOTINS, NO_ATOMIC_Warrow)
           IF (NOMP_SHARED.EQ.1) THEN
             CALL DMUMPS_ARROW_TRY_PROGRESS_COMM (IOMP,
     &        MYID, NPROCS, NBRECORDS, 
     &        MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &        COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &        IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &        IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &        RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &        EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &        )
           ELSE
              IF (.NOT.TH_InvolvedinTreatRecv) THEN
               CALL MUMPS_USLEEP(200)
              ELSE
               CALL MUMPS_USLEEP(20)
              ENDIF
           ENDIF
         GOTO 100
       ENDIF
       BUFSENDI(IREQ*2,IACT_P,ISLAVE)   = ISEND
       BUFSENDI(IREQ*2+1,IACT_P,ISLAVE) = JSEND
       BUFSENDR(IREQ,IACT_P,ISLAVE )    = VAL
       IF (IREQ.EQ.NBRECORDS) THEN
!$OMP    ATOMIC WRITE
            BUFSENDI(1,IACT_P,ISLAVE) = NBRECORDS
!$OMP    END ATOMIC 
       ENDIF
      ELSE
!$OMP  ATOMIC CAPTURE
       BUFSEND_POSRESERVED(IACT_P,ISLAVE) = 
     &         BUFSEND_POSRESERVED(IACT_P, ISLAVE) + 1
       IREQ = BUFSEND_POSRESERVED(IACT_P, ISLAVE)
!$OMP  END ATOMIC
       IF (IREQ > huge(NBRECORDS)-NOMP_SHARED-2) THEN
!$OMP   ATOMIC UPDATE
        BUFSEND_POSRESERVED(IACT_P, ISLAVE) =  min(NBRECORDS+1, 
     &                     BUFSEND_POSRESERVED(IACT_P, ISLAVE) )
!$OMP   END ATOMIC
       ENDIF
       IF (IREQ > NBRECORDS ) THEN
          IF (NOMP_SHARED.EQ.1) THEN
            CALL DMUMPS_ARROW_TRY_PROGRESS_COMM (IOMP,
     &       MYID, NPROCS, NBRECORDS, 
     &       MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &       COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &       IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &       IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &       RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &       EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &       )
          ENDIF
          IF (TH_InvolvedinTreatRecv) 
     &    CALL DMUMPS_ARROW_TRY_TREAT_RECV_BUF (IOMP,
     &             BUFRECVI, BUFRECVR, RECV_BUF_STATUS, 
     &             NPROCS, NBRECORDS, N, IW4,
     &             KEEP, LOCAL_M, LOCAL_N, root, PTR_ROOT, S, LA,
     &             NB_END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF,
     &             PTRAW, 
     &             PTR8ARR, NINCOLARR, NINROWARR, PTRDEBARR,
     &             PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR, 
     &             NOMP_SHARED, EARLYT3ROOTINS, NO_ATOMIC_Warrow)
          IF (NOMP_SHARED.EQ.1) THEN
             CALL DMUMPS_ARROW_TRY_PROGRESS_COMM (IOMP,
     &        MYID, NPROCS, NBRECORDS, 
     &        MPI_InvolvedinSend, MPI_InvolvedinRecv, 
     &        COMM, ISENDREQI, ISENDREQR, BUFSENDI, BUFSENDR, 
     &        IACT, SEND_ACTIVE, BUFSEND_POSRESERVED, 
     &        IRECVREQI, IRECVREQR, BUFRECVI, BUFRECVR,
     &        RECV_BUF_STATUS, NB_END_MSG_2_RECV, 
     &        EndNZloc, MPI_End_Send, End_TreatRecvBuf
     &        )
          ELSE 
            IF (.NOT.TH_InvolvedinTreatRecv) THEN
             CALL MUMPS_USLEEP(200)
            ELSE
             CALL MUMPS_USLEEP(20)
            ENDIF
          ENDIF
         GOTO 100
         ENDIF
         BUFSENDI(IREQ*2,IACT_P,ISLAVE)   = ISEND
         BUFSENDI(IREQ*2+1,IACT_P,ISLAVE) = JSEND
         BUFSENDR(IREQ,IACT_P,ISLAVE )    = VAL
!$OMP    ATOMIC UPDATE
         BUFSENDI(1,IACT_P,ISLAVE) = BUFSENDI(1,IACT_P,ISLAVE) + 1
!$OMP    END ATOMIC 
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_DIST_FILL_SEND_BUFFER
      END SUBROUTINE DMUMPS_FAC_ARROW_ARRANGE
#endif
