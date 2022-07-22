C
C  This file is part of MUMPS 5.2.1, released
C  on Fri Jun 14 14:46:05 UTC 2019
C
C
C  Copyright 1991-2019 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license:
C  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
C
      SUBROUTINE DMUMPS_BUILD_MAPPING
     & ( N, MAPPING, NNZ, IRN, JCN, PROCNODE, STEP,
     &   SLAVEF, PERM, FILS,
     &   RG2L, KEEP,KEEP8, MBLOCK, NBLOCK, NPROW, NPCOL )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N, SLAVEF, MBLOCK, NBLOCK, NPROW, NPCOL
      iNTEGER(8) :: NNZ
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER IRN( NNZ ), JCN( NNZ ) 
      INTEGER MAPPING( NNZ ), STEP( N )
      INTEGER PROCNODE( KEEP(28) ), PERM( N ), FILS( N ), RG2L( N )
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE
      INTEGER K4, IOLD, JOLD, INEW, JNEW, ISEND, JSEND, IARR, INODE
      INTEGER(8) :: K8
      INTEGER TYPE_NODE, DEST
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID
      INODE = KEEP(38)
      K4 = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = K4
        INODE = FILS( INODE )
        K4 = K4 + 1
      END DO
      DO K8 = 1_8, NNZ
        IOLD = IRN( K8 )
        JOLD = JCN( K8 )
        IF ( IOLD .GT. N .OR. IOLD .LT. 1 .OR.
     &       JOLD .GT. N .OR. JOLD .LT. 1 ) THEN
           MAPPING( K8 ) = -1
           CYCLE
        END IF
        IF ( IOLD .eq. JOLD ) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM( IOLD )
          JNEW = PERM( JOLD )
          IF ( INEW .LT. JNEW ) THEN
            ISEND = IOLD
            IF ( KEEP(50) .ne. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          END IF
        END IF
        IARR = abs( ISEND )
        TYPE_NODE = MUMPS_TYPENODE( PROCNODE(abs(STEP(IARR))),
     &                              SLAVEF )
        IF ( TYPE_NODE .eq. 1 .or. TYPE_NODE .eq. 2 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = MUMPS_PROCNODE( PROCNODE(abs(STEP(IARR))),
     &                             SLAVEF ) + 1
          ELSE
            DEST = MUMPS_PROCNODE( PROCNODE(abs(STEP(IARR))),
     &                             SLAVEF )
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L( JSEND )
            JPOSROOT = RG2L( IARR  )
          ELSE
            IPOSROOT = RG2L( IARR  )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/MBLOCK, NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/NBLOCK, NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * NPCOL + JCOL_GRID
          END IF
        END IF
        MAPPING( K8 ) = DEST
      END DO
      RETURN
      END SUBROUTINE DMUMPS_BUILD_MAPPING
      SUBROUTINE DMUMPS_REDISTRIBUTION(
     & N, NZ_loc8, id,
     & DBLARR, LDBLARR, INTARR, LINTARR,
     & PTRAIW, PTRARW, KEEP,KEEP8, MYID, COMM, NBRECORDS,
     &
     & A, LA, root, PROCNODE_STEPS, SLAVEF, PERM, STEP,
     & ICNTL, INFO, NSEND8, NLOCAL8,
     & ISTEP_TO_INIV2, CANDIDATES
     & )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      INTEGER N
      INTEGER(8) :: NZ_loc8
      TYPE (DMUMPS_STRUC) :: id
      INTEGER(8) :: LDBLARR, LINTARR
      DOUBLE PRECISION DBLARR( LDBLARR )
      INTEGER INTARR( LINTARR )
      INTEGER(8), INTENT(IN) :: PTRAIW( N ), PTRARW( N )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID, COMM, NBRECORDS
      INTEGER(8) :: LA
      INTEGER SLAVEF
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      INTEGER CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      DOUBLE PRECISION A( LA )
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER PROCNODE_STEPS(KEEP(28)), PERM( N ), STEP( N )
      INTEGER INFO( 80 ), ICNTL(60)
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mpif.h'
      INTEGER :: IERR, MSGSOU
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IW4
      INTEGER END_MSG_2_RECV
      INTEGER I
      INTEGER(8) :: I18, IA8
      INTEGER(8) :: K8
      INTEGER TYPE_NODE, DEST
      INTEGER IOLD, JOLD, IARR, ISEND, JSEND, INEW, JNEW
      INTEGER allocok,  TYPESPLIT, T4MASTER, INIV2
      LOGICAL T4_MASTER_CONCERNED
      DOUBLE PRECISION VAL
      INTEGER(8) :: PTR_ROOT
      INTEGER LOCAL_M, LOCAL_N, ARROW_ROOT
      INTEGER IROW_GRID, JCOL_GRID, IPOSROOT, JPOSROOT
      INTEGER MP,LP
      INTEGER KPROBE, FREQPROBE
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: BUFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: BUFR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFRECI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: BUFRECR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IACT, IREQI, IREQR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: SEND_ACTIVE
      LOGICAL FLAG
      INTEGER(8), INTENT(OUT) :: NSEND8, NLOCAL8
      INTEGER MASTER_NODE, ISTEP
      NSEND8  = 0_8
      NLOCAL8 = 0_8
      LP = ICNTL(1)
      MP = ICNTL(2)
      END_MSG_2_RECV = SLAVEF
      ALLOCATE( IACT(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IACT in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQI(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQI in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( IREQR(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating IREQR in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( SEND_ACTIVE(SLAVEF), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating SEND_ACTIVE in matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = SLAVEF
        GOTO 20
      END IF
      ALLOCATE( BUFI( NBRECORDS * 2 + 1, 2, SLAVEF ), stat=allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating int buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = ( NBRECORDS * 2 + 1 ) * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFR( NBRECORDS, 2, SLAVEF), stat = allocok)
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &     '** Error allocating real buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * SLAVEF * 2
        GOTO 20
      END IF
      ALLOCATE( BUFRECI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS * 2 + 1
        GOTO 20
      END IF
      ALLOCATE( BUFRECR( NBRECORDS ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        IF ( LP > 0 ) THEN
          WRITE(LP,*)
     &    '** Error allocating int recv buffer for matrix distribution'
        END IF
        INFO(1) = -13
        INFO(2) = NBRECORDS
        GOTO 20
      END IF
      ALLOCATE( IW4( N, 2 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        WRITE(LP,*) '** Error allocating IW4 for matrix distribution'
        INFO(1) = -13
        INFO(2) = N * 2
      END IF
 20   CONTINUE
      CALL MUMPS_PROPINFO( ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT. 0 ) GOTO 100
      ARROW_ROOT = 0
      DO I = 1, N
          I18 = PTRAIW( I )
          IA8 = PTRARW( I )
          IF ( IA8 .GT. 0_8 ) THEN
            DBLARR( IA8 ) = ZERO
            IW4( I, 1 ) = INTARR( I18 )
            IW4( I, 2 ) = -INTARR( I18 + 1_8 )
            INTARR( I18 + 2_8 ) = I
          END IF
      END DO
      IF ( KEEP(38) .NE. 0 ) THEN
          IF (KEEP(60)==0) THEN
          LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &               root%MYROW, 0, root%NPROW )
          LOCAL_M = max( 1, LOCAL_M )
          LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &               root%MYCOL, 0, root%NPCOL )
          PTR_ROOT = LA - int(LOCAL_M,8) * int(LOCAL_N,8) + 1_8
          IF ( PTR_ROOT .LE. LA ) THEN
            A( PTR_ROOT:LA ) = ZERO
          END IF
          ELSE
            DO I = 1, root%SCHUR_NLOC
              root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     &        (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC) = ZERO
            ENDDO
          ENDIF
      END IF
      DO I = 1, SLAVEF
        BUFI( 1, 1, I ) = 0
      END DO
      DO I = 1, SLAVEF
        BUFI( 1, 2, I ) = 0
      END DO
      DO I = 1, SLAVEF
        SEND_ACTIVE( I ) = .FALSE.
        IACT( I ) = 1
      END DO
      KPROBE = 0
      FREQPROBE = max(1,NBRECORDS/10)
      DO K8 = 1_8, NZ_loc8
        KPROBE = KPROBE + 1
        IF ( KPROBE .eq. FREQPROBE ) THEN
          KPROBE = 0
          CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                     FLAG, STATUS, IERR )
          IF ( FLAG ) THEN
            MSGSOU = STATUS( MPI_SOURCE )
            CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, 
     &                 MPI_INTEGER,
     &                 MSGSOU, ARR_INT, COMM, STATUS, IERR )
            CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
            CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &             BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &             A, LA,
     &             END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &             ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR
     &             )
          END IF
        END IF
        IOLD = id%IRN_loc(K8)
        JOLD = id%JCN_loc(K8)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) CYCLE
        VAL = id%A_loc(K8)
        IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) THEN
          VAL = VAL * id%ROWSCA(IOLD)*id%COLSCA(JOLD)
        ENDIF
        IF (IOLD.EQ.JOLD) THEN
          ISEND = IOLD
          JSEND = JOLD
        ELSE
          INEW = PERM(IOLD)
          JNEW = PERM(JOLD)
          IF (INEW.LT.JNEW) THEN
            ISEND = IOLD
            IF ( KEEP(50) .NE. 0 ) ISEND = -IOLD
            JSEND = JOLD
          ELSE
            ISEND = -JOLD
            JSEND = IOLD
          ENDIF
        ENDIF
        IARR = abs( ISEND )
        ISTEP = abs(STEP(IARR))
        TYPE_NODE = MUMPS_TYPENODE(   PROCNODE_STEPS(ISTEP),
     &                                SLAVEF )
        MASTER_NODE= MUMPS_PROCNODE(  PROCNODE_STEPS(ISTEP),
     &                                SLAVEF )
        TYPESPLIT  = MUMPS_TYPESPLIT( PROCNODE_STEPS(ISTEP),
     &                               SLAVEF )
        T4_MASTER_CONCERNED = .FALSE.
        T4MASTER            = -9999
        IF (TYPE_NODE.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
          T4_MASTER_CONCERNED = .TRUE.
          T4MASTER=CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
         ENDIF
        ENDIF
        IF ( TYPE_NODE .eq. 1 ) THEN
          DEST = MASTER_NODE
        ELSE IF ( TYPE_NODE .eq. 2 ) THEN
          IF ( ISEND .LT. 0 ) THEN
            DEST = -1
          ELSE
            DEST = MASTER_NODE
          END IF
        ELSE
          IF ( ISEND < 0 ) THEN
            IPOSROOT = root%RG2L_ROW(JSEND)
            JPOSROOT = root%RG2L_ROW(IARR )
          ELSE
            IPOSROOT = root%RG2L_ROW(IARR )
            JPOSROOT = root%RG2L_ROW(JSEND)
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          DEST = IROW_GRID * root%NPCOL + JCOL_GRID
        END IF
        if (DEST .eq. -1) then
          NLOCAL8 = NLOCAL8 + 1_8
          NSEND8  = NSEND8 + int(SLAVEF -1,8)
        else
          if (DEST .eq.MYID ) then
            NLOCAL8 = NLOCAL8 + 1_8
          else
            NSEND8 = NSEND8 + 1_8
          endif
        end if
        IF ( DEST.EQ.-1) THEN
         DO I=1, CANDIDATES(SLAVEF+1,ISTEP_TO_INIV2(ISTEP))
            DEST=CANDIDATES(I,ISTEP_TO_INIV2(ISTEP))
            CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     &   root, KEEP,KEEP8 )
         ENDDO
         DEST=MASTER_NODE
         CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &   root, KEEP,KEEP8 )
         IF (T4_MASTER_CONCERNED) THEN
          DEST = T4MASTER
          CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &    BUFI, BUFR, BUFRECI, BUFRECR,
     &    NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &    SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &    N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &    PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &    root, KEEP,KEEP8 )
         ENDIF
        ELSE
         CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &   BUFI, BUFR, BUFRECI, BUFRECR,
     &   NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &   SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &   N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &   PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1), 
     &   root, KEEP,KEEP8 )
         IF (T4_MASTER_CONCERNED) THEN
          DEST = T4MASTER
          CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &    BUFI, BUFR, BUFRECI, BUFRECR,
     &    NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &    SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &    N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &    PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4(1,1),
     &    root, KEEP,KEEP8 )
         ENDIF
        ENDIF
      END DO
      DEST = -2
        CALL DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &  BUFI, BUFR, BUFRECI, BUFRECR,
     &  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR,
     &  N, PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, 
     &  IW4(1,1), root, KEEP,KEEP8 )
      DO WHILE ( END_MSG_2_RECV .NE. 0 )
        CALL MPI_RECV( BUFRECI(1), NBRECORDS * 2 + 1, MPI_INTEGER,
     &                 MPI_ANY_SOURCE, ARR_INT, COMM, STATUS, IERR )
        MSGSOU = STATUS( MPI_SOURCE )
        CALL MPI_RECV( BUFRECR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                 MSGSOU, ARR_REAL, COMM, STATUS, IERR )
        CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &           BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &           KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &           A, LA,
     &           END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &           ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &           INTARR, LINTARR, DBLARR, LDBLARR
     &           )
      END DO
      DO I = 1, SLAVEF
        IF ( SEND_ACTIVE( I ) ) THEN
          CALL MPI_WAIT( IREQI( I ), STATUS, IERR )
          CALL MPI_WAIT( IREQR( I ), STATUS, IERR )
        END IF
      END DO
      KEEP(49) = ARROW_ROOT
 100  CONTINUE
      IF (ALLOCATED(IW4))     DEALLOCATE( IW4 )
      IF (ALLOCATED(BUFI))    DEALLOCATE( BUFI )
      IF (ALLOCATED(BUFR))    DEALLOCATE( BUFR )
      IF (ALLOCATED(BUFRECI)) DEALLOCATE( BUFRECI )
      IF (ALLOCATED(BUFRECR)) DEALLOCATE( BUFRECR )
      IF (ALLOCATED(IACT))    DEALLOCATE( IACT )
      IF (ALLOCATED(IREQI))   DEALLOCATE( IREQI )
      IF (ALLOCATED(IREQR))   DEALLOCATE( IREQR )
      IF (ALLOCATED(SEND_ACTIVE)) DEALLOCATE( SEND_ACTIVE )
      RETURN
      END SUBROUTINE DMUMPS_REDISTRIBUTION
      SUBROUTINE DMUMPS_DIST_FILL_BUFFER( DEST, ISEND, JSEND, VAL,
     &  BUFI, BUFR, BUFRECI, BUFRECR,
     &  NBRECORDS, SLAVEF, COMM, MYID, IACT, IREQI, IREQR,
     &  SEND_ACTIVE, INTARR, LINTARR, DBLARR, LDBLARR, N,
     &  PTRAIW, PTRARW, PERM, STEP, ARROW_ROOT, END_MSG_2_RECV,
     &  PROCNODE_STEPS, A, LA, PTR_ROOT, LOCAL_M, LOCAL_N, IW4, root,
     &  KEEP,KEEP8 )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER ISEND, JSEND, DEST, NBRECORDS, SLAVEF, COMM, MYID, N
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER ARROW_ROOT, END_MSG_2_RECV, LOCAL_M, LOCAL_N
      INTEGER(8) :: LINTARR, LDBLARR
      INTEGER(8) :: LA, PTR_ROOT
      INTEGER BUFI( NBRECORDS * 2 + 1, 2, SLAVEF )
      INTEGER BUFRECI( NBRECORDS * 2 + 1 )
      INTEGER IREQI(SLAVEF), IREQR(SLAVEF), IACT(SLAVEF)
      INTEGER IW4( N, 2 )
      INTEGER(8) PTRAIW( N ), PTRARW( N )
      INTEGER PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER INTARR( LINTARR )
      DOUBLE PRECISION DBLARR( LDBLARR ), A( LA )
      LOGICAL SEND_ACTIVE(SLAVEF)
      DOUBLE PRECISION BUFR( NBRECORDS, 2, SLAVEF )
      DOUBLE PRECISION BUFRECR( NBRECORDS )
      DOUBLE PRECISION VAL
      INTEGER ISLAVE, IBEG, IEND, NBREC, IREQ
      INTEGER TAILLE_SEND_I, TAILLE_SEND_R, MSGSOU
      LOGICAL FLAG, SEND_LOCAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      IF ( DEST .eq. -2 ) THEN
        IBEG = 1
        IEND = SLAVEF
      ELSE
        IBEG = DEST + 1
        IEND = DEST + 1
      END IF
      SEND_LOCAL = .FALSE.
      DO ISLAVE = IBEG, IEND
        NBREC = BUFI(1,IACT(ISLAVE),ISLAVE)
        IF ( DEST .eq. -2 ) THEN
          BUFI(1,IACT(ISLAVE),ISLAVE) = - NBREC
        END IF
        IF ( DEST .eq. -2 .or. NBREC + 1 > NBRECORDS ) THEN
          DO WHILE ( SEND_ACTIVE( ISLAVE ) )
            CALL MPI_TEST( IREQR( ISLAVE ), FLAG, STATUS, IERR )
            IF ( .NOT. FLAG ) THEN
                CALL MPI_IPROBE( MPI_ANY_SOURCE, ARR_INT, COMM,
     &                           FLAG, STATUS, IERR )
                IF ( FLAG ) THEN
                  MSGSOU = STATUS(MPI_SOURCE)
                  CALL MPI_RECV( BUFRECI(1), 2*NBRECORDS+1,
     &                  MPI_INTEGER, MSGSOU, ARR_INT, COMM,
     &                  STATUS, IERR )
                  CALL MPI_RECV( BUFRECR(1), NBRECORDS,
     &                  MPI_DOUBLE_PRECISION, MSGSOU,
     &                  ARR_REAL, COMM, STATUS, IERR )
                  CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &              BUFRECI, BUFRECR, NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
                END IF
            ELSE
                CALL MPI_WAIT( IREQI( ISLAVE ), STATUS, IERR )
                SEND_ACTIVE( ISLAVE ) = .FALSE.
            END IF
          END DO
          IF ( ISLAVE - 1 .ne. MYID ) THEN
            TAILLE_SEND_I = NBREC * 2 + 1
            TAILLE_SEND_R = NBREC
            CALL MPI_ISEND( BUFI(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_I,
     &          MPI_INTEGER, ISLAVE - 1, ARR_INT, COMM,
     &          IREQI( ISLAVE ), IERR )
            CALL MPI_ISEND( BUFR(1, IACT(ISLAVE), ISLAVE ),
     &          TAILLE_SEND_R,
     &          MPI_DOUBLE_PRECISION, ISLAVE - 1, ARR_REAL, COMM,
     &          IREQR( ISLAVE ), IERR )
            SEND_ACTIVE( ISLAVE ) = .TRUE.
          ELSE
            SEND_LOCAL = .TRUE.
          END IF
          IACT( ISLAVE ) = 3 - IACT( ISLAVE )
          BUFI( 1, IACT( ISLAVE ), ISLAVE ) = 0
        END IF
        IF ( DEST .ne. -2 ) THEN
          IREQ = BUFI(1,IACT(ISLAVE),ISLAVE) + 1
          BUFI(1,IACT(ISLAVE),ISLAVE) = IREQ
          BUFI(IREQ*2,IACT(ISLAVE),ISLAVE)  = ISEND
          BUFI(IREQ*2+1,IACT(ISLAVE),ISLAVE) = JSEND
          BUFR(IREQ,IACT(ISLAVE),ISLAVE )    = VAL
        END IF
      END DO
      IF ( SEND_LOCAL ) THEN
            ISLAVE = MYID + 1
            CALL DMUMPS_DIST_TREAT_RECV_BUF(
     &              BUFI(1,3-IACT(ISLAVE),ISLAVE),
     &              BUFR(1,3-IACT(ISLAVE),ISLAVE),
     &              NBRECORDS, N, IW4(1,1),
     &              KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT,
     &              A, LA,
     &              END_MSG_2_RECV, MYID, PROCNODE_STEPS, SLAVEF,
     &              ARROW_ROOT, PTRAIW, PTRARW, PERM, STEP,
     &              INTARR, LINTARR, DBLARR, LDBLARR
     &              )
      END IF
      RETURN
      END SUBROUTINE DMUMPS_DIST_FILL_BUFFER
      SUBROUTINE DMUMPS_DIST_TREAT_RECV_BUF
     &           ( BUFI, BUFR, NBRECORDS, N, IW4,
     &             KEEP,KEEP8, LOCAL_M, LOCAL_N, root, PTR_ROOT, A, LA,
     &             END_MSG_2_RECV, MYID, PROCNODE_STEPS,
     &             SLAVEF, ARROW_ROOT,
     &             PTRAIW, PTRARW, PERM, STEP,
     &             INTARR, LINTARR, DBLARR, LDBLARR )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER NBRECORDS, N, ARROW_ROOT, MYID, SLAVEF
      INTEGER BUFI( NBRECORDS * 2 + 1 )
      DOUBLE PRECISION BUFR( NBRECORDS )
      INTEGER IW4( N, 2 )
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER END_MSG_2_RECV
      INTEGER(8) :: PTRAIW( N ), PTRARW( N )
      INTEGER :: PERM( N ), STEP( N )
      INTEGER PROCNODE_STEPS( KEEP(28) )
      INTEGER(8), INTENT(IN) :: LINTARR, LDBLARR
      INTEGER INTARR( LINTARR )
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: PTR_ROOT, LA
      DOUBLE PRECISION A( LA ), DBLARR( LDBLARR )
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      INTEGER IREC, NB_REC, NODE_TYPE, IPROC
      INTEGER IPOSROOT, JPOSROOT, IROW_GRID, JCOL_GRID,
     &        ILOCROOT, JLOCROOT
      INTEGER(8) :: IA8, IS18, IIW8, IS8, IAS8
      INTEGER ISHIFT, IARR, JARR
      INTEGER TAILLE
      DOUBLE PRECISION VAL
      NB_REC = BUFI( 1 )
      IF ( NB_REC .LE. 0 ) THEN
        END_MSG_2_RECV = END_MSG_2_RECV - 1
        NB_REC = - NB_REC
      END IF
      IF ( NB_REC .eq. 0 ) GOTO 100
      DO IREC = 1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        NODE_TYPE = MUMPS_TYPENODE( 
     &              PROCNODE_STEPS(abs(STEP(abs( IARR )))),
     &              SLAVEF )
        IF ( NODE_TYPE .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( IROW_GRID .NE. root%MYROW .OR.
     &         JCOL_GRID .NE. root%MYCOL ) THEN
            WRITE(*,*) MYID,':INTERNAL Error: recvd root arrowhead '
            WRITE(*,*) MYID,':not belonging to me. IARR,JARR=',IARR,JARR
            WRITE(*,*) MYID,':IROW_GRID,JCOL_GRID=',IROW_GRID,JCOL_GRID
            WRITE(*,*) MYID,':MYROW, MYCOL=', root%MYROW, root%MYCOL
            WRITE(*,*) MYID,':IPOSROOT,JPOSROOT=', IPOSROOT, JPOSROOT
            CALL MUMPS_ABORT()
          END IF
          ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
          JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
          IF (KEEP(60)==0) THEN
            A( PTR_ROOT + int(JLOCROOT-1,8) * int(LOCAL_M,8)
     &        + int(ILOCROOT-1,8)) =  A( PTR_ROOT
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
          IA8 = PTRARW(IARR)
          DBLARR(IA8) = DBLARR(IA8) + VAL
         ELSE
          IS18         = PTRAIW(IARR)
          ISHIFT       = INTARR(IS18) + IW4(IARR,2)
          IW4(IARR,2)  = IW4(IARR,2) - 1
          IIW8         = IS18 + ISHIFT + 2
          INTARR(IIW8) = JARR
          IS8          = PTRARW(IARR)
          IAS8         = IS8 + ISHIFT
          DBLARR(IAS8) = VAL
         ENDIF
        ELSE
           IARR = -IARR
           IS8          = PTRAIW(IARR)+IW4(IARR,1)+2
           INTARR(IS8)  = JARR
           IAS8         = PTRARW(IARR)+IW4(IARR,1)
           IW4(IARR,1)  = IW4(IARR,1) - 1
           DBLARR(IAS8) = VAL
           IPROC = MUMPS_PROCNODE( PROCNODE_STEPS(abs(STEP(IARR))),
     &                             SLAVEF )
           IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0)
     &          .AND.
     &          IW4(IARR,1) .EQ. 0 .AND.
     &          IPROC .EQ. MYID
     &          .AND. STEP(IARR) > 0 ) THEN
             TAILLE = INTARR( PTRAIW(IARR) )
             CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &            INTARR( PTRAIW(IARR) + 3 ),
     &            DBLARR( PTRARW(IARR) + 1 ),
     &            TAILLE, 1, TAILLE )
           END IF
        ENDIF
      ENDDO
 100  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_DIST_TREAT_RECV_BUF
