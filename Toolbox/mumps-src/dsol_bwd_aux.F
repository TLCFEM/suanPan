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
      RECURSIVE SUBROUTINE DMUMPS_BACKSLV_RECV_AND_TREAT(
     &     BLOQ, FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP, FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP, KEEP8, DKEEP,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &     NRHS, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &     , TO_PROCESS, SIZE_TO_PROCESS
     &     , FROM_PP
     &    )
      IMPLICIT NONE
      LOGICAL BLOQ, FLAG
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER(8), intent(in) :: LWC
      DOUBLE PRECISION W( LWC )
      INTEGER POSIWCB
      INTEGER IIPOOL, LPOOL
      INTEGER IPOOL( LPOOL )
      INTEGER LPANEL_POS
      INTEGER PANEL_POS( LPANEL_POS )
      INTEGER NBFINF, INFO(80), KEEP(500)
      INTEGER(8) :: POSWCB, PLEFTW
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER PROCNODE_STEPS( KEEP(28) ), FRERE( KEEP(28) )
      INTEGER PTRICB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER(8) :: PTRACB(KEEP(28))
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER PTRIST(KEEP(28)), IW( LIW )
      INTEGER (8) :: PTRFAC(KEEP(28))
      DOUBLE PRECISION A( LA ), W2( KEEP(133) )
      INTEGER NRHS
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP_BWD(N)
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP(NRHS,LRHSCOMP)
#else
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
#endif
      INTEGER SIZE_TO_PROCESS
      LOGICAL TO_PROCESS(SIZE_TO_PROCESS)
      LOGICAL, intent(in) :: FROM_PP
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MSGSOU, MSGTAG, MSGLEN
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      FLAG = .FALSE.
      IF ( BLOQ ) THEN
           CALL MPI_PROBE( MPI_ANY_SOURCE, MPI_ANY_TAG,
     &                     COMM, STATUS, IERR )
        FLAG = .TRUE.
      ELSE
        CALL MPI_IPROBE( MPI_ANY_SOURCE, MPI_ANY_TAG, COMM,
     &                   FLAG, STATUS, IERR )
      END IF
      IF (FLAG) THEN
         KEEP(266)=KEEP(266)-1
         MSGSOU=STATUS(MPI_SOURCE)
         MSGTAG=STATUS(MPI_TAG)
         CALL MPI_GET_COUNT( STATUS, MPI_PACKED, MSGLEN, IERR )
         IF ( MSGLEN .GT. LBUFR_BYTES ) THEN
           INFO(1) = -20
           INFO(2) = MSGLEN
           CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
         ELSE
           CALL MPI_RECV(BUFR, LBUFR_BYTES, MPI_PACKED, MSGSOU,
     &                   MSGTAG, COMM, STATUS, IERR)
           CALL DMUMPS_BACKSLV_TRAITER_MESSAGE( MSGTAG, MSGSOU,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &                FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     &                KEEP, KEEP8, DKEEP,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                NRHS, MTYPE, 
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &               , TO_PROCESS, SIZE_TO_PROCESS
     &               , FROM_PP
     &          )
         END IF
      END IF
      RETURN
      END SUBROUTINE DMUMPS_BACKSLV_RECV_AND_TREAT
      RECURSIVE SUBROUTINE DMUMPS_BACKSLV_TRAITER_MESSAGE(
     &                MSGTAG, MSGSOU,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &                FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     &                KEEP, KEEP8, DKEEP,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                NRHS, MTYPE, 
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &               , TO_PROCESS, SIZE_TO_PROCESS
     &               , FROM_PP
     &           )
      USE DMUMPS_OOC
      USE DMUMPS_SOL_LR, ONLY: DMUMPS_SOL_SLAVE_LR_U,
     &                         DMUMPS_SOL_BWD_LR_SU
      USE DMUMPS_BUF
      IMPLICIT NONE
      INTEGER MSGTAG, MSGSOU
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR( LBUFR )
      INTEGER MYID, SLAVEF, COMM
      INTEGER N, LIWW
      INTEGER IWCB( LIWW )
      INTEGER(8), intent(in) :: LWC
      DOUBLE PRECISION W( LWC )
      INTEGER POSIWCB
      INTEGER IIPOOL, LPOOL, LPANEL_POS
      INTEGER IPOOL( LPOOL )
      INTEGER PANEL_POS( LPANEL_POS )
      INTEGER NBFINF, INFO(80), KEEP(500)
      INTEGER(8) :: POSWCB, PLEFTW
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER PTRICB(KEEP(28)), STEP( N ), FILS( N )
      INTEGER(8) :: PTRACB(KEEP(28))
      INTEGER FRERE(KEEP(28))
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER LIW
      INTEGER(8) :: LA
      INTEGER IW( LIW ), PTRIST( KEEP(28) )
      INTEGER(8) :: PTRFAC(KEEP(28))
      DOUBLE PRECISION A( LA ), W2( KEEP(133) )
      INTEGER NRHS
      INTEGER MYLEAFE, MTYPE
      INTEGER LRHSCOMP, POSINRHSCOMP_BWD(N)
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP(NRHS,LRHSCOMP)
#else
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
#endif
      INTEGER TMP_NBPANELS, I_PIVRPTR, I_PIVR
      LOGICAL MUST_BE_PERMUTED
      INTEGER  SIZE_TO_PROCESS
      LOGICAL TO_PROCESS(SIZE_TO_PROCESS), NO_CHILDREN
      LOGICAL, intent(in) :: FROM_PP
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED
      LOGICAL  OOCWRITE_COMPATIBLE_WITH_BLR 
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER POSITION, IF, INODE, IERR, LONG, DUMMY(1)
      INTEGER    :: LIELL, K
      INTEGER(8) :: APOS, IST
      INTEGER NPIV, NROW_L, IPOS, NROW_RECU
      INTEGER(8) :: IFR8
      INTEGER I, JJ, IN, PROCDEST, J1, J2, LDA
      INTEGER NSLAVES, NELIM, J, POSINDICES, INODEPOS,
     &     IPOSINRHSCOMP, IPOSINRHSCOMP_TMP, IPOSINRHSCOMP_PANEL
      INTEGER JBDEB, JBFIN, NRHS_B, allocok
      INTEGER(8) :: P_UPDATE, P_SOL_MAS         
      INTEGER :: IWHDLR, MTYPE_SLAVE, LDA_SLAVE 
      LOGICAL FLAG
      DOUBLE PRECISION ZERO, ALPHA, ONE
      PARAMETER (ZERO=0.0D0, ONE = 1.0D0, ALPHA=-1.0D0)
      INCLUDE 'mumps_headers.h'
      INTEGER POOL_FIRST_POS, TMP
      LOGICAL,DIMENSION(:),ALLOCATABLE :: DEJA_SEND
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE, dtrsv, dtrsm, dgemv, dgemm
      INTEGER :: NCB
      INTEGER(8) :: APOSDEB, NBENTRIES_ALLPANELS
      INTEGER(8) :: PTWCB, PTWCB_PANEL
      INTEGER LDAJ, NBJ, LIWFAC,
     &        NBJLAST, NPIV_LAST, PANEL_SIZE,
     &        NCB_PANEL, TYPEF
      LOGICAL TWOBYTWO
      INTEGER BEG_PANEL
      INTEGER IPANEL, NPANELS
      ALLOCATE(DEJA_SEND( 0:SLAVEF-1 ), stat=allocok)
      if(allocok.ne.0) then
         INFO(1)=-13
         INFO(2)=SLAVEF
         WRITE(6,*) MYID,' Allocation error of DEJA_SEND '
     &        //'in bwd solve COMPSO'
         GOTO 260
      END IF
      DUMMY(1)=0
      IF (MSGTAG .EQ. FEUILLE) THEN
          NBFINF = NBFINF - 1
      ELSE IF (MSGTAG .EQ. NOEUD) THEN
          POSITION = 0
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        INODE, 1, MPI_INTEGER,
     &        COMM, IERR)
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBDEB, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &        JBFIN, 1, MPI_INTEGER, COMM, IERR )
          CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &        LONG, 1, MPI_INTEGER,
     &        COMM, IERR)
         NRHS_B = JBFIN-JBDEB+1
          IF (   POSIWCB - LONG .LT. 0
     &      .OR. POSWCB - PLEFTW + 1_8 .LT. LONG ) THEN
            CALL DMUMPS_COMPSO(N, KEEP(28), IWCB,
     &      LIWW, W, LWC,
     &      POSWCB, POSIWCB, PTRICB, PTRACB)
            IF (POSIWCB - LONG .LT. 0) THEN
              INFO(1)=-14
              INFO(2)=-POSIWCB + LONG
              WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
              GOTO 260
            END IF
            IF ( POSWCB - PLEFTW + 1_8 .LT. LONG ) THEN
              INFO(1) = -11
              CALL MUMPS_SET_IERROR(LONG + PLEFTW - POSWCB - 1_8,
     &                             INFO(2))
              WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
              GOTO 260
            END IF
          ENDIF
          POSIWCB = POSIWCB - LONG
          POSWCB = POSWCB - LONG
          IF (LONG .GT. 0) THEN
            CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          IWCB(POSIWCB + 1), 
     &          LONG, MPI_INTEGER, COMM, IERR)
            DO K=JBDEB,JBFIN
             CALL MPI_UNPACK(BUFR, LBUFR_BYTES, POSITION,
     &          W(POSWCB + 1), LONG, 
     &          MPI_DOUBLE_PRECISION, COMM, IERR)
             DO JJ=0, LONG-1
              IPOSINRHSCOMP = abs(POSINRHSCOMP_BWD(IWCB(POSIWCB+1+JJ)))
              IF ( (IPOSINRHSCOMP.EQ.0) .OR.
     &           ( IPOSINRHSCOMP.GT.N ) ) CYCLE  
#if defined(RHSCOMP_BYROWS)
              RHSCOMP(K,IPOSINRHSCOMP) = W(POSWCB+1+JJ)
#else
              RHSCOMP(IPOSINRHSCOMP,K) = W(POSWCB+1+JJ)
#endif
             ENDDO
            ENDDO
            POSIWCB = POSIWCB + LONG
            POSWCB = POSWCB + LONG
          ENDIF
          POOL_FIRST_POS = IIPOOL
          IF ( KEEP(237).GT. 0 ) THEN
             IF (.NOT.TO_PROCESS(STEP(INODE))) 
     &            GOTO 1010
          ENDIF
             IPOOL( IIPOOL ) = INODE
             IIPOOL = IIPOOL + 1
 1010     CONTINUE
          IF = FRERE( STEP(INODE) )
          DO WHILE ( IF .GT. 0 )
             IF ( MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),
     &            SLAVEF) .eq. MYID ) THEN
                IF ( KEEP(237).GT. 0 ) THEN
                   IF (.NOT.TO_PROCESS(STEP(IF))) THEN
                      IF = FRERE(STEP(IF))
                      CYCLE
                   ENDIF
                ENDIF
                   IPOOL( IIPOOL ) = IF
                   IIPOOL = IIPOOL + 1
             END IF
             IF = FRERE( STEP( IF ) )
          END DO
             DO I=1,(IIPOOL-POOL_FIRST_POS)/2
                TMP=IPOOL(POOL_FIRST_POS+I-1)
                IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
                IPOOL(IIPOOL-I)=TMP
             ENDDO      
      ELSE IF ( MSGTAG .EQ. BACKSLV_MASTER2SLAVE ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   INODE, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NROW_RECU, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   JBDEB, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   JBFIN, 1, MPI_INTEGER, COMM, IERR )
        NRHS_B = JBFIN-JBDEB+1
        LR_ACTIVATED   = (IW(PTRIST(STEP(INODE))+XXLR).GT.0)
        COMPRESS_PANEL = (IW(PTRIST(STEP(INODE))+XXLR).GE.2)
        OOCWRITE_COMPATIBLE_WITH_BLR =
     &     ( .NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(485).EQ.0) 
     &     )
        IPOS   = PTRIST( STEP(INODE) ) + KEEP(IXSZ)
        NPIV   = - IW( IPOS     )
        NROW_L =   IW( IPOS + 1 )
        IF ( NROW_L .NE. NROW_RECU ) THEN
          WRITE(*,*) 'Error1 : NROW L/RECU=',NROW_L, NROW_RECU
          CALL MUMPS_ABORT()
        END IF
        LONG = NROW_L + NPIV
        IF ( POSWCB - LONG*NRHS_B .LT. PLEFTW - 1_8 ) THEN
           CALL DMUMPS_COMPSO( N, KEEP(28), IWCB,
     &          LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
           IF ( POSWCB - LONG*NRHS_B .LT. PLEFTW - 1_8 ) THEN
             INFO(1) = -11
             CALL MUMPS_SET_IERROR(LONG * NRHS_B- POSWCB,INFO(2))
             WRITE(6,*) MYID,' Internal error in bwd solve COMPSO'
             GOTO 260
           END IF
        END IF
        P_UPDATE  = PLEFTW
        P_SOL_MAS = PLEFTW + NPIV * NRHS_B
        PLEFTW    = P_SOL_MAS + NROW_L * NRHS_B
        DO K=JBDEB, JBFIN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   W( P_SOL_MAS+(K-JBDEB)*NROW_L),NROW_L,
     &                   MPI_DOUBLE_PRECISION,
     &                   COMM, IERR )
        ENDDO
        IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
           CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &     INODE,PTRFAC,KEEP,A,LA,STEP,
     &     KEEP8,N,MUST_BE_PERMUTED,IERR)           
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
        ENDIF
        APOS   =   PTRFAC( STEP(INODE))
        IF ( IW(PTRIST(STEP(INODE))+XXLR) .GE. 2 .AND.
     &       KEEP(485) .EQ. 1 ) THEN
          IWHDLR = IW(PTRIST(STEP(INODE))+XXF)
          MTYPE_SLAVE = 0 
          W(P_UPDATE:P_UPDATE+NPIV*NRHS_B-1)=ZERO
          CALL DMUMPS_SOL_SLAVE_LR_U(INODE, IWHDLR, -9999,
     &         W, LWC,
     &         NROW_L, NPIV, 
     &         P_SOL_MAS, P_UPDATE,
     &         JBDEB, JBFIN,
     &         MTYPE_SLAVE, KEEP,
     &         INFO(1), INFO(2) )
        ELSE
          IF (KEEP(201) .EQ. 1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &    THEN 
            MTYPE_SLAVE = 1
            LDA_SLAVE = NROW_L
          ELSE
            MTYPE_SLAVE = 0
            LDA_SLAVE = NPIV
          ENDIF
          CALL DMUMPS_SOLVE_GEMM_UPDATE(
     &           A, LA, APOS, NROW_L,
     &           LDA_SLAVE,    
     &           NPIV,
     &           NRHS_B, W, LWC,
     &           P_SOL_MAS, NROW_L,
     &           P_UPDATE, NPIV,
     &           MTYPE_SLAVE, KEEP, ZERO)
        ENDIF
        IF (KEEP(201) .EQ. 1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &  THEN 
         CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &          A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 260
         ENDIF
        ENDIF
        PLEFTW = PLEFTW - NROW_L * NRHS_B
 100    CONTINUE
        CALL DMUMPS_BUF_SEND_BACKVEC( NRHS_B, INODE, 
     &                               W(P_UPDATE),
     &                               NPIV, NPIV,
     &                                MSGSOU, 
     &                                BACKSLV_UPDATERHS,
     &                                JBDEB, JBFIN,
     &                                KEEP, COMM, IERR )
        IF ( IERR .EQ. -1 ) THEN
          CALL DMUMPS_BACKSLV_RECV_AND_TREAT(
     &     .FALSE., FLAG,
     &     BUFR, LBUFR, LBUFR_BYTES,
     &     MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &     FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     &     KEEP, KEEP8, DKEEP,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE,
     &     NRHS, MTYPE,
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &      , TO_PROCESS, SIZE_TO_PROCESS
     &      , FROM_PP
     &      )
          IF ( INFO( 1 ) .LT. 0 ) GOTO 270
          GOTO 100
       ELSE IF ( IERR .EQ. -2 ) THEN
          INFO( 1 ) = -17
          INFO( 2 ) = NRHS_B * NPIV * KEEP(35) + 4 * KEEP(34)
          GOTO 260
       ELSE IF ( IERR .EQ. -3 ) THEN
          INFO( 1 ) = -20
          INFO( 2 ) = NRHS_B * NPIV * KEEP(35) + 4 * KEEP(34)
          GOTO 260
        END IF
        PLEFTW = PLEFTW - NPIV * NRHS_B
      ELSE IF ( MSGTAG .EQ. BACKSLV_UPDATERHS ) THEN
        POSITION = 0
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   INODE, 1, MPI_INTEGER, COMM, IERR )
        LR_ACTIVATED   = (IW(PTRIST(STEP(INODE))+XXLR).GT.0)
        COMPRESS_PANEL = (IW(PTRIST(STEP(INODE))+XXLR).GE.2)
        OOCWRITE_COMPATIBLE_WITH_BLR =
     &     (.NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(485).EQ.0) 
     &     )
        IPOS  = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
        LIELL = IW(IPOS-2)+IW(IPOS+1)
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   NPIV, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   JBDEB, 1, MPI_INTEGER, COMM, IERR )
        CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                   JBFIN, 1, MPI_INTEGER, COMM, IERR )
        NRHS_B = JBFIN-JBDEB+1
        NELIM = IW(IPOS-1)
        IPOS = IPOS + 1
        NPIV = IW(IPOS)
        IPOS = IPOS + 1
        NSLAVES = IW( IPOS + 1 )
        IPOS = IPOS + 1 + NSLAVES
        INODEPOS = PTRIST(STEP(INODE)) + KEEP(IXSZ) + 4
        IF ( KEEP(50) .eq. 0 ) THEN
          LDA = LIELL
        ELSE
          LDA = NPIV
        ENDIF
        IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
             J1 = IPOS + LIELL + 1
             J2 = IPOS + NPIV + LIELL
        ELSE
             J1 = IPOS + 1
             J2 = IPOS + NPIV
        ENDIF
        IPOSINRHSCOMP =  POSINRHSCOMP_BWD(IW(J1)) 
        DO K=JBDEB, JBFIN
          CALL MPI_UNPACK( BUFR, LBUFR_BYTES, POSITION,
     &                     W2, NPIV, MPI_DOUBLE_PRECISION,
     &                     COMM, IERR )
          I = 1
          IF ( (KEEP(253).NE.0) .AND.
     &         (IW(PTRIST(STEP(INODE))+XXS).EQ.C_FINI+NSLAVES) 
     &       ) THEN
          DO JJ = J1,J2   
#if defined(RHSCOMP_BYROWS)
            RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) = W2(I)
#else
            RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = W2(I)
#endif
            I = I+1
          ENDDO
         ELSE
          DO JJ = J1,J2   
#if defined(RHSCOMP_BYROWS)
            RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) = 
     &      RHSCOMP(K,IPOSINRHSCOMP+JJ-J1) + W2(I)
#else
            RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = 
     &      RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) + W2(I)
#endif
            I = I+1
          ENDDO
         ENDIF
        ENDDO  
        IW(PTRIST(STEP(INODE))+XXS) = 
     &      IW(PTRIST(STEP(INODE))+XXS) - 1
        IF ( IW(PTRIST(STEP(INODE))+XXS).EQ.C_FINI ) THEN
          IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &    THEN
             CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &            INODE,PTRFAC,KEEP,A,LA,STEP,
     &            KEEP8,N,MUST_BE_PERMUTED,IERR)
             IF(IERR.LT.0)THEN
                INFO(1)=IERR
                INFO(2)=0
                GOTO 260
             ENDIF
             IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1) THEN
               CALL DMUMPS_OOC_PP_CHECK_PERM_FREED(
     &              IW(IPOS+1+2*LIELL),
     &              MUST_BE_PERMUTED )
             ENDIF
          ENDIF  
          APOS = PTRFAC(IW(INODEPOS))
          IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) 
     &    THEN 
             LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
             TYPEF = TYPEF_L
             NROW_L   = NPIV+NELIM 
             PANEL_SIZE = DMUMPS_OOC_PANEL_SIZE(NROW_L)
             IF (PANEL_SIZE.LT.0) THEN
               WRITE(6,*) ' Internal error in bwd solve PANEL_SIZE=',
     &         PANEL_SIZE
               CALL MUMPS_ABORT()
             ENDIF
          ENDIF 
           IF ( POSIWCB - 2 .LT. 0 .or.
     &         POSWCB - LIELL*NRHS_B .LT. PLEFTW - 1_8 ) THEN
            CALL DMUMPS_COMPSO( N, KEEP(28), IWCB, 
     &          LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
            IF ( POSWCB - LIELL*NRHS_B .LT. PLEFTW - 1_8 ) THEN
              INFO( 1 ) = -11
              CALL MUMPS_SET_IERROR( LIELL*NRHS_B - POSWCB-PLEFTW+1_8,
     &                              INFO(2) )
              GOTO 260
            END IF
            IF ( POSIWCB - 2 .LT. 0 ) THEN
              INFO( 1 ) = -14
              INFO( 2 ) = 2 - POSIWCB
              GO TO 260
            END IF
           END IF
           POSIWCB = POSIWCB - 2
           POSWCB  = POSWCB - LIELL*NRHS_B
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1_8
           IWCB( PTRICB(STEP( INODE ))     ) = LIELL*NRHS_B
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IPOS = PTRIST(STEP(INODE)) + KEEP(IXSZ) + 5 + NSLAVES
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
             POSINDICES = IPOS + LIELL + 1
           ELSE
             POSINDICES = IPOS + 1
           END IF
           IPOSINRHSCOMP =  POSINRHSCOMP_BWD(IW(J1)) 
           IFR8 = PTRACB(STEP( INODE ))
           IF (KEEP(350).EQ.0) THEN
             DO K=JBDEB, JBFIN
               DO JJ = J1, J2
                 W(IFR8+JJ-J1+(K-JBDEB)*LIELL) = 
#if defined(RHSCOMP_BYROWS)
     &           RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else
     &           RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif
               END DO
             END DO
           ELSE IF (KEEP(350).eq.1.OR.KEEP(350).EQ.2) THEN
           ELSE
             WRITE(*,*) "Internal error DMUMPS_BACKSLV_TRAITER_MESSAGE"
             CALL MUMPS_ABORT()
           ENDIF
           IFR8 = PTRACB(STEP(INODE))+NPIV-1
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
             J1 = IPOS + LIELL + NPIV + 1
             J2 = IPOS + 2 * LIELL
           ELSE
             J1 = IPOS + NPIV + 1
             J2 = IPOS + LIELL
           END IF
           IF (KEEP(350).EQ.0) THEN
             DO JJ = J1, J2-KEEP(253)   
               J = IW(JJ)
               IFR8 = IFR8 + 1
               IPOSINRHSCOMP_TMP =  abs(POSINRHSCOMP_BWD(J))
               DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                 W(IFR8+(K-JBDEB)*LIELL) = RHSCOMP(K,IPOSINRHSCOMP_TMP)
#else
                 W(IFR8+(K-JBDEB)*LIELL) = RHSCOMP(IPOSINRHSCOMP_TMP,K)
#endif
               ENDDO
             ENDDO
           ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
             CALL DMUMPS_SOL_BWD_GTHR( JBDEB, JBFIN, J1, J2,
     &       RHSCOMP, NRHS, LRHSCOMP,
     &       W(PTRACB(STEP(INODE))), LIELL, NPIV+1,
     &       IW, LIW, KEEP, N, POSINRHSCOMP_BWD )
             IFR8 = IFR8 + J2-KEEP(253)-J1+1
          ELSE
             WRITE(*,*) "Internal error DMUMPS_BACKSLV_TRAITER_MESSAGE"
             CALL MUMPS_ABORT()
           ENDIF
       IF ( KEEP(201).EQ.1 .AND. OOCWRITE_COMPATIBLE_WITH_BLR .AND.
     &    (( NELIM .GT. 0 ).OR. (MTYPE.NE.1 )))  THEN
          J = NPIV / PANEL_SIZE  
          TWOBYTWO = KEEP(50).EQ.2 .AND. KEEP(105).GT.0
          IF (TWOBYTWO) THEN
            CALL DMUMPS_BUILD_PANEL_POS(PANEL_SIZE, PANEL_POS,
     &           LPANEL_POS, IW(IPOS+1+LIELL), NPIV, NPANELS,
     &           NROW_L, NBENTRIES_ALLPANELS)
          ELSE
            IF (NPIV.EQ.J*PANEL_SIZE) THEN
              NPIV_LAST = NPIV
              NBJLAST   = PANEL_SIZE
              NPANELS   = J
            ELSE
              NPIV_LAST = (J+1)* PANEL_SIZE
              NBJLAST   = NPIV-J*PANEL_SIZE
              NPANELS   = J+1
            ENDIF
            NBENTRIES_ALLPANELS =
     &  int(NROW_L,8) * int(NPIV,8) 
     &  - int( ( J * ( J - 1 ) ) /2,8 ) 
     &    * int(PANEL_SIZE,8) * int(PANEL_SIZE,8) 
     &  - int(J,8)                       
     &    * int(mod(NPIV, PANEL_SIZE),8) 
     &    * int(PANEL_SIZE,8)    
            JJ=NPIV_LAST
          ENDIF
          APOSDEB = APOS + NBENTRIES_ALLPANELS 
          DO IPANEL=NPANELS,1,-1
            IF (TWOBYTWO) THEN
              NBJ = PANEL_POS(IPANEL+1)-PANEL_POS(IPANEL)
              BEG_PANEL = PANEL_POS(IPANEL)
            ELSE
              IF (JJ.EQ.NPIV_LAST) THEN
                NBJ = NBJLAST
              ELSE
                NBJ = PANEL_SIZE
              ENDIF
              BEG_PANEL = JJ- PANEL_SIZE+1
            ENDIF
            LDAJ    = NROW_L-BEG_PANEL+1 
            APOSDEB = APOSDEB - int(NBJ,8)*int(LDAJ,8)
            PTWCB   = PTRACB(STEP(INODE))
            PTWCB_PANEL =  PTRACB(STEP(INODE)) + int(BEG_PANEL - 1,8)
            IPOSINRHSCOMP_PANEL = IPOSINRHSCOMP + BEG_PANEL - 1
            NCB_PANEL   = LDAJ - NBJ
            NCB     = NROW_L - NPIV
            IF (KEEP(50).NE.1 .AND.MUST_BE_PERMUTED) THEN
              CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &        I_PIVRPTR, I_PIVR, IPOS + 1 + 2 * LIELL, IW, LIW)
              CALL DMUMPS_PERMUTE_PANEL(
     &        IW(I_PIVR + IW(I_PIVRPTR+IPANEL-1)-IW(I_PIVRPTR)),
     &        NPIV-IW(I_PIVRPTR+IPANEL-1)+1,
     &        IW(I_PIVRPTR+IPANEL-1)-1,
     &        A(APOSDEB),
     &        LDAJ, NBJ, BEG_PANEL-1)
            ENDIF
#if defined(MUMPS_USE_BLAS2)
            IF ( NRHS_B == 1 ) THEN
              IF (NCB_PANEL.NE.0) THEN
                IF (KEEP(350).EQ.0) THEN
                  CALL dgemv( 'T', NCB_PANEL, NBJ, ALPHA, 
     &                A( APOSDEB + int(NBJ,8) ), LDAJ,
     &                W( PTWCB_PANEL + int(NBJ,8) ),
     &                1, ONE,
     &                W(PTWCB_PANEL), 1 )
                ELSE
                  IF (NCB_PANEL - NCB.NE. 0) THEN
                    CALL dgemv( 'T', NCB_PANEL-NCB, NBJ, ALPHA, 
     &              A( APOSDEB + int(NBJ,8) ), LDAJ,
#                   if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL+NBJ), 
     &              1, ONE,
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1 ) 
#                   else
     &              RHSCOMP(IPOSINRHSCOMP_PANEL+NBJ,JBDEB),
     &              1, ONE,
     &              RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1 )
#                   endif
                  ENDIF
                  IF (NCB .NE. 0) THEN
                    CALL dgemv( 'T', NCB, NBJ, ALPHA,
     &              A( APOSDEB + int(LDAJ-NCB,8) ), LDAJ,
     &              W( PTWCB  + NPIV ),
     &              1, ONE,
#                   if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1 )
#                   else
     &              RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1 )
#                   endif
                  ENDIF
                ENDIF
              ENDIF
              IF (KEEP(350).eq.0) THEN
                CALL dtrsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
     &              W(PTWCB_PANEL), 1)
              ELSE
                CALL dtrsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
#if defined(RHSCOMP_BYROWS)
     &            RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1)
#else
     &            RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1)
#endif
              ENDIF
            ELSE
#endif
              IF (NCB_PANEL.NE.0) THEN
                IF (KEEP(350).eq.0) THEN
                  CALL dgemm( 'T', 'N', NBJ, NRHS_B, NCB_PANEL, ALPHA,
     &              A(APOSDEB +int(NBJ,8)), LDAJ,
     &              W(PTWCB_PANEL+int(NBJ,8)),LIELL,
     &              ONE, W(PTWCB_PANEL),LIELL)
                ELSE
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*)
     &            "Internal error in DMUMPS_BACKSLV_TRAITER_MESSAGE"
                  CALL MUMPS_ABORT()
#else
                  IF (NCB_PANEL - NCB .NE. 0) THEN
                  CALL dgemm( 'T', 'N', NBJ, NRHS_B,
     &                                         NCB_PANEL-NCB, ALPHA,
     &              A(APOSDEB +int(NBJ,8)), LDAJ,
     &              RHSCOMP(IPOSINRHSCOMP_PANEL+NBJ,JBDEB), LRHSCOMP,
     &              ONE, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), LRHSCOMP)
                  ENDIF
                  IF (NCB .NE. 0) THEN
                  CALL dgemm( 'T', 'N', NBJ, NRHS_B, NCB, ALPHA,
     &              A(APOSDEB +int(LDAJ-NCB,8)), LDAJ,
     &              W( PTWCB+NPIV ), LIELL,
     &              ONE, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB),LRHSCOMP)
                  ENDIF
#endif
                ENDIF
              ENDIF
              IF (KEEP(350).eq.0) THEN
                CALL dtrsm('L','L','T','U',NBJ, NRHS_B, ONE, 
     &          A(APOSDEB),
     &          LDAJ, W(PTWCB_PANEL), LIELL)
              ELSE
#if defined(RHSCOMP_BYROWS)
                WRITE(*,*)
     &          "Internal error in DMUMPS_BACKSLV_TRAITER_MESSAGE"
                CALL MUMPS_ABORT()
#else
                CALL dtrsm('L','L','T','U',NBJ, NRHS_B, ONE,
     &          A(APOSDEB), 
     &          LDAJ, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), LRHSCOMP)
#endif
              ENDIF
#if defined(MUMPS_USE_BLAS2)
            ENDIF
#endif
            IF (.NOT. TWOBYTWO) JJ=BEG_PANEL-1
          ENDDO 
        GOTO 1234  
       ENDIF 
       IF ( IW(PTRIST(STEP(INODE))+XXLR) .GE. 2
     &      .AND. KEEP(485) .EQ. 1 ) THEN
          IWHDLR = IW(PTRIST(STEP(INODE))+XXF)
          CALL  DMUMPS_SOL_BWD_LR_SU 
     &     ( INODE, IWHDLR, NPIV, NSLAVES,
     &     LIELL, W, LWC, NRHS_B, PTRACB(STEP(INODE)), 
     &     RHSCOMP, LRHSCOMP, NRHS,
     &     IPOSINRHSCOMP, JBDEB, 
     &     MTYPE, KEEP,
     &     INFO(1), INFO(1) )
       ELSE
          IF (NELIM .GT.0) THEN
            IF ( KEEP(50) .eq. 0 ) THEN
                IST = APOS + int(NPIV,8) * int(LIELL,8)
            ELSE
                IST = APOS + int(NPIV,8) * int(NPIV,8)
            END IF
            IF ( NRHS_B == 1 ) THEN
              IF (KEEP(350).EQ.0) THEN
                CALL dgemv( 'N', NPIV, NELIM, ALPHA,
     &              A( IST ), NPIV,
     &              W( NPIV + PTRACB(STEP(INODE)) ),
     &              1, ONE,
     &              W(PTRACB(STEP(INODE))), 1 )
              ELSE
                CALL dgemv( 'N', NPIV, NELIM, ALPHA,
     &              A( IST ), NPIV,
     &              W( NPIV + PTRACB(STEP(INODE)) ),
     &              1, ONE,
#if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP), 1)
#else
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB), 1)
#endif
              ENDIF
            ELSE
              IF (KEEP(350).EQ.0) THEN
                CALL dgemm( 'N', 'N', NPIV, NRHS_B, NELIM, ALPHA,
     &                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     &                ONE, W(PTRACB(STEP(INODE))),LIELL)
              ELSE
#if defined(RHSCOMP_BYROWS)
                WRITE(*,*)
     &          "Internal error in DMUMPS_BACKSLV_TRAITER_MESSAGE"
                CALL MUMPS_ABORT()
#else
                CALL dgemm( 'N', 'N', NPIV, NRHS_B, NELIM, ALPHA,
     &                A(IST), NPIV, W(NPIV+PTRACB(STEP(INODE))),LIELL,
     &                ONE, RHSCOMP(IPOSINRHSCOMP,JBDEB), LRHSCOMP)
#endif
              ENDIF
            END IF
          ENDIF 
#if defined(MUMPS_USE_BLAS2)
          IF ( NRHS_B == 1 ) THEN
            IF (KEEP(350).EQ.0) THEN
              CALL dtrsv( 'U', 'N', 'U', NPIV, A(APOS), LDA,
     &                  W(PTRACB(STEP(INODE))),1)
            ELSE
              CALL dtrsv( 'U', 'N', 'U', NPIV, A(APOS), LDA,
#if defined(RHSCOMP_BYROWS)
     &                RHSCOMP(JBDEB,IPOSINRHSCOMP), 1)
#else
     &                RHSCOMP(IPOSINRHSCOMP,JBDEB), 1)
#endif
            ENDIF
          ELSE
#endif
             IF (KEEP(350).EQ.0) THEN
               CALL dtrsm( 'L','U', 'N', 'U', NPIV, NRHS_B, ONE,
     &                   A(APOS), LDA,
     &                   W(PTRACB(STEP(INODE))),LIELL)
             ELSE
#if defined(RHSCOMP_BYROWS)
                WRITE(*,*)
     &          "Internal error in DMUMPS_BACKSLV_TRAITER_MESSAGE"
                CALL MUMPS_ABORT()
#else
                CALL dtrsm( 'L','U', 'N', 'U', NPIV, NRHS_B, ONE,
     &                    A(APOS), LDA,
     &                    RHSCOMP(IPOSINRHSCOMP,JBDEB), LRHSCOMP)
#endif
             ENDIF
#if defined(MUMPS_USE_BLAS2)
          END IF
#endif
       ENDIF
 1234     CONTINUE   
          IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
           CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &          A,LA,.TRUE.,IERR)
           IF(IERR.LT.0)THEN
              INFO(1)=IERR
              INFO(2)=0
              GOTO 260
           ENDIF
          ENDIF
          IPOS =   PTRIST(STEP(INODE)) +  KEEP(IXSZ) + 6 + NSLAVES
          IPOSINRHSCOMP     = POSINRHSCOMP_BWD(IW(IPOS))
          IF (KEEP(350).EQ.0) THEN
            IPOSINRHSCOMP_TMP = IPOSINRHSCOMP
            DO I = 1, NPIV
              DO K=JBDEB,JBFIN
#if defined(RHSCOMP_BYROWS)
                RHSCOMP( K, IPOSINRHSCOMP_TMP ) = 
     &           W( PTRACB(STEP(INODE))+I-1 + (K-JBDEB)*LIELL )
#else
                RHSCOMP( IPOSINRHSCOMP_TMP , K ) = 
     &           W( PTRACB(STEP(INODE))+I-1 + (K-JBDEB)*LIELL )
#endif
              ENDDO
              IPOSINRHSCOMP_TMP =  IPOSINRHSCOMP_TMP + 1
            END DO
          ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
          ELSE
            WRITE(*,*)"Internal error in DMUMPS_BACKSLV_TRAITER_MESSAGE"
            CALL MUMPS_ABORT()
          ENDIF
          IN = INODE
  200     IN = FILS(IN)
          IF (IN .GT. 0) GOTO 200
          IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
              CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &                       FEUILLE, SLAVEF, KEEP )
              NBFINF = NBFINF - 1
            ENDIF
            IWCB( PTRICB(STEP(INODE)) + 1 ) = 0
            CALL DMUMPS_FREETOPSO(N, KEEP(28),
     &          IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
            GOTO 270
          ENDIF  
          DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
          END DO
          IN = -IN
          IF ( KEEP(237).GT.0 ) THEN 
            NO_CHILDREN = .TRUE. 
          ELSE
            NO_CHILDREN = .FALSE.
          ENDIF
          DO WHILE (IN.GT.0) 
            IF ( KEEP(237).GT.0 ) THEN
               IF (.NOT.TO_PROCESS(STEP(IN))) THEN
                  IN = FRERE(STEP(IN))
                  CYCLE
               ELSE
                 NO_CHILDREN = .FALSE.
               ENDIF
            ENDIF
           POOL_FIRST_POS  = IIPOOL
            IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IN)),
     &          SLAVEF) .EQ. MYID) THEN
                  IPOOL(IIPOOL ) = IN
                  IIPOOL = IIPOOL + 1
            ELSE
              PROCDEST = MUMPS_PROCNODE( PROCNODE_STEPS(STEP(IN)), 
     &                   SLAVEF )
              IF ( .NOT. DEJA_SEND( PROCDEST ) ) THEN
 110            CONTINUE
                CALL DMUMPS_BUF_SEND_VCB( NRHS_B, IN, 0, 0,
     &          LIELL, LIELL-KEEP(253),
     &          IW( POSINDICES ) ,
     &          W( PTRACB(STEP(INODE)) ), JBDEB, JBFIN,
     &          RHSCOMP(1, 1), NRHS, LRHSCOMP,
     &          IPOSINRHSCOMP, NPIV, KEEP,
     &          PROCDEST, NOEUD, COMM, IERR )
                IF ( IERR .EQ. -1 ) THEN
                  CALL DMUMPS_BACKSLV_RECV_AND_TREAT(
     &            .FALSE., FLAG,
     &            BUFR, LBUFR, LBUFR_BYTES,
     &            MYID, SLAVEF, COMM,
     &            N, IWCB, LIWW, POSIWCB,
     &            W, LWC, POSWCB,
     &            IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &            IPOOL, LPOOL, PANEL_POS, LPANEL_POS, STEP,
     &            FRERE, FILS, PROCNODE_STEPS, PLEFTW,
     &            KEEP, KEEP8, DKEEP,
     &            PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &            NRHS, MTYPE, 
     &            RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &            , TO_PROCESS, SIZE_TO_PROCESS
     &            , FROM_PP
     &            )
                  IF ( INFO( 1 ) .LT. 0 ) GOTO 270
                  GOTO 110
               ELSE IF ( IERR .eq. -2 ) THEN
                  INFO(1) = -17
                  INFO(2) = LIELL * NRHS_B * KEEP(35) +
     &                    ( LIELL + 4 ) * KEEP(34)
                  GOTO 260
               ELSE IF ( IERR .eq. -3 ) THEN
                  INFO(1) = -20
                  INFO(2) = LIELL * NRHS_B * KEEP(35) +
     &                    ( LIELL + 4 ) * KEEP(34)
                  GOTO 260
                END IF
                DEJA_SEND( PROCDEST ) = .TRUE.
              END IF
            END IF
            IN = FRERE( STEP( IN ) )
          END DO
          IF (NO_CHILDREN) THEN
                   MYLEAFE = MYLEAFE - 1
                   IF (MYLEAFE .EQ. 0) THEN
                      CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, 
     &                     COMM, FEUILLE, SLAVEF, KEEP )
                      NBFINF = NBFINF - 1
                   ENDIF
                   IWCB( PTRICB(STEP(INODE)) + 1 ) = 0
                   CALL DMUMPS_FREETOPSO( N, KEEP(28),
     &                  IWCB, LIWW, W, LWC,
     &                  POSWCB, POSIWCB, PTRICB, PTRACB)
                   GOTO 270
           ENDIF
          DO I=1,(IIPOOL-POOL_FIRST_POS)/2
           TMP=IPOOL(POOL_FIRST_POS+I-1)
           IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
           IPOOL(IIPOOL-I)=TMP
          ENDDO 
          IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
          CALL DMUMPS_FREETOPSO( N, KEEP(28),
     &          IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
        END IF   
      ELSE IF (MSGTAG.EQ.TERREUR) THEN
          INFO(1) = -001
          INFO(2) = MSGSOU
          GO TO 270
       ELSE IF ( (MSGTAG.EQ.UPDATE_LOAD).OR.
     &      (MSGTAG.EQ.TAG_DUMMY) ) THEN
          GO TO 270
       ELSE
          INFO(1) = -100
          INFO(2) = MSGTAG
          GOTO 260
      ENDIF
      GO TO 270
 260  CONTINUE
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
 270  CONTINUE
      DEALLOCATE(DEJA_SEND)
      RETURN
      CALL DMUMPS_BDC_ERROR( MYID, SLAVEF, COMM, KEEP )
      RETURN
      END SUBROUTINE DMUMPS_BACKSLV_TRAITER_MESSAGE
      SUBROUTINE DMUMPS_BUILD_PANEL_POS(PANEL_SIZE, PANEL_POS,
     &                           LEN_PANEL_POS, INDICES, NPIV,
     &                           NPANELS, NFRONT_OR_NASS,
     &                           NBENTRIES_ALLPANELS)
      IMPLICIT NONE
      INTEGER, intent (in)   :: PANEL_SIZE, NPIV
      INTEGER, intent (in)   :: INDICES(NPIV)
      INTEGER, intent (in)   :: LEN_PANEL_POS
      INTEGER, intent (out)  :: NPANELS
      INTEGER, intent (out)  :: PANEL_POS(LEN_PANEL_POS)
      INTEGER, intent (in)   :: NFRONT_OR_NASS
      INTEGER(8), intent(out):: NBENTRIES_ALLPANELS
      INTEGER NPANELS_MAX, I, NBeff
      INTEGER(8) :: NBENTRIES_THISPANEL
      NBENTRIES_ALLPANELS = 0_8
      NPANELS_MAX = (NPIV+PANEL_SIZE-1)/PANEL_SIZE
      IF (LEN_PANEL_POS .LT. NPANELS_MAX + 1) THEN
        WRITE(*,*) "Error 1 in DMUMPS_BUILD_PANEL_POS",
     &              LEN_PANEL_POS,NPANELS_MAX
        CALL MUMPS_ABORT()
      ENDIF
      I = 1
      NPANELS = 0
      IF (I .GT. NPIV) RETURN 
 10   CONTINUE
      NPANELS = NPANELS + 1
      PANEL_POS(NPANELS) = I
      NBeff = min(PANEL_SIZE, NPIV-I+1)
      IF ( INDICES(I+NBeff-1) < 0) THEN
        NBeff=NBeff+1
      ENDIF
      NBENTRIES_THISPANEL = int(NFRONT_OR_NASS-I+1,8) * int(NBeff,8)
      NBENTRIES_ALLPANELS = NBENTRIES_ALLPANELS + NBENTRIES_THISPANEL
      I=I+NBeff
      IF ( I .LE. NPIV ) GOTO 10
      PANEL_POS(NPANELS+1)=NPIV+1
      RETURN
      END SUBROUTINE DMUMPS_BUILD_PANEL_POS
