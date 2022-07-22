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
      SUBROUTINE DMUMPS_SOL_S(N, A, LA, IW, LIW, W, LWC,
     &    NRHS, 
     &    RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD,
     &    PTRICB, PTRACB, IWCB, LIWW, W2, 
     &    NE_STEPS, NA, LNA, STEP,
     &    FRERE, DAD, FILS, IPOOL, LPOOL, PTRIST, PTRFAC, 
     &    MYLEAF, ICNTL, INFO, 
     &    PROCNODE_STEPS,
     &    SLAVEF, COMM,MYID, BUFR, LBUFR, LBUFR_BYTES,
     &    KEEP,KEEP8, DKEEP, RHS_ROOT, LRHS_ROOT, MTYPE, 
     &
     &    ISTEP_TO_INIV2, TAB_POS_IN_PERE, PANEL_POS, LPANEL_POS
     &    , TO_PROCESS, SIZE_TO_PROCESS
     &    , RHS_BOUNDS, LRHS_BOUNDS, DO_NBSPARSE, FROM_PP
     &    )
      USE DMUMPS_OOC
      USE DMUMPS_BUF
      USE DMUMPS_SOL_LR, only : DMUMPS_SOL_BWD_LR_SU
      IMPLICIT NONE
      INTEGER MTYPE
      INTEGER(8) :: LA
      INTEGER(8), intent(in) :: LWC
      INTEGER N,LIW,LIWW,LPOOL,LNA
      INTEGER SLAVEF,MYLEAF,COMM,MYID
      INTEGER LPANEL_POS
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION, INTENT(INOUT) :: DKEEP(230)
      INTEGER PROCNODE_STEPS(KEEP(28))
      INTEGER NA(LNA),NE_STEPS(KEEP(28))
      INTEGER IPOOL(LPOOL)
      INTEGER PANEL_POS(LPANEL_POS)
      INTEGER ICNTL(60), INFO(80)
      INTEGER PTRIST(KEEP(28)),
     &        PTRICB(KEEP(28))
      INTEGER(8) :: PTRACB(KEEP(28))
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER NRHS
      DOUBLE PRECISION A(LA), W(LWC)
      DOUBLE PRECISION W2(KEEP(133))
      INTEGER IW(LIW),IWCB(LIWW)
      INTEGER STEP(N), FRERE(KEEP(28)),DAD(KEEP(28)),FILS(N)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER BUFR(LBUFR)
      INTEGER ISTEP_TO_INIV2(KEEP(71)), 
     &        TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER LRHSCOMP, POSINRHSCOMP_BWD(N)
#if defined(RHSCOMP_BYROWS)
      DOUBLE PRECISION RHSCOMP(NRHS,LRHSCOMP)
#else
      DOUBLE PRECISION RHSCOMP(LRHSCOMP,NRHS)
#endif
      INTEGER(8), intent(in) :: LRHS_ROOT
      DOUBLE PRECISION RHS_ROOT( LRHS_ROOT )
      INTEGER, intent(in)           :: SIZE_TO_PROCESS
      LOGICAL, intent(in)           :: TO_PROCESS(SIZE_TO_PROCESS)
      LOGICAL, intent(in) :: DO_NBSPARSE
      INTEGER, intent(in) :: LRHS_BOUNDS
      INTEGER, intent(in) :: RHS_BOUNDS(LRHS_BOUNDS)
      LOGICAL, intent(in) :: FROM_PP
      LOGICAL COMPRESS_PANEL, LR_ACTIVATED
      LOGICAL  OOCWRITE_COMPATIBLE_WITH_BLR 
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      LOGICAL FLAG
      INTEGER POSIWCB,K
      INTEGER(8) :: APOS, IST
      INTEGER(8) :: IFR
      INTEGER NPIV
      INTEGER IPOS,LIELL,NELIM,JJ,I
      INTEGER J1,J2,J,NCB,NBFINF
      INTEGER NBLEAF,INODE,NBROOT,NROOT,NBFILS
      INTEGER IN,IF,LONG,POOL_FIRST_POS,TMP
      INTEGER III,IIPOOL,MYLEAFE
      INTEGER NSLAVES
      INTEGER JBDEB, JBFIN, NRHS_B
      DOUBLE PRECISION ALPHA,ONE,ZERO
      PARAMETER (ZERO=0.0D0, ONE = 1.0D0, ALPHA=-1.0D0)
      LOGICAL BLOQ,DEBUT
      INTEGER PROCDEST, DEST
      INTEGER POSINDICES, IPOSINRHSCOMP, IPOSINRHSCOMP_PANEL,
     &        IPOSINRHSCOMP_TMP
      INTEGER DUMMY(1)
      INTEGER(8) :: POSWCB, PLEFTW, PTWCB
      INTEGER Offset, EffectiveSize, ISLAVE, FirstIndex
      LOGICAL LTLEVEL2, IN_SUBTREE
      INTEGER TYPENODE
      INCLUDE 'mumps_headers.h'
      LOGICAL BLOCK_SEQUENCE
      INTEGER TMP_NBPANELS, I_PIVRPTR, I_PIVR
      LOGICAL MUST_BE_PERMUTED
      LOGICAL NO_CHILDREN
      LOGICAL Exploit_Sparsity, AM1 
      LOGICAL, DIMENSION(:), ALLOCATABLE :: DEJA_SEND
      INTEGER :: allocok
      INTEGER(8) :: APOSDEB, NBENTRIES_ALLPANELS
      INTEGER(8) :: PTWCB_PANEL
      INTEGER LDAJ, NBJ, LIWFAC,
     &        NBJLAST, NPIV_LAST, PANEL_SIZE,
     &        NCB_PANEL, TYPEF
      INTEGER BEG_PANEL
      LOGICAL TWOBYTWO
      INTEGER NPANELS, IPANEL
      INTEGER :: IWHDLR
      LOGICAL MUMPS_IN_OR_ROOT_SSARBR
      INTEGER MUMPS_TYPENODE
      EXTERNAL dgemv, dtrsv, dtrsm, dgemm,
     &         MUMPS_TYPENODE, 
     &         MUMPS_IN_OR_ROOT_SSARBR
      DUMMY(1)=0
      KEEP(266)=0
      ALLOCATE(DEJA_SEND( 0:SLAVEF-1 ), stat=allocok)
      if(allocok.ne.0) then
         WRITE(6,*) ' Allocation error of DEJA_SEND in '
     &        //'routine DMUMPS_SOL_S '
         INFO(1)=-13
         INFO(2)=SLAVEF
      endif
      CALL MUMPS_PROPINFO(ICNTL, INFO, COMM, MYID )
      IF ( INFO(1) .LT.0 ) GOTO 340
      PLEFTW = 1_8
      POSIWCB = LIWW
      POSWCB = LWC
      NROOT = 0
      NBLEAF = NA(1)
      NBROOT = NA(2)
      DO I = NBROOT, 1, -1
        INODE = NA(NBLEAF+I+2)
        IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &      SLAVEF) .EQ. MYID) THEN
          NROOT = NROOT + 1
          IPOOL(NROOT) = INODE
        ENDIF
      END DO
      III = 1
      IIPOOL = NROOT + 1
      BLOCK_SEQUENCE = .FALSE.
      Exploit_Sparsity = .FALSE.
      AM1 = .FALSE.
      IF (KEEP(235).NE.0) Exploit_Sparsity = .TRUE.
      IF (KEEP(237).NE.0) AM1 = .TRUE.
      NO_CHILDREN = .FALSE.
      IF (Exploit_Sparsity .OR. AM1) MYLEAF = -1
      IF (MYLEAF .EQ. -1) THEN
        MYLEAF = 0
        DO I=1, NBLEAF
          INODE=NA(I+2)
          IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF) .EQ. MYID) THEN
            MYLEAF = MYLEAF + 1
          ENDIF
        ENDDO
      ENDIF
      MYLEAFE=MYLEAF
      NBFINF = SLAVEF
      IF (MYLEAFE .EQ. 0) THEN
        CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM, FEUILLE,
     &                  SLAVEF, KEEP)
        NBFINF = NBFINF - 1
        IF (NBFINF .EQ. 0) THEN
          GOTO 340
        ENDIF
      ENDIF
 50   CONTINUE
      BLOQ = ( (  III .EQ. IIPOOL  )
     &     )
      CALL DMUMPS_BACKSLV_RECV_AND_TREAT( BLOQ, FLAG, BUFR, LBUFR,
     &     LBUFR_BYTES, MYID, SLAVEF, COMM,
     &     N, IWCB, LIWW, POSIWCB,
     &     W, LWC, POSWCB,
     &     IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &     IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &     STEP,  FRERE, FILS, PROCNODE_STEPS,
     &     PLEFTW, KEEP,KEEP8, DKEEP,
     &     PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &     NRHS, MTYPE, 
     &     RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &     , TO_PROCESS, SIZE_TO_PROCESS
     &     , FROM_PP
     &     )
      IF ( INFO(1) .LT. 0 ) GOTO 340
      IF ( .NOT. FLAG ) THEN
        IF (III .NE. IIPOOL) THEN
          INODE = IPOOL(IIPOOL-1)
          IIPOOL = IIPOOL - 1
          GO TO 60
        ENDIF
      END IF                    
      IF ( NBFINF .eq. 0 ) GOTO 340
      GOTO 50
   60 CONTINUE
      IF (DO_NBSPARSE) THEN
        JBDEB= RHS_BOUNDS(2*STEP(INODE)-1)
        JBFIN= RHS_BOUNDS(2*STEP(INODE))
        NRHS_B = JBFIN-JBDEB+1
      ELSE
        JBDEB = 1
        JBFIN = NRHS
        NRHS_B = NRHS
      ENDIF
      IF ( INODE .EQ. KEEP( 38 ) .OR. INODE .EQ. KEEP( 20 ) ) THEN
         IPOS = PTRIST(STEP(INODE))+KEEP(IXSZ)
         NPIV  = IW(IPOS+3)
         LIELL = IW(IPOS) + NPIV  
         IPOS =  PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ)
         IF ( MTYPE .EQ. 1 .AND. KEEP(50) .EQ. 0) THEN
            J1   = IPOS + LIELL + 1
            J2   = IPOS + LIELL + NPIV
         ELSE
            J1   = IPOS + 1
            J2   = IPOS + NPIV
         END IF
         IFR  = 0_8
         IPOSINRHSCOMP = POSINRHSCOMP_BWD(IW(J1))  
         IF (KEEP(350).EQ.0) THEN
           IPOSINRHSCOMP_TMP = IPOSINRHSCOMP
           DO JJ = J1, J2
              IFR = IFR + 1_8
              DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
                RHSCOMP(K,IPOSINRHSCOMP_TMP) = RHS_ROOT(IFR+NPIV*(K-1))
#else
                RHSCOMP(IPOSINRHSCOMP_TMP,K) = RHS_ROOT(IFR+NPIV*(K-1))
#endif
              END DO
              IPOSINRHSCOMP_TMP =  IPOSINRHSCOMP_TMP + 1  
           END DO
         ELSE
           CALL DMUMPS_SOL_CPY_FS2RHSCOMP(JBDEB, JBFIN, J2-J1+1,
     &          KEEP, RHSCOMP, NRHS, LRHSCOMP, IPOSINRHSCOMP,
     &          RHS_ROOT(1+NPIV*(JBDEB-1)), NPIV, 1)
         ENDIF
         IN = INODE
 270     IN = FILS(IN)
         IF (IN .GT. 0) GOTO 270
         IF (IN .EQ. 0) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF, KEEP )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
            ENDIF
            GOTO 50
         ENDIF
         IF   = -IN
         LONG = NPIV
         NBFILS = NE_STEPS(STEP(INODE))
         IF ( AM1 ) THEN
            I = NBFILS
            NBFILS = 0
            DO WHILE (I.GT.0)
               IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
               IF = FRERE(STEP(IF))
               I = I -1
            ENDDO
            IF (NBFILS.EQ.0) THEN
               NO_CHILDREN = .TRUE.
            ELSE
               NO_CHILDREN = .FALSE.
            ENDIF
            IF = -IN
         ENDIF
         DEBUT = .TRUE.
         DO I = 0, SLAVEF - 1
            DEJA_SEND( I ) = .FALSE.
         END DO
         POOL_FIRST_POS=IIPOOL
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1030          IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1030
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
            IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),SLAVEF)
     &           .EQ. MYID) THEN
                  IPOOL(IIPOOL) = IF
                  IIPOOL = IIPOOL + 1
            ELSE
               PROCDEST = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),
     &              SLAVEF)
               IF (.NOT. DEJA_SEND( PROCDEST ))  THEN
 600              CONTINUE
                  CALL DMUMPS_BUF_SEND_VCB( NRHS_B, IF, 0, 0,
     &                 LONG, LONG, IW( J1 ),
     &                 RHS_ROOT( 1+NPIV*(JBDEB-1) ), 
     &                 JBDEB, JBFIN,
     &                 RHSCOMP(1, 1), NRHS, LRHSCOMP,
     &                 IPOSINRHSCOMP, NPIV,
     &                 KEEP, PROCDEST,
     &                 NOEUD, COMM, IERR )
                  IF ( IERR .EQ. -1 ) THEN
                     CALL DMUMPS_BACKSLV_RECV_AND_TREAT(
     &                    .FALSE., FLAG,
     &                    BUFR, LBUFR, LBUFR_BYTES,
     &                    MYID, SLAVEF, COMM,
     &                    N, IWCB, LIWW, POSIWCB,
     &                    W, LWC, POSWCB,
     &                    IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                    IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                    STEP, FRERE, FILS, PROCNODE_STEPS,
     &                    PLEFTW, KEEP,KEEP8, DKEEP,
     &                    PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                    NRHS, MTYPE,
     &                    RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                    , TO_PROCESS, SIZE_TO_PROCESS
     &                    , FROM_PP
     &                    )
                     IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                     GOTO 600
                  ELSE IF ( IERR .EQ. -2 ) THEN
                     INFO( 1 ) = -17
                     INFO( 2 ) = NRHS_B * LONG * KEEP(35) +
     &                    ( LONG + 4 ) * KEEP(34)
                     GOTO 330
                  ELSE IF ( IERR .EQ. -3 ) THEN
                     INFO( 1 ) = -20
                     INFO( 2 ) = NRHS_B * LONG * KEEP(35) +
     &                    ( LONG + 4 ) * KEEP(34)
                     GOTO 330
                  END IF
                  DEJA_SEND( PROCDEST ) = .TRUE.
               END IF
               IF ( IERR .NE. 0 ) THEN
                  CALL MUMPS_ABORT()
               ENDIF
            ENDIF
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND.NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF, KEEP )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
            IF (IIPOOL.NE.POOL_FIRST_POS) THEN
               DO I=1,(IIPOOL-POOL_FIRST_POS)/2
                  TMP=IPOOL(POOL_FIRST_POS+I-1)
                  IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
                  IPOOL(IIPOOL-I)=TMP
               ENDDO
            ENDIF
         GOTO 50
      END IF
      IN_SUBTREE = MUMPS_IN_OR_ROOT_SSARBR( 
     &               PROCNODE_STEPS(STEP(INODE)), SLAVEF ) 
      TYPENODE = MUMPS_TYPENODE(PROCNODE_STEPS(STEP(INODE)),
     &         SLAVEF)
      LTLEVEL2= ( 
     &   (TYPENODE .eq.2 ) .AND.
     &   (MTYPE.NE.1)   )
      NPIV = IW(PTRIST(STEP(INODE))+2+KEEP(IXSZ)+1)
      IF ((NPIV.NE.0).AND.(LTLEVEL2)) THEN
            IPOS  = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
            LIELL = IW(IPOS-2)+IW(IPOS+1)
            NELIM = IW(IPOS-1)
            IPOS  = IPOS + 1
            NPIV  = IW(IPOS)
            NCB   = LIELL - NPIV - NELIM
            IPOS  = IPOS + 2
            NSLAVES = IW( IPOS )
            Offset = 0  
            IPOS = IPOS + NSLAVES   
            IW(PTRIST(STEP(INODE))+XXS)= C_FINI+NSLAVES
           IF ( POSIWCB - 2 .LT. 0 .or.
     &          POSWCB-int(NCB,8)*int(NRHS_B,8) .LT. PLEFTW-1_8 ) THEN
             CALL DMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
             IF ( POSWCB-int(NCB,8)*int(NRHS_B,8) .LT. PLEFTW-1_8 ) THEN
               INFO( 1 ) = -11
               CALL MUMPS_SET_IERROR(NCB * NRHS_B - POSWCB-PLEFTW+1_8,
     &                              INFO(2))
               GOTO 330
             END IF
             IF ( POSIWCB - 2 .LT. 0 ) THEN
               INFO( 1 ) = -14
               INFO( 2 ) = 2 - POSIWCB
               GO TO 330
             END IF
           END IF
           POSIWCB = POSIWCB - 2
           POSWCB  = POSWCB - int(NCB,8)*int(NRHS_B,8)
           PTRICB(STEP( INODE )) = POSIWCB + 1
           PTRACB(STEP( INODE )) = POSWCB  + 1_8
           IWCB( PTRICB(STEP( INODE ))     ) = NCB*NRHS_B
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
           IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
              POSINDICES = IPOS + LIELL + 1
           ELSE
              POSINDICES = IPOS + 1
           END IF
           IF ( NCB.EQ.0 ) THEN
             write(6,*) ' Internal Error type 2 node with no CB '
             CALL MUMPS_ABORT()
           ENDIF
           IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
               J1 = IPOS + LIELL + NPIV + NELIM +1
               J2 = IPOS + 2 * LIELL
           ELSE
               J1 = IPOS + NPIV + NELIM +1
               J2 = IPOS + LIELL
           END IF
           IFR = PTRACB(STEP( INODE )) - 1_8
           IF (KEEP(350).EQ.0) THEN
             DO JJ = J1, J2 - KEEP(253)
               J = IW(JJ)
               IFR = IFR + 1_8
               IPOSINRHSCOMP_TMP =  abs(POSINRHSCOMP_BWD(J))
               DO K=JBDEB, JBFIN
                 W(IFR+int(K-JBDEB,8)*int(NCB,8)) =
#if defined(RHSCOMP_BYROWS)
     &           RHSCOMP(K,IPOSINRHSCOMP_TMP)
#else
     &           RHSCOMP(IPOSINRHSCOMP_TMP,K)
#endif
               ENDDO
             ENDDO
           ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
             CALL DMUMPS_SOL_BWD_GTHR( JBDEB, JBFIN, J1, J2,
     &       RHSCOMP, NRHS, LRHSCOMP,
     &       W(PTRACB(STEP(INODE))), NCB, 1,
     &       IW, LIW, KEEP, N, POSINRHSCOMP_BWD )
             IFR = IFR + int(J2-KEEP(253)-J1+1,8)
          ELSE
             WRITE(*,*) "Internal error DMUMPS_SOL_S"
             CALL MUMPS_ABORT()
           END IF
           IF (KEEP(252).NE.0) THEN
             DO JJ = J2-KEEP(253)+1, J2
              IFR = IFR + 1_8
              DO K=JBDEB, JBFIN
               IF (K.EQ.JJ-J2+KEEP(253)) THEN
                 W(IFR+int(K-JBDEB,8)*int(NCB,8)) = ALPHA
               ELSE
                 W(IFR+int(K-JBDEB,8)*int(NCB,8)) = ZERO
               ENDIF
              ENDDO
             ENDDO
           ENDIF
           DO ISLAVE = 1, NSLAVES
              CALL MUMPS_BLOC2_GET_SLAVE_INFO( 
     &                KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &                ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &                ISLAVE, NCB, 
     &                NSLAVES, 
     &                EffectiveSize,
     &                FirstIndex )
 500         CONTINUE
             DEST = IW( PTRIST(STEP(INODE))+5+ISLAVE+KEEP(IXSZ))
             CALL DMUMPS_BUF_SEND_BACKVEC(NRHS_B, INODE,
     &             W(Offset+PTRACB(STEP(INODE))), 
     &             EffectiveSize, 
     &             NCB, DEST,
     &             BACKSLV_MASTER2SLAVE, JBDEB, JBFIN,
     &             KEEP, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                 CALL DMUMPS_BACKSLV_RECV_AND_TREAT(
     &                .FALSE., FLAG,
     &                BUFR, LBUFR, LBUFR_BYTES,
     &                MYID, SLAVEF, COMM,
     &                N, IWCB, LIWW, POSIWCB,
     &                W, LWC, POSWCB,
     &                IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &                IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &                STEP, FRERE, FILS,
     &                PROCNODE_STEPS, PLEFTW, KEEP,KEEP8, DKEEP,
     &                PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &                NRHS, MTYPE,
     &                RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                , FROM_PP
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 500
             ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NRHS_B * EffectiveSize * KEEP(35) +
     &               2 * KEEP(34)
                GOTO 330
             ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NRHS_B * EffectiveSize * KEEP(35) +
     &                            2 * KEEP(34)
                GOTO 330
              END IF
              Offset = Offset + EffectiveSize
           END DO
           IWCB( PTRICB(STEP( INODE )) + 1 ) = 0
           CALL DMUMPS_FREETOPSO(N, KEEP(28), IWCB, LIWW, W, LWC,
     &             POSWCB,POSIWCB,PTRICB,PTRACB)
           GOTO 50
      ENDIF   
      LR_ACTIVATED = (IW(PTRIST(STEP(INODE))+XXLR).GT.0)
      COMPRESS_PANEL = (IW(PTRIST(STEP(INODE))+XXLR).GE.2)
      OOCWRITE_COMPATIBLE_WITH_BLR =
     &     (.NOT.LR_ACTIVATED.OR.(.NOT.COMPRESS_PANEL).OR.
     &            (KEEP(485).EQ.0) 
     &     )
      IPOS = PTRIST(STEP(INODE)) + 2 + KEEP(IXSZ)
      LIELL = IW(IPOS-2)+IW(IPOS+1)
      NELIM = IW(IPOS-1)
      IPOS = IPOS + 1
      NPIV = IW(IPOS)
      NCB   = LIELL - NPIV
      IPOS = IPOS + 1
      IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
         CALL DMUMPS_SOLVE_GET_OOC_NODE(
     &        INODE,PTRFAC,KEEP,A,LA,STEP,
     &        KEEP8,N,MUST_BE_PERMUTED,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF                     
      APOS = PTRFAC(IW(IPOS))
      NSLAVES = IW( PTRIST(STEP(INODE)) + 5 + KEEP(IXSZ) )
      IPOS = IPOS + 1 + NSLAVES
      IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
           LIWFAC =  IW(PTRIST(STEP(INODE))+XXI)
           IF (MTYPE.NE.1) THEN
            TYPEF = TYPEF_L
           ELSE
            TYPEF = TYPEF_U
           ENDIF
           PANEL_SIZE =  DMUMPS_OOC_PANEL_SIZE( LIELL )
           IF (KEEP(50).NE.1) THEN
             CALL DMUMPS_OOC_PP_CHECK_PERM_FREED(
     &                   IW(IPOS+1+2*LIELL),
     &                   MUST_BE_PERMUTED )
           ENDIF
      ENDIF  
      LONG = 0
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0 ) THEN
        J1 = IPOS + LIELL + 1
        J2 = IPOS + NPIV + LIELL
      ELSE
        J1 = IPOS + 1
        J2 = IPOS + NPIV
      ENDIF
      IF (IN_SUBTREE) THEN
        PTWCB = PLEFTW
        IF ( POSWCB .LT. int(LIELL,8)*int(NRHS_B,8) ) THEN
          CALL DMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &                POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB .LT. int(LIELL,8)*int(NRHS_B,8) ) THEN
            INFO(1) = -11
            CALL MUMPS_SET_IERROR(int(LIELL,8)*int(NRHS_B,8)-POSWCB,
     &                 INFO(2))
            GOTO 330
          END IF
        END IF
      ELSE
        IF ( POSIWCB - 2 .LT. 0 .or.
     &     POSWCB-int(LIELL,8)*int(NRHS_B,8) .LT. PLEFTW-1_8 ) THEN
          CALL DMUMPS_COMPSO( N, KEEP(28), IWCB, LIWW, W, LWC,
     &          POSWCB, POSIWCB, PTRICB, PTRACB)
          IF ( POSWCB-int(LIELL,8)*int(NRHS_B,8) .LT. PLEFTW-1_8 ) THEN
            INFO( 1 ) = -11
            CALL MUMPS_SET_IERROR( int(LIELL,8)*int(NRHS_B,8)-
     &                            POSWCB-PLEFTW+1_8,
     &                            INFO(2) )
            GOTO 330
          END IF
          IF ( POSIWCB - 2 .LT. 0 ) THEN
            INFO( 1 ) = -14
            INFO( 2 ) = 2 - POSIWCB
            GO TO 330
          END IF
        END IF
        POSIWCB = POSIWCB - 2
        POSWCB  = POSWCB - int(LIELL,8)*int(NRHS_B,8)
        PTRICB(STEP( INODE )) = POSIWCB + 1
        PTRACB(STEP( INODE )) = POSWCB  + 1_8
        IWCB( PTRICB(STEP( INODE ))     ) = LIELL*NRHS_B
        IWCB( PTRICB(STEP( INODE )) + 1 ) = 1  
        IF ( MTYPE.EQ.1 .AND. KEEP(50).EQ.0 ) THEN
           POSINDICES = IPOS + LIELL + 1
        ELSE
           POSINDICES = IPOS + 1
        END IF
        PTWCB = PTRACB(STEP( INODE )) 
      ENDIF
      IF (J2.GE.J1) THEN
        IPOSINRHSCOMP = POSINRHSCOMP_BWD(IW(J1)) 
      ELSE
        IPOSINRHSCOMP = -99999
      ENDIF
      IF (J2.GE.J1) THEN
        DO K=JBDEB, JBFIN
          IF (KEEP(350).eq.0) THEN
            IF (KEEP(252).NE.0) THEN
              DO JJ = J1, J2
                W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = ZERO
              ENDDO
            ELSE
              DO JJ = J1, J2
#if defined(RHSCOMP_BYROWS)
                W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = 
     &               RHSCOMP(K,IPOSINRHSCOMP+JJ-J1)
#else          
                W(PTWCB+JJ-J1+(K-JBDEB)*LIELL) = 
     &               RHSCOMP(IPOSINRHSCOMP+JJ-J1,K)
#endif        
              ENDDO
            ENDIF
          ELSE IF (KEEP(350).eq.1.OR.KEEP(350).EQ.2) THEN
            IF (KEEP(252).NE.0) THEN
              DO JJ = J1, J2
                RHSCOMP(IPOSINRHSCOMP+JJ-J1,K) = ZERO
              ENDDO
            ENDIF
         ELSE
            WRITE(*,*) "Internal error DMUMPS_SOL_BWD"
            CALL MUMPS_ABORT()
          ENDIF
        END DO
      ENDIF
      IFR   = PTWCB + int(NPIV - 1,8)
      IF ( LIELL .GT. NPIV ) THEN
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
            IFR = IFR + 1_8
            IPOSINRHSCOMP_TMP = abs(POSINRHSCOMP_BWD(J))
            DO K=JBDEB, JBFIN
              W(IFR+int(K-JBDEB,8)*int(LIELL,8)) =
#if defined(RHSCOMP_BYROWS)
     &        RHSCOMP(K,IPOSINRHSCOMP_TMP)
#else
     &        RHSCOMP(IPOSINRHSCOMP_TMP,K)
#endif
            ENDDO
          ENDDO
        ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
          CALL DMUMPS_SOL_BWD_GTHR( JBDEB, JBFIN, J1, J2,
     &       RHSCOMP, NRHS, LRHSCOMP,
     &       W(PTWCB), LIELL, NPIV+1,
     &       IW, LIW, KEEP, N, POSINRHSCOMP_BWD )
          IFR = IFR + int(J2-KEEP(253)-J1+1,8)
       ELSE
          WRITE(*,*) "Internal error DMUMPS_SOL_S"
          CALL MUMPS_ABORT()
        ENDIF
        IF (KEEP(252).NE.0) THEN
          DO JJ = J2-KEEP(253)+1, J2
           IFR = IFR + 1_8
           DO K=JBDEB, JBFIN
            IF (K.EQ.JJ-J2+KEEP(253)) THEN
              W(IFR+int(K-JBDEB,8)*int(LIELL,8)) = ALPHA   
            ELSE
              W(IFR+int(K-JBDEB,8)*int(LIELL,8)) = ZERO
            ENDIF
           ENDDO
          ENDDO
        ENDIF
        NCB = LIELL - NPIV
        IF (NPIV .EQ. 0) GOTO 160
      ENDIF
      IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
       J = NPIV / PANEL_SIZE 
       TWOBYTWO = KEEP(50).EQ.2 .AND.
     & ((TYPENODE.EQ.1.AND.KEEP(103).GT.0) .OR.
     &  (TYPENODE.EQ.2.AND.KEEP(105).GT.0))
       IF (TWOBYTWO) THEN 
         CALL DMUMPS_BUILD_PANEL_POS(PANEL_SIZE, PANEL_POS, LPANEL_POS,
     &        IW(IPOS+1+LIELL), NPIV, NPANELS, LIELL,
     &        NBENTRIES_ALLPANELS)
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
     &  int(LIELL,8) * int(NPIV,8) 
     &  - int( ( J * ( J - 1 ) ) /2,8 ) 
     &    * int(PANEL_SIZE,8) * int(PANEL_SIZE,8) 
     &  - int(J,8)                       
     &    * int(mod(NPIV, PANEL_SIZE),8) 
     &    * int(PANEL_SIZE,8)    
         JJ=NPIV_LAST
       ENDIF
       APOSDEB = APOS + NBENTRIES_ALLPANELS 
       DO IPANEL = NPANELS, 1, -1
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
            LDAJ    = LIELL-BEG_PANEL+1 
            APOSDEB = APOSDEB - int(NBJ,8)*int(LDAJ,8)
            PTWCB_PANEL = PTWCB + int(BEG_PANEL - 1,8)
            IPOSINRHSCOMP_PANEL = IPOSINRHSCOMP + BEG_PANEL - 1
            NCB_PANEL   = LDAJ - NBJ
            IF (KEEP(50).NE.1 .AND. MUST_BE_PERMUTED) THEN
              CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF, TMP_NBPANELS,
     &        I_PIVRPTR, I_PIVR, IPOS + 1 + 2 * LIELL, IW, LIW)
              IF (NPIV.EQ.(IW(I_PIVRPTR)-1)) THEN
                MUST_BE_PERMUTED=.FALSE. 
              ELSE
               CALL DMUMPS_PERMUTE_PANEL(
     &         IW(I_PIVR + IW(I_PIVRPTR+IPANEL-1)-IW(I_PIVRPTR)),
     &         NPIV-IW(I_PIVRPTR+IPANEL-1)+1,
     &         IW(I_PIVRPTR+IPANEL-1)-1,
     &         A(APOSDEB),
     &         LDAJ, NBJ, BEG_PANEL-1)
              ENDIF
            ENDIF
#if defined(MUMPS_USE_BLAS2)
            IF ( NRHS_B == 1 ) THEN
              IF (NCB_PANEL.NE.0) THEN
                IF (KEEP(350).EQ.0) THEN
                  CALL dgemv( 'T', NCB_PANEL, NBJ, ALPHA, 
     &                A( APOSDEB + int(NBJ,8) ), LDAJ,
     &                W( PTWCB_PANEL+int(NBJ,8) ),
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
     &              W( PTWCB  + int(NPIV,8) ),
     &              1, ONE,
#                   if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1 )
#                   else
     &              RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1 )
#                   endif
                  ENDIF
                ENDIF
              ENDIF
              IF (MTYPE.NE.1) THEN
                IF (KEEP(350).eq.0) THEN
                  CALL dtrsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
     &                W(PTWCB_PANEL), 1)
                ELSE
                  CALL dtrsv('L','T','U', NBJ, A(APOSDEB), LDAJ,
#if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1)
#else
     &              RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1)
#endif
                ENDIF
              ELSE
                IF (KEEP(350).eq.0) THEN
                  CALL dtrsv('L','T','N', NBJ, A(APOSDEB), LDAJ,
     &                W(PTWCB_PANEL), 1)
                ELSE
                  CALL dtrsv('L','T','N', NBJ, A(APOSDEB), LDAJ,
#if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP_PANEL), 1)
#else
     &              RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), 1)
#endif
                ENDIF
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
                  WRITE(*,*) "Internal error in DMUMPS_SOL_S"
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
     &              W( PTWCB+int(NPIV,8) ), LIELL,
     &              ONE, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB),LRHSCOMP)
                  ENDIF
#endif
                ENDIF
              ENDIF
              IF (MTYPE.NE.1) THEN
                IF (KEEP(350).eq.0) THEN
                  CALL dtrsm('L','L','T','U',NBJ, NRHS_B, ONE, 
     &            A(APOSDEB),
     &            LDAJ, W(PTWCB_PANEL), LIELL)
                ELSE
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*) "Internal error in DMUMPS_SOL_S"
                  CALL MUMPS_ABORT()
#else
                  CALL dtrsm('L','L','T','U',NBJ, NRHS_B, ONE,
     &            A(APOSDEB), 
     &            LDAJ, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), LRHSCOMP)
#endif
                ENDIF
              ELSE
                IF (KEEP(350).eq.0) THEN
                  CALL dtrsm('L','L','T','N',NBJ, NRHS_B, ONE,
     &            A(APOSDEB), 
     &            LDAJ, W(PTWCB_PANEL), LIELL)
                ELSE
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*) "Internal error in DMUMPS_SOL_S"
                  CALL MUMPS_ABORT()
#else
                  CALL dtrsm('L','L','T','N',NBJ, NRHS_B, ONE,
     &            A(APOSDEB),
     &            LDAJ, RHSCOMP(IPOSINRHSCOMP_PANEL,JBDEB), LRHSCOMP)
#endif
                ENDIF
              ENDIF
#if defined(MUMPS_USE_BLAS2)
            ENDIF
#endif
            IF (.NOT. TWOBYTWO) JJ=BEG_PANEL-1 
       ENDDO 
      ELSE
        IF ( IW(PTRIST(STEP(INODE))+XXLR) .GE. 2
     &        .AND. KEEP(485) .EQ. 1 ) THEN
           IWHDLR = IW(PTRIST(STEP(INODE))+XXF)
          CALL DMUMPS_SOL_BWD_LR_SU (
     &           INODE, IWHDLR, NPIV, NSLAVES, 
     &           LIELL, W, LWC, NRHS_B, PTWCB,
     &           RHSCOMP, LRHSCOMP, NRHS,
     &           IPOSINRHSCOMP, JBDEB, 
     &           MTYPE, KEEP,
     &           INFO(1), INFO(2) )
          IF (INFO(1).LT.0) GOTO 330
        ELSE
       IF ( LIELL .GT. NPIV ) THEN
        IF ( MTYPE .eq. 1 ) THEN
          IST = APOS + int(NPIV,8)
#if defined(MUMPS_USE_BLAS2)
          IF (NRHS_B == 1) THEN
            IF (KEEP(350).EQ.0) THEN
              CALL dgemv( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     &              W(PTWCB+int(NPIV,8)), 1,
     &              ONE,
     &              W(PTWCB), 1 )
            ELSE
              CALL dgemv( 'T', NCB, NPIV, ALPHA, A(IST), LIELL,
     &              W(PTWCB+int(NPIV,8)), 1,
     &              ONE,
#if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP), 1 )
#else
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB), 1 )
#endif
            ENDIF
          ELSE
#endif
            IF (KEEP(350).EQ.0) THEN
              CALL dgemm('T','N', NPIV, NRHS_B, NCB, ALPHA,
     &              A(IST), LIELL,
     &              W(PTWCB+int(NPIV,8)), LIELL, ONE,
     &              W(PTWCB), LIELL)
            ELSE
#if defined(RHSCOMP_BYROWS)
              WRITE(*,*) "Internal error in DMUMPS_SOL_S"
              CALL MUMPS_ABORT()
#else
              CALL dgemm('T','N', NPIV, NRHS_B, NCB, ALPHA,
     &              A(IST), LIELL,
     &              W(PTWCB+int(NPIV,8)), LIELL, ONE,
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB), LRHSCOMP)
#endif
            ENDIF
#if defined(MUMPS_USE_BLAS2)
          ENDIF
#endif
        ELSE
          IF ( KEEP(50) .eq. 0 ) THEN
            IST = APOS + int(NPIV,8) * int(LIELL,8)
          ELSE
            IST = APOS + int(NPIV,8) * int(NPIV,8)
          END IF
#if defined(MUMPS_USE_BLAS2)
          IF ( NRHS_B == 1 ) THEN
              IF (KEEP(350).EQ.0) THEN
                CALL dgemv( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     &                W( PTWCB+int(NPIV,8) ),
     &                1, ONE,
     &                W(PTWCB), 1 )
              ELSE
                CALL dgemv( 'N', NPIV, NCB, ALPHA, A( IST ), NPIV,
     &                W( PTWCB + int(NPIV,8) ),
     &                1, ONE,
#if defined(RHSCOMP_BYROWS)
     &                RHSCOMP(JBDEB,IPOSINRHSCOMP), 1 )
#else
     &                RHSCOMP(IPOSINRHSCOMP,JBDEB), 1 )
#endif
              ENDIF
          ELSE
#endif
                IF (KEEP(350).EQ.0) THEN
                  CALL dgemm( 'N', 'N', NPIV, NRHS_B, NCB, ALPHA,
     &                 A(IST), NPIV, W(PTWCB+int(NPIV,8)),LIELL,
     &                 ONE, W(PTWCB),LIELL)
                ELSE
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*) "Internal error in DMUMPS_SOL_S"
                  CALL MUMPS_ABORT()
#else
                  CALL dgemm( 'N', 'N', NPIV, NRHS_B, NCB, ALPHA,
     &                 A(IST), NPIV, W(PTWCB+int(NPIV,8)),LIELL,
     &                 ONE, RHSCOMP(IPOSINRHSCOMP,JBDEB),LRHSCOMP)
#endif
                ENDIF
#if defined(MUMPS_USE_BLAS2)
          END IF
#endif
        END IF 
       ENDIF  
       IF ( MTYPE .eq. 1 ) THEN
#if defined(MUMPS_USE_BLAS2)
         IF ( NRHS_B == 1 ) THEN
           IF (KEEP(350).EQ.0) THEN
             CALL dtrsv('L', 'T', 'N', NPIV, A(APOS), LIELL,
     &                W(PTWCB), 1)
           ELSE
             CALL dtrsv('L', 'T', 'N', NPIV, A(APOS), LIELL,
#if defined(RHSCOMP_BYROWS)
     &                RHSCOMP(JBDEB,IPOSINRHSCOMP), 1)
#else
     &                RHSCOMP(IPOSINRHSCOMP,JBDEB), 1)
#endif
           ENDIF
         ELSE
#endif
           IF (KEEP(350).EQ.0) THEN
             CALL dtrsm('L','L','T','N', NPIV, NRHS_B, ONE, A(APOS),
     &                LIELL, W(PTWCB), LIELL)
           ELSE
#if defined(RHSCOMP_BYROWS)
                  WRITE(*,*) "Internal error in DMUMPS_SOL_S"
                  CALL MUMPS_ABORT()
#else
             CALL dtrsm('L','L','T','N', NPIV, NRHS_B, ONE, A(APOS),
     &              LIELL, RHSCOMP(IPOSINRHSCOMP,JBDEB), LRHSCOMP)
#endif
           ENDIF
#if defined(MUMPS_USE_BLAS2)
         ENDIF
#endif
       ELSE
         IF ( KEEP(50) .EQ. 0 ) THEN
           LDAJ=LIELL
         ELSE
           LDAJ=NPIV
         ENDIF
#if defined(MUMPS_USE_BLAS2)
         IF ( NRHS_B == 1 ) THEN
            IF (KEEP(350).EQ.0) THEN
              CALL dtrsv('U','N','U', NPIV, A(APOS), LDAJ,
     &              W(PTWCB), 1)
            ELSE
              CALL dtrsv('U','N','U', NPIV, A(APOS), LDAJ,
#if defined(RHSCOMP_BYROWS)
     &              RHSCOMP(JBDEB,IPOSINRHSCOMP), 1)
#else
     &              RHSCOMP(IPOSINRHSCOMP,JBDEB), 1)
#endif
            ENDIF
         ELSE
#endif
            IF (KEEP(350).EQ.0) THEN
              CALL dtrsm('L','U','N','U', NPIV, NRHS_B, ONE, A(APOS),
     &                   LDAJ,W(PTWCB),LIELL)
            ELSE
#if defined(RHSCOMP_BYROWS)
              WRITE(*,*) "Internal error in DMUMPS_SOL_S"
              CALL MUMPS_ABORT()
#else
              CALL dtrsm('L','U','N','U', NPIV, NRHS_B, ONE, A(APOS),
     &                LDAJ, RHSCOMP(IPOSINRHSCOMP,JBDEB), LRHSCOMP)
#endif
            ENDIF
#if defined(MUMPS_USE_BLAS2)
         END IF
#endif
       END IF 
      ENDIF  
      ENDIF 
      IF ( MTYPE .EQ. 1 .AND. KEEP(50).EQ.0) THEN
        J1 = IPOS + LIELL + 1
      ELSE
        J1 = IPOS + 1
      END IF
      IPOSINRHSCOMP =  POSINRHSCOMP_BWD(IW(J1)) 
      IF (KEEP(350).EQ.0) THEN
        IPOSINRHSCOMP_TMP = IPOSINRHSCOMP
        DO 150 I = 1, NPIV
          DO K=JBDEB, JBFIN
#if defined(RHSCOMP_BYROWS)
            RHSCOMP(K,IPOSINRHSCOMP_TMP) = W(PTWCB+I-1+(K-JBDEB)*LIELL)
#else
            RHSCOMP(IPOSINRHSCOMP_TMP, K) = W(PTWCB+I-1+(K-JBDEB)*LIELL)
#endif
          ENDDO
          IPOSINRHSCOMP_TMP =  IPOSINRHSCOMP_TMP + 1
  150   CONTINUE
      ELSE IF (KEEP(350).EQ.1.OR.KEEP(350).EQ.2) THEN
      ELSE
        WRITE(*,*)"Internal error in DMUMPS_SOL_S"
        CALL MUMPS_ABORT()
      ENDIF
  160 CONTINUE
      IF (KEEP(201).GT.0.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
         CALL DMUMPS_FREE_FACTORS_FOR_SOLVE(INODE,PTRFAC,KEEP(28),
     &        A,LA,.TRUE.,IERR)
         IF(IERR.LT.0)THEN
            INFO(1)=IERR
            INFO(2)=0
            GOTO 330
         ENDIF
      ENDIF
      IN = INODE
  170 IN = FILS(IN)
      IF (IN .GT. 0) GOTO 170
      IF (IN .EQ. 0) THEN
        MYLEAFE = MYLEAFE - 1
        IF (.NOT. IN_SUBTREE ) THEN
          IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
          CALL DMUMPS_FREETOPSO(N, KEEP(28), IWCB, LIWW, 
     &    W, LWC,
     &    POSWCB,POSIWCB,PTRICB,PTRACB)
        ENDIF
        IF (MYLEAFE .EQ. 0) THEN
          CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &                     FEUILLE, SLAVEF, KEEP )
          NBFINF = NBFINF - 1
          IF (NBFINF .EQ. 0) GOTO 340
        ENDIF
        GOTO 50
      ENDIF
      IF = -IN
      NBFILS = NE_STEPS(STEP(INODE))
      IF (AM1) THEN
         I = NBFILS
         NBFILS = 0
         DO WHILE (I.GT.0)
            IF ( TO_PROCESS(STEP(IF)) ) NBFILS = NBFILS+1
            IF = FRERE(STEP(IF))
            I = I -1
         ENDDO
         IF (NBFILS.EQ.0) THEN
            NO_CHILDREN = .TRUE.
         ELSE
            NO_CHILDREN = .FALSE.
         ENDIF
         IF = -IN
      ENDIF
      IF (IN_SUBTREE) THEN
         DO I = 1, NBFILS
            IF ( AM1 ) THEN
 1010          IF ( .NOT.TO_PROCESS(STEP(IF)) )  THEN
                  IF = FRERE(STEP(IF))
                  GOTO 1010
               ENDIF
               NO_CHILDREN = .FALSE.
            ENDIF
               IPOOL((IIPOOL-I+1)+NBFILS-I) = IF
               IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
         ENDDO
         IF (AM1 .AND. NO_CHILDREN) THEN
            MYLEAFE = MYLEAFE - 1
            IF (MYLEAFE .EQ. 0) THEN
               CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &              FEUILLE, SLAVEF, KEEP )
               NBFINF = NBFINF - 1
               IF (NBFINF .EQ. 0) GOTO 340
               GOTO 50
            ENDIF
         ENDIF
      ELSE
        DEBUT = .TRUE.
        DO I = 0, SLAVEF - 1
          DEJA_SEND( I ) = .FALSE.
        END DO
        POOL_FIRST_POS=IIPOOL
        DO 190 I = 1, NBFILS
           IF ( AM1 ) THEN
1020      IF ( .NOT.TO_PROCESS(STEP(IF)) ) THEN
                 IF = FRERE(STEP(IF))
                 GOTO 1020
              ENDIF
              NO_CHILDREN = .FALSE.
           ENDIF
          IF (MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),
     &      SLAVEF) .EQ. MYID) THEN
                IPOOL(IIPOOL) = IF
                IIPOOL = IIPOOL + 1
            IF = FRERE(STEP(IF))
          ELSE
            PROCDEST = MUMPS_PROCNODE(PROCNODE_STEPS(STEP(IF)),SLAVEF)
            IF (.not. DEJA_SEND( PROCDEST ))  THEN
 400          CONTINUE
              CALL DMUMPS_BUF_SEND_VCB( NRHS_B, IF, 0, 0, LIELL,
     &          LIELL - KEEP(253),
     &          IW( POSINDICES ), 
     &          W ( PTRACB(STEP( INODE )) ), 
     &          JBDEB, JBFIN, 
     &          RHSCOMP(1, 1), NRHS, LRHSCOMP,
     &          IPOSINRHSCOMP, NPIV,
     &          KEEP, PROCDEST,     NOEUD, COMM, IERR )
              IF ( IERR .EQ. -1 ) THEN
                CALL DMUMPS_BACKSLV_RECV_AND_TREAT(
     &          .FALSE., FLAG,
     &          BUFR, LBUFR, LBUFR_BYTES,
     &          MYID, SLAVEF, COMM,
     &          N, IWCB, LIWW, POSIWCB,
     &          W, LWC, POSWCB,
     &          IIPOOL, NBFINF, PTRICB, PTRACB, INFO,
     &          IPOOL, LPOOL, PANEL_POS, LPANEL_POS,
     &          STEP, FRERE, FILS, PROCNODE_STEPS,
     &          PLEFTW, KEEP,KEEP8, DKEEP,
     &          PTRIST, PTRFAC, IW, LIW, A, LA, W2, MYLEAFE, 
     &          NRHS, MTYPE, 
     &          RHSCOMP, LRHSCOMP, POSINRHSCOMP_BWD
     &                , TO_PROCESS, SIZE_TO_PROCESS
     &                , FROM_PP
     &                )
                IF ( INFO( 1 ) .LT. 0 ) GOTO 340
                GOTO 400
             ELSE IF ( IERR .EQ. -2 ) THEN
                INFO( 1 ) = -17
                INFO( 2 ) = NRHS_B * LIELL * KEEP(35) + 4 * KEEP(34)
                GOTO 330
             ELSE IF ( IERR .EQ. -3 ) THEN
                INFO( 1 ) = -20
                INFO( 2 ) = NRHS_B * LIELL * KEEP(35) + 4 * KEEP(34)
                GOTO 330
              END IF
              DEJA_SEND( PROCDEST ) = .TRUE.
            END IF
            IF = FRERE(STEP(IF))
          ENDIF
  190   CONTINUE
        IF (AM1 .AND. NO_CHILDREN) THEN
           MYLEAFE = MYLEAFE - 1
           IF (MYLEAFE .EQ. 0) THEN
              CALL DMUMPS_MCAST2( DUMMY, 1, MPI_INTEGER, MYID, COMM,
     &             FEUILLE, SLAVEF, KEEP )
              NBFINF = NBFINF - 1
              IF (NBFINF .EQ. 0) GOTO 340
              GOTO 50
           ENDIF
        ENDIF
           DO I=1,(IIPOOL-POOL_FIRST_POS)/2
              TMP=IPOOL(POOL_FIRST_POS+I-1)
              IPOOL(POOL_FIRST_POS+I-1)=IPOOL(IIPOOL-I)
              IPOOL(IIPOOL-I)=TMP
           ENDDO 
        IWCB(PTRICB(STEP(INODE))+1) = IWCB(PTRICB(STEP(INODE))+1)-1
        CALL DMUMPS_FREETOPSO(N, KEEP(28), IWCB, LIWW, 
     &     W, LWC,
     &     POSWCB,POSIWCB,PTRICB,PTRACB)
      ENDIF
      GOTO 50
  330 CONTINUE
      CALL DMUMPS_MCAST2(DUMMY, 1, MPI_INTEGER, MYID, COMM, TERREUR,
     & SLAVEF, KEEP)
  340 CONTINUE
      CALL DMUMPS_CLEAN_PENDING(INFO(1), KEEP, BUFR, LBUFR,LBUFR_BYTES,
     &     COMM, DUMMY(1),  
     &     SLAVEF, .TRUE., .FALSE.) 
      IF (ALLOCATED(DEJA_SEND)) DEALLOCATE(DEJA_SEND)
      RETURN
      END SUBROUTINE DMUMPS_SOL_S
