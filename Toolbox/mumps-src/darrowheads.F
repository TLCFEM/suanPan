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
      SUBROUTINE DMUMPS_ANA_DIST_ARROWHEADS( MYID, SLAVEF, N,
     &           PROCNODE, STEP, PTRAIW, PTRARW, ISTEP_TO_INIV2,
     &           I_AM_CAND,
     &           KEEP, KEEP8, ICNTL, id )
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
      TYPE (DMUMPS_STRUC) :: id
      INTEGER MYID, N, SLAVEF
      INTEGER KEEP( 500 ), ICNTL( 60 )
      INTEGER(8) KEEP8(150)
      INTEGER PROCNODE( KEEP(28) ), STEP( N )
      INTEGER(8), INTENT(INOUT) :: PTRAIW( N ), PTRARW( N )
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      LOGICAL I_AM_CAND(max(1,KEEP(56)))
      LOGICAL I_AM_SLAVE
      LOGICAL I_AM_CAND_LOC
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE, MUMPS_TYPESPLIT
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE, MUMPS_TYPESPLIT
      INTEGER ISTEP, I, NCOL, NROW, allocok
      INTEGER TYPE_PARALL, ITYPE, IRANK, INIV2, TYPESPLIT 
      LOGICAL T4_MASTER_CONCERNED
      INTEGER(8) :: IPTRI, IPTRR
      TYPE_PARALL = KEEP(46)
      I_AM_SLAVE = (KEEP(46).EQ.1 .OR. MYID.NE.0)
      KEEP8(26) = 0_8
      KEEP8(27) = 0_8
      DO I = 1, N
        ISTEP=abs(STEP(I))
        ITYPE = MUMPS_TYPENODE( PROCNODE(ISTEP), SLAVEF )
        IRANK = MUMPS_PROCNODE( PROCNODE(ISTEP), SLAVEF )
        I_AM_CAND_LOC = .FALSE.
        TYPESPLIT = MUMPS_TYPESPLIT ( PROCNODE(ISTEP), SLAVEF )
        T4_MASTER_CONCERNED = .FALSE.
        IF (ITYPE.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF (I_AM_SLAVE)  THEN 
           I_AM_CAND_LOC = I_AM_CAND(INIV2)
          IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
           IF ( TYPE_PARALL .eq. 0 ) THEN
            T4_MASTER_CONCERNED = 
     &     ( id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
     &       .EQ.MYID-1 )
           ELSE
            T4_MASTER_CONCERNED = 
     &     ( id%CANDIDATES (id%CANDIDATES(SLAVEF+1, INIV2)+1,INIV2 ) 
     &       .EQ.MYID )
           ENDIF
          ENDIF
         ENDIF
        ENDIF
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK = IRANK + 1
        END IF
        IF (
     &       ( (ITYPE .EQ. 1.OR.ITYPE.EQ.2) .AND.
     &            IRANK .EQ. MYID ) 
     &       .OR.
     &       ( T4_MASTER_CONCERNED ) 
     &     ) THEN
          KEEP8(26) = KEEP8(26) + 1_8 + PTRAIW(I)+PTRARW(I)
          KEEP8(27) = KEEP8(27) + 3_8 + PTRAIW(I)+PTRARW(I)
        ELSE IF ( ITYPE .EQ. 3 ) THEN
        ELSE IF ( ITYPE .EQ. 2 .AND. I_AM_CAND_LOC ) THEN
           PTRARW( I ) = 0_8
           KEEP8(26) = KEEP8(26) + 1_8 + PTRAIW(I)+PTRARW(I)
           KEEP8(27) = KEEP8(27) + 3_8 + PTRAIW(I)+PTRARW(I)
        END IF
      END DO
      IF ( associated( id%INTARR ) ) THEN
        DEALLOCATE( id%INTARR )
        NULLIFY( id%INTARR )
      END IF
      IF ( KEEP8(27) > 0 ) THEN
      ALLOCATE( id%INTARR( KEEP8(27) ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        id%INFO(1) = -7
        CALL  MUMPS_SET_IERROR(KEEP8(27),id%INFO(2))
        RETURN
      END IF
      ELSE
      ALLOCATE( id%INTARR( 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        id%INFO(1) = -7
        id%INFO(2) = 1
        RETURN
      END IF
      END IF
      IPTRI = 1_8
      IPTRR = 1_8
      DO I = 1, N
        ISTEP = abs(STEP(I))
        ITYPE = MUMPS_TYPENODE( PROCNODE(ISTEP), SLAVEF )
        IRANK = MUMPS_PROCNODE( PROCNODE(ISTEP), SLAVEF )
        TYPESPLIT = MUMPS_TYPESPLIT ( PROCNODE(ISTEP), SLAVEF )
        I_AM_CAND_LOC = .FALSE.
        T4_MASTER_CONCERNED = .FALSE.
        IF (ITYPE.EQ.2) THEN
          INIV2         = ISTEP_TO_INIV2(ISTEP)
          IF (I_AM_SLAVE)  THEN
           I_AM_CAND_LOC = I_AM_CAND(INIV2)
           IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
            IF ( TYPE_PARALL .eq. 0 ) THEN
             T4_MASTER_CONCERNED = 
     &       (id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
     &         .EQ.MYID-1 )
            ELSE
              T4_MASTER_CONCERNED = 
     &        (id%CANDIDATES (id%CANDIDATES(SLAVEF+1,INIV2)+1,INIV2) 
     &         .EQ.MYID )
            ENDIF
           ENDIF
          ENDIF
        ENDIF
        IF ( TYPE_PARALL .eq. 0 ) THEN
          IRANK =IRANK + 1
        END IF
        IF (
     &      ( ITYPE .eq. 2 .and.
     &        IRANK .eq. MYID )
     & .or.
     &      ( ITYPE .eq. 1 .and.
     &        IRANK .eq. MYID )
     & .or.
     &      ( T4_MASTER_CONCERNED )
     &     )  THEN
          NCOL = int(PTRAIW( I ))
          NROW = int(PTRARW( I ))
          id%INTARR( IPTRI     ) = NCOL
          id%INTARR( IPTRI + 1 ) = -NROW
          id%INTARR( IPTRI + 2 ) = I
          PTRAIW( I ) = IPTRI
          PTRARW( I ) = IPTRR
          IPTRI = IPTRI + int(NCOL + NROW + 3,8)
          IPTRR = IPTRR + int(NCOL + NROW + 1,8)
        ELSE IF ( ITYPE .eq. 2  .AND. I_AM_CAND_LOC ) THEN
           NCOL = int(PTRAIW( I ))
           NROW = 0
           id%INTARR( IPTRI     ) = NCOL
           id%INTARR( IPTRI + 1 ) = -NROW
           id%INTARR( IPTRI + 2 ) = I
           PTRAIW( I ) = IPTRI
           PTRARW( I ) = IPTRR
           IPTRI = IPTRI + int(NCOL + NROW + 3, 8)
           IPTRR = IPTRR + int(NCOL + NROW + 1, 8)
        ELSE
          PTRAIW(I) = 0_8
          PTRARW(I) = 0_8
        END IF
      END DO
      IF ( IPTRI - 1_8 .NE. KEEP8(27) ) THEN
        WRITE(*,*) 'Error 1 in ana_arrowheads',  
     &      ' IPTRI - 1, KEEP8(27)=', IPTRI - 1, KEEP8(27)
        CALL MUMPS_ABORT()
      END IF
      IF ( IPTRR - 1_8 .NE. KEEP8(26) ) THEN
        WRITE(*,*) 'Error 2 in ana_arrowheads'
        CALL MUMPS_ABORT()
      END IF
      RETURN
      END SUBROUTINE DMUMPS_ANA_DIST_ARROWHEADS
      SUBROUTINE DMUMPS_FACTO_SEND_ARROWHEADS(N, NZ, ASPK, 
     &   IRN, ICN, PERM,
     &   LSCAL,COLSCA,ROWSCA,
     &   MYID, SLAVEF, PROCNODE_STEPS, NBRECORDS,
     &   LP, COMM, root, KEEP, KEEP8, FILS, RG2L,
     &   INTARR, LINTARR, DBLARR, LDBLARR, PTRAIW, PTRARW, FRERE_STEPS,
     &   STEP, A, LA, ISTEP_TO_INIV2, I_AM_CAND, CANDIDATES )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER    :: N, COMM, NBRECORDS
      INTEGER(8), INTENT(IN) :: NZ
      INTEGER KEEP( 500 )
      INTEGER(8) KEEP8(150)
      DOUBLE PRECISION ASPK(NZ)
      DOUBLE PRECISION COLSCA(*), ROWSCA(*)
      INTEGER IRN(NZ), ICN(NZ) 
      INTEGER PERM(N), PROCNODE_STEPS(KEEP(28))
      INTEGER RG2L( N ), FILS( N )
      INTEGER ISTEP_TO_INIV2(KEEP(71))
      LOGICAL I_AM_CAND(max(1,KEEP(56)))
      INTEGER LP, SLAVEF, MYID
      INTEGER CANDIDATES(SLAVEF+1, max(1,KEEP(56)))
      LOGICAL LSCAL
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER(8), INTENT(IN)    :: LA
      INTEGER(8), INTENT(INOUT) :: PTRAIW( N ), PTRARW( N ) 
      INTEGER    :: FRERE_STEPS( KEEP(28) )
      INTEGER    :: STEP(N)
      INTEGER(8) :: LINTARR, LDBLARR
      INTEGER    :: INTARR( LINTARR )
      DOUBLE PRECISION    :: DBLARR( LDBLARR )
      DOUBLE PRECISION    :: A( LA )
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: BUFI
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BUFR
      INTEGER MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      EXTERNAL MUMPS_PROCNODE, MUMPS_TYPENODE, numroc, 
     &        MUMPS_TYPESPLIT
      DOUBLE PRECISION VAL
      INTEGER IOLD,JOLD,INEW,JNEW,ISEND,JSEND,DEST,I,IARR
      INTEGER IPOSROOT, JPOSROOT
      INTEGER IROW_GRID, JCOL_GRID
      INTEGER INODE, ISTEP
      INTEGER NBUFS
      INTEGER ARROW_ROOT, TAILLE
      INTEGER LOCAL_M, LOCAL_N
      INTEGER(8) :: PTR_ROOT
      INTEGER TYPENODE_TMP, MASTER_NODE
      LOGICAL I_AM_CAND_LOC, I_AM_SLAVE
      INTEGER JARR, ILOCROOT, JLOCROOT
      INTEGER allocok, INIV2, TYPESPLIT, T4MASTER
      INTEGER(8) ::  I1, IA, IIW, IS1, IS, IAS, ISHIFT, K
      INTEGER NCAND
      LOGICAL T4_MASTER_CONCERNED
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      ARROW_ROOT = 0
      I_AM_SLAVE=(MYID.NE.0.OR.KEEP(46).EQ.1)
      IF ( KEEP(46) .eq. 0 ) THEN
        NBUFS = SLAVEF
      ELSE
        NBUFS = SLAVEF - 1
        ALLOCATE( IW4( N, 2 ), stat = allocok )
        IF ( allocok .GT. 0 ) THEN
          WRITE(*,*) 'Error allocating IW4'
          CALL MUMPS_ABORT()
        END IF
        DO I = 1, N
          I1 = PTRAIW( I )
          IA = PTRARW( I )
          IF ( IA .GT. 0 ) THEN
            DBLARR( IA ) = ZERO
            IW4( I, 1 ) = INTARR( I1 )       
            IW4( I, 2 ) = -INTARR( I1 + 1 )  
            INTARR( I1 + 2 ) = I
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
              root%SCHUR_POINTER(int(I-1,8)*int(root%SCHUR_LLD,8)+1_8:
     &        int(I-1,8)*int(root%SCHUR_LLD,8)+int(root%SCHUR_MLOC,8))=
     &        ZERO
            ENDDO
          ENDIF
        END IF
      END IF
      IF (NBUFS.GT.0) THEN
       ALLOCATE( BUFI(NBRECORDS*2+1,NBUFS),stat=allocok )
       IF ( allocok .GT. 0 ) THEN
        WRITE(*,*) 'Error allocating BUFI'
        CALL MUMPS_ABORT()
       END IF
       ALLOCATE( BUFR( NBRECORDS, NBUFS ), stat=allocok )
       IF ( allocok .GT. 0 ) THEN
         WRITE(*,*) 'Error allocating BUFR'
         CALL MUMPS_ABORT()
       END IF
       DO I = 1, NBUFS
        BUFI( 1, I ) = 0
       ENDDO
      ENDIF
      INODE = KEEP(38)
      I     = 1
      DO WHILE ( INODE .GT. 0 )
        RG2L( INODE ) = I
        INODE = FILS( INODE )
        I = I + 1
      END DO
      DO 120 K=1,NZ
        IOLD = IRN(K)
        JOLD = ICN(K)
        IF ( (IOLD.GT.N).OR.(JOLD.GT.N).OR.(IOLD.LT.1)
     &                 .OR.(JOLD.LT.1) ) THEN
           GOTO 120
        END IF
        IF (LSCAL) THEN
          VAL = ASPK(K)*ROWSCA(IOLD)*COLSCA(JOLD)
        ELSE
          VAL = ASPK(K)
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
        IARR  = abs( ISEND )
        ISTEP = abs( STEP(IARR) )
        TYPENODE_TMP = MUMPS_TYPENODE( PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF ) 
        MASTER_NODE  = MUMPS_PROCNODE( PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF )
        I_AM_CAND_LOC          = .FALSE.
        TYPESPLIT  = MUMPS_TYPESPLIT(  PROCNODE_STEPS(ISTEP),
     &                                 SLAVEF )
        T4_MASTER_CONCERNED = .FALSE.
        T4MASTER               = -9999
        IF (TYPENODE_TMP.EQ.2) THEN
         INIV2         = ISTEP_TO_INIV2(ISTEP)
         IF (I_AM_SLAVE)  I_AM_CAND_LOC = I_AM_CAND(INIV2)
         IF ( (TYPESPLIT.EQ.5).OR.(TYPESPLIT.EQ.6)) THEN
          T4_MASTER_CONCERNED = .TRUE.
          T4MASTER=CANDIDATES(CANDIDATES(SLAVEF+1,INIV2)+1,INIV2)
          IF ( KEEP(46) .eq. 0 ) THEN
           T4MASTER=T4MASTER+1
          ENDIF
         ENDIF
        ENDIF
        IF ( TYPENODE_TMP .EQ. 1 ) THEN
          IF ( KEEP(46) .eq. 0 ) THEN
            DEST = MASTER_NODE + 1
          ELSE
            DEST = MASTER_NODE
          END IF
        ELSE IF ( TYPENODE_TMP .EQ. 2 ) THEN
          IF ( ISEND .LT. 0  ) THEN
            DEST = -1
          ELSE
            IF ( KEEP( 46 ) .eq. 0 ) THEN
              DEST = MASTER_NODE + 1
            ELSE 
              DEST = MASTER_NODE
            END IF
          END IF
        ELSE
          IF ( ISEND .LT. 0 ) THEN
            IPOSROOT = RG2L(JSEND)
            JPOSROOT = RG2L(IARR)
          ELSE
            IPOSROOT = RG2L( IARR )
            JPOSROOT = RG2L( JSEND )
          END IF
          IROW_GRID = mod( ( IPOSROOT - 1 )/root%MBLOCK, root%NPROW )
          JCOL_GRID = mod( ( JPOSROOT - 1 )/root%NBLOCK, root%NPCOL )
          IF ( KEEP( 46 ) .eq. 0 ) THEN
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID + 1
          ELSE
            DEST = IROW_GRID * root%NPCOL + JCOL_GRID
          END IF
        END IF
        IF ( DEST .eq. 0 .or.
     &      ( DEST .eq. -1 .and. KEEP( 46 ) .eq. 1 .AND.
     &       ( I_AM_CAND_LOC .OR. MASTER_NODE .EQ. 0 ) )
     &     .or. 
     &      ( T4MASTER.EQ.0 )
     &     ) THEN
          IARR = ISEND  
          JARR = JSEND
          IF ( TYPENODE_TMP .eq. 3 ) THEN
            ARROW_ROOT = ARROW_ROOT + 1
            IF ( IROW_GRID .EQ. root%MYROW .AND.
     &         JCOL_GRID .EQ. root%MYCOL ) THEN
              ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
              JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
             IF (KEEP(60)==0) THEN
               A( PTR_ROOT
     &           + int(JLOCROOT - 1,8) * int(LOCAL_M,8) 
     &           + int(ILOCROOT - 1,8) )
     &         =  A( PTR_ROOT
     &           + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &           + int(ILOCROOT - 1,8) )
     &         + VAL
             ELSE
               root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                           * int(root%SCHUR_LLD,8)
     &                           + int(ILOCROOT,8) )
     &          = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                           *    int(root%SCHUR_LLD,8)
     &                           +    int(ILOCROOT,8))
     &          + VAL
             ENDIF
            ELSE
              WRITE(*,*) MYID,':INTERNAL Error: root arrowhead '
              WRITE(*,*) MYID,':is not belonging to me. IARR,JARR='
     &        ,IARR,JARR
              CALL MUMPS_ABORT()
            END IF
          ELSE IF ( IARR .GE. 0 ) THEN
            IF ( IARR .eq. JARR ) THEN
              IA = PTRARW( IARR )
              DBLARR( IA ) = DBLARR( IA ) + VAL
            ELSE
              IS1 =  PTRAIW(IARR)
              ISHIFT      = int(INTARR(IS1) + IW4(IARR,2),8)
              IW4(IARR,2) = IW4(IARR,2) - 1
              IIW         = IS1 + ISHIFT + 2_8
              INTARR(IIW)     = JARR
              IS          = PTRARW(IARR)
              IAS         = IS + ISHIFT
              DBLARR(IAS) = VAL
            END IF
          ELSE
            IARR = -IARR
            ISHIFT      = int(PTRAIW(IARR)+IW4(IARR,1)+2,8)
            INTARR(ISHIFT)  = JARR
            IAS         = PTRARW(IARR)+int(IW4(IARR,1),8)
            IW4(IARR,1) = IW4(IARR,1) - 1
            DBLARR(IAS)      = VAL
            IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0 .OR. .TRUE.)
     &           .AND.  IW4(IARR,1) .EQ. 0 .AND.
     &           STEP( IARR) > 0 ) THEN
              IF (MUMPS_PROCNODE( PROCNODE_STEPS(abs(STEP(IARR))),
     &                            SLAVEF ) == MYID) THEN
                TAILLE = INTARR( PTRAIW(IARR) )
                CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &             INTARR( PTRAIW(IARR) + 3 ),
     &             DBLARR( PTRARW(IARR) + 1 ),
     &             TAILLE, 1, TAILLE )
              END IF
            END IF
          ENDIF
        END IF
        IF ( DEST.EQ. -1 ) THEN
         INIV2 = ISTEP_TO_INIV2(ISTEP)
         NCAND = CANDIDATES(SLAVEF+1,INIV2)
         IF (KEEP(79).GT.0) THEN
          DO I=1, SLAVEF
           DEST=CANDIDATES(I,INIV2)
           IF (KEEP(46).EQ.0.AND.(DEST.GE.0)) DEST=DEST+1
           IF (DEST.LT.0) EXIT 
           IF (I.EQ.NCAND+1) CYCLE 
           IF (DEST.NE.0)
     &     CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &     DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
          ENDDO
         ELSE
          DO I=1, NCAND
           DEST=CANDIDATES(I,INIV2)
           IF (KEEP(46).EQ.0) DEST=DEST+1
           IF (DEST.NE.0)
     &     CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &     DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
          ENDDO
         ENDIF
         DEST = MASTER_NODE
         IF (KEEP(46).EQ.0) DEST=DEST+1
         IF ( DEST .NE. 0 ) THEN
           CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &     DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
         ENDIF
         IF ((T4_MASTER_CONCERNED).AND.(T4MASTER.GT.0)) THEN 
          CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &     T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &     LP, COMM, KEEP(46))
         ENDIF 
        ELSE IF ( DEST .GT. 0 ) THEN
         CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &    DEST, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
         IF ( T4MASTER.GT.0 ) THEN
          CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &    T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
         ENDIF
        ELSE IF ( T4MASTER.GT.0 ) THEN
         CALL DMUMPS_ARROW_FILL_SEND_BUF( ISEND, JSEND, VAL,
     &    T4MASTER, BUFI, BUFR, NBRECORDS, NBUFS, 
     &    LP, COMM, KEEP(46))
        END IF
  120 CONTINUE
      KEEP(49) = ARROW_ROOT
      IF (NBUFS.GT.0) THEN
       CALL DMUMPS_ARROW_FINISH_SEND_BUF(
     &   BUFI, BUFR, NBRECORDS, NBUFS,
     &   LP, COMM, KEEP( 46 ) )
      ENDIF
      IF ( KEEP( 46 ) .NE. 0 ) DEALLOCATE( IW4 )
      IF (NBUFS.GT.0) THEN
        DEALLOCATE( BUFI )
        DEALLOCATE( BUFR )
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FACTO_SEND_ARROWHEADS
      SUBROUTINE DMUMPS_ARROW_FILL_SEND_BUF(ISEND, JSEND, VAL,
     &   DEST, BUFI, BUFR, NBRECORDS, NBUFS, LP, COMM,
     &   TYPE_PARALL )
      IMPLICIT NONE
      INTEGER ISEND, JSEND, DEST, NBUFS, NBRECORDS, TYPE_PARALL
      INTEGER BUFI( NBRECORDS * 2 + 1, NBUFS )
      DOUBLE PRECISION BUFR( NBRECORDS, NBUFS )
      INTEGER COMM
      INTEGER LP
      DOUBLE PRECISION VAL
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER IERR
      INTEGER TAILLE_SENDI, TAILLE_SENDR, IREQ
         IF (BUFI(1,DEST)+1.GT.NBRECORDS) THEN
          TAILLE_SENDI = BUFI(1,DEST) * 2 + 1
          TAILLE_SENDR = BUFI(1,DEST)
          CALL MPI_SEND(BUFI(1,DEST),TAILLE_SENDI,
     &                   MPI_INTEGER,
     &                   DEST, ARROWHEAD, COMM, IERR )
          CALL MPI_SEND( BUFR(1,DEST), TAILLE_SENDR,
     &                   MPI_DOUBLE_PRECISION, DEST,
     &                   ARROWHEAD, COMM, IERR )
          BUFI(1,DEST) = 0
         ENDIF
         IREQ = BUFI(1,DEST) + 1
         BUFI(1,DEST) = IREQ
         BUFI( IREQ * 2, DEST )     = ISEND
         BUFI( IREQ * 2 + 1, DEST ) = JSEND
         BUFR( IREQ, DEST )         = VAL
      RETURN
      END SUBROUTINE DMUMPS_ARROW_FILL_SEND_BUF
      SUBROUTINE DMUMPS_ARROW_FINISH_SEND_BUF(
     &   BUFI, BUFR, NBRECORDS, NBUFS, LP, COMM,
     &   TYPE_PARALL )
      IMPLICIT NONE
      INTEGER NBUFS, NBRECORDS, TYPE_PARALL
      INTEGER BUFI( NBRECORDS * 2 + 1, NBUFS )
      DOUBLE PRECISION BUFR( NBRECORDS, NBUFS )
      INTEGER COMM
      INTEGER LP
      INTEGER ISLAVE, TAILLE_SENDI, TAILLE_SENDR, IERR
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
        DO ISLAVE = 1,NBUFS 
          TAILLE_SENDI = BUFI(1,ISLAVE) * 2 + 1
          TAILLE_SENDR = BUFI(1,ISLAVE)
          BUFI(1,ISLAVE) = - BUFI(1,ISLAVE)
          CALL MPI_SEND(BUFI(1,ISLAVE),TAILLE_SENDI,
     &                   MPI_INTEGER,
     &                   ISLAVE, ARROWHEAD, COMM, IERR )
          IF ( TAILLE_SENDR .NE. 0 ) THEN
            CALL MPI_SEND( BUFR(1,ISLAVE), TAILLE_SENDR,
     &                     MPI_DOUBLE_PRECISION, ISLAVE,
     &                     ARROWHEAD, COMM, IERR )
          END IF
        ENDDO
      RETURN
      END SUBROUTINE DMUMPS_ARROW_FINISH_SEND_BUF
      RECURSIVE SUBROUTINE DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM, 
     &            INTLIST, DBLLIST, TAILLE, LO, HI )
      IMPLICIT NONE
      INTEGER N, TAILLE
      INTEGER PERM( N )
      INTEGER INTLIST( TAILLE )
      DOUBLE PRECISION DBLLIST( TAILLE )
      INTEGER LO, HI
      INTEGER I,J
      INTEGER ISWAP, PIVOT
      DOUBLE PRECISION dswap
      I = LO
      J = HI
      PIVOT = PERM(INTLIST((I+J)/2))
 10   IF (PERM(INTLIST(I)) < PIVOT) THEN
        I=I+1
        GOTO 10
      ENDIF
 20   IF (PERM(INTLIST(J)) > PIVOT) THEN
        J=J-1
        GOTO 20
      ENDIF
      IF (I < J) THEN
        ISWAP = INTLIST(I)
        INTLIST(I) = INTLIST(J)
        INTLIST(J)=ISWAP
        dswap = DBLLIST(I)
        DBLLIST(I) = DBLLIST(J)
        DBLLIST(J) = dswap
      ENDIF
      IF ( I <= J) THEN
        I = I+1
        J = J-1
      ENDIF
      IF ( I <= J ) GOTO 10
      IF ( LO < J ) CALL DMUMPS_QUICK_SORT_ARROWHEADS(N, PERM,
     &              INTLIST, DBLLIST, TAILLE, LO, J)
      IF ( I < HI ) CALL DMUMPS_QUICK_SORT_ARROWHEADS(N, PERM,
     &              INTLIST, DBLLIST, TAILLE, I, HI)
      RETURN
      END SUBROUTINE DMUMPS_QUICK_SORT_ARROWHEADS
      SUBROUTINE DMUMPS_FACTO_RECV_ARROWHD2(  N,
     &    DBLARR, LDBLARR, INTARR, LINTARR, PTRAIW, PTRARW, 
     &    KEEP, KEEP8, MYID,  COMM, NBRECORDS,
     &    A, LA, root,
     &    PROCNODE_STEPS,
     &    SLAVEF, PERM, FRERE_STEPS, STEP, INFO1, INFO2
     &   )
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INTEGER N, MYID, COMM
      INTEGER(8), INTENT(IN) :: LDBLARR, LINTARR
      INTEGER INTARR(LINTARR) 
      INTEGER(8), INTENT(IN) :: PTRAIW(N), PTRARW(N) 
      INTEGER   KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), intent(IN) :: LA
      INTEGER PROCNODE_STEPS( KEEP(28) ), PERM( N )
      INTEGER SLAVEF, NBRECORDS
      DOUBLE PRECISION A( LA )
      INTEGER INFO1, INFO2
      DOUBLE PRECISION DBLARR(LDBLARR)
      INTEGER FRERE_STEPS( KEEP(28) ), STEP(N)
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER, POINTER, DIMENSION(:) :: BUFI
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: BUFR
      INTEGER, POINTER, DIMENSION(:,:) :: IW4
      LOGICAL FINI 
      INTEGER IREC, NB_REC, IARR, JARR, I, allocok
      INTEGER(8) :: I18, IA8, IS18, IIW8, IS8, IAS8
      INTEGER ISHIFT
      INTEGER LOCAL_M, LOCAL_N, ILOCROOT, JLOCROOT, 
     &        IPOSROOT, JPOSROOT, TAILLE,
     &        IPROC
      INTEGER(8) :: PTR_ROOT
      INTEGER ARROW_ROOT, TYPE_PARALL
      INTEGER MUMPS_TYPENODE, MUMPS_PROCNODE
      EXTERNAL MUMPS_TYPENODE, MUMPS_PROCNODE
      DOUBLE PRECISION VAL
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER MASTER
      PARAMETER(MASTER=0)
      INTEGER :: IERR
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER numroc
      EXTERNAL numroc
      TYPE_PARALL = KEEP(46)
      ARROW_ROOT=0
      ALLOCATE( BUFI( NBRECORDS * 2 + 1 ), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS * 2 + 1
        WRITE(*,*) MYID,': Could not allocate BUFI: goto 500'
        GOTO 500
      END IF
      ALLOCATE( BUFR( NBRECORDS )        , stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = NBRECORDS
        WRITE(*,*) MYID,': Could not allocate BUFR: goto 500'
        GOTO 500
      END IF
      ALLOCATE( IW4(N,2), stat = allocok )
      IF ( allocok .GT. 0 ) THEN
        INFO1 = -13
        INFO2 = 2 * N
        WRITE(*,*) MYID,': Could not allocate IW4: goto 500'
        GOTO 500
      END IF
      IF ( KEEP(38).NE.0) THEN
        IF (KEEP(60)==0) THEN
         LOCAL_M = numroc( root%ROOT_SIZE, root%MBLOCK,
     &             root%MYROW, 0, root%NPROW )
         LOCAL_M = max( 1, LOCAL_M )
         LOCAL_N = numroc( root%ROOT_SIZE, root%NBLOCK,
     &             root%MYCOL, 0, root%NPCOL )
         PTR_ROOT = LA - int(LOCAL_M,8) * int(LOCAL_N,8) + 1_8
         IF ( PTR_ROOT .LE. LA ) THEN
           A( PTR_ROOT:LA ) = ZERO
         END IF
        ELSE
         DO I=1, root%SCHUR_NLOC
           root%SCHUR_POINTER((I-1)*root%SCHUR_LLD+1:
     &     (I-1)*root%SCHUR_LLD+root%SCHUR_MLOC)=ZERO
         ENDDO
        ENDIF
      END IF
      FINI = .FALSE.
      DO I=1,N
       I18 = PTRAIW(I)
       IA8 = PTRARW(I)
       IF (IA8.GT.0_8) THEN
        DBLARR(IA8) = ZERO
        IW4(I,1) = INTARR(I18)
        IW4(I,2) = -INTARR(I18+1_8)
        INTARR(I18+2)=I
       ENDIF
      ENDDO
      DO WHILE (.NOT.FINI) 
       CALL MPI_RECV( BUFI(1), 2*NBRECORDS+1, 
     &                MPI_INTEGER, MASTER, 
     &                ARROWHEAD,
     &                COMM, STATUS, IERR )
       NB_REC = BUFI(1)
       IF (NB_REC.LE.0) THEN
         FINI = .TRUE.
         NB_REC = -NB_REC 
       ENDIF
       IF (NB_REC.EQ.0) EXIT
       CALL MPI_RECV( BUFR(1), NBRECORDS, MPI_DOUBLE_PRECISION,
     &                  MASTER, ARROWHEAD,
     &                COMM, STATUS, IERR )
       DO IREC=1, NB_REC
        IARR = BUFI( IREC * 2 )
        JARR = BUFI( IREC * 2 + 1 )
        VAL  = BUFR( IREC )
        IF ( MUMPS_TYPENODE( PROCNODE_STEPS(abs(STEP(abs(IARR)))),
     &                       SLAVEF ) .eq. 3 ) THEN
          ARROW_ROOT = ARROW_ROOT + 1
          IF ( IARR .GT. 0 ) THEN
            IPOSROOT = root%RG2L_ROW( IARR )
            JPOSROOT = root%RG2L_COL( JARR )
          ELSE
            IPOSROOT = root%RG2L_ROW( JARR )
            JPOSROOT = root%RG2L_COL( -IARR )
          END IF
            ILOCROOT = root%MBLOCK * ( ( IPOSROOT - 1 ) /
     &                 ( root%MBLOCK * root%NPROW ) )
     &               + mod( IPOSROOT - 1, root%MBLOCK ) + 1
            JLOCROOT = root%NBLOCK * ( ( JPOSROOT - 1 ) /
     &                 ( root%NBLOCK * root%NPCOL ) )
     &               + mod( JPOSROOT - 1, root%NBLOCK ) + 1
           IF (KEEP(60)==0) THEN
             A( PTR_ROOT + int(JLOCROOT - 1,8) * int(LOCAL_M,8)
     &                   + int(ILOCROOT - 1,8) )
     &       =  A( PTR_ROOT + int(JLOCROOT - 1,8)
     &                      * int(LOCAL_M,8)
     &                      + int(ILOCROOT - 1,8))
     &        + VAL
           ELSE
             root%SCHUR_POINTER( int(JLOCROOT-1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8) )
     &       = root%SCHUR_POINTER( int(JLOCROOT - 1,8)
     &                         * int(root%SCHUR_LLD,8)
     &                         + int(ILOCROOT,8))
     &       + VAL
           ENDIF
        ELSE IF (IARR.GE.0) THEN
         IF (IARR.EQ.JARR) THEN
          IA8 = PTRARW(IARR)
          DBLARR(IA8) = DBLARR(IA8) + VAL
         ELSE
          IS18 =  PTRAIW(IARR)
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
           IF ( (KEEP(50) .NE. 0 .OR. KEEP(234).NE.0 .OR. .TRUE.)
     &          .AND.  IW4(IARR,1) .EQ. 0 
     &          .AND. STEP(IARR) > 0 ) THEN
              IPROC = MUMPS_PROCNODE( PROCNODE_STEPS(abs(STEP(IARR))),
     &                                SLAVEF )
              IF ( TYPE_PARALL .eq. 0 ) THEN
                IPROC = IPROC + 1
              END IF 
              IF (IPROC .EQ. MYID) THEN
                TAILLE = INTARR( PTRAIW(IARR) )
                CALL DMUMPS_QUICK_SORT_ARROWHEADS( N, PERM,
     &            INTARR( PTRAIW(IARR) + 3 ),
     &            DBLARR( PTRARW(IARR) + 1 ),
     &            TAILLE, 1, TAILLE )
              END IF
           END IF
        ENDIF
       ENDDO
      END DO
      DEALLOCATE( BUFI )
      DEALLOCATE( BUFR )
      DEALLOCATE( IW4 )
 500  CONTINUE
      KEEP(49) = ARROW_ROOT
      RETURN 
      END SUBROUTINE DMUMPS_FACTO_RECV_ARROWHD2
