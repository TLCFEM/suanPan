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
      SUBROUTINE DMUMPS_FAC_B(N, NSTEPS,
     & A, LA, IW, LIW, SYM_PERM, NA, LNA,
     & NE_STEPS, NFSIZ, FILS,
     & STEP, FRERE, DAD, CAND, 
     & ISTEP_TO_INIV2, TAB_POS_IN_PERE, 
     & PTRAR, LDPTRAR,
     & PTRIST, PTLUST_S, PTRFAC, IW1, IW2, ITLOC, RHS_MUMPS,
     & POOL, LPOOL,  
     & CNTL1, ICNTL, INFO, RINFO, KEEP,KEEP8,PROCNODE_STEPS,
     & SLAVEF,
     & COMM_NODES, MYID, MYID_NODES,
     & BUFR,LBUFR,LBUFR_BYTES,INTARR,DBLARR,
     & root, NELT, FRTPTR, FRTELT, COMM_LOAD,
     & ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     & MEM_DISTRIB,
     & DKEEP,PIVNUL_LIST,LPN_LIST
     &       ,LRGROUPS
     &     )
      USE DMUMPS_LOAD 
      USE DMUMPS_FAC_PAR_M
      USE DMUMPS_STRUC_DEF, ONLY : DMUMPS_ROOT_STRUC
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      TYPE (DMUMPS_ROOT_STRUC) :: root
      INTEGER(8) :: LA
      INTEGER N,NSTEPS,LIW,LPOOL,SLAVEF,COMM_NODES
      INTEGER MYID, MYID_NODES,LNA
      DOUBLE PRECISION A(LA)
      DOUBLE PRECISION RINFO(40)
      INTEGER LBUFR, LBUFR_BYTES
      INTEGER, DIMENSION(0:SLAVEF - 1) :: MEM_DISTRIB
      INTEGER BUFR( LBUFR )
      INTEGER NELT, LDPTRAR
      INTEGER FRTPTR(*), FRTELT(*)
      INTEGER LRGROUPS(N)
      DOUBLE PRECISION CNTL1
      INTEGER   ICNTL(60)
      INTEGER   INFO(80), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER   IW(LIW), SYM_PERM(N), NA(LNA),
     &          NE_STEPS(KEEP(28)), FILS(N),
     &          FRERE(KEEP(28)), NFSIZ(KEEP(28)), 
     &          DAD(KEEP(28))
      INTEGER   CAND(SLAVEF+1, max(1,KEEP(56)))
      INTEGER   STEP(N)
      INTEGER(8), INTENT(IN) :: PTRAR(LDPTRAR,2)
      INTEGER(8) :: PTRFAC(KEEP(28))
      INTEGER   PTRIST(KEEP(28)), PTLUST_S(KEEP(28))
      INTEGER   IW1(3*KEEP(28)), ITLOC(N+KEEP(253)), POOL(LPOOL)
      DOUBLE PRECISION :: RHS_MUMPS(KEEP(255))
      INTEGER(8) :: IW2(2*KEEP(28))
      INTEGER   PROCNODE_STEPS(KEEP(28))
      INTEGER   COMM_LOAD, ASS_IRECV
      INTEGER   ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      DOUBLE PRECISION   DBLARR(KEEP8(26))
      INTEGER   INTARR(KEEP8(27))
      DOUBLE PRECISION SEUIL, SEUIL_LDLT_NIV2
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
       INTEGER MUMPS_PROCNODE
       EXTERNAL MUMPS_PROCNODE
      DOUBLE PRECISION UULOC
      INTEGER LP, MPRINT
      INTEGER NSTK,PTRAST, NBPROCFILS
      INTEGER PIMASTER, PAMASTER
      LOGICAL PROK
      DOUBLE PRECISION ZERO, ONE
      DATA ZERO /0.0D0/
      DATA ONE /1.0D0/
      INTRINSIC int,real,log
      INTEGER IERR
      INTEGER NTOTPV, NTOTPVTOT, NMAXNPIV
      INTEGER(8) :: POSFAC, LRLU, IPTRLU, LRLUS
      INTEGER IWPOS, LEAF, NBROOT, NROOT
      KEEP(41)=0
      KEEP(42)=0
      NSTEPS   = 0
      LP     = ICNTL(1)
      MPRINT = ICNTL(2)
      PROK   = ((MPRINT.GT.0).AND.(ICNTL(4).GE.2))
      UULOC = CNTL1
      IF (UULOC.GT.ONE)   UULOC=ONE
      IF (UULOC.LT.ZERO)  UULOC=ZERO
      IF (KEEP(50).NE.0.AND.UULOC.GT.0.5D0) THEN
        UULOC = 0.5D0
      ENDIF
      PIMASTER   = 1
      NSTK       = PIMASTER + KEEP(28)
      NBPROCFILS = NSTK + KEEP(28)
      PTRAST = 1
      PAMASTER = 1 + KEEP(28)
      IF (KEEP(4).LE.0) KEEP(4)=32
      IF (KEEP(5).LE.0) KEEP(5)=16
      IF (KEEP(5).GT.KEEP(4)) KEEP(5) = KEEP(4)
      IF (KEEP(6).LE.0) KEEP(6)=24
      IF (KEEP(3).LE.KEEP(4)) KEEP(3)=KEEP(4)*2
      IF (KEEP(6).GT.KEEP(3)) KEEP(6) = KEEP(3)
      POSFAC = 1_8
      IWPOS  = 1
      LRLU = LA
      LRLUS = LRLU
      KEEP8(67) = LRLUS
      KEEP8(68) = 0_8
      KEEP8(69) = 0_8
      KEEP8(70) = 0_8
      KEEP8(71) = 0_8
      KEEP8(72) = 0_8
      KEEP8(73) = 0_8
      KEEP8(74) = 0_8
      IPTRLU = LRLU
      NTOTPV   = 0
      NMAXNPIV = 0
      IW1(NSTK:NSTK+KEEP(28)-1) = NE_STEPS(1:KEEP(28))
      CALL MUMPS_INIT_NROOT_DIST(N, NBROOT, NROOT,
     &                     MYID_NODES,
     &                     SLAVEF, NA, LNA,
     &                     KEEP, STEP,
     &                     PROCNODE_STEPS)
      CALL MUMPS_INIT_POOL_DIST(N, LEAF,
     &                     MYID_NODES,
     &                     SLAVEF, NA, LNA,
     &                     KEEP,KEEP8, STEP,
     &                     PROCNODE_STEPS,
     &                     POOL, LPOOL)
      CALL DMUMPS_INIT_POOL_LAST3(POOL, LPOOL, LEAF)     
      CALL DMUMPS_LOAD_INIT_SBTR_STRUCT(POOL, LPOOL,KEEP,KEEP8)
      IF ( KEEP( 38 ) .NE. 0 ) THEN
        NBROOT = NBROOT + root%NPROW * root%NPCOL - 1
      END IF
      IF ( root%yes )  THEN 
         IF ( MUMPS_PROCNODE( PROCNODE_STEPS(STEP(KEEP(38))), SLAVEF )
     &         .NE. MYID_NODES ) THEN
             NROOT = NROOT + 1
         END IF
      END IF
      CALL DMUMPS_FAC_PAR(N,IW,LIW,A,LA,IW1(NSTK),IW1(NBPROCFILS),
     &     NFSIZ,FILS,STEP,FRERE, DAD, CAND,
     &     ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &     INFO(11), NTOTPV, NMAXNPIV, PTRIST,IW2(PTRAST),
     &     IW1(PIMASTER), IW2(PAMASTER), PTRAR(1,2), 
     &     PTRAR(1,1),          
     &     ITLOC, RHS_MUMPS,
     &     POOL, LPOOL,
     &     RINFO, POSFAC,IWPOS,LRLU,IPTRLU, 
     &     LRLUS, LEAF, NROOT, NBROOT,
     &     UULOC,ICNTL,PTLUST_S,PTRFAC,NSTEPS,INFO,
     &     KEEP,KEEP8, PROCNODE_STEPS,SLAVEF,MYID,COMM_NODES,
     &     MYID_NODES, BUFR,LBUFR, LBUFR_BYTES,
     &     INTARR, DBLARR, root, SYM_PERM,
     &     NELT, FRTPTR, FRTELT, LDPTRAR, 
     &     COMM_LOAD, ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2,
     &     MEM_DISTRIB,NE_STEPS,
     &     DKEEP(1),PIVNUL_LIST(1),LPN_LIST
     &     ,LRGROUPS(1)
     &     )
      POSFAC = POSFAC -1_8
      IWPOS = IWPOS -1
      IF (KEEP(201).LE.0) THEN
        IF (KEEP(201) .EQ. -1 .AND. INFO(1) .LT. 0) THEN
          POSFAC = 0_8
        ENDIF
        KEEP8(31) = POSFAC
         RINFO(6)  = ZERO
      ELSE
         RINFO(6)  = dble(KEEP8(31))/1D6
      ENDIF
      KEEP8(48) = KEEP8(31)+KEEP8(71)
      KEEP(32) = IWPOS
      CALL MUMPS_SETI8TOI4(KEEP8(48), INFO(9))
      INFO(10) = KEEP(32)
      KEEP8(67) = LA - KEEP8(67)
      KEEP(89)  = NTOTPV
      KEEP(246) = NMAXNPIV
      INFO(23) = KEEP(89)
      CALL MPI_ALLREDUCE(NTOTPV, NTOTPVTOT, 1, MPI_INTEGER, MPI_SUM, 
     &                COMM_NODES, IERR)
      IF ( ( (INFO(1).EQ.-10 .OR. INFO(1).EQ.-40)
     &       .AND. (NTOTPVTOT.EQ.N) )
     &     .OR. ( NTOTPVTOT.GT.N ) ) THEN
       write(*,*) ' Error 1 NTOTPVTOT=', NTOTPVTOT,N
       CALL MUMPS_ABORT()
      ENDIF
      IF ( (KEEP(19).NE.0 ) .AND. (NTOTPVTOT.NE.N) .AND. 
     &     (INFO(1).GE.0) )  THEN
       write(*,*) ' Error 2 NTOTPVTOT=', NTOTPVTOT 
       CALL MUMPS_ABORT()
      ENDIF
      IF ( (INFO(1) .GE. 0 ) 
     &      .AND. (NTOTPVTOT.NE.N) ) THEN
         INFO(1) = -10
      ENDIF
      IF (INFO(1).EQ.-10) THEN
         INFO(2) = NTOTPVTOT
      ENDIF
      IF (PROK) THEN
        WRITE (MPRINT,99980) INFO(1), INFO(2),
     &       KEEP(28), KEEP8(31), INFO(10), INFO(11)
        IF(KEEP(50) .EQ. 0) THEN
          WRITE(MPRINT,99982) INFO(12)
        ENDIF
        IF (KEEP(50) .NE. 0) THEN
          WRITE(MPRINT,99984) INFO(12)
        ENDIF
        WRITE (MPRINT, 99986)
     &       INFO(13), INFO(14), INFO(25), RINFO(2), RINFO(3)
      ENDIF
      RETURN
99980 FORMAT (/' LEAVING FACTORIZATION PHASE WITH ...'/
     &      ' INFO (1)                                      =',I15/
     &      '  --- (2)                                      =',I15/
     &      '           Number of nodes in the tree         =',I15/
     &      ' INFO (9)  Real space for factors              =',I15/
     &      '  --- (10) Integer space for factors           =',I15/
     &      '  --- (11) Maximum size of frontal matrices    =',I15)
99982 FORMAT ('  --- (12) Number of off diagonal pivots       =',I15)
99984 FORMAT ('  --- (12) Number of negative pivots           =',I15)
99986 FORMAT ('  --- (13) Number of delayed pivots            =',I15/
     &      '  --- (14) Number of memory compresses         =',I15/
     &      '  --- (25) Number of entries in factors        =',I15/
     &  ' RINFO(2)  Operations during node assembly     =',1PD10.3/
     &  ' -----(3)  Operations during node elimination  =',1PD10.3)
      END SUBROUTINE DMUMPS_FAC_B
