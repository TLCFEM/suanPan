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
      SUBROUTINE MUMPS_SYMQAMD_NEW
     &                ( JOB, THRESH, NDENSE, 
     &                 N, TOTEL,
     &                 IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                 PERM, COMPLEM_LIST, SIZE_COMPLEM_LIST, 
     &                 AGG6 ) 
      IMPLICIT NONE
      INTEGER, INTENT(IN)    :: N, SIZE_COMPLEM_LIST
      INTEGER, INTENT(IN)    :: TOTEL
      INTEGER(8), INTENT(IN) :: IWLEN
      INTEGER, INTENT(IN)    :: THRESH
      LOGICAL, INTENT(IN)  ::  AGG6
      INTEGER, INTENT (IN) :: COMPLEM_LIST(max(1,SIZE_COMPLEM_LIST))
      INTEGER, INTENT(INOUT) :: JOB
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
      INTEGER, INTENT(INOUT)    :: PERM(N)
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(INOUT) :: NV(N)
      INTEGER, INTENT(OUT)   :: LAST(N) 
      INTEGER, INTENT(INOUT) :: ELEN(N)
      INTEGER, INTENT(OUT) :: NDENSE(N), DEGREE(N), 
     &                        HEAD(N), NEXT(N), W(N)
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
      INTEGER :: FDEG, ThresMin, ThresPrev, IBEGSchur, 
     &        ThresMinINIT
      LOGICAL :: AGG6_loc
      INTEGER :: THD_AGG
      LOGICAL :: SchurON, COMPRESS
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) MAXMEM, MEM, NEWMEM
      INTEGER :: MAXINT_N
      INTEGER(8) :: HASH, HMOD 
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &              PN, PSRC, PLN, PELN
      INTRINSIC max, min, mod
        IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0_8
           IF (NV(1).LT.0) NV(1) = 1
           RETURN
        ENDIF
        AGG6_loc = AGG6
        THD_AGG = max(128, min(TOTEL/2048, 1024))
        IF ( SIZE_COMPLEM_LIST < 0 .OR. SIZE_COMPLEM_LIST > N ) THEN
          WRITE(*,*) "Internal MUMPS_SYMQAMD_NEW", SIZE_COMPLEM_LIST,N
          CALL MUMPS_ABORT()
        ENDIF
        IF (JOB.EQ.2) THEN
          SchurON = .FALSE.        
        ENDIF
        THRESM    = THRESH  
        IF (JOB.NE.2) THEN
          SchurON   = (SIZE_COMPLEM_LIST > 0)
          IF ((JOB.EQ.1) .AND. (.NOT.SchurON) .AND. (N .GT. 0)) THEN
          ENDIF
          IBEGSchur = N-SIZE_COMPLEM_LIST+1
          IF (THRESM.GT.N) THRESM = N
          IF (THRESM.LT.0) THRESM = 0
          IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0_8
                ENDIF
             ENDIF
           ENDDO
          ENDIF
        ENDIF
        IF (SchurON) THEN
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
        ELSE
             THRESM    = max(int(31*N/32),THRESM)
             THRESM    = max(THRESM,1)
             ThresMin  = max( 3*THRESM / 4, 1)
             ThresPrev = THRESM
        ENDIF
        ThresMinINIT = ThresMin/4
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
          THRESM = N
       ENDIF
      ENDIF
      IF (JOB.EQ.2) THEN
      ENDIF
      PERMeqN = 0
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      NEL   = 0
      WFLG   = 2
      MAXINT_N=huge(WFLG)-TOTEL
      MINDEG = 1
      NCMPA  = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0  
      MEM  = PFREE - 1
      MAXMEM = MEM
      DO I = 1, N
        NDENSE(I)= 0
        LAST (I) = 0
        HEAD (I) = 0
        NEXT (I) = 0
        W (I) = 1
      ENDDO
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF (.NOT.COMPRESS) THEN
        IF (JOB.EQ.2) THEN
          DO I = 1,SIZE_COMPLEM_LIST
               X       = COMPLEM_LIST(I)
               ELEN(X) = -I       
               NV(X)   = LEN(X)+1 
               DEGREE(X) = LEN(X)
               DMAX = max(DMAX, LEN(X))
          ENDDO
          NEL = NEL + SIZE_COMPLEM_LIST  
          DO I=1, N
            IF (ELEN(I).LT.0) CYCLE
            DEGREE (I) = LEN (I)
            NV(I)      = 1
          ENDDO
        ELSE
          DO I=1, N
            ELEN(I)    = 0
            DEGREE (I) = LEN (I)
            NV(I)      = 1
          ENDDO
        ENDIF
      ELSE
        IF (JOB.EQ.2) THEN
          DO I = 1,SIZE_COMPLEM_LIST
               X       = COMPLEM_LIST(I)
               ELEN(X) = -I       
               NV(X)=1
               DO P=PE(X), PE(X)+int(LEN(X)-1,8)
                NV(X) = NV(X) + NV(IW(P))
               ENDDO 
               DEGREE(X) = NV(X)-1
               DMAX      = max(DMAX,DEGREE(X))
          ENDDO
          NEL = NEL + SIZE_COMPLEM_LIST  
          DO I=1,N
            IF (ELEN(I).LT.0) CYCLE
             DEGREE (I) = LEN (I) 
          ENDDO
        ELSE
          DO I=1, N
            ELEN (I)   = 0
            DEGREE(I) = 0
            DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
               DEGREE(I) = DEGREE(I) + NV(IW(P))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      DO 20 I = 1, N
        IF (ELEN(I).LT.0) CYCLE   
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
            NBD = NBD+NV(I)
            IF (FDEG.NE.N+1) THEN
             DEGREE(I) = DEGREE(I)+TOTEL+2
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
             NBED = NBED+NV(I)
             DEGREE(I) = TOTEL+1
             DEG = N
             IF (LASTD.EQ.0) THEN
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
   20 CONTINUE
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
   30 IF (NEL .LT. TOTEL) THEN
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.TOTEL)  THEN
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          MINDEG = 1
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.TOTEL+1) THEN
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + int(LEN(ME) -1,8)
            PLN       = P1
            PELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0_8) THEN
                X = E
  53            X = int(-PE(X))
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0 ) GOTO 53
                E = X
              ENDIF
              IF (ELEN(E).LT.0) THEN
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(PLN) = IW(PELN)
               IW(PELN) = E
               PLN  = PLN  + 1_8
               PELN = PELN + 1_8
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+LEN(E)-1
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(PLN)=E
               PLN = PLN+1_8
              ENDIF
  55        CONTINUE
            WFLG     = WFLG + 1
            LEN(ME)  = int(PLN-P1)
            ELEN(ME) = int(PELN-P1)
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(TOTEL+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
            IF (THRESM.LT.N) THEN
             ThresMin  = max(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = min(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = max(
     &         THRESM + int(sqrt(dble(ThresMin)))+ ThresMinINIT ,
     &         ThresPrev)
             THRESM    = min(THRESM,N) 
             ThresMin  = min(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
            GOTO 30
          ENDIF
          IF (DEGREE(ME).EQ.TOTEL+1) THEN
             IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in MUMPS_SYMQAMD_NEW ',
     &                ' quasi dense rows remains'
          CALL MUMPS_ABORT()
           ENDIF
           IF (JOB.EQ.1) THEN
            DO I = 1,SIZE_COMPLEM_LIST
             X       = COMPLEM_LIST(I)
             ELEN(X) = -(N-SIZE_COMPLEM_LIST+I)   
             NV(X)   = 1                
             PE(X)   = 0_8              
            ENDDO
            GOTO 265
           ENDIF
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0_8) .AND. (ELEN(X).LT.0)) THEN
             PE(X) = int(-COMPLEM_LIST(1),8)
            ELSEIF (DEGREE(X).EQ.TOTEL+1) THEN
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0_8
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal error 3 detected in MUMPS_SYMQAMD_NEW:'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           IF (ME.NE. COMPLEM_LIST(1)) THEN
             DO I=1, SIZE_COMPLEM_LIST
               PE(COMPLEM_LIST(I)) = int(-COMPLEM_LIST(1),8)
             ENDDO
             PE(COMPLEM_LIST(1)) = 0_8
             NV( COMPLEM_LIST(1))= NV(ME)
             NV(ME)               = 0
             ELEN( COMPLEM_LIST(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + int(LEN (ME) - 1,8)
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
              DEGME = DEGME + NVI
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).LE.TOTEL) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
          NEWMEM = 0
        ELSE
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
              E = ME
              PJ = P
              LN = SLENME
            ELSE
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
                IF (PFREE .GT. IWLEN) THEN
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
                  KNT1_UPDATED = 0
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
                  KNT2_UPDATED = 0
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + int(LENJ - 1,8)
                      PSRC = PSRC + int(LENJ - 1,8)
                    ENDIF
                    GO TO 80
                  ENDIF
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
                DEGME = DEGME + NVI
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
                IF (DEGREE(I).LE.TOTEL) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
        AGG6_loc = (AGG6 .OR. (DEGREE(ME) .LT. THD_AGG))
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 180
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
          DO 160 P = P1, P2
            E = IW (P)
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (.NOT. AGG6_loc .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (AGG6_loc .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
                PE (E) = int(-ME,8)
                W (E)  = 0
             ELSE IF (AGG6_loc .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
            ENDIF
  160     CONTINUE
          ELEN (I) = int(PN - P1 + 1_8)
          P3 = PN
          DO 170 P = P2 + 1, P1 + LEN (I) - 1
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
              IF (DEGREE(J).LE.TOTEL) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     &     .OR.
     &         (AGG6_loc.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     &       )
     &    THEN
            PE (I) = int(-ME, 8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
            IW (PN) = IW (P3)
            IW (P3) = IW (P1)
            IW (P1) = ME
            LEN (I) = int(PN - P1 + 1)
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.TOTEL) ) THEN
            HASH = int(LAST (I),kind=8)
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
              I = -J
              HEAD (HASH) = 0
            ELSE
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
              LN = LEN (I)
              ELN = ELEN (I)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
              JLAST = I
              J = NEXT (I)
  220         CONTINUE
              IF (J .NE. 0) THEN
                IF (LEN (J) .NE. LN) GO TO 240
                IF (ELEN (J) .NE. ELN) GO TO 240
                DO 230 P = PE (J) + 1, PE (J) + int(LN - 1,8)
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
                IF (PERM(J).GT.PERM(X)) THEN
                  PE (J) = int(-X,8)
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                  PE (X) = int(-J,8)
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
  240           CONTINUE
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
        IF ( (THRESM .GT. 0).AND.(THRESM.LT.N) ) THEN 
          THRESM = max(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
            NV (I) = NVI
            IF (DEGREE(I).LE.TOTEL) THEN
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
               IDENSE = .TRUE.
               DEGREE(I) = DEGREE(I)+TOTEL+2
             ENDIF
             IF (IDENSE) THEN
               P1 = PE(I)
               P2 = P1 + int(ELEN(I) - 1, 8)
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
               NBD = NBD+NVI
               FDEG = N
               DEG = N
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
             ENDIF
            ENDIF
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
            MINDEG = min (MINDEG, FDEG)
            ENDIF
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
        NV (ME) = NVPIV + DEGME
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
      GO TO 30
      ENDIF
  265 CONTINUE
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
          J = int(-PE (I))
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
            K = -ELEN (E)
            J = I
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J)= int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      DO 300 I = 1, N
        K = abs (ELEN (I))
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
        IF (PERMeqN.GT.0) PERM(PERMeqN) = N
      ENDIF
      PFREE = MAXMEM
      RETURN
      END SUBROUTINE MUMPS_SYMQAMD_NEW
      SUBROUTINE MUMPS_WRAP_GINP94 
     &       ( N, IPE, IW, LIW8,
     &         PERM, SizeOfBlocks,
     &         KEEP60, LISTVAR_SCHUR, SIZE_SCHUR, KEEP378,
     &   COLCOUNT, PARENT, 
     &   PORDER, IWTMP1, IWTMP2, IWTMP3, IWTMP4, 
     &   IWTMP5,
     &   INFO )  
      IMPLICIT NONE
      INTEGER, INTENT(IN)      :: N
      INTEGER(8), INTENT(IN)   :: LIW8
      INTEGER(8), INTENT(IN)   :: IPE(N+1)
      INTEGER,    INTENT(IN)   :: SizeOfBlocks(N)
      INTEGER, INTENT(INOUT)   :: PERM(N)
      INTEGER, INTENT(IN)      :: IW(LIW8)
      INTEGER, INTENT(OUT)     :: COLCOUNT(N)
      INTEGER, INTENT(OUT)     :: PARENT(N)
      INTEGER, INTENT(IN)      :: KEEP60, SIZE_SCHUR
      INTEGER, INTENT(IN)      :: LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)      :: KEEP378
      INTEGER, INTENT(INOUT)   :: INFO(2)
      INTEGER, INTENT(OUT):: PORDER(N), IWTMP1(N), IWTMP2(N) 
      INTEGER, INTENT(OUT):: IWTMP3(N), IWTMP4(N), IWTMP5(N)
      INTEGER :: I, KEEP378_loc, SIZE_SCHUR_EFF
      LOGICAL :: SizeOfBlocks_provided 
      SizeOfBlocks_provided = (SizeOfBlocks(1).NE.-1)
      IF (KEEP378.NE.0) KEEP378_loc=1 
      DO I=1, N
        IWTMP1(PERM(I)) = I
      END DO
      CALL MUMPS_GINP94_ELIM_TREE (
     &            N, IPE, IW, LIW8, IWTMP1, PERM, PARENT, 
     &            IWTMP2, INFO)
      IF (INFO(1).LT.0) RETURN
      CALL MUMPS_GINP94_POSTORDER(PARENT, N, PORDER, 
     &     IWTMP1, IWTMP2, IWTMP3, 
     &     INFO)
      IF (INFO(1).LT.0) RETURN
      IF (KEEP60.NE.0) THEN
           SIZE_SCHUR_EFF = SIZE_SCHUR
      ELSE
           SIZE_SCHUR_EFF = 0
      ENDIF
      CALL MUMPS_GINP94_COLCOUNTS(
     &         N, LIW8, IPE, IW, PARENT, PORDER, COLCOUNT, 
     &         SizeOfBlocks_provided, SizeOfBlocks, KEEP378_loc,
     &         SIZE_SCHUR_EFF, PERM,
     &         IWTMP1, IWTMP2, IWTMP3, IWTMP4, IWTMP5,
     &         INFO)
      IF (INFO(1).LT.0) RETURN
      IF (KEEP60.NE.0) THEN
       CALL MUMPS_GINP94_POSTPROCESS_SCHUR (
     &      N, PARENT, COLCOUNT, PERM,
     &      LISTVAR_SCHUR, SIZE_SCHUR )
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_WRAP_GINP94
        SUBROUTINE MUMPS_GINP94_ELIM_TREE (
     &              n, iptr, jcn, ljcn, invperm, perm, parent,
     &              ancestor, info)
        IMPLICIT NONE
        INTEGER(8), INTENT(IN)   :: ljcn
        integer                  :: n
        INTEGER(8), INTENT(IN)   :: iptr(n+1)
        integer, INTENT(IN)      :: jcn(ljcn), invperm(n), perm(n)
        integer, INTENT(OUT)     :: parent(n)
        integer, INTENT(INOUT)   :: INFO(2)
        integer, INTENT(OUT)     :: ancestor(n)
        integer                  :: jpos, i, j, k
        integer(8)               :: iidx8
       ancestor=0
       parent  =0
       do jpos = 1, n
         j = invperm(jpos)
         do iidx8 = iptr(j), iptr(j+1)-1
           i = jcn(iidx8)
           if (perm(i).ge.jpos) cycle 
           k = i
           call add_node(n, ancestor, parent, j, k)
         end do 
       end do 
       return
       contains
       subroutine add_node(n, ancestor, parent, j, i)
         implicit none
         integer, intent(in):: n
         integer      :: parent(n)
         integer      :: ancestor(n)
         integer      :: i, j
         integer      :: r, k
         if(i.eq.0) return
         k = i
         do
            r = ancestor(k)
            if (r .eq. j) then
               return
            end if
            ancestor(k) = j
            if(r .eq. 0) then
               parent(k) = j
               return
            end if
            k = r
         end do
         end subroutine add_node
        END SUBROUTINE MUMPS_GINP94_ELIM_TREE
        SUBROUTINE MUMPS_GINP94_POSTORDER(parent, n, porder, 
     &        son, brother, stack, 
     &        INFO
     &         )
        IMPLICIT NONE
        integer, intent(in)   :: n
        integer, intent(in)   :: parent(n)
        integer, intent(out)  :: porder(n) 
        integer, intent(inout):: INFO(2)
        integer, intent(out) :: son(n), brother(n), stack(n)
        integer              :: i, father, br, head, hp, pp, t
        son = 0
           do i=n, 1, -1
              father = parent(i)
              if (father .ne. 0) then
                 br          = son(father)
                 brother(i)  = br
                 son(father) = i
              end if
           end do
        head = 0
        hp   = 0
        pp   = 1
        do t=1, n
           if (parent(t) .ne. 0) cycle
           hp        = hp+1
           stack(hp) = t
           head      = t
           do
              if(son(head) .eq. 0) then
                 porder(pp) = head
                 pp = pp+1
                 hp = hp-1
                 if (parent(head) .ne. 0) then
                    son(parent(head)) = brother(head)
                 end if
                 if (hp .eq. 0) then
                    exit
                 end if
                 head = stack(hp)
              else
                 hp = hp+1
                 stack(hp) = son(head)
                 head = son(head)
              end if
           end do
        end do
        RETURN
        END SUBROUTINE MUMPS_GINP94_POSTORDER
        SUBROUTINE MUMPS_GINP94_COLCOUNTS( 
     &            n, ljcn, iptr, jcn, parent, porder, cc, 
     &            SizeOfBlocks_provided, SizeOfBlocks, KEEP378,
     &            SIZE_SCHUR_EFF, PERM,
     &            fst_desc, iporder, prev_p,  prev_nbr, setpath,
     &            INFO)
        IMPLICIT NONE
        integer, intent(in)    :: n 
        integer(8), intent(in) :: ljcn
        integer(8), intent(in) :: iptr(n+1)
        integer, intent(in)    ::  jcn(ljcn)
        integer, intent(inout) :: parent(n), porder(n) 
        integer, intent(in)    ::  SizeOfBlocks(n)
        logical, intent(in)    ::  SizeOfBlocks_provided
        integer, intent(in)    ::  KEEP378, SIZE_SCHUR_EFF, PERM(n)
        integer, intent(out)   :: cc(n)
         integer, intent(inout):: INFO(2)
        integer, intent(out)  :: fst_desc(n), iporder(n), prev_p(n) 
        integer, intent(out)  :: prev_nbr(n), setpath(n)
        integer               :: i, curr, fd, j, jidx, k
        integer(8)            :: iidx8
        integer               :: f, ref, p_leaf, q, jj
        integer               :: FIRSTinSchur, pi, pj
        logical               :: SCHUR_ON
        do j=1, n
           iporder(porder(j)) = j
        end do
        SCHUR_ON     = (SIZE_SCHUR_EFF.GT.0)
        FIRSTinSchur = N-SIZE_SCHUR_EFF+1
        cc = 0
        fst_desc=-1
        do i=1, n
           curr = porder(i)
           fd   = curr
           if(fst_desc(curr) .eq. -1) then
              if (SizeOfBlocks_provided) then
                  cc(curr) = SizeOfBlocks(curr)
              else
                  cc(curr) = 1
              endif
           end if
           do
              if (fst_desc(curr) .gt. 0) exit
              fst_desc(curr) = fd
              if (parent(curr) .eq. 0) exit
              curr = parent(curr)
           end do
        end do
        do j=1, n
           setpath(j)=j
        end do
        prev_p   = 0
        prev_nbr = 0
        do jidx=1, n
           j = abs(porder(jidx))
           if(parent(j) .ne. 0) then
              if (KEEP378.eq.1) then
               if (cc(parent(j)) .lt. 0) then
                  porder(iporder(parent(j)))= -parent(j)
               endif
              endif
              if (SizeOfBlocks_provided) then
                cc(parent(j)) = cc(parent(j)) - SizeOfBlocks(j)
              else
                cc(parent(j)) = cc(parent(j))-1
              endif
           endif
           do iidx8=iptr(j), iptr(j+1)-1
              i = jcn(iidx8)
              if (iporder(i).le.jidx) cycle
              if(prev_nbr(i) .eq. 0) then
                 ref = 0
              else
                 ref = iporder(prev_nbr(i))
              end if
              if(iporder(fst_desc(j)) .gt. ref) then
                 if (KEEP378.eq.1) then
                     porder(iporder(j))= -j
                 endif
                 if (SizeOfBlocks_provided) then
                   cc(j) = cc(j) + SizeOfBlocks(i)
                 else
                   cc(j) = cc(j) + 1
                 endif
                 p_leaf = prev_p(i)
                 if (p_leaf .ne. 0) then
                    q = setfind(setpath, p_leaf)
                    if (SizeOfBlocks_provided) then
                      cc(q) = cc(q) - SizeOfBlocks(i)
                    else
                      cc(q) = cc(q) - 1
                    endif
                 end if
                 prev_p(i) = j
              end if
              prev_nbr(i) = j
           end do
           if (parent(j).ne.0) setpath(j)=parent(j)
        end do
        do jj=1, n-1
           j=abs(porder(jj))
           if(parent(j) .ne. 0) cc(parent(j)) = cc(parent(j)) + cc(j)
        end do
        if (KEEP378.eq.1) then
          i=1
          do while (i .lt. n)
             porder(i) = abs(porder(i))
             pi        = porder(i)
             if (SCHUR_ON) then
               if (PERM(pi).GE.FIRSTinSchur) THEN
                i= i+1
                cycle
               endif
             endif
             j = i+1
             pj= porder(j)
             if (SCHUR_ON) then
              if (PERM(abs(pj)).GE.FIRSTinSchur) THEN
                i= j + 1
                cycle
              endif
             endif
             if (parent(pi).ne.0) then
              do while (pj.gt.0)
                j = j+1
                if (parent(abs(porder(j-1))).eq.0) exit
                if (j.gt.n) exit
                pj = porder(j)
                if (SCHUR_ON) then
                 if (PERM(abs(pj)).GE.FIRSTinSchur) exit
                endif
              end do
             endif
             parent(porder(i)) = parent(porder(j-1))
             do k=i+1, j-1
                parent(porder(k)) = -porder(i)
                cc(porder(k)) = 0
            end do
            i = j
          end do
          porder(n) = abs(porder(n))
        do i=1,n-1
          f = abs(parent(i))
          if (f.eq.0) cycle
          if (cc(f).eq.0) then
            parent(i) = parent(f)
          endif
        end do
        endif
        do i=1,n
          f = parent(i)
          if (parent(i).gt.0) then 
            parent(i) = -parent(i)
          endif
        end do
        return
        contains
        function setfind(setpath, p_leaf)
          implicit none
          integer :: setpath(:), p_leaf, setfind
          integer :: q, c, tmp
          q=p_leaf
          do while (setpath(q) .ne.q)
             q = setpath(q)
          end do
          c = p_leaf
          do while (c .ne.q)
             tmp = setpath(c)
             setpath(c) = q
             c = tmp
          end do
          setfind = q
          return
        end function setfind
        END SUBROUTINE MUMPS_GINP94_COLCOUNTS
        SUBROUTINE MUMPS_GINP94_POSTPROCESS_SCHUR (
     &      N, PARENT, COLCOUNT, PERM,
     &      LISTVAR_SCHUR, SIZE_SCHUR 
     &      )
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, SIZE_SCHUR 
        INTEGER, INTENT(IN) :: PERM(N), LISTVAR_SCHUR(SIZE_SCHUR)
        INTEGER, INTENT(INOUT) :: PARENT(N), COLCOUNT(N)
        INTEGER I, FIRSTinSchur, PrincipalVarSchur
        FIRSTinSchur      = N-SIZE_SCHUR+1
        PrincipalVarSchur = LISTVAR_SCHUR(1)
        DO I=1, N
         IF (I.EQ.PrincipalVarSchur) THEN
           IF ( PARENT(I) .NE. 0 ) THEN
             PARENT(I)   = 0 
           ENDIF
             COLCOUNT(I) = SIZE_SCHUR
         ELSE IF (PERM(I).GE.FIRSTinSchur) THEN
          PARENT(I)    = -PrincipalVarSchur
            COLCOUNT (I) = 0
         ELSE IF (PARENT(I) .NE. 0) THEN
           IF (PERM(-PARENT(I)).GE.FIRSTinSchur) THEN
            PARENT(I) =  -PrincipalVarSchur
           ENDIF
         ENDIF
        ENDDO
        RETURN
        END SUBROUTINE MUMPS_GINP94_POSTPROCESS_SCHUR 
