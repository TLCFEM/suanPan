      SUBROUTINE AMD(N,PE,IW,LEN,IWLEN,PFREE,NV,NEXT,LAST,HEAD,ELEN,
     +               DEGREE,NCMPA,W)
      IMPLICIT NONE
      INTEGER N,IWLEN,PFREE,NCMPA,IW(IWLEN),PE(N),DEGREE(N),NV(N),
     +        NEXT(N),LAST(N),HEAD(N),ELEN(N),W(N),LEN(N)
      INTEGER DEG,DEGME,DEXT,DMAX,E,ELENME,ELN,HASH,HMOD,I,ILAST,INEXT,
     +        J,JLAST,JNEXT,K,KNT1,KNT2,KNT3,LENJ,LN,MAXMEM,ME,MEM,
     +        MINDEG,NEL,NEWMEM,NLEFT,NVI,NVJ,NVPIV,SLENME,WE,WFLG,WNVI,
     +        X
      INTEGER P,P1,P2,P3,PDST,PEND,PJ,PME,PME1,PME2,PN,PSRC
      INTRINSIC MAX,MIN,MOD
      WFLG=2
      MINDEG=1
      NCMPA=0
      NEL=0
      HMOD=MAX(1,N-1)
      DMAX=0
      MEM=PFREE-1
      MAXMEM=MEM
      ME=0
      DO I=1,N
          LAST(I)=0
          HEAD(I)=0
          NV(I)=1
          W(I)=1
          ELEN(I)=0
          DEGREE(I)=LEN(I)
      ENDDO
      DO I=1,N
          DEG=DEGREE(I)
          IF (DEG.GT.0) THEN
              INEXT=HEAD(DEG)
              IF (INEXT.NE.0) LAST(INEXT)=I
              NEXT(I)=INEXT
              HEAD(DEG)=I
          ELSE
              NEL=NEL+1
              ELEN(I)=-NEL
              PE(I)=0
              W(I)=0
          ENDIF
      ENDDO
  100 IF (NEL.LT.N) THEN
          DO DEG=MINDEG,N
              ME=HEAD(DEG)
              IF (ME.GT.0) GOTO 150
          ENDDO
  150     MINDEG=DEG
          INEXT=NEXT(ME)
          IF (INEXT.NE.0) LAST(INEXT)=0
          HEAD(DEG)=INEXT
          ELENME=ELEN(ME)
          ELEN(ME)=-(NEL+1)
          NVPIV=NV(ME)
          NEL=NEL+NVPIV
          NV(ME)=-NVPIV
          DEGME=0
          IF (ELENME.EQ.0) THEN
              PME1=PE(ME)
              PME2=PME1-1
              DO P=PME1,PME1+LEN(ME)-1
                  I=IW(P)
                  NVI=NV(I)
                  IF (NVI.GT.0) THEN
                      DEGME=DEGME+NVI
                      NV(I)=-NVI
                      PME2=PME2+1
                      IW(PME2)=I
                      ILAST=LAST(I)
                      INEXT=NEXT(I)
                      IF (INEXT.NE.0) LAST(INEXT)=ILAST
                      IF (ILAST.NE.0) THEN
                          NEXT(ILAST)=INEXT
                      ELSE
                          HEAD(DEGREE(I))=INEXT
                      ENDIF
                  ENDIF
              ENDDO
              NEWMEM=0
          ELSE
              P=PE(ME)
              PME1=PFREE
              SLENME=LEN(ME)-ELENME
              DO KNT1=1,ELENME+1
                  IF (KNT1.GT.ELENME) THEN
                      E=ME
                      PJ=P
                      LN=SLENME
                  ELSE
                      E=IW(P)
                      P=P+1
                      PJ=PE(E)
                      LN=LEN(E)
                  ENDIF
                  DO KNT2=1,LN
                      I=IW(PJ)
                      PJ=PJ+1
                      NVI=NV(I)
                      IF (NVI.GT.0) THEN
                          IF (PFREE.GT.IWLEN) THEN
                              PE(ME)=P
                              LEN(ME)=LEN(ME)-KNT1
                              IF (LEN(ME).EQ.0) PE(ME)=0
                              PE(E)=PJ
                              LEN(E)=LN-KNT2
                              IF (LEN(E).EQ.0) PE(E)=0
                              NCMPA=NCMPA+1
                              DO J=1,N
                                  PN=PE(J)
                                  IF (PN.GT.0) THEN
                                      PE(J)=IW(PN)
                                      IW(PN)=-J
                                  ENDIF
                              ENDDO
                              PDST=1
                              PSRC=1
                              PEND=PME1-1
  152                         IF (PSRC.LE.PEND) THEN
                                  J=-IW(PSRC)
                                  PSRC=PSRC+1
                                  IF (J.GT.0) THEN
                                      IW(PDST)=PE(J)
                                      PE(J)=PDST
                                      PDST=PDST+1
                                      LENJ=LEN(J)
                                      DO KNT3=0,LENJ-2
                                         IW(PDST+KNT3)=IW(PSRC+KNT3)
                                      ENDDO
                                      PDST=PDST+LENJ-1
                                      PSRC=PSRC+LENJ-1
                                  ENDIF
                                  GOTO 152
                              ENDIF
                              P1=PDST
                              DO PSRC=PME1,PFREE-1
                                  IW(PDST)=IW(PSRC)
                                  PDST=PDST+1
                              ENDDO
                              PME1=P1
                              PFREE=PDST
                              PJ=PE(E)
                              P=PE(ME)
                          ENDIF
                          DEGME=DEGME+NVI
                          NV(I)=-NVI
                          IW(PFREE)=I
                          PFREE=PFREE+1
                          ILAST=LAST(I)
                          INEXT=NEXT(I)
                          IF (INEXT.NE.0) LAST(INEXT)=ILAST
                          IF (ILAST.NE.0) THEN
                              NEXT(ILAST)=INEXT
                          ELSE
                              HEAD(DEGREE(I))=INEXT
                          ENDIF
                      ENDIF
                  ENDDO
                  IF (E.NE.ME) THEN
                      PE(E)=-ME
                      W(E)=0
                  ENDIF
              ENDDO
              PME2=PFREE-1
              NEWMEM=PFREE-PME1
              MEM=MEM+NEWMEM
              MAXMEM=MAX(MAXMEM,MEM)
          ENDIF
          DEGREE(ME)=DEGME
          PE(ME)=PME1
          LEN(ME)=PME2-PME1+1
          IF (WFLG+N.LE.WFLG) THEN
              DO X=1,N
                  IF (W(X).NE.0) W(X)=1
              ENDDO
              WFLG=2
          ENDIF
          DO PME=PME1,PME2
              I=IW(PME)
              ELN=ELEN(I)
              IF (ELN.GT.0) THEN
                  NVI=-NV(I)
                  WNVI=WFLG-NVI
                  DO P=PE(I),PE(I)+ELN-1
                      E=IW(P)
                      WE=W(E)
                      IF (WE.GE.WFLG) THEN
                          WE=WE-NVI
                      ELSEIF (WE.NE.0) THEN
                          WE=DEGREE(E)+WNVI
                      ENDIF
                      W(E)=WE
                  ENDDO
              ENDIF
          ENDDO
          DO PME=PME1,PME2
              I=IW(PME)
              P1=PE(I)
              P2=P1+ELEN(I)-1
              PN=P1
              HASH=0
              DEG=0
              DO P=P1,P2
                  E=IW(P)
                  DEXT=W(E)-WFLG
                  IF (DEXT.GT.0) THEN
                      DEG=DEG+DEXT
                      IW(PN)=E
                      PN=PN+1
                      HASH=HASH+E
                  ELSEIF (DEXT.EQ.0) THEN
                      PE(E)=-ME
                      W(E)=0
                  ENDIF
              ENDDO
              ELEN(I)=PN-P1+1
              P3=PN
              DO P=P2+1,P1+LEN(I)-1
                  J=IW(P)
                  NVJ=NV(J)
                  IF (NVJ.GT.0) THEN
                      DEG=DEG+NVJ
                      IW(PN)=J
                      PN=PN+1
                      HASH=HASH+J
                  ENDIF
              ENDDO
              IF (DEG.EQ.0) THEN
                  PE(I)=-ME
                  NVI=-NV(I)
                  DEGME=DEGME-NVI
                  NVPIV=NVPIV+NVI
                  NEL=NEL+NVI
                  NV(I)=0
                  ELEN(I)=0
              ELSE
                  DEGREE(I)=MIN(DEGREE(I),DEG)
                  IW(PN)=IW(P3)
                  IW(P3)=IW(P1)
                  IW(P1)=ME
                  LEN(I)=PN-P1+1
                  HASH=MOD(HASH,HMOD)+1
                  J=HEAD(HASH)
                  IF (J.LE.0) THEN
                      NEXT(I)=-J
                      HEAD(HASH)=-I
                  ELSE
                      NEXT(I)=LAST(J)
                      LAST(J)=I
                  ENDIF
                  LAST(I)=HASH
              ENDIF
          ENDDO
          DEGREE(ME)=DEGME
          DMAX=MAX(DMAX,DEGME)
          WFLG=WFLG+DMAX
          IF (WFLG+N.LE.WFLG) THEN
              DO X=1,N
                  IF (W(X).NE.0) W(X)=1
              ENDDO
              WFLG=2
          ENDIF
          DO PME=PME1,PME2
              I=IW(PME)
              IF (NV(I).LT.0) THEN
                  HASH=LAST(I)
                  J=HEAD(HASH)
                  IF (J.NE.0) THEN
                      IF (J.LT.0) THEN
                          I=-J
                          HEAD(HASH)=0
                      ELSE
                          I=LAST(J)
                          LAST(J)=0
                      ENDIF
                      IF (I.NE.0) THEN
  154                     IF (NEXT(I).NE.0) THEN
                              LN=LEN(I)
                              ELN=ELEN(I)
                              DO P=PE(I)+1,PE(I)+LN-1
                                  W(IW(P))=WFLG
                              ENDDO
                              JLAST=I
                              J=NEXT(I)
  156                         IF (J.NE.0) THEN
                                  IF (LEN(J).EQ.LN) THEN
                                      IF (ELEN(J).EQ.ELN) THEN
                                         DO P=PE(J)+1,PE(J)+LN-1
                                         IF (W(IW(P)).NE.WFLG) GOTO 158
                                         ENDDO
                                         PE(J)=-I
                                         NV(I)=NV(I)+NV(J)
                                         NV(J)=0
                                         ELEN(J)=0
                                         J=NEXT(J)
                                         NEXT(JLAST)=J
                                         GOTO 156
                                      ENDIF
                                  ENDIF
  158                             JLAST=J
                                  J=NEXT(J)
                                  GOTO 156
                              ENDIF
                              WFLG=WFLG+1
                              I=NEXT(I)
                              IF (I.NE.0) GOTO 154
                          ENDIF
                      ENDIF
                  ENDIF
              ENDIF
          ENDDO
          P=PME1
          NLEFT=N-NEL
          DO PME=PME1,PME2
              I=IW(PME)
              NVI=-NV(I)
              IF (NVI.GT.0) THEN
                  NV(I)=NVI
                  DEG=MIN(DEGREE(I)+DEGME-NVI,NLEFT-NVI)
                  INEXT=HEAD(DEG)
                  IF (INEXT.NE.0) LAST(INEXT)=I
                  NEXT(I)=INEXT
                  LAST(I)=0
                  HEAD(DEG)=I
                  MINDEG=MIN(MINDEG,DEG)
                  DEGREE(I)=DEG
                  IW(P)=I
                  P=P+1
              ENDIF
          ENDDO
          NV(ME)=NVPIV+DEGME
          LEN(ME)=P-PME1
          IF (LEN(ME).EQ.0) THEN
              PE(ME)=0
              W(ME)=0
          ENDIF
          IF (NEWMEM.NE.0) THEN
              PFREE=P
              MEM=MEM-NEWMEM+LEN(ME)
          ENDIF
          GOTO 100
      ENDIF
      DO I=1,N
          IF (ELEN(I).EQ.0) THEN
              J=-PE(I)
  160         IF (ELEN(J).GE.0) THEN
                  J=-PE(J)
                  GOTO 160
              ENDIF
              E=J
              K=-ELEN(E)
              J=I
  180         IF (ELEN(J).GE.0) THEN
                  JNEXT=-PE(J)
                  PE(J)=-E
                  IF (ELEN(J).EQ.0) THEN
                      ELEN(J)=K
                      K=K+1
                  ENDIF
                  J=JNEXT
                  GOTO 180
              ENDIF
              ELEN(E)=-K
          ENDIF
      ENDDO
      DO I=1,N
          K=ABS(ELEN(I))
          LAST(K)=I
          ELEN(I)=K
      ENDDO
      PFREE=MAXMEM
      END
