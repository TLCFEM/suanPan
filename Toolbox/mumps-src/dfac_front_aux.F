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
      MODULE DMUMPS_FAC_FRONT_AUX_M
      CONTAINS
      SUBROUTINE DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &   INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,KEEP, DKEEP,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U,MAXFROMN,IS_MAXFROMN_AVAIL, 
     &     Inextpiv, OOCWRITE_COMPATIBLE_WITH_BLR
     &)
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER NFRONT,NASS,LIW,INOPV
      INTEGER(8) :: LA
      INTEGER KEEP(500)
      DOUBLE PRECISION    DKEEP(230)
      DOUBLE PRECISION UU, SEUIL
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION, intent(in) :: MAXFROMN
      LOGICAL, intent(inout) :: IS_MAXFROMN_AVAIL
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOCWRITE_COMPATIBLE_WITH_BLR
      DOUBLE PRECISION AMROW
      DOUBLE PRECISION RMAX
      DOUBLE PRECISION  SWOP
      INTEGER(8) :: APOS, POSELT
      INTEGER(8) :: J1, J2, J3_8, JJ, IDIAG
      INTEGER(8) :: J1_ini
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER NOFFW,NPIV,IPIV,IPIV_SHIFT
      INTEGER J, J3
      INTEGER NPIVP1,JMAX,ISW,ISWPS1
      INTEGER ISWPS2,KSW,XSIZE
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER ISHIFT, K206
      INTEGER DMUMPS_IXAMAX
      INCLUDE 'mumps_headers.h'
      INTRINSIC max
      DOUBLE PRECISION, PARAMETER :: RZERO = 0.0D0
#if defined(_OPENMP)
      INTEGER :: NOMP, CHUNK, K360
      K360 = KEEP(360)
      NOMP    = OMP_GET_MAX_THREADS()
#endif
        NFRONT8 = int(NFRONT,8)
        INOPV   = 0
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        NPIVP1  = NPIV + 1
        K206    = KEEP(206)
        IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1.AND.
     &      OOCWRITE_COMPATIBLE_WITH_BLR) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)
     &              +KEEP(IXSZ),
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        ISHIFT = 0   
        IF (K206.GE.1) THEN
          IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.NASS) THEN
            ISHIFT = Inextpiv - NPIVP1
          ENDIF
          IF (ISHIFT.GT.0.AND.IS_MAXFROMN_AVAIL) THEN
            IPIV = NPIVP1
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)).GT.max(UU*MAXFROMN,SEUIL,
     &                               tiny(MAXFROMN))) THEN
              ISHIFT = 0
            ENDIF
          ENDIF
          IF ( ISHIFT .GT. 0) THEN
            IS_MAXFROMN_AVAIL = .FALSE.
          ENDIF
        ENDIF  
          DO 460 IPIV_SHIFT=NPIVP1+ISHIFT,NASS+ISHIFT
            IF (IPIV_SHIFT .LE. NASS) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV=IPIV_SHIFT-NASS-1+NPIVP1
            ENDIF
            APOS = POSELT + NFRONT8*int(NPIV,8) + int(IPIV-1,8)
            JMAX = 1
            AMROW = RZERO
            J1    = APOS
            J3    = NASS -NPIV
            JMAX  = DMUMPS_IXAMAX(J3,A(J1),NFRONT,KEEP(360))
            JJ    = J1 + int(JMAX-1,8)*NFRONT8
            AMROW = abs(A(JJ))
            RMAX  = AMROW
            J1    = APOS +  int(NASS-NPIV,8) * NFRONT8
            J3 = NFRONT - NASS - KEEP(253)
            IF (IS_MAXFROMN_AVAIL) THEN
              RMAX = max(MAXFROMN,RMAX)
              IS_MAXFROMN_AVAIL = .FALSE.
            ELSE
              IF (J3.EQ.0) GOTO 370
              IF (KEEP(351).EQ.1) THEN
                J1_ini = J1
!$              CHUNK = max(K360/2,(J3+NOMP-1)/NOMP)
!$OMP  PARALLEL DO schedule(static, CHUNK)
!$OMP& FIRSTPRIVATE(J1_ini,NFRONT8,J3)
!$OMP& REDUCTION(max:RMAX) IF (J3.GE.K360)
                DO J=1,J3
                  RMAX = max(abs(A(J1_ini + int(J-1,8) * NFRONT8)),
     &                       RMAX)
                END DO
!$OMP  END PARALLEL DO
              ELSE
                DO J=1,J3
                  RMAX = max(abs(A(J1)), RMAX)
                  J1 = J1 + NFRONT8
                END DO
              ENDIF
            END IF
  370       IF (RMAX.LE.tiny(RMAX)) GO TO 460
            IDIAG = APOS + int(IPIV - NPIVP1,8)*NFRONT8
            IF (abs(A(IDIAG)).GT.max(UU*RMAX,SEUIL,tiny(RMAX))) THEN
              JMAX = IPIV - NPIV
              GO TO 380
            ENDIF
            IF (AMROW.LE.max(UU*RMAX,SEUIL,tiny(RMAX))) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
            IF (K206.GE.1) THEN
              Inextpiv = IPIV + 1 
            ENDIF
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(
     &             A(APOS + int(JMAX - 1,8) * NFRONT8 ),
     &             DKEEP(6),
     &             KEEP(259) )
            ENDIF
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1   = POSELT + int(NPIV,8)
            J3_8 = POSELT + int(IPIV-1,8)
            DO J= 1,NFRONT
              SWOP  = A(J1)
              A(J1) = A(J3_8)
              A(J3_8) = SWOP
              J1 = J1 + NFRONT8
              J3_8 = J3_8 + NFRONT8
            END DO
            ISWPS1 = IOLDPS + 5 + NPIVP1 + NFRONT + XSIZE
            ISWPS2 = IOLDPS + 5 + IPIV + NFRONT + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8) * NFRONT8
            J2 = POSELT + int(NPIV + JMAX - 1,8) * NFRONT8
            DO KSW=1,NFRONT
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + 1_8
              J2 = J2 + 1_8
            END DO
            ISWPS1 = IOLDPS + 5 + NPIV + 1 + XSIZE
            ISWPS2 = IOLDPS + 5 + NPIV + JMAX + XSIZE
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420
  460    CONTINUE
       INOPV = 1
       GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIV+JMAX,
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, IPIV,
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
 430  CONTINUE
      IS_MAXFROMN_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_FAC_H
      SUBROUTINE DMUMPS_FAC_M(IBEG_BLOCK,
     &     NFRONT,NASS,N,INODE,IW,LIW,A,LA,
     &     IOLDPS,POSELT,IFINB,LKJIB,LKJIT,XSIZE)
      IMPLICIT NONE
      INTEGER NFRONT,NASS,N,LIW,INODE,IFINB,LKJIB,IBEG_BLOCK
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS
      INTEGER LKJIT, XSIZE
      DOUBLE PRECISION ONE, ALPHA
      INTEGER NPIV,JROW2
      INTEGER NEL2,NPIVP1,KROW,NEL
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1
        IFINB  = 0
        IF (IW(IOLDPS+3+XSIZE).LE.0) THEN
          IF (NASS.LT.LKJIT) THEN
           IW(IOLDPS+3+XSIZE) = NASS
          ELSE
           IW(IOLDPS+3+XSIZE) = min0(NASS,LKJIB)
          ENDIF
        ENDIF
        JROW2 = IW(IOLDPS+3+XSIZE)
        NEL2   = JROW2 - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (JROW2.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
          IW(IOLDPS+3+XSIZE) = min0(JROW2+LKJIB,NASS)
          IBEG_BLOCK = NPIVP1+1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_M
      SUBROUTINE DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &       IOLDPS,POSELT,IFINB,XSIZE,
     &       KEEP,MAXFROMN,IS_MAXFROMN_AVAIL)
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mumps_headers.h'
      INTEGER NFRONT,NASS,LIW,IFINB
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW)
      DOUBLE PRECISION    ALPHA,VALPIV
      INTEGER(8) :: APOS, POSELT, UUPOS, LPOS, IRWPOS
      INTEGER(8) :: NFRONT8
      INTEGER IOLDPS,NPIV,XSIZE
      INTEGER, intent(in) :: KEEP(500)
      DOUBLE PRECISION, intent(inout) :: MAXFROMN
      LOGICAL, intent(inout) :: IS_MAXFROMN_AVAIL
      INTEGER NEL,IROW,NEL2,JCOL, NCB
      INTEGER NPIVP1
      DOUBLE PRECISION, PARAMETER :: ONE = 1.0D0
#if defined(_OPENMP)
      LOGICAL:: OMP_FLAG
      INTEGER:: NOMP, K360, CHUNK
      NOMP = OMP_GET_MAX_THREADS()
      K360 = KEEP(360)
#endif
        NFRONT8=int(NFRONT,8)
        NPIV   = IW(IOLDPS+1+XSIZE)
        NPIVP1 = NPIV + 1
        NEL    = NFRONT - NPIVP1 
        NEL2   = NASS - NPIVP1
        NCB    = NFRONT - NASS - KEEP(253)
        IFINB  = 0
        IF (NPIVP1.EQ.NASS) IFINB = 1
        APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
        VALPIV = ONE/A(APOS)
#if defined(_OPENMP)
        OMP_FLAG = .FALSE.
        CHUNK = max(NEL,1)
        IF (NOMP.GT.1) THEN
          IF (NEL.LT.K360) THEN
              IF (NEL*NEL2.GE.KEEP(361)) THEN
                OMP_FLAG = .TRUE.
                CHUNK = max(20, (NEL+NOMP-1)/NOMP)
              ENDIF
          ELSE
            OMP_FLAG = .TRUE.
            CHUNK = max(K360/2, (NEL+NOMP-1)/NOMP)
          ENDIF
        ENDIF
#endif
        IF (KEEP(351).EQ.2) THEN
          MAXFROMN = 0.0D0
          IF (NEL2 > 0) THEN
            IS_MAXFROMN_AVAIL = .TRUE.
          ENDIF
!$OMP PARALLEL DO schedule(static, CHUNK)
!$OMP& PRIVATE(LPOS, UUPOS, IRWPOS, ALPHA, JCOL)
!$OMP& FIRSTPRIVATE(APOS,NFRONT8,VALPIV,NEL,NEL2)
!$OMP& REDUCTION(max:MAXFROMN) IF(OMP_FLAG)
          DO IROW = 1, NEL
            LPOS = APOS + NFRONT8*int(IROW,8)
            A(LPOS) = A(LPOS)*VALPIV
            ALPHA   = -A(LPOS)
            IRWPOS  = LPOS + 1_8
            UUPOS  = APOS + 1_8
            IF (NEL2 > 0) THEN
              A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
              MAXFROMN=max(MAXFROMN, abs(A(IRWPOS)))
              IRWPOS = IRWPOS+1_8
              UUPOS  = UUPOS+1_8
              DO JCOL = 2, NEL2
                A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
                IRWPOS = IRWPOS+1_8
                UUPOS  = UUPOS+1_8
              ENDDO
            ENDIF
          END DO
!$OMP END PARALLEL DO
        ELSE
!$OMP PARALLEL DO schedule(static, CHUNK)
!$OMP& FIRSTPRIVATE(APOS,NFRONT8,VALPIV,NEL,NEL2)
!$OMP& PRIVATE(LPOS, UUPOS, IRWPOS, ALPHA, JCOL) IF(OMP_FLAG)
          DO IROW = 1, NEL
            LPOS = APOS + NFRONT8*int(IROW,8)
            A(LPOS) = A(LPOS)*VALPIV
            ALPHA   = -A(LPOS)
            IRWPOS  = LPOS + 1_8
            UUPOS  = APOS + 1_8
            DO JCOL = 1, NEL2
              A(IRWPOS) = A(IRWPOS) + ALPHA*A(UUPOS)
              IRWPOS = IRWPOS+1_8
              UUPOS  = UUPOS+1_8
            ENDDO
          ENDDO
!$OMP END PARALLEL DO
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_N
      SUBROUTINE DMUMPS_FAC_P(A,LA,NFRONT,
     &       NPIV,NASS,POSELT,CALL_UTRSM
     &      )
      IMPLICIT NONE
      INTEGER(8) :: LA,POSELT
      DOUBLE PRECISION    A(LA)
      INTEGER NFRONT, NPIV, NASS
      LOGICAL CALL_UTRSM
      INTEGER(8) :: LPOS, LPOS1, LPOS2, UPOS
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = POSELT + int(NASS,8)*int(NFRONT,8)
         CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(POSELT),NFRONT,
     &        A(LPOS2),NFRONT)
       IF (CALL_UTRSM) THEN
         UPOS  = POSELT + int(NASS,8)
           CALL dtrsm('R', 'U', 'N', 'U', NEL1, NPIV, ONE,
     &                 A(POSELT), NFRONT, A(UPOS), NFRONT)
       ENDIF
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = POSELT + int(NPIV,8)
         CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &        NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
      RETURN
      END SUBROUTINE DMUMPS_FAC_P
      SUBROUTINE DMUMPS_FAC_P_PANEL(A,LAFAC,NFRONT,
     &      NPIV,NASS, IW, LIWFAC,
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
      USE DMUMPS_OOC   
      IMPLICIT NONE
      INTEGER NFRONT, NPIV, NASS
      INTEGER(8) :: LAFAC
      INTEGER  LIWFAC, TYPEFile, MYID, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten, STRAT
      DOUBLE PRECISION  A(LAFAC)
      INTEGER  IW(LIWFAC)
      INTEGER(8) KEEP8(150)
      TYPE(IO_BLOCK) :: MonBloc 
      INTEGER(8) :: LPOS2,LPOS1,LPOS
      INTEGER NEL1,NEL11
      DOUBLE PRECISION ALPHA, ONE
      LOGICAL LAST_CALL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        LPOS2  = 1_8 + int(NASS,8) * int(NFRONT,8)
        CALL dtrsm('L','L','N','N',NPIV,NEL1,ONE,A(1),NFRONT,
     &              A(LPOS2),NFRONT)
        LAST_CALL=.FALSE.
           CALL DMUMPS_OOC_IO_LU_PANEL
     &          ( STRAT, TYPEFile, 
     &           A, LAFAC, MonBloc,
     &           LNextPiv2beWritten, UNextPiv2beWritten,
     &           IW, LIWFAC, 
     &           MYID, KEEP8(31), IFLAG_OOC,LAST_CALL )
        LPOS   = LPOS2 + int(NPIV,8)
        LPOS1  = int(1 + NPIV,8)
        CALL dgemm('N','N',NEL11,NEL1,NPIV,ALPHA,A(LPOS1),
     &          NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_FAC_P_PANEL
      SUBROUTINE DMUMPS_FAC_T(A,LA,NPIVB,NFRONT,
     &                             NPIV,NASS,POSELT)
      IMPLICIT NONE
      INTEGER NPIVB,NASS
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER(8) :: APOS, POSELT
      INTEGER NFRONT, NPIV, NASSL
      INTEGER(8) :: LPOS, LPOS1, LPOS2
      INTEGER NEL1, NEL11, NPIVE
      DOUBLE PRECISION    ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NEL1   = NFRONT - NASS
        NEL11  = NFRONT - NPIV
        NPIVE  = NPIV - NPIVB
        NASSL  = NASS - NPIVB
        APOS   = POSELT + int(NPIVB,8)*int(NFRONT,8)
     &                  + int(NPIVB,8)
        LPOS2  = APOS + int(NASSL,8)
        CALL dtrsm('R','U','N','U',NEL1,NPIVE,ONE,A(APOS),NFRONT,
     &              A(LPOS2),NFRONT)
        LPOS   = LPOS2 + int(NFRONT,8)*int(NPIVE,8)
        LPOS1  = APOS  + int(NFRONT,8)*int(NPIVE,8)
        CALL dgemm('N','N',NEL1,NEL11,NPIVE,ALPHA,A(LPOS2),
     &          NFRONT,A(LPOS1),NFRONT,ONE,A(LPOS),NFRONT)
        RETURN
        END SUBROUTINE DMUMPS_FAC_T
      SUBROUTINE DMUMPS_FAC_SQ(IBEG_BLOCK, IEND_BLOCK, NPIV,
     &    NFRONT, LAST_ROW, LAST_COL, A, LA, POSELT, 
     &    FIRST_COL, CALL_LTRSM, CALL_UTRSM, CALL_GEMM, 
     &    WITH_COMM_THREAD )
!$    USE OMP_LIB
!$    USE DMUMPS_BUF
      IMPLICIT NONE
      INTEGER, intent(in)     :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(in)     :: NPIV, NFRONT, LAST_ROW, LAST_COL
      INTEGER, intent(in)     :: FIRST_COL
      INTEGER(8), intent(in)  :: LA
      DOUBLE PRECISION, intent(inout)  :: A(LA)
      INTEGER(8), intent(in)  :: POSELT 
      LOGICAL, intent(in)     :: CALL_LTRSM, CALL_UTRSM, CALL_GEMM
      LOGICAL, intent(in) :: WITH_COMM_THREAD
      INTEGER(8) :: NFRONT8, LPOSN, LPOS2N
      INTEGER(8) :: LPOS, LPOS1, LPOS2, UPOS, POSELT_LOCAL
      INTEGER NELIM, LKJIW, NEL1, NEL11, UTRSM_NCOLS
      DOUBLE PRECISION ALPHA, ONE
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
!$    INTEGER :: NOMP
!$    LOGICAL :: TRSM_GEMM_FINISHED
!$    LOGICAL :: SAVE_NESTED, SAVE_DYNAMIC
      NFRONT8= int(NFRONT,8)
      NELIM  = IEND_BLOCK - NPIV
      NEL1   = LAST_ROW - IEND_BLOCK
      IF ( NEL1 < 0 ) THEN
        WRITE(*,*)
     &  "Internal error 1 in DMUMPS_FAC_SQ,IEND_BLOCK>LAST_ROW",
     &  IEND_BLOCK, LAST_ROW
        CALL MUMPS_ABORT()
      ENDIF
      LKJIW  = NPIV - IBEG_BLOCK + 1
      NEL11  = LAST_COL - NPIV
      LPOS2  = POSELT + int(IEND_BLOCK,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
      UTRSM_NCOLS = LAST_COL - FIRST_COL
      UPOS   = POSELT + int(IBEG_BLOCK-1,8)*NFRONT8 + int(FIRST_COL,8)
      POSELT_LOCAL = POSELT + int(IBEG_BLOCK-1,8)*NFRONT8 
     &                      + int(IBEG_BLOCK-1,8)
      IF ((NEL1.NE.0).AND.(LKJIW.NE.0)) THEN
        IF (WITH_COMM_THREAD .EQV. .FALSE.) THEN
           IF (CALL_LTRSM) THEN
             CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &                A(POSELT_LOCAL),NFRONT,
     &                A(LPOS2),NFRONT)
           ENDIF
           IF (CALL_UTRSM) THEN
             CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                A(POSELT_LOCAL),NFRONT,
     &                A(UPOS),NFRONT)
             LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
             LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
             CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &               LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &               NFRONT,ONE,A(LPOSN),NFRONT)
           ENDIF
           IF (CALL_GEMM) THEN
            LPOS   = LPOS2 + int(LKJIW,8)
            LPOS1  = POSELT_LOCAL + int(LKJIW,8)
            CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &           NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
           END IF
        ELSE
!$        NOMP = OMP_GET_MAX_THREADS()
!$        CALL OMP_SET_NUM_THREADS(2)
!$        SAVE_NESTED = OMP_GET_NESTED()
!$        SAVE_DYNAMIC = OMP_GET_DYNAMIC()
!$        CALL OMP_SET_NESTED(.TRUE.)
!$        CALL OMP_SET_DYNAMIC(.FALSE.)
!$        TRSM_GEMM_FINISHED = .FALSE.
!$OMP     PARALLEL SHARED(TRSM_GEMM_FINISHED)
!$        IF (OMP_GET_THREAD_NUM() .EQ. 1) THEN
!$          CALL OMP_SET_NUM_THREADS(NOMP)
          IF (CALL_LTRSM) THEN
            CALL dtrsm('L','L','N','N',LKJIW,NEL1,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(LPOS2),NFRONT)
          ENDIF
          IF (CALL_UTRSM) THEN
            CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(UPOS),NFRONT)
            LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
            LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
            CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &                LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &                NFRONT,ONE,A(LPOSN),NFRONT)
          ENDIF
          IF (CALL_GEMM) THEN
            LPOS   = LPOS2 + int(LKJIW,8)
            LPOS1  = POSELT_LOCAL + int(LKJIW,8)
            CALL dgemm('N','N',NEL11,NEL1,LKJIW,ALPHA,A(LPOS1),
     &             NFRONT,A(LPOS2),NFRONT,ONE,A(LPOS),NFRONT)
          END IF
!$           TRSM_GEMM_FINISHED = .TRUE.
!$        ELSE
!$           DO WHILE (.NOT. TRSM_GEMM_FINISHED)
!$             CALL DMUMPS_BUF_TEST()
!$             CALL MUMPS_USLEEP(10000)
!$           END DO
!$        END IF
!$OMP     END PARALLEL
!$        CALL OMP_SET_NESTED(SAVE_NESTED)
!$        CALL OMP_SET_DYNAMIC(SAVE_DYNAMIC)
!$        CALL OMP_SET_NUM_THREADS(NOMP)
        ENDIF
      ELSE
        IF (CALL_UTRSM.AND.UTRSM_NCOLS.NE.0) THEN
          CALL dtrsm('R','U','N','U',UTRSM_NCOLS,LKJIW,ONE,
     &                 A(POSELT_LOCAL),NFRONT,
     &                 A(UPOS),NFRONT)
          LPOS2N = POSELT + int(NPIV,8)*NFRONT8 + int(IBEG_BLOCK-1,8)
          LPOSN  = POSELT + int(NPIV,8)*NFRONT8 + int(FIRST_COL,8)
          CALL dgemm('N','N',UTRSM_NCOLS,NELIM,
     &              LKJIW,ALPHA,A(UPOS),NFRONT,A(LPOS2N),
     &              NFRONT,ONE,A(LPOSN),NFRONT)
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ
      SUBROUTINE DMUMPS_FAC_MQ(IBEG_BLOCK,IEND_BLOCK,
     &     NFRONT, NASS, NPIV, LAST_COL, A, LA, POSELT, IFINB)
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK, IEND_BLOCK, NFRONT, 
     &                          NASS, NPIV, LAST_COL
      INTEGER, intent(out)   ::  IFINB
      INTEGER(8), intent(in) :: LA, POSELT
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION    :: VALPIV
      INTEGER(8) :: APOS,  UUPOS, LPOS
      INTEGER(8) :: NFRONT8
      DOUBLE PRECISION    :: ONE, ALPHA
      INTEGER    :: NEL2,NPIVP1,KROW,NEL
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
        NFRONT8= int(NFRONT,8)
        NPIVP1 = NPIV + 1
        NEL  = LAST_COL - NPIVP1
        IFINB  = 0
        NEL2   = IEND_BLOCK - NPIVP1
        IF (NEL2.EQ.0) THEN
         IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
         ELSE
          IFINB        = 1
         ENDIF
        ELSE
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + NFRONT8
         DO 541 KROW = 1,NEL2
             A(LPOS) = A(LPOS)*VALPIV
             LPOS    = LPOS + NFRONT8
 541     CONTINUE
         LPOS   = APOS + NFRONT8
         UUPOS  = APOS + 1_8
#if defined(MUMPS_USE_BLAS2)
         CALL dger(NEL,NEL2,ALPHA,A(UUPOS),1,A(LPOS),NFRONT,
     &              A(LPOS+1_8),NFRONT)
#else
         CALL dgemm('N','N',NEL,NEL2,1,ALPHA,A(UUPOS),NEL,
     &               A(LPOS),NFRONT,ONE,A(LPOS+1_8),NFRONT)
#endif
        ENDIF
        RETURN
        END SUBROUTINE DMUMPS_FAC_MQ
      SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS( INODE, NFRONT, NASS, 
     &     CALL_UTRSM, A, LA, LAFAC, POSELT, IW, LIW, IOLDPS,
     &     MonBloc, MYID, NOFFW, LIWFAC,
     &     PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &     LNextPiv2beWritten, UNextPiv2beWritten, 
     &     PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &     
     &     XSIZE, SEUIL, UU, DKEEP, KEEP8, KEEP, IFLAG, 
     &     OOCWRITE_COMPATIBLE_WITH_BLR)
      USE DMUMPS_OOC 
      IMPLICIT NONE
      INTEGER, intent(in)    :: INODE, NFRONT, NASS,
     &                          LIW, MYID, XSIZE, IOLDPS, LIWFAC
      INTEGER(8), intent(in) :: LA, POSELT
      INTEGER, intent(inout) :: NOFFW,
     &                   PP_FIRST2SWAP_L, PP_FIRST2SWAP_U,
     &                   LNextPiv2beWritten, UNextPiv2beWritten,
     &                   PP_LastPIVRPTRFilled_L, PP_LastPIVRPTRFilled_U,
     &                   IFLAG
      LOGICAL, intent(in)    :: CALL_UTRSM
      INTEGER, intent(inout) :: IW(LIW)
      DOUBLE PRECISION, intent(inout) :: A(LA)
      DOUBLE PRECISION, intent(in)       :: SEUIL, UU, DKEEP(230)
      INTEGER, intent(in)    :: KEEP( 500 ) 
      INTEGER(8), intent(inout) :: LAFAC
      INTEGER(8)             :: KEEP8(150)
      TYPE(IO_BLOCK), intent(inout)  :: MonBloc
      LOGICAL, intent(in)    :: OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER  :: NPIV, NEL1, STRAT, TYPEFile, IFLAG_OOC, 
     &            IBEG_BLOCK, IFINB, INOPV
      INTEGER Inextpiv
      DOUBLE PRECISION :: MAXFROMN
      LOGICAL :: IS_MAXFROMN_AVAIL
      NPIV   = IW(IOLDPS+1+XSIZE)
      NEL1   = NFRONT - NASS
      IF (KEEP(206).GE.1) THEN
        Inextpiv = 1   
      ELSE 
        Inextpiv = 0   
      ENDIF
      IF ((NPIV.GT.0).AND.(NEL1.GT.0)) THEN
        IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN 
         STRAT          = STRAT_TRY_WRITE
         TYPEFile     = TYPEF_BOTH_LU
         MonBloc%LastPiv= NPIV
         CALL DMUMPS_FAC_P_PANEL(A(POSELT), LAFAC, NFRONT, 
     &      NPIV, NASS, IW(IOLDPS), LIWFAC, 
     &      MonBloc, TYPEFile, MYID, KEEP8,
     &      STRAT, IFLAG_OOC,
     &      LNextPiv2beWritten, UNextPiv2beWritten)
          IF (IFLAG_OOC < 0 ) IFLAG=IFLAG_OOC
        ELSE
          CALL DMUMPS_FAC_P(A,LA,NFRONT, NPIV, NASS, POSELT, 
     &                      CALL_UTRSM
     &                     )
        ENDIF
       ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        IBEG_BLOCK = NPIV
        IF (NASS.EQ.NPIV) GOTO 500
        IS_MAXFROMN_AVAIL = .FALSE.
 120    CALL DMUMPS_FAC_H(NFRONT,NASS,IW,LIW,A,LA,
     &     INOPV,NOFFW,IOLDPS,POSELT,UU,SEUIL,
     &     KEEP, DKEEP,
     &     PP_FIRST2SWAP_L,  MonBloc%LastPanelWritten_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U,  MonBloc%LastPanelWritten_U,
     &     PP_LastPIVRPTRFilled_U, MAXFROMN, IS_MAXFROMN_AVAIL,
     &     Inextpiv, OOCWRITE_COMPATIBLE_WITH_BLR
     &     )
        IF (INOPV.NE.1) THEN
         CALL DMUMPS_FAC_N(NFRONT,NASS,IW,LIW,A,LA,
     &                 IOLDPS,POSELT,IFINB,XSIZE,
     &                 KEEP, MAXFROMN, IS_MAXFROMN_AVAIL)
         IW(IOLDPS+1+XSIZE) = IW(IOLDPS+1+XSIZE) + 1
         IF (IFINB.EQ.0) GOTO 120
        ENDIF
        NPIV   = IW(IOLDPS+1+XSIZE)
        NEL1   = NFRONT - NASS
        IF ((NPIV.LE.IBEG_BLOCK).OR.(NEL1.EQ.0)) GO TO 500
        CALL DMUMPS_FAC_T(A,LA,IBEG_BLOCK,
     &                NFRONT,NPIV,NASS,POSELT)
 500  CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_FR_UPDATE_CBROWS
      SUBROUTINE DMUMPS_FAC_I(NFRONT,NASS,LAST_ROW,
     &    IBEG_BLOCK, IEND_BLOCK,
     &    N,INODE,IW,LIW,A,LA,
     &    INOPV,NOFFW,IFLAG,IOLDPS,POSELT,UU,SEUIL,KEEP,KEEP8,
     &     DKEEP,PIVNUL_LIST,LPN_LIST,
     &
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &     PP_LastPIVRPTRFilled_L,
     &     PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &     PP_LastPIVRPTRFilled_U,
     &     PIVOT_OPTION, LR_ACTIVATED, IEND_BLR, Inextpiv,
     &     OOCWRITE_COMPATIBLE_WITH_BLR,
     &     TIPIV    
     &     )
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER, intent(in)    :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(inout), OPTIONAL :: TIPIV(:)
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in)    :: NFRONT,NASS,N,LIW,INODE,LAST_ROW
      INTEGER, intent(inout) :: IFLAG,INOPV,NOFFW
      DOUBLE PRECISION, intent(in)       :: UU, SEUIL
      INTEGER, intent(inout) :: IW(LIW)
      INTEGER, intent(in)    :: IOLDPS
      INTEGER(8), intent(in) :: POSELT
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, intent(in)    :: LPN_LIST
      INTEGER, intent(inout) :: PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION    DKEEP(230)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk_L,
     &        PP_LastPIVRPTRFilled_L,
     &        PP_FIRST2SWAP_U, PP_LastPanelonDisk_U,
     &        PP_LastPIVRPTRFilled_U
      INTEGER, intent(in)    :: PIVOT_OPTION, IEND_BLR
      LOGICAL, intent(in)    :: LR_ACTIVATED
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOCWRITE_COMPATIBLE_WITH_BLR
      INCLUDE 'mumps_headers.h'
      DOUBLE PRECISION SWOP
      INTEGER XSIZE
      INTEGER(8) :: APOS, IDIAG
      INTEGER(8) :: J1, J2, JJ, J3
      INTEGER(8) :: NFRONT8
      INTEGER ILOC
      DOUBLE PRECISION ZERO
      PARAMETER( ZERO = 0.0D0 )
      DOUBLE PRECISION RZERO, RMAX, AMROW
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      INTEGER NPIV,IPIV
      INTEGER NPIVP1,JMAX,J,ISW,ISWPS1
      INTEGER ISWPS2,KSW, HF
      INTEGER DMUMPS_IXAMAX
      INTEGER :: ISHIFT, K206
      INTEGER :: IPIV_SHIFT,IPIV_END
      INTRINSIC max
      DATA RZERO /0.0D0/
#if defined(_OPENMP)
      INTEGER :: NOMP,CHUNK,K361
#endif
      INTEGER I_PIVRPTR_L, I_PIVR_L, NBPANELS_L
      INTEGER I_PIVRPTR_U, I_PIVR_U, NBPANELS_U
#if defined(_OPENMP)
        NOMP    = OMP_GET_MAX_THREADS()
        K361    = KEEP(361)
#endif
        PIVNUL  = DKEEP(1)
        FIXA    = DKEEP(2)
        CSEUIL  = SEUIL
        NFRONT8 = int(NFRONT,8)
        K206    = KEEP(206)
        XSIZE   = KEEP(IXSZ)
        NPIV    = IW(IOLDPS+1+XSIZE)
        HF = 6 + IW(IOLDPS+5+XSIZE)+XSIZE
        NPIVP1  = NPIV + 1
        IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR_L, I_PIVR_L, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
          CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_U, NBPANELS_U, 
     &       I_PIVRPTR_U, I_PIVR_U, 
     &       IOLDPS+2*NFRONT+6+IW(IOLDPS+5+XSIZE)+XSIZE,
     &       IW, LIW)
        ENDIF
        IF ( present(TIPIV) ) THEN
          ILOC    = NPIVP1 - IBEG_BLOCK + 1
          TIPIV(ILOC) = ILOC
        ENDIF
        IF (INOPV .EQ. -1) THEN
           APOS = POSELT + NFRONT8*int(NPIVP1-1,8) + int(NPIV,8)
           IDIAG = APOS
           IF(abs(A(APOS)).LT.SEUIL) THEN
              IF (dble(A(APOS)) .GE. RZERO) THEN
                 A(APOS) = CSEUIL
              ELSE
                 A(APOS) = -CSEUIL
              ENDIF
              KEEP(98) = KEEP(98)+1
           ELSE IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DKEEP(6), KEEP(259))
           ENDIF
           IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
             IF (KEEP(251).EQ.0) THEN 
               CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
             ENDIF
             CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIVP1, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
           ENDIF
           GO TO 420
        ENDIF
        INOPV   = 0
      ISHIFT   = 0            
      IPIV_END = IEND_BLOCK   
      IF (K206.GE.1) THEN
        IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.IEND_BLOCK) THEN
          ISHIFT = Inextpiv - NPIVP1
        ENDIF
        IF ( K206.EQ.1
     &      .OR.  (K206 .GT.1 .AND. IEND_BLOCK.EQ.IEND_BLR) ) THEN
          IPIV_END = IEND_BLOCK + ISHIFT
        ENDIF
      ENDIF
          DO 460 IPIV_SHIFT = NPIVP1+ISHIFT, IPIV_END
            IF (IPIV_SHIFT .LE. IEND_BLOCK) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV = IPIV_SHIFT-IEND_BLOCK-1+NPIVP1
              IF (IBEG_BLOCK.EQ.NPIVP1) THEN
                EXIT
              ENDIF
            ENDIF
            APOS = POSELT + NFRONT8*int(IPIV-1,8) + int(NPIV,8)
            JMAX = 1
            IF (PIVOT_OPTION.GT.0.AND.UU.GT.RZERO) GO TO 340
            IF (A(APOS).EQ.ZERO) GO TO 630 
            GO TO 380 
  340       CONTINUE
            AMROW = RZERO
            J1 = APOS
            IF (PIVOT_OPTION.EQ.1 .OR. 
     &         (LR_ACTIVATED.AND.KEEP(480).GE.2)) THEN
              J = IEND_BLR - NPIV
            ELSE
              J = NASS - NPIV
            ENDIF
            J2 = J1 + J - 1_8
            JMAX  = DMUMPS_IXAMAX(J,A(J1),1,KEEP(361))
            JJ    = J1 + int(JMAX - 1,8)
            AMROW = abs(A(JJ))
            RMAX = AMROW
            IF (PIVOT_OPTION.GE.2) THEN
              J1 = J2 + 1_8
              IF (PIVOT_OPTION.GE.3) THEN
                J2 = APOS +int(- NPIV + NFRONT - 1 - KEEP(253),8)
              ELSE
                J2 = APOS +int(- NPIV + NASS - 1 - KEEP(253),8)
              ENDIF
              IF (J2.LT.J1) GO TO 370
              IF (KEEP(351).EQ.1) THEN
!$              CHUNK = max(K361/2,(int(J2-J1)+NOMP-1)/NOMP)
!$OMP  PARALLEL DO schedule(static, CHUNK) PRIVATE(JJ)
!$OMP& FIRSTPRIVATE(J1,J2)
!$OMP& REDUCTION(max:RMAX) IF ((J2-J1).GE.K361)
                DO JJ=J1,J2
                  RMAX = max(abs(A(JJ)),RMAX)
                ENDDO
!$OMP  END PARALLEL DO
              ELSE
                DO 360 JJ=J1,J2
                  RMAX = max(abs(A(JJ)),RMAX)
  360           CONTINUE
              ENDIF
            ENDIF
  370       IDIAG = APOS + int(IPIV - NPIVP1,8)
            IF ( RMAX .LE. PIVNUL ) THEN
               IF (NFRONT - KEEP(253) .EQ. NASS) THEN 
                 IF (IEND_BLOCK.NE.NASS ) THEN 
                   GOTO 460
                 ENDIF
                 J1=POSELT+int(IPIV-1,8)+int(NPIV,8)*NFRONT8
                 J2=POSELT+int(IPIV-1,8)+int(LAST_ROW-1,8)*NFRONT8
               ELSE
                 J1=POSELT+int(IPIV-1,8)
                 J2=POSELT+int(IPIV-1,8)+int(LAST_ROW-1,8)*NFRONT8
               ENDIF
               DO JJ=J1, J2, NFRONT8
                 IF ( abs(A(JJ)) .GT. PIVNUL ) THEN
                   GOTO 460
                 END IF
               ENDDO
               KEEP(109) = KEEP(109)+1
               ISW = IOLDPS+HF+
     &                      IW(IOLDPS+1+XSIZE)+IPIV-NPIVP1
               PIVNUL_LIST(KEEP(109)) = IW(ISW)
               IF(dble(FIXA).GT.RZERO) THEN
                  IF(dble(A(IDIAG)) .GE. RZERO) THEN
                     A(IDIAG) = FIXA
                  ELSE
                     A(IDIAG) = -FIXA
                  ENDIF
               ELSE
                 J1 = APOS
                 J2 = APOS +int(- NPIV + NFRONT - 1 - KEEP(253),8)
                 DO JJ=J1,J2
                   A(JJ) = ZERO
                 ENDDO
                 A(IDIAG) = -FIXA
               ENDIF
               JMAX = IPIV - NPIV
               GOTO 385   
            ENDIF
            IF (abs(A(IDIAG)) .GE. UU*RMAX .AND.
     &          abs(A(IDIAG)) .GT. max(SEUIL,tiny(RMAX))) THEN
               JMAX = IPIV - NPIV
               GO TO 380
            ENDIF
            IF ( .NOT. (AMROW .GE. UU*RMAX .AND.
     &                  AMROW .GT. max(SEUIL,tiny(RMAX))) ) GO TO 460
            NOFFW = NOFFW + 1
  380       CONTINUE
          IF (K206.GE.1) THEN
             Inextpiv = IPIV + 1
          ENDIF
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER( A(APOS+int(JMAX-1,8)),
     &                                 DKEEP(6),
     &                                 KEEP(259))
            ENDIF
  385       CONTINUE
            IF (IPIV.EQ.NPIVP1) GO TO 400
            KEEP(260)=-KEEP(260)
            J1 = POSELT + int(NPIV,8)*NFRONT8
            J2 = J1 + NFRONT8 - 1_8
            J3 = POSELT + int(IPIV-1,8)*NFRONT8
            DO 390 JJ=J1,J2
              SWOP = A(JJ)
              A(JJ) = A(J3)
              A(J3) = SWOP
              J3 = J3 + 1_8
  390       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NPIVP1
            ISWPS2 = IOLDPS + HF - 1 + IPIV
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
  400       IF (JMAX.EQ.1) GO TO 420
            KEEP(260)=-KEEP(260)
            IF ( present(TIPIV) ) THEN
              TIPIV(ILOC) = ILOC + JMAX - 1
            ENDIF
            J1 = POSELT + int(NPIV,8)
            J2 = POSELT + int(NPIV + JMAX - 1,8)
            DO 410 KSW=1,LAST_ROW
              SWOP = A(J1)
              A(J1) = A(J2)
              A(J2) = SWOP
              J1 = J1 + NFRONT8
              J2 = J2 + NFRONT8
  410       CONTINUE
            ISWPS1 = IOLDPS + HF - 1 + NFRONT + NPIV + 1
            ISWPS2 = IOLDPS + HF - 1 + NFRONT + NPIV + JMAX
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            GO TO 420  
  460     CONTINUE
          IF (K206 .GE. 1) THEN
            Inextpiv=IEND_BLOCK+1
          ENDIF
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 430
  630 CONTINUE
      IFLAG = -10
      GOTO 430
  420 CONTINUE
              IF (KEEP(201).EQ.1.AND.OOCWRITE_COMPATIBLE_WITH_BLR) THEN
                IF (KEEP(251).EQ.0) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_L), 
     &               NBPANELS_L,
     &               IW(I_PIVR_L), NASS, NPIVP1, IPIV, 
     &               PP_LastPanelonDisk_L,
     &               PP_LastPIVRPTRFilled_L)
                ENDIF
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR_U), 
     &               NBPANELS_U,
     &               IW(I_PIVR_U), NASS, NPIVP1, NPIV+JMAX, 
     &               PP_LastPanelonDisk_U,
     &               PP_LastPIVRPTRFilled_U)
              ENDIF
  430 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_I
      SUBROUTINE DMUMPS_FAC_I_LDLT 
     &   (NFRONT,NASS,INODE,IBEG_BLOCK,IEND_BLOCK,
     &    IW,LIW, A,LA, INOPV,
     &    NNEG,
     &    IFLAG,IOLDPS,POSELT,UU, SEUIL,KEEP,KEEP8,PIVSIZ,
     &     DKEEP,PIVNUL_LIST,LPN_LIST, XSIZE,
     &     PP_FIRST2SWAP_L, PP_LastPanelonDisk,
     &     PP_LastPIVRPTRIndexFilled,MAXFROMM,IS_MAXFROMM_AVAIL, 
     &     PIVOT_OPTION, IEND_BLR, Inextpiv, 
     &     OOCWRITE_COMPATIBLE_WITH_BLR)
!$    USE OMP_LIB
      USE MUMPS_OOC_COMMON 
      IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER NFRONT,NASS,LIW,INODE,IFLAG,INOPV,
     &        IOLDPS, NNEG
      INTEGER, intent(in) :: IBEG_BLOCK, IEND_BLOCK
      INTEGER, intent(in)    :: PIVOT_OPTION,IEND_BLR
      INTEGER, intent(inout) :: Inextpiv
      LOGICAL, intent(in)    :: OOCWRITE_COMPATIBLE_WITH_BLR
      INTEGER PIVSIZ,LPIV, XSIZE
      DOUBLE PRECISION A(LA) 
      DOUBLE PRECISION UU, UULOC, SEUIL
      INTEGER IW(LIW)
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER LPN_LIST
      INTEGER PIVNUL_LIST(LPN_LIST)
      DOUBLE PRECISION DKEEP(230)
      INTEGER PP_FIRST2SWAP_L, PP_LastPanelonDisk
      INTEGER PP_LastPIVRPTRIndexFilled
      DOUBLE PRECISION, intent(in) :: MAXFROMM
      LOGICAL, intent(inout) :: IS_MAXFROMM_AVAIL
      include 'mpif.h'
      INTEGER (8) :: POSPV1,POSPV2,OFFDAG,APOSJ
      INTEGER JMAX, LIM
      DOUBLE PRECISION RMAX,AMAX,TMAX
      DOUBLE PRECISION MAXPIV
      DOUBLE PRECISION PIVNUL
      DOUBLE PRECISION FIXA, CSEUIL
      DOUBLE PRECISION PIVOT,DETPIV
      INCLUDE 'mumps_headers.h'
      INTEGER :: J
      INTEGER(8) :: APOS, J1, J2, JJ, NFRONT8, KK, J1_ini, JJ_ini
      INTEGER    :: LDA
      INTEGER(8) :: LDA8
      INTEGER NPIV,IPIV
      INTEGER NPIVP1,K 
      INTEGER :: ISHIFT, K206, IPIV_SHIFT, IPIV_END
      INTRINSIC max
      DOUBLE PRECISION ZERO, ONE
      PARAMETER( ZERO = 0.0D0 )
      PARAMETER( ONE = 1.0D0 )
      DOUBLE PRECISION RZERO,RONE
      PARAMETER(RZERO=0.0D0, RONE=1.0D0)
#if defined(_OPENMP)
      LOGICAL :: OMP_FLAG
      INTEGER :: NOMP, CHUNK, J1_end
#endif
      INTEGER I_PIVRPTR, I_PIVR, NBPANELS_L
!$    NOMP = OMP_GET_MAX_THREADS()
      PIVNUL = DKEEP(1)
      FIXA   = DKEEP(2)
      CSEUIL = SEUIL
      LDA     = NFRONT
      LDA8    = int(LDA,8)
      NFRONT8 = int(NFRONT,8)
      K206    = KEEP(206)
      IF (KEEP(201).EQ.1 .AND. KEEP(50).NE.1 .AND.
     &   OOCWRITE_COMPATIBLE_WITH_BLR ) THEN
             CALL DMUMPS_GET_OOC_PERM_PTR(TYPEF_L, NBPANELS_L, 
     &       I_PIVRPTR, I_PIVR, IOLDPS+2*NFRONT+6+KEEP(IXSZ),
     &       IW, LIW)
      ENDIF
      UULOC = UU
      PIVSIZ = 1
      NPIV    = IW(IOLDPS+1+XSIZE)
      NPIVP1  = NPIV + 1
      IF(INOPV .EQ. -1) THEN
         APOS = POSELT + (LDA8+1_8) * int(NPIV,8)
         POSPV1 = APOS
         IF(abs(A(APOS)).LT.SEUIL) THEN
            IF(dble(A(APOS)) .GE. RZERO) THEN
               A(APOS) = CSEUIL
            ELSE
               A(APOS) = -CSEUIL
               NNEG = NNEG+1
            ENDIF
            KEEP(98) = KEEP(98)+1
         ELSE IF (KEEP(258) .NE. 0) THEN
            CALL DMUMPS_UPDATEDETER( A(APOS), DKEEP(6), KEEP(259) )
         ENDIF
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1 .AND.
     &            OOCWRITE_COMPATIBLE_WITH_BLR) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, NPIVP1,
     &               PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
         GO TO 420
      ENDIF
      INOPV   = 0
      ISHIFT = 0              
      IPIV_END = IEND_BLOCK   
      IF (K206.GE.1) THEN
        IF (Inextpiv.GT.NPIVP1.AND.Inextpiv.LE.IEND_BLOCK) THEN
          ISHIFT = Inextpiv - NPIVP1
        ENDIF
        IF ( K206.EQ.1
     &      .OR.  (K206 .GT.1 .AND. IEND_BLOCK.EQ.IEND_BLR) ) THEN
          IPIV_END = IEND_BLOCK + ISHIFT
        ENDIF
        IF (ISHIFT.GT.0.AND.IS_MAXFROMM_AVAIL) THEN
           IPIV = NPIVP1
           APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
           POSPV1 = APOS + int(IPIV - NPIVP1,8)
           PIVOT = A(POSPV1)
           IF ( MAXFROMM .GT. PIVNUL ) THEN
               IF ( abs(PIVOT) .GE. UULOC*MAXFROMM
     &         .AND. abs(PIVOT) .GT. max(SEUIL,tiny(MAXFROMM)) ) THEN
                 ISHIFT = 0
               ENDIF
           ENDIF
        ENDIF  
        IF ( ISHIFT .GT. 0) THEN
          IS_MAXFROMM_AVAIL = .FALSE.
        ENDIF
      ENDIF  
       DO 460 IPIV_SHIFT = NPIVP1+ISHIFT, IPIV_END
            IF (IPIV_SHIFT .LE. IEND_BLOCK) THEN
              IPIV=IPIV_SHIFT
            ELSE
              IPIV = IPIV_SHIFT-IEND_BLOCK-1+NPIVP1
              IF (IBEG_BLOCK.EQ.NPIVP1) THEN
                EXIT
              ENDIF
            ENDIF
         APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIV,8)
         POSPV1 = APOS + int(IPIV - NPIVP1,8)
         PIVOT = A(POSPV1)
         IF (UULOC.EQ.RZERO.OR.PIVOT_OPTION.EQ.0) THEN 
            IF (abs(A(APOS)).EQ.RZERO) GO TO 630
            IF (A(APOS).LT.RZERO) NNEG = NNEG+1
            IF (KEEP(258) .NE. 0) THEN
              CALL DMUMPS_UPDATEDETER(A(APOS), DKEEP(6), KEEP(259))
            ENDIF
            GO TO 420
         ENDIF
         IF ( IS_MAXFROMM_AVAIL ) THEN
            IF ( MAXFROMM .GT. PIVNUL ) THEN
               IF ( abs(PIVOT) .GE. UULOC*MAXFROMM
     &         .AND. abs(PIVOT) .GT. max(SEUIL,tiny(MAXFROMM)) ) THEN
                 IF (PIVOT .LT. RZERO) NNEG = NNEG+1
                 IF (KEEP(258) .NE. 0) THEN
                   CALL DMUMPS_UPDATEDETER(PIVOT, DKEEP(6), KEEP(259))
                 ENDIF
                 GOTO 415
               ENDIF
            ENDIF
            IS_MAXFROMM_AVAIL = .FALSE.
         ENDIF
         AMAX = -RONE
         JMAX = 0
         J1 = APOS
         J2 = POSPV1 - 1_8
         DO JJ=J1,J2
            IF(abs(A(JJ)) .GT. AMAX) THEN
               AMAX = abs(A(JJ))
               JMAX = IPIV - int(POSPV1-JJ)
            ENDIF
         ENDDO
         J1 = POSPV1 + LDA8
         DO J=1, IEND_BLOCK - IPIV
            IF(abs(A(J1)) .GT. AMAX) THEN
               AMAX = abs(A(J1))
               JMAX = IPIV + J
            ENDIF
            J1 = J1 + LDA8
         ENDDO
           RMAX = RZERO
           IF (PIVOT_OPTION.EQ.3) THEN
             LIM = NFRONT
           ELSEIF (PIVOT_OPTION.EQ.2) THEN
             LIM = NASS
           ELSEIF (PIVOT_OPTION.EQ.1) THEN
             LIM = IEND_BLR
           ELSE
             write(*,*) 'Internal error in FAC_I_LDLT: PIVOT_OPTION=',
     &                  PIVOT_OPTION
             CALL MUMPS_ABORT()
           ENDIF
           J1_ini = J1
#if defined(_OPENMP)
           J1_end = LIM - KEEP(253) - IEND_BLOCK
           CHUNK = max(J1_end,1)
           IF ( J1_end.GE.KEEP(360)) THEN
             OMP_FLAG = .TRUE.
             CHUNK = max(KEEP(360)/2,(J1_end+NOMP-1)/NOMP)
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
#endif
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) PRIVATE(J1)
!$OMP& REDUCTION(max:RMAX) IF(OMP_FLAG)
           DO J=1, LIM - KEEP(253) - IEND_BLOCK
              J1 = J1_ini + int(J-1,8) * LDA8
              RMAX = max(abs(A(J1)),RMAX)
           ENDDO
!$OMP END PARALLEL DO
         IF (max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL) THEN
            KEEP(109) = KEEP(109)+1
            PIVNUL_LIST(KEEP(109)) = -1
            IF(dble(FIXA).GT.RZERO) THEN
               IF(dble(PIVOT) .GE. RZERO) THEN
                  A(POSPV1) = FIXA
               ELSE
                  A(POSPV1) = -FIXA
               ENDIF
            ELSE
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  A(JJ) = ZERO
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1, IEND_BLOCK - IPIV
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               DO J=1,NFRONT - IEND_BLOCK
                  A(J1) = ZERO
                  J1 = J1 + LDA8
               ENDDO
               A(POSPV1) = ONE
            ENDIF
            PIVOT = A(POSPV1)
            GO TO 415
         ENDIF
         IF ( abs(PIVOT).GE.UULOC*max(RMAX,AMAX)
     &        .AND. abs(PIVOT) .GT. max(SEUIL,tiny(RMAX)) ) THEN
               IF (PIVOT .LT. ZERO) NNEG = NNEG+1
               IF (KEEP(258) .NE.0 ) THEN
                 CALL DMUMPS_UPDATEDETER(PIVOT, DKEEP(6), KEEP(259))
               ENDIF
               GO TO 415
         END IF
         IF (NPIVP1.EQ.IEND_BLOCK) THEN
           GOTO 460
         ELSE IF (JMAX.EQ.0) THEN
           GOTO 460
         ENDIF
         IF (max(abs(PIVOT),RMAX,AMAX).LE.tiny(RMAX)) THEN
           GOTO 460
         ENDIF
         IF (
     &   (KEEP(19).NE.0).AND.(max(AMAX,RMAX,abs(PIVOT)).LE.SEUIL)
     &      )
     &       THEN
           GO TO 460
         ENDIF
         IF (RMAX.LT.AMAX) THEN
               J1 = APOS
               J2 = POSPV1 - 1_8
               DO JJ=J1,J2
                  IF(int(POSPV1-JJ) .NE. IPIV-JMAX) THEN
                     RMAX = max(RMAX,abs(A(JJ)))
                  ENDIF
               ENDDO
               J1 = POSPV1 + LDA8
               DO J=1,NASS-IPIV
                  IF(IPIV+J .NE. JMAX) THEN
                     RMAX = max(abs(A(J1)),RMAX)
                  ENDIF
                  J1 = J1 + LDA8
               ENDDO
           ENDIF
           IF (PIVOT_OPTION.EQ.3) THEN
             LIM = NFRONT
           ELSEIF (PIVOT_OPTION.EQ.2) THEN
             LIM = NASS
           ELSEIF (PIVOT_OPTION.EQ.1) THEN
             LIM = IEND_BLR
           ELSE
             write(*,*) 'Internal error in FAC_I_LDLT: PIVOT_OPTION=',
     &                  PIVOT_OPTION
             CALL MUMPS_ABORT()
           ENDIF
           APOSJ = POSELT + int(JMAX-1,8)*LDA8 + int(NPIV,8)
           POSPV2 = APOSJ + int(JMAX - NPIVP1,8)
           IF (IPIV.LT.JMAX) THEN
              OFFDAG = APOSJ + int(IPIV - NPIVP1,8)
           ELSE
              OFFDAG = APOS + int(JMAX - NPIVP1,8)
           END IF
           TMAX = RZERO
#if defined(_OPENMP)
           J1_end = LIM-JMAX-KEEP(253)
           CHUNK = max(J1_end,1)
           IF (J1_end.GE.KEEP(360)) THEN
             OMP_FLAG = .TRUE.
             CHUNK = max(KEEP(360)/2,(J1_end+NOMP-1)/NOMP)
           ELSE
             OMP_FLAG = .FALSE.
           ENDIF
#endif
           IF (JMAX .LT. IPIV) THEN
              JJ_ini = POSPV2
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) IF (OMP_FLAG)
!$OMP& PRIVATE(JJ) REDUCTION(max:TMAX)
              DO K = 1, LIM - JMAX - KEEP(253)
                 JJ = JJ_ini+ int(K,8)*NFRONT8
                 IF (JMAX+K.NE.IPIV) THEN
                    TMAX=max(TMAX,abs(A(JJ)))
                 ENDIF
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2-1_8
                 TMAX = max(TMAX,abs(A(KK)))
              ENDDO
           ELSE
              JJ_ini = POSPV2
!$OMP PARALLEL DO SCHEDULE(static, CHUNK) PRIVATE(JJ) 
!$OMP& REDUCTION(max:TMAX) IF(OMP_FLAG)
              DO K = 1, LIM-JMAX-KEEP(253)
                 JJ = JJ_ini + int(K,8)*NFRONT8
                 TMAX=max(TMAX,abs(A(JJ)))
              ENDDO
!$OMP END PARALLEL DO
              DO KK =  APOSJ, POSPV2 - 1_8
                 IF (KK.NE.OFFDAG) THEN
                    TMAX = max(TMAX,abs(A(KK)))
                 ENDIF
              ENDDO
           ENDIF
           DETPIV = A(POSPV1)*A(POSPV2) - A(OFFDAG)**2
           IF (SEUIL.GT.RZERO) THEN
                IF (sqrt(abs(DETPIV)) .LE. SEUIL ) THEN
                  GOTO 460
                ENDIF
           ENDIF
           MAXPIV = max(abs(A(POSPV1)),abs(A(POSPV2)))
           IF (MAXPIV.EQ.RZERO) MAXPIV = RONE
           IF ((abs(A(POSPV2))*RMAX+AMAX*TMAX)*UULOC.GT.
     &          abs(DETPIV) .OR. abs(DETPIV) .EQ. RZERO) THEN
             GO TO 460
           ENDIF
           IF ((abs(A(POSPV1))*TMAX+AMAX*RMAX)*UULOC.GT.
     &          abs(DETPIV) .OR. abs(DETPIV) .EQ. RZERO) THEN
             GO TO 460
           ENDIF
           IF (KEEP(258) .NE.0 ) THEN
             CALL DMUMPS_UPDATEDETER(DETPIV, DKEEP(6), KEEP(259))
           ENDIF
           PIVSIZ = 2
           KEEP(103) = KEEP(103)+1
           IF(DETPIV .LT. RZERO) THEN
             NNEG = NNEG+1
           ELSE IF(A(POSPV2) .LT. RZERO) THEN
             NNEG = NNEG+2
           ENDIF
 415       CONTINUE
          IF (K206.GE.1) THEN
             Inextpiv = max(NPIVP1+PIVSIZ, IPIV+1)
          ENDIF
           DO K=1,PIVSIZ
              IF (PIVSIZ .EQ. 2) THEN
                IF (K==1) THEN
                  LPIV = min(IPIV,JMAX)
                ELSE
                  LPIV   = max(IPIV,JMAX)
                ENDIF
              ELSE
                LPIV = IPIV
              ENDIF
              IF (LPIV.EQ.NPIVP1) GOTO 416 
              CALL DMUMPS_SWAP_LDLT( A, LA, IW, LIW,
     &             IOLDPS, NPIVP1, LPIV, POSELT, NASS,
     &             LDA, NFRONT, 1, KEEP(219), KEEP(50),
     &             KEEP(IXSZ), -9999)
 416          CONTINUE
              IF (KEEP(201).EQ.1.AND.KEEP(50).NE.1 .AND.
     &            OOCWRITE_COMPATIBLE_WITH_BLR ) THEN
                CALL DMUMPS_STORE_PERMINFO( IW(I_PIVRPTR), NBPANELS_L,
     &               IW(I_PIVR), NASS, NPIVP1, LPIV, PP_LastPanelonDisk,
     &               PP_LastPIVRPTRIndexFilled)
              ENDIF
              NPIVP1 = NPIVP1 + 1
           ENDDO
           IF(PIVSIZ .EQ. 2) THEN
              A(POSELT+(LDA8+1_8)*int(NPIV,8)+1_8) = DETPIV
           ENDIF
           GOTO 420
  460   CONTINUE
          IF (K206 .GE. 1) THEN
            Inextpiv=IEND_BLOCK+1
          ENDIF
      IF (IEND_BLOCK.EQ.NASS) THEN
       INOPV = 1
      ELSE
       INOPV = 2
      ENDIF
      GO TO 420
  630 CONTINUE
      PIVSIZ = 0
      IFLAG = -10
  420 CONTINUE
      IS_MAXFROMM_AVAIL = .FALSE.
      RETURN
      END SUBROUTINE DMUMPS_FAC_I_LDLT
      SUBROUTINE DMUMPS_FAC_MQ_LDLT(IEND_BLOCK,
     &     NFRONT,NASS,NPIV,INODE,
     &     A,LA,LDA, 
     &     POSELT,IFINB,PIVSIZ,
     &     MAXFROMM, IS_MAXFROMM_AVAIL, IS_MAX_USEFUL,
     &     KEEP253, PIVOT_OPTION, IEND_BLR
     &     )
      IMPLICIT NONE
      INTEGER, intent(out):: IFINB
      INTEGER, intent(in) :: INODE, NFRONT, NASS, NPIV
      INTEGER, intent(in) :: IEND_BLOCK
      INTEGER, intent(in) :: LDA
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in)    :: PIVOT_OPTION, IEND_BLR
      INTEGER(8) :: POSELT
      DOUBLE PRECISION, intent(out) :: MAXFROMM
      LOGICAL, intent(out) :: IS_MAXFROMM_AVAIL
      LOGICAL, intent(in) :: IS_MAX_USEFUL
      INTEGER, INTENT(in) :: KEEP253
      DOUBLE PRECISION    VALPIV
      DOUBLE PRECISION :: MAXFROMMTMP
      INTEGER  NCB1
      INTEGER(8) :: NFRONT8
      INTEGER(8) :: LDA8
      INTEGER(8) :: K1POS
      INTEGER NEL2, NEL, LIM
      DOUBLE PRECISION ONE, ZERO
      DOUBLE PRECISION A11,A22,A12
      INTEGER(8) :: APOS, LPOS, LPOS1, LPOS2
      INTEGER(8) :: POSPV1, POSPV2
      INTEGER PIVSIZ,NPIV_NEW,J2,I
      INTEGER(8) :: OFFDAG, OFFDAG_OLD, IBEG, IEND
      INTEGER(8) :: JJ, K1, K2, IROW
      DOUBLE PRECISION SWOP,DETPIV,MULT1,MULT2
      INCLUDE 'mumps_headers.h'
      PARAMETER(ONE  = 1.0D0,
     &          ZERO = 0.0D0)
      LDA8     = int(LDA,8)
      NFRONT8  = int(NFRONT,8)
      NPIV_NEW = NPIV + PIVSIZ
      NEL      = NFRONT - NPIV_NEW
      IFINB    = 0
      IS_MAXFROMM_AVAIL = .FALSE.
      NEL2   = IEND_BLOCK - NPIV_NEW
      IF (NEL2.EQ.0) THEN
        IF (IEND_BLOCK.EQ.NASS) THEN
          IFINB        = -1
        ELSE
          IFINB        = 1
        ENDIF
      ENDIF
      IF(PIVSIZ .EQ. 1) THEN
         APOS   = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         VALPIV = ONE/A(APOS)
         LPOS   = APOS + LDA8
         MAXFROMM = 0.0D00
         IF (NEL2 > 0) THEN
           IF (.NOT. IS_MAX_USEFUL) THEN
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               DO JJ=1_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ELSE
             IS_MAXFROMM_AVAIL = .TRUE.
             DO I=1, NEL2
               K1POS = LPOS + int(I-1,8)*LDA8
               A(APOS+int(I,8))=A(K1POS)
               A(K1POS) = A(K1POS) * VALPIV
               A(K1POS+1_8)=A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMM=max( MAXFROMM,abs(A(K1POS+1_8)) )
               DO JJ = 2_8, int(I,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDDO
           ENDIF
         ENDIF
         IF (PIVOT_OPTION.EQ.3) THEN
           LIM = NFRONT
         ELSEIF (PIVOT_OPTION.EQ.2) THEN
           LIM = NASS
         ELSE
           LIM = IEND_BLR
         ENDIF
         NCB1 = LIM - IEND_BLOCK
         IF (.NOT. IS_MAX_USEFUL) THEN
!$OMP      PARALLEL DO PRIVATE(JJ,K1POS) IF (NCB1 > 300)
           DO I=NEL2+1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
!$OMP      END PARALLEL DO
         ELSE
           MAXFROMMTMP=0.0D0
!$OMP      PARALLEL DO PRIVATE(JJ,K1POS)
!$OMP&     REDUCTION(max:MAXFROMMTMP) IF (NCB1 - KEEP253 > 300)
           DO I=NEL2+1, NEL2 + NCB1 - KEEP253
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             IF (NEL2 > 0) THEN
               A(K1POS+1_8) = A(K1POS+1_8) - A(K1POS) * A(APOS+1_8)
               MAXFROMMTMP=max(MAXFROMMTMP, abs(A(K1POS+1_8)))
               DO JJ = 2_8, int(NEL2,8)
                 A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
               ENDDO
             ENDIF
           ENDDO
!$OMP      END PARALLEL DO
           DO I = NEL2 + NCB1 - KEEP253 + 1, NEL2 + NCB1
             K1POS = LPOS+ int(I-1,8)*LDA8
             A(APOS+int(I,8))=A(K1POS)
             A(K1POS) = A(K1POS) * VALPIV
             DO JJ = 1_8, int(NEL2,8)
               A(K1POS+JJ)=A(K1POS+JJ) - A(K1POS) * A(APOS+JJ)
             ENDDO
           ENDDO
           MAXFROMM=max(MAXFROMM, MAXFROMMTMP)
         ENDIF
      ELSE
         IF (PIVOT_OPTION.EQ.3) THEN
           LIM = NFRONT
         ELSEIF (PIVOT_OPTION.EQ.2) THEN
           LIM = NASS
         ELSE
           LIM = IEND_BLR
         ENDIF
         POSPV1 = POSELT + int(NPIV,8)*(NFRONT8 + 1_8)
         POSPV2 = POSPV1 + NFRONT8 + 1_8
         OFFDAG_OLD = POSPV2 - 1_8
         OFFDAG = POSPV1 + 1_8
         SWOP = A(POSPV2)
         DETPIV = A(OFFDAG)
          A22 = A(POSPV1)/DETPIV   
          A11 =  SWOP/DETPIV       
          A12 = -A(OFFDAG_OLD)/DETPIV   
          A(OFFDAG)     = A(OFFDAG_OLD)  
          A(OFFDAG_OLD) = ZERO
         LPOS1   = POSPV2 + LDA8 - 1_8
         LPOS2   = LPOS1 + 1_8
         CALL dcopy(LIM-NPIV_NEW, A(LPOS1), LDA, A(POSPV1+2_8), 1)
         CALL dcopy(LIM-NPIV_NEW, A(LPOS2), LDA, A(POSPV2+1_8), 1)
         JJ = POSPV2 + NFRONT8-1_8  
         IBEG = JJ + 2_8
         IEND = IBEG
         DO J2 = 1,NEL2
           K1 = JJ
           K2 = JJ+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG, IEND
              A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
              K1 = K1 + 1_8
              K2 = K2 + 1_8
           ENDDO
           A( JJ       ) = -MULT1
           A( JJ + 1_8 ) = -MULT2
           IBEG = IBEG + NFRONT8
           IEND = IEND + NFRONT8 + 1_8
           JJ = JJ+NFRONT8
         ENDDO
         IEND = IEND-1_8
         DO J2 = IEND_BLOCK+1,LIM
           K1 = JJ
           K2 = JJ+1_8
           MULT1 = - (A11*A(K1)+A12*A(K2))
           MULT2 = - (A12*A(K1)+A22*A(K2))
           K1 = POSPV1 + 2_8
           K2 = POSPV2 + 1_8
           DO IROW = IBEG, IEND
               A(IROW) = A(IROW) + MULT1*A(K1) + MULT2*A(K2)
               K1 = K1 + 1_8
               K2 = K2 + 1_8
           ENDDO
           A( JJ       ) = -MULT1
           A( JJ + 1_8 ) = -MULT2
           IBEG = IBEG + NFRONT8
           IEND = IEND + NFRONT8
           JJ   = JJ   + NFRONT8
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_FAC_MQ_LDLT
      SUBROUTINE DMUMPS_FAC_SQ_LDLT(IBEG_BLOCK,IEND_BLOCK,NPIV,
     &    NFRONT,NASS,LAST_VAR,INODE,A,LA,
     &    LDA,
     &    POSELT,
     &    KEEP,KEEP8,
     &    PIVOT_OPTION, CALL_TRSM, UPDATE_SCHUR)
      IMPLICIT NONE
      INTEGER, intent(in) :: NPIV
      INTEGER, intent(in) :: NFRONT, NASS, IBEG_BLOCK, IEND_BLOCK
      INTEGER(8), intent(in) :: LA
      DOUBLE PRECISION, intent(inout) :: A(LA)
      INTEGER, intent(in) :: INODE
      INTEGER, intent(in) :: LAST_VAR 
      INTEGER :: KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8), intent(in) :: POSELT
      INTEGER, intent(in) :: LDA
      INTEGER, intent(in) :: PIVOT_OPTION
      LOGICAL, intent(in) :: CALL_TRSM, UPDATE_SCHUR
      INTEGER(8) :: LDA8
      INTEGER NPIV_BLOCK, NEL1, I, II
      INTEGER(8) :: LPOS,UPOS,APOS
      INTEGER IROW
      INTEGER Block
      INTEGER BLSIZE, ELSIZE
      DOUBLE PRECISION ONE, ALPHA, VALPIV
      INCLUDE 'mumps_headers.h'
      PARAMETER (ONE=1.0D0, ALPHA=-1.0D0)
      LDA8 = int(LDA,8)
      ELSIZE = IEND_BLOCK - IBEG_BLOCK +1
      NEL1 = LAST_VAR - IEND_BLOCK
      NPIV_BLOCK  = NPIV - IBEG_BLOCK + 1
      IF (NPIV_BLOCK.EQ.0) GO TO 500
      IF (NEL1.NE.0) THEN
        IF (PIVOT_OPTION.LE.1.AND.CALL_TRSM) THEN
          APOS = POSELT + LDA8*int(IBEG_BLOCK-1,8) + int(IBEG_BLOCK-1,8)
          LPOS = POSELT + LDA8*int(IEND_BLOCK,8) + int(IBEG_BLOCK-1,8)
          UPOS = POSELT + LDA8*int(IBEG_BLOCK-1,8) + int(IEND_BLOCK,8)
          CALL dtrsm('L', 'U', 'T', 'U', ELSIZE, NEL1, ONE, 
     &              A(APOS), LDA, A(LPOS), LDA)
!$OMP PARALLEL PRIVATE(VALPIV,I,II)          
         DO I = 1, ELSIZE
          VALPIV = ONE/A(POSELT+(LDA8+1_8)*int(IBEG_BLOCK+I-2,8))
!$OMP DO 
          DO II = 1,NEL1
            A(UPOS+int(I-1,8)*LDA8 + int(II-1,8)) = 
     &          A(LPOS+int(I-1,8) + int(II-1,8)*LDA8)
            A(LPOS+int(I-1,8) + int(II-1,8)*LDA8) = 
     &          A(LPOS+int(I-1,8) + int(II-1,8)*LDA8)*VALPIV
          ENDDO
!$OMP END DO NOWAIT
         ENDDO
!$OMP END PARALLEL
        ENDIF
        IF (UPDATE_SCHUR) THEN
#if defined(GEMMT_AVAILABLE)
         IF ( KEEP(421).EQ. -1) THEN
           LPOS = POSELT + LDA8*int(IEND_BLOCK,8) + int(IBEG_BLOCK-1,8)
           UPOS = POSELT + LDA8*int(IBEG_BLOCK-1,8) + int(IEND_BLOCK,8)
           APOS = POSELT + LDA8*int(IEND_BLOCK,8) + int(IEND_BLOCK,8)
           CALL dgemmt( 'U','N','N', NEL1,
     &                NPIV_BLOCK,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA, ONE, A( APOS ), LDA )
         ELSE
#endif
          IF ( LAST_VAR - IEND_BLOCK > KEEP(7) ) THEN
           BLSIZE = KEEP(8)
          ELSE
           BLSIZE = LAST_VAR - IEND_BLOCK
          END IF
          IF ( NASS - IEND_BLOCK .GT. 0 ) THEN
#if defined(SAK_BYROW)
           DO IROW = IEND_BLOCK+1, LAST_VAR, BLSIZE
            Block = min( BLSIZE, NASS - IROW + 1 )
            LPOS = POSELT + int(IROW  - 1,8) * LDA8 +
     &                      int(IBEG_BLOCK - 1,8)
            UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                      int(IROW - 1,8)
            APOS = POSELT + int(IROW  - 1,8) * LDA8 + 
     &             int(IEND_BLOCK,8)
            CALL dgemm( 'N','N', IROW + Block - IEND_BLOCK - 1, 
     &             Block, NPIV_BLOCK,
     &                 ALPHA, A( UPOS ), LDA,
     &                 A( LPOS ), LDA, ONE, A( APOS ), LDA )
           ENDDO
#else
           DO IROW = IEND_BLOCK+1, LAST_VAR, BLSIZE
            Block = min( BLSIZE, LAST_VAR - IROW + 1 )
            LPOS = POSELT + int( IROW - 1,8) * LDA8 +
     &                      int(IBEG_BLOCK - 1,8)
            UPOS = POSELT + int(IBEG_BLOCK - 1,8) * LDA8 +
     &                      int( IROW - 1,8)
            APOS = POSELT + int( IROW - 1,8) * LDA8 + int( IROW - 1,8)
            CALL dgemm( 'N','N', Block, LAST_VAR - IROW + 1, NPIV_BLOCK,
     &                 ALPHA, A( UPOS ), LDA,
     &                 A( LPOS ), LDA, ONE, A( APOS ), LDA )
           END DO
#endif
          END IF
#if defined(GEMMT_AVAILABLE)
         END IF 
#endif
         LPOS = POSELT + int(LAST_VAR,8)*LDA8 + int(IBEG_BLOCK - 1,8)
         UPOS = POSELT + int(IBEG_BLOCK-1,8) * LDA8 +
     &                 int(IEND_BLOCK,8)
         APOS = POSELT + int(LAST_VAR,8)*LDA8 + int(IEND_BLOCK,8)
         IF (PIVOT_OPTION.EQ.3) THEN
           CALL dgemm('N', 'N', NEL1, NFRONT-LAST_VAR, NPIV_BLOCK, 
     &              ALPHA, A(UPOS), LDA, A(LPOS), LDA, ONE, 
     &              A(APOS), LDA)
         ELSEIF (PIVOT_OPTION.EQ.2.AND.(NASS.GT. LAST_VAR)) THEN
           CALL dgemm('N', 'N', NEL1, NASS-LAST_VAR, NPIV_BLOCK, 
     &              ALPHA, A(UPOS), LDA, A(LPOS), LDA, ONE, 
     &              A(APOS), LDA)
         ENDIF
        ENDIF
      ENDIF
  500 CONTINUE
      RETURN
      END SUBROUTINE DMUMPS_FAC_SQ_LDLT
        SUBROUTINE DMUMPS_SWAP_LDLT( A, LA, IW, LIW,
     &                       IOLDPS, NPIVP1, IPIV, POSELT, NASS,
     &                       LDA, NFRONT, LEVEL, K219, K50, XSIZE,
     &                       IBEG_BLOCK_TO_SEND )
        IMPLICIT NONE
      INTEGER(8) :: POSELT, LA
      INTEGER LIW, IOLDPS, NPIVP1, IPIV
      INTEGER LDA, NFRONT, NASS, LEVEL, K219, K50, XSIZE
      DOUBLE PRECISION A( LA )
      INTEGER IW( LIW )
      INTEGER, INTENT(IN) :: IBEG_BLOCK_TO_SEND
      INCLUDE 'mumps_headers.h'
      INTEGER :: LASTROW2SWAP, IBEG
      INTEGER ISW, ISWPS1, ISWPS2, HF
      INTEGER(8) :: IDIAG, APOS
      INTEGER(8) :: LDA8
      DOUBLE PRECISION SWOP
            LDA8 = int(LDA,8)
            APOS = POSELT + LDA8*int(IPIV-1,8) + int(NPIVP1-1,8)
            IDIAG = APOS + int(IPIV - NPIVP1,8)
            HF = 6 + IW( IOLDPS + 5 + XSIZE) + XSIZE
            ISWPS1 = IOLDPS + HF + NPIVP1 - 1
            ISWPS2 = IOLDPS + HF + IPIV - 1
            ISW = IW(ISWPS1)
            IW(ISWPS1) = IW(ISWPS2)
            IW(ISWPS2) = ISW
            ISW = IW(ISWPS1+NFRONT)
            IW(ISWPS1+NFRONT) = IW(ISWPS2+NFRONT)
            IW(ISWPS2+NFRONT) = ISW
            IF ( LEVEL .eq. 2 ) THEN
              IBEG = IBEG_BLOCK_TO_SEND
              CALL dswap( NPIVP1 - 1 - IBEG + 1,
     &            A( POSELT + int(NPIVP1-1,8) +
     &                     int(IBEG-1,8) * LDA8), LDA,
     &            A( POSELT + int(IPIV-1,8)   +
     &                     int(IBEG-1,8) * LDA8), LDA )
            END IF
            CALL dswap( NPIVP1-1,
     &           A( POSELT+int(NPIVP1-1,8) * LDA8 ), 1,
     &           A( POSELT + int(IPIV-1,8) * LDA8 ), 1 )
             CALL dswap( IPIV - NPIVP1 - 1,
     &           A( POSELT+int(NPIVP1,8) * LDA8 + int(NPIVP1-1,8) ),
     &           LDA, A( APOS + 1_8 ), 1 )
            SWOP = A(IDIAG)
            A(IDIAG) = A( POSELT+int(NPIVP1-1,8)*LDA8+int(NPIVP1-1,8) )
            A( POSELT + int(NPIVP1-1,8)*LDA8 + int(NPIVP1-1,8) ) = SWOP
            IF (LEVEL .EQ. 1) THEN
              LASTROW2SWAP = NFRONT
            ELSE
              LASTROW2SWAP = NASS
            ENDIF
            CALL dswap( LASTROW2SWAP - IPIV,
     &      A( APOS  + LDA8 ), LDA,
     &      A( IDIAG + LDA8 ), LDA )
            IF (K219.NE.0 .AND.K50.EQ.2) THEN
             IF ( LEVEL .eq. 2) THEN
              APOS                 = POSELT+LDA8*LDA8-1_8
              SWOP                 = A(APOS+int(NPIVP1,8))
              A(APOS+int(NPIVP1,8))= A(APOS+int(IPIV,8))
              A(APOS+int(IPIV,8))  = SWOP
             ENDIF
            ENDIF
        RETURN
        END SUBROUTINE DMUMPS_SWAP_LDLT
      SUBROUTINE DMUMPS_FAC_T_LDLT_COPY2U_SCALEL( IROWMAX, IROWMIN,
     &                  SIZECOPY, LDA, NPIV, LIW, IW, OFFSET_IW,
     &                  LA, A, POSELT, A_LPOS, A_UPOS )
!$    USE OMP_LIB
      INTEGER, INTENT(IN) :: IROWMAX, IROWMIN
      INTEGER, INTENT(IN) :: SIZECOPY
      INTEGER, INTENT(IN) :: LDA, NPIV
      INTEGER, INTENT(IN) :: LIW
      INTEGER, INTENT(IN) :: IW(LIW)
      INTEGER, INTENT(IN) :: OFFSET_IW
      INTEGER(8), INTENT(IN) :: LA
      DOUBLE PRECISION, INTENT(INOUT) :: A(LA)
      INTEGER(8), INTENT(IN) :: POSELT, A_LPOS, A_UPOS
      INTEGER(8) :: LPOS, UPOS
      INTEGER(8) :: DPOS, POSPV1, POSPV2, OFFDAG
      INTEGER(8) :: LDA8
      INTEGER :: IROWEND, IROW, Block2
      INTEGER :: I, J
      DOUBLE PRECISION :: MULT1, MULT2, A11, DETPIV, A22, A12
      INTEGER :: BLSIZECOPY
      DOUBLE PRECISION :: ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER(8) :: LPOSI, UPOSI
      LOGICAL :: PIVOT_2X2
!$    LOGICAL :: OMP_FLAG
!$    INTEGER :: NOMP, CHUNK
      LDA8 = int(LDA,8)
      IF (SIZECOPY.NE.0) THEN
        BLSIZECOPY = SIZECOPY
      ELSE
        BLSIZECOPY = 250
      ENDIF
!$    NOMP = OMP_GET_MAX_THREADS()
!$    OMP_FLAG = .FALSE.
!$    CHUNK = (64/4)
!$    IF (NOMP .GT. 1 .AND. NPIV .GE. 4*CHUNK) THEN
!$      OMP_FLAG = .TRUE.
!$      CHUNK = max(2*CHUNK, NPIV/NOMP)
!$    ENDIF
      DO IROWEND = IROWMAX, IROWMIN, -BLSIZECOPY
        Block2 = min(BLSIZECOPY, IROWEND)
        IROW = IROWEND - Block2 + 1
        LPOS = A_LPOS + int(IROW-1,8)*LDA8
        UPOS = A_UPOS + int(IROW-1,8)
!$OMP  PARALLEL DO PRIVATE(PIVOT_2X2, A11, DPOS,
!$OMP&   POSPV1, POSPV2, OFFDAG, A22, A12, DETPIV, J, MULT1, MULT2
!$OMP&   , LPOSI, UPOSI
!$OMP&   ) FIRSTPRIVATE(Block2, LDA, LDA8, LPOS, UPOS, POSELT)
!$OMP& SCHEDULE(STATIC,CHUNK) IF(OMP_FLAG)
        DO I=1, NPIV
          PIVOT_2X2 = .FALSE.
          IF(IW(OFFSET_IW+I-1) .LE. 0) THEN
            PIVOT_2X2 = .TRUE.
          ELSE
            IF (I .GT. 1) THEN
              IF (IW(OFFSET_IW+I-2) .LE. 0) THEN
                cycle
              ENDIF
            ENDIF
          ENDIF
          DPOS = POSELT + LDA8*int(I-1,8) + int(I-1,8)
          IF(.not. PIVOT_2X2) THEN
            A11 = ONE/A(DPOS)
            LPOSI = LPOS+int(I-1,8)
            UPOSI = UPOS+int(I-1,8)*LDA8
            DO J = 1, Block2
              A(UPOSI+int(J-1,8)) = A(LPOSI+int(J-1,8)*LDA8)
            END DO
            DO J = 1, Block2
              A(LPOSI+int(J-1,8)*LDA8) = A(LPOSI+int(J-1,8)*LDA8)*A11
            END DO
          ELSE
            CALL dcopy(Block2, A(LPOS+int(I-1,8)),
     &                 LDA, A(UPOS+int(I-1,8)*LDA8), 1)
            CALL dcopy(Block2, A(LPOS+int(I,8)),
     &                 LDA, A(UPOS+int(I,8)*LDA8), 1)
            POSPV1 = DPOS
            POSPV2 = DPOS + int(LDA+1,8)
            OFFDAG = POSPV1+1_8
            A11 = A(POSPV1)
            A22 = A(POSPV2)
            A12 = A(OFFDAG)
            DETPIV = A11*A22 - A12**2
            A22 = A11/DETPIV
            A11 = A(POSPV2)/DETPIV
            A12 = -A12/DETPIV
            DO J = 1,Block2
              MULT1 = A11*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A12*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              MULT2 = A12*A(LPOS+int(J-1,8)*LDA8+int(I-1,8))
     &              + A22*A(LPOS+int(J-1,8)*LDA8+int(I,8))
              A(LPOS+int(J-1,8)*LDA8+int(I-1,8)) = MULT1
              A(LPOS+int(J-1,8)*LDA8+int(I,8))   = MULT2
            ENDDO
          ENDIF
        ENDDO
!$OMP END PARALLEL DO
      ENDDO
      END SUBROUTINE DMUMPS_FAC_T_LDLT_COPY2U_SCALEL
      SUBROUTINE DMUMPS_FAC_T_LDLT(NFRONT,NASS,
     &    IW,LIW,A,LA,
     &    LDA,
     &    IOLDPS,POSELT,KEEP,KEEP8,
     &    POSTPONE_COL_UPDATE, ETATASS,
     &    TYPEFile, LAFAC, MonBloc, NextPiv2beWritten,
     &    LIWFAC, MYID, IFLAG, OFFSET_IW
     &    )
      USE DMUMPS_OOC
      IMPLICIT NONE
      INTEGER NFRONT, NASS,LIW
      INTEGER(8) :: LA
      DOUBLE PRECISION    A(LA)
      INTEGER IW(LIW) 
      INTEGER KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER(8) :: POSELT
      INTEGER LDA
      INTEGER IOLDPS, ETATASS
      LOGICAL POSTPONE_COL_UPDATE
      INTEGER(8) :: LAFAC
      INTEGER TYPEFile, NextPiv2beWritten
      INTEGER LIWFAC, MYID, IFLAG
      TYPE(IO_BLOCK):: MonBloc
      INTEGER IDUMMY
      LOGICAL LAST_CALL
      INTEGER :: OFFSET_IW
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: UPOS, APOS, LPOS
      INTEGER(8) :: LDA8
      INTEGER BLSIZE, BLSIZE2, Block, IROW, NPIV, IROWEND
      INTEGER I2, I2END, Block2
      DOUBLE PRECISION  ONE, ALPHA, BETA, ZERO
      PARAMETER (ONE = 1.0D0, ALPHA=-1.0D0)
      PARAMETER (ZERO=0.0D0)
      LDA8 = int(LDA,8)
      IF (ETATASS.EQ.1) THEN
        BETA = ZERO
      ELSE
        BETA = ONE
      ENDIF
      IF ( NFRONT - NASS > KEEP(58) ) THEN
        IF ( NFRONT - NASS > KEEP(57) ) THEN
          BLSIZE = KEEP(58)
        ELSE
          BLSIZE = (NFRONT - NASS)/2
        END IF
      ELSE
        BLSIZE = NFRONT - NASS
      END IF
      BLSIZE2 = KEEP(218)
      NPIV = IW( IOLDPS + 1 + KEEP(IXSZ))
      IF ( NFRONT - NASS .GT. 0 ) THEN
       IF ( POSTPONE_COL_UPDATE ) THEN
         CALL dtrsm( 'L', 'U', 'T', 'U',
     &               NPIV, NFRONT-NPIV, ONE,
     &               A( POSELT ), LDA,
     &               A( POSELT + LDA8 * int(NPIV,8) ), LDA )
       ENDIF
#if defined(GEMMT_AVAILABLE)
       IF ( KEEP(421).EQ. -1) THEN
         LPOS = POSELT + int(NASS,8)*LDA8
         UPOS = POSELT + int(NASS,8)
         APOS = POSELT + int(NASS,8)*LDA8 + int(NASS,8)
         IF (POSTPONE_COL_UPDATE) THEN
           CALL DMUMPS_FAC_T_LDLT_COPY2U_SCALEL( NFRONT - NASS, 1,
     &                 KEEP(424), NFRONT, NPIV,
     &                 LIW, IW, OFFSET_IW, LA, A, POSELT, LPOS, UPOS )
         ENDIF
         CALL dgemmt('U', 'N', 'N', NFRONT-NASS, NPIV,
     &                ALPHA, A( UPOS ), LDA,
     &                A( LPOS ), LDA,
     &                BETA,
     &                A( APOS ), LDA )
       ELSE
#endif
       DO IROWEND = NFRONT - NASS, 1, -BLSIZE
        Block = min( BLSIZE, IROWEND )
        IROW  = IROWEND - Block + 1
        LPOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8
        APOS = POSELT + int(NASS,8)*LDA8 + int(IROW-1,8) * LDA8 +
     &                  int(NASS + IROW - 1,8)
        UPOS = POSELT + int(NASS,8)
        IF (.NOT. POSTPONE_COL_UPDATE) THEN
          UPOS = POSELT + int(NASS + IROW - 1,8)
        ENDIF
        IF (POSTPONE_COL_UPDATE) THEN
          CALL DMUMPS_FAC_T_LDLT_COPY2U_SCALEL( Block, 1,
     &                KEEP(424), NFRONT, NPIV,
     &                LIW, IW, OFFSET_IW, LA, A, POSELT, LPOS, UPOS )
        ENDIF
        DO I2END = Block, 1, -BLSIZE2
          Block2 = min(BLSIZE2, I2END)
          I2 = I2END - Block2+1
          CALL dgemm('N', 'N', Block2, Block-I2+1, NPIV, ALPHA,
     &               A(UPOS+int(I2-1,8)), LDA,
     &               A(LPOS+int(I2-1,8)*LDA8), LDA,
     &               BETA,
     &               A(APOS + int(I2-1,8) + int(I2-1,8)*LDA8), LDA)
          IF (KEEP(201).EQ.1) THEN
            IF (NextPiv2beWritten.LE.NPIV) THEN
              LAST_CALL=.FALSE.
              CALL DMUMPS_OOC_IO_LU_PANEL(
     &        STRAT_TRY_WRITE, TYPEFile,
     &        A(POSELT), LAFAC, MonBloc,
     &        NextPiv2beWritten, IDUMMY,
     &        IW(IOLDPS), LIWFAC, MYID,
     &        KEEP8(31),
     &        IFLAG,LAST_CALL )
              IF (IFLAG .LT. 0 ) RETURN
            ENDIF
          ENDIF
        ENDDO
        IF ( NFRONT - NASS - IROW + 1 - Block > 0 ) THEN
        CALL dgemm( 'N', 'N', Block, NFRONT-NASS-Block-IROW+1, NPIV,
     &              ALPHA,  A( UPOS ), LDA,
     &              A( LPOS + LDA8 * int(Block,8) ), LDA,
     &              BETA,
     &              A( APOS + LDA8 * int(Block,8) ), LDA )
        ENDIF
       END DO
#if defined(GEMMT_AVAILABLE)
      END IF
#endif
      END IF
      RETURN
      END SUBROUTINE DMUMPS_FAC_T_LDLT
      SUBROUTINE DMUMPS_STORE_PERMINFO( PIVRPTR, NBPANELS, PIVR, NASS,
     &                                  K, P, LastPanelonDisk,
     &                                  LastPIVRPTRIndexFilled )
      IMPLICIT NONE
      INTEGER, intent(in) :: NBPANELS, NASS, K, P
      INTEGER, intent(inout) :: PIVRPTR(NBPANELS), PIVR(NASS)
      INTEGER LastPanelonDisk, LastPIVRPTRIndexFilled
      INTEGER I
      IF ( LastPanelonDisk+1 > NBPANELS ) THEN
           WRITE(*,*) "INTERNAL ERROR IN DMUMPS_STORE_PERMINFO!"
           WRITE(*,*) "NASS=",NASS,"PIVRPTR=",PIVRPTR(1:NBPANELS)
           WRITE(*,*) "K=",K, "P=",P, "LastPanelonDisk=",LastPanelonDisk
           WRITE(*,*) "LastPIVRPTRIndexFilled=", LastPIVRPTRIndexFilled
           CALL MUMPS_ABORT()
      ENDIF
      PIVRPTR(LastPanelonDisk+1) = K + 1
      IF (LastPanelonDisk.NE.0) THEN
        PIVR(K - PIVRPTR(1) + 1) = P
        DO I = LastPIVRPTRIndexFilled + 1, LastPanelonDisk
          PIVRPTR(I)=PIVRPTR(LastPIVRPTRIndexFilled)
        ENDDO
      ENDIF
      LastPIVRPTRIndexFilled = LastPanelonDisk + 1
      RETURN
      END SUBROUTINE DMUMPS_STORE_PERMINFO
      END MODULE DMUMPS_FAC_FRONT_AUX_M
