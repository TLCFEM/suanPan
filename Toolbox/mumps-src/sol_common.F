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
      SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS ( ISTEP, KEEP,
     &           NPIV, LIELL, IPOS, IW, LIW, PTRIST, STEP, N )
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: ISTEP, LIW, KEEP(500), N
      INTEGER, INTENT(IN)  :: IW( LIW )
      INTEGER, INTENT(IN)  :: STEP( N ), PTRIST( KEEP(28) )
      INTEGER, INTENT(OUT) :: NPIV, LIELL, IPOS
      INCLUDE 'mumps_headers.h'
      INTEGER :: SROOT
      IF (KEEP(38) .NE. 0) THEN
        SROOT = STEP(KEEP(38))
      ELSE IF (KEEP(20) .NE. 0) THEN
        SROOT = STEP(KEEP(20))
      ELSE
        SROOT = 0
      ENDIF
      IPOS = PTRIST(ISTEP)
      IF (IPOS .LE. 0) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_SOL_GET_NPIV_LIELL_IPOS",
     &  ISTEP 
        CALL MUMPS_ABORT()
      ENDIF
      NPIV = IW(IPOS+3+KEEP(IXSZ))
      IF ( ISTEP.EQ.SROOT ) THEN
             IPOS = PTRIST(ISTEP)
             LIELL = IW(IPOS+3+KEEP(IXSZ))
             NPIV = LIELL
             IPOS= PTRIST(ISTEP)+5+KEEP(IXSZ)
      ELSE
             IPOS = PTRIST(ISTEP) + 2+ KEEP(IXSZ)
             LIELL = IW(IPOS-2)+IW(IPOS+1)
             IPOS= IPOS+1
             NPIV = IW(IPOS)
             IPOS= IPOS+1
             IPOS= IPOS+1+IW( PTRIST(ISTEP) + 5 +KEEP(IXSZ))
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SOL_GET_NPIV_LIELL_IPOS
