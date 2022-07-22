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
      SUBROUTINE DMUMPS_FREE_BLOCK_CB_STATIC(
     &       SSARBR, MYID, N, IPOSBLOCK,
     &       IW, LIW,
     &       LRLU, LRLUS, IPTRLU,
     &       IWPOSCB, LA, KEEP, KEEP8, IN_PLACE_STATS
     &     )
      USE DMUMPS_LOAD
      IMPLICIT NONE
      INTEGER IPOSBLOCK,
     &         LIW, IWPOSCB, N
      INTEGER(8) :: LA, LRLU, LRLUS, IPTRLU
      LOGICAL IN_PLACE_STATS
      INTEGER IW( LIW ), KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER MYID
      LOGICAL SSARBR
      INTEGER SIZFI_BLOCK, SIZFI
      INTEGER IPOSSHIFT
      INTEGER(8) :: SIZFR, SIZFR_BLOCK, SIZFR_BLOCK_EFF,
     &              SIZEHOLE, MEM_INC, DYNSIZE_BLOCK
      INCLUDE 'mumps_headers.h'
      IPOSSHIFT = IPOSBLOCK + KEEP(IXSZ)
      SIZFI_BLOCK=IW(IPOSBLOCK+XXI)
      CALL MUMPS_GETI8( SIZFR_BLOCK, IW(IPOSBLOCK+XXR) )
      CALL MUMPS_GETI8( DYNSIZE_BLOCK,IW(IPOSBLOCK+XXD) )
      IF (DYNSIZE_BLOCK .GT. 0_8) THEN
        SIZFR_BLOCK_EFF = 0_8
      ELSE IF (KEEP(216).eq.3 
     &    ) THEN
        SIZFR_BLOCK_EFF = SIZFR_BLOCK
      ELSE
        CALL DMUMPS_SIZEFREEINREC( IW(IPOSBLOCK),
     &                     LIW-IPOSBLOCK+1,
     &                     SIZEHOLE, KEEP(IXSZ))
        SIZFR_BLOCK_EFF = SIZFR_BLOCK - SIZEHOLE
      ENDIF
      IF (.NOT. IN_PLACE_STATS) THEN
           LRLUS   = LRLUS + SIZFR_BLOCK_EFF
           KEEP8(69) = KEEP8(69) - SIZFR_BLOCK_EFF
      ENDIF
      IF ( IPOSBLOCK .eq. IWPOSCB + 1 ) THEN
         IPTRLU  = IPTRLU  + SIZFR_BLOCK
         IWPOSCB = IWPOSCB + SIZFI_BLOCK
         LRLU    = LRLU  + SIZFR_BLOCK
         MEM_INC = -SIZFR_BLOCK_EFF
         IF (IN_PLACE_STATS) THEN
           MEM_INC= 0_8
         ENDIF
         CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &         LA-LRLUS,0_8,MEM_INC,KEEP,KEEP8,LRLUS)
 90      IF ( IWPOSCB .eq. LIW ) GO TO 100
         IPOSSHIFT = IWPOSCB + KEEP(IXSZ)
         SIZFI = IW( IWPOSCB+1+XXI )
         CALL MUMPS_GETI8( SIZFR, IW(IWPOSCB+1+XXR) )
         IF ( IW( IWPOSCB+1+XXS ) .EQ. S_FREE ) THEN
              IPTRLU  = IPTRLU + SIZFR
              LRLU    = LRLU + SIZFR
              IWPOSCB = IWPOSCB + SIZFI
              GO TO 90
         ENDIF
 100     CONTINUE
         IW( IWPOSCB+1+XXP)=TOP_OF_STACK
      ELSE
         IW( IPOSBLOCK +XXS)=S_FREE
         CALL DMUMPS_LOAD_MEM_UPDATE(SSARBR,.FALSE.,
     &                LA-LRLUS,0_8,-SIZFR_BLOCK_EFF,KEEP,KEEP8,LRLUS)
      END IF
      RETURN
      END SUBROUTINE DMUMPS_FREE_BLOCK_CB_STATIC
