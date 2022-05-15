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
      SUBROUTINE MUMPS_SET_ORDERING(N, SYM, NPROCS, IORD,
     &                       symmetry, NBQD, AvgDens,
     &                       PROK, MP)
      IMPLICIT NONE
      INTEGER, intent(in)    :: N, NPROCS, SYM
      INTEGER, intent(in)    :: symmetry, NBQD, AvgDens, MP
      LOGICAL, intent(in)    :: PROK
      INTEGER, intent(inout) :: IORD
      INTEGER MAXQD
      PARAMETER (MAXQD=2)
      INTEGER SMALLSYM, SMALLUNS
      PARAMETER (SMALLUNS=5000, SMALLSYM=10000)
#if ! defined(metis) && ! defined(parmetis) && ! defined(metis4) && ! defined(parmetis3)
      IF ( IORD .EQ. 5 ) THEN
        IF (PROK) WRITE(MP,*)
     &  'WARNING: METIS not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
#if ! defined(pord)
      IF ( IORD .EQ. 4 ) THEN
        IF (PROK) WRITE(MP,*)
     &  'WARNING: PORD not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
#if ! defined(scotch) && !  defined(ptscotch)
      IF ( IORD .EQ. 3 ) THEN
        IF (PROK) WRITE(MP,*)
     &  'WARNING: SCOTCH not available. Ordering set to default.'
        IORD = 7
      END IF
#endif
      IF (IORD.EQ.7) THEN
        IF (SYM.NE.0) THEN
          IF ( N.LE.SMALLSYM ) THEN 
             IF (NBQD.GE.MAXQD) THEN
               IORD = 6         
             ELSE
               IORD = 2         
             ENDIF
          ELSE
#if  defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
             IORD = 5
#else
#  if defined(scotch) || defined(ptscotch)
             IORD = 3
#  else
#    if defined(pord)
               IORD = 4
#    else
               IORD = 6
#    endif
#  endif
#endif
          ENDIF
        ELSE
          IF ( N.LE.SMALLUNS ) THEN
            IF (NBQD.GE.MAXQD) THEN
              IORD = 6  
            ELSE
              IORD = 2   
            ENDIF
          ELSE
#if  defined(metis) || defined(parmetis) || defined(metis4) || defined(parmetis3)
            IORD = 5
#else
#  if defined(scotch) || defined(ptscotch)
            IORD = 3
#  else
#    if defined(pord)
              IORD = 4
#    else
              IORD = 6
#    endif
#  endif
#endif
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_SET_ORDERING
