C
C  This file is part of MUMPS 5.7.1, released
C  on Thu May  2 10:15:09 UTC 2024
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
      INTEGER FUNCTION SMUMPS_IXAMAX(N,X,INCX,GRAIN)
      IMPLICIT NONE
      REAL, intent(in) :: X(*)
      INTEGER, intent(in) :: INCX,N
      INTEGER, intent(in) :: GRAIN
      INTEGER isamax
      SMUMPS_IXAMAX = isamax(N,X,INCX)
      RETURN
      END FUNCTION SMUMPS_IXAMAX
