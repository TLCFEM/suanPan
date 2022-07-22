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
      INTEGER FUNCTION DMUMPS_IXAMAX(N,X,INCX,GRAIN)
      IMPLICIT NONE
      DOUBLE PRECISION, intent(in) :: X(*)
      INTEGER, intent(in) :: INCX,N
      INTEGER, intent(in) :: GRAIN
      INTEGER idamax
      DMUMPS_IXAMAX = idamax(N,X,INCX)
      RETURN
      END FUNCTION DMUMPS_IXAMAX
