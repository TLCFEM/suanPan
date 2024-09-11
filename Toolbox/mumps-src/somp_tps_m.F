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
      MODULE SMUMPS_TPS_M
      TYPE SMUMPS_TPS_T
          REAL, DIMENSION(:), POINTER :: A
      END TYPE SMUMPS_TPS_T
      END MODULE SMUMPS_TPS_M
      SUBROUTINE SMUMPS_TPS_M_RETURN()
      RETURN
      END SUBROUTINE SMUMPS_TPS_M_RETURN
