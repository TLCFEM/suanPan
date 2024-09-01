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
C     Define constants for possible modes:
C
C       memory_save  = compute the size of the save
C                      file and of the structure
C       save         = save the instance
C       restore      = restore the instance
C       restore_ooc  = restore the ooc part of the
C                      instance
C       fake_restore = extract from the saved file
C                      the size of the save file and of
C                      the structure
C
      INTEGER, PARAMETER :: memory_save_mode=1
      INTEGER, PARAMETER :: save_mode=2
      INTEGER, PARAMETER :: restore_mode=3
      INTEGER, PARAMETER :: restore_ooc_mode=4
      INTEGER, PARAMETER :: fake_restore_mode=5
