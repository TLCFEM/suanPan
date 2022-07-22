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
C -----------------------------------------
C  This file contains the definition
C  of all tags.
C -----------------------------------------
C
C -----------------------------------------
C  Tag for arrowheads distribution
C -----------------------------------------
      INTEGER ARROWHEAD, ARR_INT, ARR_REAL, ELT_INT, ELT_REAL
      PARAMETER ( ARROWHEAD = 20,
     &            ARR_INT = 29,
     &            ARR_REAL = 30,
     &            ELT_INT = 31,
     &            ELT_REAL = 32 )
C ----------------------------------------------------
C   Tags for collecting distributed integer info
C   for analysis in case of initial distributed matrix
C ----------------------------------------------------
      INTEGER COLLECT_NZ, COLLECT_IRN, COLLECT_JCN
      PARAMETER( COLLECT_NZ  = 35,
     &           COLLECT_IRN = 36,
     &           COLLECT_JCN = 37 )
C -----------------------------------------
C   Tags for factorization
C -----------------------------------------
      INTEGER RACINE,
     &        NOEUD,
     &        TERREUR,
     &        MAITRE_DESC_BANDE,
     &        MAITRE2,
     &        BLOC_FACTO_RELAY,
     &        CONTRIB_TYPE2,
     &        MAPLIG,
     &        FACTOR,
     &        BLOC_FACTO
      PARAMETER ( RACINE            = 2,
     &            NOEUD             = 3,
     &            MAITRE_DESC_BANDE = 4,
     &            MAITRE2           = 5,
     &            BLOC_FACTO_RELAY  = 6,
     &            CONTRIB_TYPE2     = 7,
     &            MAPLIG            = 8,
     &            FACTOR            = 9,
     &            BLOC_FACTO        = 10,
     &            TERREUR           = 99 )
C -----------------------------------------
C   Tags for assembly of root (in facto)
C -----------------------------------------
      INTEGER ROOT_NELIM_INDICES,
     &        ROOT_CONT_STATIC,
     &        ROOT_NON_ELIM_CB,
     &        ROOT_2SLAVE,
     &        ROOT_2SON
       PARAMETER( ROOT_NELIM_INDICES = 15,
     &        ROOT_CONT_STATIC       = 16,
     &        ROOT_NON_ELIM_CB       = 17,
     &        ROOT_2SLAVE            = 18,
     &        ROOT_2SON              = 19 )
C -----------------------------------------
C   Tags for solve
C -----------------------------------------
      INTEGER RACINE_SOLVE,
     &        ContVec,
     &        Master2Slave,
     &        GatherSol,
     &        ScatterRhsI,
     &        ScatterRhsR,
     &        DistRhsI,
     &        DistRhsR
      PARAMETER( RACINE_SOLVE = 14,
     &           ContVec      = 11,
     &           Master2Slave = 12,
     &           GatherSol    = 13,
     &           ScatterRhsI  = 54,
     &           ScatterRhsR  = 55,
     &           DistRhsI     = 51,
     &           DistRhsR     = 52)
      INTEGER, PARAMETER :: DIST_RHS_INT    = 56
      INTEGER, PARAMETER :: DIST_RHS_SCALAR = 57
C -----------------------------------------
C   Tags for backsolve
C -----------------------------------------
      INTEGER FEUILLE,
     &        BACKSLV_UPDATERHS,
     &        BACKSLV_MASTER2SLAVE
      PARAMETER( FEUILLE = 21,
     &           BACKSLV_UPDATERHS = 22,
     &           BACKSLV_MASTER2SLAVE = 23 )
C ------------------------
C   Tag for symmetrization
C ------------------------
      INTEGER SYMMETRIZE
      PARAMETER ( SYMMETRIZE = 24 )
C ----------------------------
C   Tags specific to symmetric
C ----------------------------
      INTEGER BLOC_FACTO_SYM,
     &        BLOC_FACTO_SYM_SLAVE, END_NIV2_LDLT,
     &        END_NIV2
      PARAMETER ( BLOC_FACTO_SYM = 25,
     &            BLOC_FACTO_SYM_SLAVE = 26, 
     &            END_NIV2_LDLT = 33,
     &            END_NIV2 = 34 )
C -------------------------------------
C   Tags specific to dynamic scheduling
C -------------------------------------
      INTEGER UPDATE_LOAD
      PARAMETER ( UPDATE_LOAD = 27 )
C   To send deficientcy
      INTEGER DEFIC_TAG
      PARAMETER(  DEFIC_TAG = 28 )
C   To send Schur
      INTEGER TAG_SCHUR
      PARAMETER( TAG_SCHUR = 38 )
C   To clean up IRECV
      INTEGER TAG_DUMMY
      PARAMETER( TAG_DUMMY = 39 )
C   To send zero pivot indices
      INTEGER ZERO_PIV
      PARAMETER( ZERO_PIV = 40 )
C   To send Singular values (if defined(try_null_space))
      INTEGER TAG_ROOT1, TAG_ROOT2
      PARAMETER( TAG_ROOT1 = 41 )
      PARAMETER( TAG_ROOT2 = 42 )
C
C   Note: tags 100-160 are reserved for
C         the parallel scaling routine
C
