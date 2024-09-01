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
      SUBROUTINE SMUMPS_FAC_DRIVER( id)
      USE SMUMPS_BUF
      USE SMUMPS_LOAD
      USE SMUMPS_OOC
      USE SMUMPS_STRUC_DEF
      USE SMUMPS_FACSOL_L0OMP_M, ONLY: SMUMPS_FREE_L0_OMP_FACTORS,
     &  SMUMPS_INIT_L0_OMP_FACTORS
      USE SMUMPS_LR_STATS
      USE SMUMPS_LR_DATA_M, only: SMUMPS_BLR_INIT_MODULE, 
     &                            SMUMPS_BLR_END_MODULE
     &                          , SMUMPS_BLR_MOD_TO_STRUC
      USE MUMPS_PIVNUL_MOD, only: PIVNUL_LIST_STRUCT_T
      USE MUMPS_FRONT_DATA_MGT_M
#if ! defined(NO_FDM_DESCBAND)
      USE MUMPS_FAC_DESCBAND_DATA_M
#endif
#if ! defined(NO_FDM_MAPROW)
      USE MUMPS_FAC_MAPROW_DATA_M
#endif
!$    USE OMP_LIB
C     Derived datatype to pass pointers with implicit interfaces
      USE SMUMPS_FAC_S_IS_POINTERS_M, ONLY : S_IS_POINTERS_T
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C  Performs scaling, sorting in arrowhead, then
C  distributes the matrix, and perform
C  factorization.
C
C
      INTERFACE
      SUBROUTINE SMUMPS_ANORMINF(id, ANORMINF, LSCAL, EFF_SIZE_SCHUR)
      USE SMUMPS_STRUC_DEF
      TYPE (SMUMPS_STRUC), TARGET :: id
      REAL, INTENT(OUT) :: ANORMINF
      LOGICAL, INTENT(IN) :: LSCAL
      INTEGER, INTENT(IN) :: EFF_SIZE_SCHUR
      END SUBROUTINE SMUMPS_ANORMINF
      END INTERFACE
C
C  Parameters
C  ==========
C
      TYPE(SMUMPS_STRUC), TARGET :: id
C
C  MPI
C  ===
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
C
C  Local variables
C  ===============
C
      INCLUDE 'mumps_headers.h'
      INTEGER(8) :: NSEND8, NSEND_TOT8
      INTEGER(8) :: NLOCAL8, NLOCAL_TOT8
      INTEGER(4) :: I4
      INTEGER :: LDPTRAR, NELT_arg, NBRECORDS
      INTEGER :: ITMP, JTMP
      INTEGER :: KEEP464COPY, KEEP465COPY
      REAL    :: RATIOK465
      INTEGER(8) :: KEEP826_SAVE
      INTEGER(8) :: K67, K68, K70, K74, K75
      INTEGER(8) ITMP8
      INTEGER  MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
      INTEGER MP, LP, MPG, allocok
      LOGICAL PROK, PROKG, LSCAL, LPOK, COMPUTE_ANORMINF
C     Reception buffer
      INTEGER    :: SMUMPS_LBUFR, SMUMPS_LBUFR_BYTES
      INTEGER(8) :: SMUMPS_LBUFR_BYTES8 ! for intermediate computation
      INTEGER, ALLOCATABLE, DIMENSION(:) :: BUFR
C     Size of send buffers (in bytes)
      INTEGER    :: SMUMPS_LBUF, SMUMPS_LBUF_INT
      INTEGER(8) :: SMUMPS_LBUF8 ! for intermediate computation
C
      INTEGER PTRIST, PTRWB, MAXELT_SIZE,
     &     ITLOC, IPOOL, LPOOL
      INTEGER IRANK, ID_ROOT
      INTEGER KKKK
      INTEGER(8) :: NZ_locMAX8
      INTEGER(8) MEMORY_MD_ARG
      INTEGER(8) MAXS_BASE8, MAXS_BASE_RELAXED8
      REAL CNTL4, AVG_FLOPS
      INTEGER MIN_PERLU, MAXIS_ESTIM
      INTEGER SUM_INFO22_THIS_NODE, MAX_SUM_INFO22_THIS_NODE
C
      TYPE (S_IS_POINTERS_T) :: S_IS_POINTERS
      INTEGER   MAXIS
      INTEGER(8) :: MAXS
      INTEGER(8) :: I8, TMPpeak 
      INTEGER    :: ICNTL49_LOC, TMP_INFOG_4
      LOGICAL    :: Compact_S_Authorized
C     For S argument to arrowhead routines:
      INTEGER(8) :: MAXS_ARG
      REAL, TARGET :: S_DUMMY_ARG(1)
      REAL, POINTER, DIMENSION(:) :: S_PTR_ARG
      TYPE(PIVNUL_LIST_STRUCT_T)     :: PIVNUL_LIST_STRUCT
C     Arrowheads
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INTARR
      REAL, POINTER, DIMENSION(:) :: DBLARR
C     (pointer to point on used-data in some cases--elt-entry)
      REAL TMPTIME
      INTEGER  NOMP
!$    INTEGER(8) :: CHUNK8
      INTEGER NB_THREADS
      DOUBLE PRECISION TIMEAVG, TIMEMAX,
     &                 FLOPAVG, FLOPMAX
      REAL TMPFLOP
      INTEGER NPIV_CRITICAL_PATH, EFF_SIZE_SCHUR
      DOUBLE PRECISION TIME, TIMEET
      REAL ZERO, ONE, MONE
      PARAMETER( ZERO = 0.0E0, ONE = 1.0E0, MONE = -1.0E0)
      REAL CZERO
      PARAMETER( CZERO = 0.0E0 )
      INTEGER PERLU, TOTAL_MBYTES, K231, K232, K233, BLR_STRAT
      INTEGER, PARAMETER :: IDUMMY = -9999
      LOGICAL, PARAMETER :: BDUMMY =.FALSE.
      INTEGER, PARAMETER :: PANEL_TABSIZE = 20
      INTEGER COLOUR, COMM_FOR_SCALING ! For Simultaneous scaling
      INTEGER(8) :: LIWK
      INTEGER(8) :: LWK, LWK_REAL, LWRKR_TH, LWRKC_TH
      INTEGER    :: NOMP_MAX
C     I_AM_SLAVE: used to determine if proc has the role of a slave
C     WK_USER_PROVIDED is set to true when WK_USER is provided by user
      LOGICAL I_AM_SLAVE, PERLU_ON, WK_USER_PROVIDED, EARLYT3ROOTINS
      LOGICAL PRINT_MAXAVG, PRINT_NODEINFO
      REAL :: ANORMINF, SEUIL, SEUIL_LDLT_NIV2, Thresh_Seuil
      REAL :: CNTL1, CNTL3, CNTL5, CNTL6, EPS
      INTEGER N, LPN_LIST,POSBUF
      INTEGER, DIMENSION (:), ALLOCATABLE :: ITMP2
      INTEGER I,K
      INTEGER(8) :: ITEMP8
      INTEGER    :: PARPIV_T1
      INTEGER FRONTWISE
C temporary variables for collecting stats from all processors
      DOUBLE PRECISION :: TMP_MRY_LU_FR
      DOUBLE PRECISION :: TMP_MRY_LU_LRGAIN
      INTEGER :: KEEP399_SAVE
      DOUBLE PRECISION :: TMP_MRY_CB_FR
      DOUBLE PRECISION :: TMP_MRY_CB_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_LRGAIN
      DOUBLE PRECISION :: TMP_FLOP_TRSM
      DOUBLE PRECISION :: TMP_FLOP_PANEL
      DOUBLE PRECISION :: TMP_FLOP_FRFRONTS
      DOUBLE PRECISION :: TMP_FLOP_TRSM_FR
      DOUBLE PRECISION :: TMP_FLOP_TRSM_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_FR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LR
      DOUBLE PRECISION :: TMP_FLOP_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_FLOP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_ACCUM_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_COMPRESS
      DOUBLE PRECISION :: TMP_FLOP_CB_DECOMPRESS
      DOUBLE PRECISION :: TMP_FLOP_FACTO_FR
      INTEGER :: TMP_CNT_NODES
      DOUBLE PRECISION :: TMP_TIME_UPDATE
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR1
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR2
      DOUBLE PRECISION :: TMP_TIME_UPDATE_LRLR3
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRLR
      DOUBLE PRECISION :: TMP_TIME_UPDATE_FRFR
      DOUBLE PRECISION :: TMP_TIME_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_MIDBLK_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_FRSWAP_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_CB_COMPRESS
      DOUBLE PRECISION :: TMP_TIME_PANEL
      DOUBLE PRECISION :: TMP_TIME_FAC_I
      DOUBLE PRECISION :: TMP_TIME_FAC_MQ
      DOUBLE PRECISION :: TMP_TIME_FAC_SQ
      DOUBLE PRECISION :: TMP_TIME_LRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRTRSM
      DOUBLE PRECISION :: TMP_TIME_FRFRONTS
      DOUBLE PRECISION :: TMP_TIME_LR_MODULE
      DOUBLE PRECISION :: TMP_TIME_DIAGCOPY
      DOUBLE PRECISION :: TMP_TIME_DECOMP
      DOUBLE PRECISION :: TMP_TIME_DECOMP_UCFS
      DOUBLE PRECISION :: TMP_TIME_LRASM_NIV1
      DOUBLE PRECISION :: TMP_TIME_LRASM_LOCASM2
      DOUBLE PRECISION :: TMP_TIME_LRASM_MAPLIG1
      DOUBLE PRECISION :: TMP_TIME_LRASM_CONTRIB2
      DOUBLE PRECISION :: TMP_TIME_FRASM_LOCASM2
      DOUBLE PRECISION :: TMP_TIME_FRASM_MAPLIG1
      DOUBLE PRECISION :: TMP_TIME_FRASM_CONTRIB2
C
C  Workspace
C
      INTEGER, DIMENSION(:), ALLOCATABLE :: IWK
      REAL, DIMENSION(:), ALLOCATABLE :: WK
      REAL, DIMENSION(:), ALLOCATABLE    :: WK_REAL
      REAL, DIMENSION(:,:), ALLOCATABLE:: WRKR_TH, 
     &                                    WRKC_TH
      INTEGER(8), DIMENSION(:), ALLOCATABLE :: IWK8
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCP
      INTEGER, DIMENSION(:), ALLOCATABLE :: BURS
      INTEGER, DIMENSION(:), ALLOCATABLE :: BUCS
      INTEGER(8) :: BUREGISTRE(12)
      INTEGER(8) :: BUINTSZ, BURESZ, NZ_loc8
      INTEGER    :: BUJOB
      INTEGER BUMAXMN, M, SCMYID, SCNPROCS
      REAL    SCONEERR, SCINFERR
C
C  Parameters arising from the structure
C  =====================================
C
*     Control parameters: see description in SMUMPSID
      REAL,DIMENSION(:),POINTER::RINFO, RINFOG
      REAL,DIMENSION(:),POINTER::    CNTL
      INTEGER,DIMENSION(:),POINTER   :: INFOG, KEEP
      INTEGER, DIMENSION(:), POINTER :: IRN_loc_PTR, JCN_loc_PTR
      REAL, DIMENSION(:), POINTER   :: COLSCA_PTR, 
     &                                 ROWSCA_PTR
      REAL, DIMENSION(:), POINTER:: A_loc_PTR
      INTEGER, TARGET :: DUMMYIRN_loc(1), DUMMYJCN_loc(1)
      REAL, TARGET    :: DUMMYSCA(1)
      REAL, TARGET :: DUMMYA_loc(1)
      INTEGER,DIMENSION(:),POINTER::ICNTL
      EXTERNAL MUMPS_GET_POOL_LENGTH
      INTEGER MUMPS_GET_POOL_LENGTH, SIZESCAL
      INTEGER(8) :: TOTAL_BYTES
      INTEGER(8) :: I8TMP, LWK_USER_SUM8
C
C  External references
C  ===================
      INTEGER, EXTERNAL :: MUMPS_NUMROC
      INTEGER:: NWORKING
      LOGICAL:: MEM_EFF_ALLOCATED
      INTEGER   :: TOTAL_MBYTES_UNDER_L0
      INTEGER(8):: TOTAL_BYTES_UNDER_L0
C  Fwd in facto:
      REAL, DIMENSION(:), POINTER :: RHS_MUMPS
      LOGICAL :: RHS_MUMPS_ALLOCATED
      INTEGER :: NB_ACTIVE_FRONTS_ESTIM
      INTEGER :: NB_FRONTS_F_ESTIM
      INTEGER :: KEEP_486_FOR_PRINT
C
C     --------------------------
C     Pointers used as shortcuts
C     --------------------------
      RINFO=>id%RINFO
      RINFOG=>id%RINFOG
      CNTL=>id%CNTL
      INFOG=>id%INFOG
      KEEP=>id%KEEP
      ICNTL=>id%ICNTL
      IF (id%KEEP8(29) .NE. 0) THEN
        IRN_loc_PTR=>id%IRN_loc
        JCN_loc_PTR=>id%JCN_loc
        A_loc_PTR=>id%A_loc
      ELSE
        IRN_loc_PTR=>DUMMYIRN_loc
        JCN_loc_PTR=>DUMMYJCN_loc
        A_loc_PTR=>DUMMYA_loc
      ENDIF
      NOMP = 1 
      N = id%N
      EPS = epsilon ( ZERO )
C     TIMINGS: reset to 0
      id%DKEEP(92)=0.0E0
      id%DKEEP(93)=0.0E0
      id%DKEEP(94)=0.0E0
      id%DKEEP(95)=0.0E0
      id%DKEEP(96)=0.0E0
      id%DKEEP(97)=0.0E0
      id%DKEEP(98)=0.0E0
      id%DKEEP(99)=0.0E0
      id%DKEEP(56)=0.0E0
C     Count of MPI messages: reset to 0
      id%KEEP(266)=0
      id%KEEP(267)=0
C     MIN/MAX pivots reset to 0
      id%DKEEP(19)=huge(0.0E0)
      id%DKEEP(20)=huge(0.0E0)
      id%DKEEP(21)=0.0E0
C     RR related
      id%KEEP(17)  = 0
      id%INFOG(28) = 0
C     Number of symmetric swaps
      id%KEEP8(80)=0_8
C     Largest increase of internal panel size
      id%KEEP(425) =0
      PRINT_MAXAVG = .NOT.(id%NSLAVES.EQ.1 .AND. KEEP(46).EQ.1)
C     Print per node informtation only in case ther are several
C     compute nodes (id%KEEP(412): #MPI procs on comupte node)
      PRINT_NODEINFO = PRINT_MAXAVG .AND. id%NPROCS .NE. id%KEEP(412)
C
C     Related to forward in facto functionality (referred to as "Fwd in facto")
      NULLIFY(RHS_MUMPS)
      NULLIFY(DBLARR)
      RHS_MUMPS_ALLOCATED = .FALSE.
C     -----------------------------------------------------------------------
C     Set WK_USER_PROVIDED to true when workspace WK_USER is provided by user
C     We can accept WK_USER to be provided on only one proc and 
C     different values of WK_USER per processor
C     
      IF (id%KEEP8(24).GT.0_8) THEN 
C                We nullify S so that later when we test
C                if (associated(S) we can free space and reallocate it).
           NULLIFY(id%S)
      ENDIF
C
C     --  KEEP8(24) can now then be reset safely
      WK_USER_PROVIDED = (id%LWK_USER.NE.0)
      IF (WK_USER_PROVIDED) THEN
          IF (id%LWK_USER.GT.0) THEN
            id%KEEP8(24) = int(id%LWK_USER,8)
          ELSE
            id%KEEP8(24) = -int(id%LWK_USER,8)* 1000000_8 
          ENDIF
      ELSE
          id%KEEP8(24) = 0_8
      ENDIF
C     Compute sum of LWK_USER provided by user 
      LWK_USER_SUM8 = 0_8
      CALL MPI_REDUCE ( id%KEEP8(24), LWK_USER_SUM8, 1, MPI_INTEGER8,
     &                  MPI_SUM, MASTER, id%COMM, IERR )
C
C     KEEP8(26) might be modified
C       (element entry format) 
C       but need be restore for 
C       future factorisation
C       with different scaling option
C 
      KEEP826_SAVE = id%KEEP8(26)
C     In case of loop on factorization with
C     different scaling options, initialize
C     DKEEP(4:5) to 0.
      id%DKEEP(4)=-1.0E0
      id%DKEEP(5)=-1.0E0
C  Mapping information used during solve. In case of several facto+solve
C  it has to be recomputed. In case of several solves with the same
C  facto, it is not recomputed.
      IF (associated(id%IPTR_WORKING)) THEN
        DEALLOCATE(id%IPTR_WORKING)
        NULLIFY(id%IPTR_WORKING)
      END IF
      IF (associated(id%WORKING)) THEN 
        DEALLOCATE(id%WORKING)
        NULLIFY(id%WORKING)
      END IF
C
C  Units for printing
C  MP: diagnostics
C  LP: errors
C
      LP  = ICNTL( 1 )
      MP  = ICNTL( 2 )
      MPG = ICNTL( 3 )
      LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
      PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
      PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
      PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
      IF ( PROK ) WRITE( MP, 130 )
      IF ( PROKG ) WRITE( MPG, 130 )
C     -------------------------------------
C     Depending on the type of parallelism,
C     the master can now (soon) potentially
C     have the role of a slave
C     -------------------------------------
      I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &             ( id%MYID .eq. MASTER .AND.
     &               KEEP(46) .eq. 1 ) )
C
C  Prepare work for out-of-core
C
      IF (id%MYID .EQ. MASTER .AND. KEEP(201) .NE. -1) THEN
C       Note that if KEEP(201)=-1, then we have decided
C       at analysis phase that factors will not be stored
C       (neither in memory nor on disk). In that case,
C       ICNTL(22) is ignored.
C       -- ICNTL(22) must be set before facto phase 
C          (=1 OOC on; =0 OOC off)
C          and cannot be changed for subsequent solve phases.
        KEEP(201)=id%ICNTL(22)
        IF (KEEP(201) .EQ. 1) THEN
#         if defined(OLD_OOC_NOPANEL)
            KEEP(201)=2
#         else
            KEEP(201)=1
#         endif
        ELSE
          id%KEEP(201)=0
        ENDIF
      ENDIF
C     ----------------------
C     Broadcast ICNTL(49)
      IF (id%MYID.EQ.MASTER) ICNTL49_LOC=id%ICNTL(49)
      CALL MPI_BCAST( ICNTL49_LOC, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     ----------------------
C
C     Broadcast KEEP options
C     defined for facto:
C     ----------------------
      CALL MPI_BCAST( KEEP(12), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(19), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(21), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(201), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(459), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      CALL MPI_BCAST( KEEP(460), 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
      IF ( KEEP(459) .GE. PANEL_TABSIZE ) THEN
        IF ( LPOK ) THEN
        WRITE(LP,'(A,I4,A,I3)') " ** WARNING ** KEEP(459)=",KEEP(459),
     &  " too large, resetting to",PANEL_TABSIZE-1
        ENDIF
        KEEP(459) = PANEL_TABSIZE - 1
      ENDIF
      PERLU = KEEP(12)
      IF (id%MYID.EQ.MASTER) THEN
C       KEEP(50)  case
C       ==============
C
C       KEEP(50)  = 0 : matrix is unsymmetric
C       KEEP(50) /= 0 : matrix is symmetric
C       KEEP(50) = 1 : Ask L L^T on the root. Matrix is PSD.
C       KEEP(50) = 2 : Ask for L U on the root
C       KEEP(50) = 3 ... L D L^T ??
C     
        CNTL1 = id%CNTL(1)
C       ---------------------------------------
C       For symmetric (non general) matrices 
C       set (directly) CNTL1 = 0.0
C       ---------------------------------------
        KEEP(17)=0
C       Automatic choice if CNTL(1)<0
C       For rank-revealing (KEEP(19).GT.0) then
C       set CNTL1=0.1 even if SYM=1 
        IF (CNTL1.LT.ZERO) THEN
C        automatic choice
         IF (KEEP(19).GT.0) THEN
          CNTL1=0.1E0
         ELSE IF (KEEP(50).EQ.1) THEN
          CNTL1=ZERO
         ELSE
          CNTL1=0.01E0
         ENDIF
        ENDIF
C       CNTL1 threshold value must be between 
C       0.0 and 1.0 (for SYM=0) and 0.5 (for SYM=1,2)
        IF (CNTL1.GT.ONE)   CNTL1=ONE
        IF (CNTL1.LT.ZERO)  CNTL1=ZERO
        IF (KEEP(50).NE.0.AND.CNTL1.GT.0.5E0) THEN
          CNTL1 = 0.5E0
        ENDIF
        PARPIV_T1 = id%KEEP(268)
        IF (PARPIV_T1.EQ.77) THEN
         PARPIV_T1 = 0 
#if defined(__ve__)
         PARPIV_T1 = -2
#endif
        ENDIF
        IF (PARPIV_T1.EQ.-3) THEN
          PARPIV_T1 = 0
        ENDIF
        IF ((PARPIV_T1.LT.-3).OR.(PARPIV_T1.GT.1)) THEN
C        out of range values
         PARPIV_T1 =0
        ENDIF
C       note that KEEP(50).EQ.1 => CNTL1=0.0
        IF (CNTL1.EQ.0.0E0.OR.(KEEP(50).eq.1)) PARPIV_T1 = 0 
C
        IF (PARPIV_T1.EQ.-2) THEN
         IF (KEEP(19).NE.0) THEN
C         switch off PARPIV_T1 if RR activated
C         but do NOT switch off PARPIV_1 with null pivot detection
          PARPIV_T1 = 0
         ENDIF
        ENDIF
        id%KEEP(269) = PARPIV_T1
      ENDIF
      CALL MPI_BCAST(CNTL1, 1, MPI_REAL,
     &             MASTER, id%COMM, IERR)
        CALL MPI_BCAST( KEEP(269), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (id%MYID.EQ.MASTER) THEN
C{
C         OMP parallelization of arrowheads
C         out of range are treated as zero
          IF (KEEP(399).LT.-1) KEEP(399)=-1
          KEEP399_SAVE = KEEP(399)
          IF (KEEP(399).EQ.-1) THEN
            IF ((KEEP(54).EQ.0).AND.(id%NPROCS.GT.1)) THEN
              KEEP(399) = 1
            ELSE
              KEEP(399) = 3
            ENDIF
          ENDIF
#if defined(PCPRET)
C         new multithreaded >=2 algo does not compile on PCPRET
          KEEP(399) = 1
#endif
C         -----------------------------------------------------
C         Decoding of ICNTL(35) for factorization: same as
C         at analysis except that we store a copy of ICNTL(35)
C         in KEEP(486) instead of KEEP(494) and need to check
C         compatibility of KEEP(486) and KEEP(494): If LR was
C         not activated during analysis, it cannot be activated 
C         at factorization.
C         ------------------------------------------------------
          id%KEEP(486) = id%ICNTL(35)
          IF (id%KEEP(486).EQ.1) THEN
C         -- Automatic BLR option setting
           id%KEEP(486)= 2
          ENDIF
          IF ( id%KEEP(486).EQ.4) id%KEEP(486)=0
          IF ((id%KEEP(486).LT.0).OR.(id%KEEP(486).GT.4)) THEN
C           Out of range values treated as 0
            id%KEEP(486) = 0
          ENDIF
          IF ((KEEP(486).NE.0).AND.(KEEP(494).EQ.0)) THEN
C           To activate BLR during factorization,
C           ICNTL(35) must have been set at analysis.
            IF (LPOK)  THEN
              WRITE(LP,'(A)') 
     &      " *** Error with BLR setting "
              WRITE(LP,'(A)') " *** BLR was not activated during ",
     &      " analysis but is requested during factorization."
            ENDIF
            id%INFO(1)=-54
            id%INFO(2)=0
            GOTO 105
          ENDIF
C         Save value of KEEP(486) before possibly
C         forcing it to 3 in case of discard factors
          KEEP_486_FOR_PRINT=KEEP(486)
          IF (KEEP(201) .EQ. -1 .AND. KEEP(486) .NE.0) THEN
            KEEP(486) = 3
          ENDIF
          KEEP464COPY  = id%ICNTL(38)
          IF (KEEP464COPY.LT.0.OR.KEEP464COPY.GT.1000) THEN
C          Out of range values treated as 1000
           KEEP464COPY = 1000
          ENDIF
          IF (id%KEEP(461).LT.1) THEN
            id%KEEP(461) = 10
          ENDIF
          IF (id%KEEP(462).LT.1) THEN
            id%KEEP(462) = 10
          ENDIF
          KEEP465COPY = id%ICNTL(39)
          IF (KEEP465COPY.LT.0.OR.(KEEP465COPY.GT.1000)) THEN
C          Out of range values treated as 1000
           KEEP465COPY = 1000
          ENDIF
          IF (id%ICNTL(36).EQ.1.OR.id%ICNTL(36).EQ.3) THEN
            IF (CNTL1.EQ.ZERO .OR. KEEP(468).LE.1) THEN
              KEEP(475) = 3
            ELSE IF ( (KEEP(269).GT.0).OR. (KEEP(269).EQ.-2)) THEN
              KEEP(475) = 2
            ELSE IF (KEEP(468).EQ.2) THEN
              KEEP(475) = 2
            ELSE
              KEEP(475) = 1
            ENDIF
          ELSE
            KEEP(475) = 0
          ENDIF
          KEEP(481)=0
          IF (id%ICNTL(36).LT.0 .OR. id%ICNTL(36).GE.2) THEN
C           Only options 1 and 2 are allowed
            KEEP(475) = 0
          ENDIF
C         K489 is set according to ICNTL(37)
          IF (id%ICNTL(37).EQ.0.OR.id%ICNTL(37).EQ.1) THEN
            KEEP(489) = id%ICNTL(37)
          ELSE
C           Other values treated as zero
            KEEP(489) = 0
          ENDIF
          IF (KEEP(79).GE.1) THEN
C          CompressCB incompatible with type4,5,6 nodes
           KEEP(489)=0
          ENDIF
C         id%KEEP(476) \in [1,100] 
          IF ((id%KEEP(476).GT.100).OR.(id%KEEP(476).LT.1)) THEN
              id%KEEP(476)=  50
          ENDIF
C         id%KEEP(477) \in [1,100] 
          IF ((id%KEEP(477).GT.100).OR.(id%KEEP(477).LT.1)) THEN
              id%KEEP(477)=  100
          ENDIF
C         id%KEEP(483) \in [1,100] 
          IF ((id%KEEP(483).GT.100).OR.(id%KEEP(483).LT.1)) THEN
              id%KEEP(483)= 50
          ENDIF
C         id%KEEP(484) \in [1,100] 
          IF ((id%KEEP(484).GT.100).OR.(id%KEEP(484).LT.1)) THEN
              id%KEEP(484)=  50
          ENDIF
C         id%KEEP(480)=0,2,3,4,5,6
          IF ((id%KEEP(480).GT.6).OR.(id%KEEP(480).LT.0)
     &                           .OR.(id%KEEP(480).EQ.1)) THEN
               id%KEEP(480)=0
          ENDIF
C         id%KEEP(473)=0 or 1
          IF ((id%KEEP(473).NE.0).AND.(id%KEEP(473).NE.1)) THEN
               id%KEEP(473)=0
          ENDIF
C         id%KEEP(474)=0,1,2,3
          IF ((id%KEEP(474).GT.3).OR.(id%KEEP(474).LT.0)) THEN
               id%KEEP(474)=0
          ENDIF
C         id%KEEP(479)>0
          IF (id%KEEP(479).LE.0) THEN
               id%KEEP(479)=1
          ENDIF
         IF (id%KEEP(474).NE.0.AND.id%KEEP(480).EQ.0) THEN
              id%KEEP(474) = 0
          ENDIF
          IF (id%KEEP(478).NE.0.AND.id%KEEP(480).LT.4) THEN
              id%KEEP(478) = 0
          ENDIF
          IF (id%KEEP(480).GE.5 .OR.
     &           (id%KEEP(480).NE.0.AND.id%KEEP(474).EQ.3)) THEN
            IF (id%KEEP(475).LT.2) THEN
C             Reset to 3 if 5 or to 4 if 6
              id%KEEP(480) = id%KEEP(480) - 2
              write(*,*) ' Resetting KEEP(480) to ', id%KEEP(480)
            ENDIF
          ENDIF
 105    CONTINUE
C}
        ENDIF  ! id%MYID .EQ. MASTER
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
        IF (id%INFO(1).LT.0) GOTO 530
        CALL MPI_BCAST( KEEP(281), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(399), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        CALL MPI_BCAST( KEEP(473), 14, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(486).NE.0) THEN
          CALL MPI_BCAST( KEEP(489), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP464COPY, 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP465COPY, 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
         ENDIF
        IF (id%MYID.EQ.MASTER) THEN
          IF (KEEP(217).GT.2.OR.KEEP(217).LT.0) THEN
            KEEP(217)=0
          ENDIF
          KEEP(214)=KEEP(217)
          IF (KEEP(214).EQ.0) THEN
            IF (KEEP(201).NE.0) THEN ! OOC or no factors
              KEEP(214)=1
            ELSE
              KEEP(214)=2
            ENDIF
            IF (KEEP(486).EQ.2) THEN
              KEEP(214)=1
            ENDIF
          ENDIF
        ENDIF
        CALL MPI_BCAST( KEEP(214), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        IF (KEEP(201).NE.0) THEN
C         -- Low Level I/O strategy
          CALL MPI_BCAST( KEEP(99), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(205), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(211), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
          CALL MPI_BCAST( KEEP(255), 1, MPI_INTEGER,
     &                  MASTER, id%COMM, IERR )
        ENDIF
C     Fwd in facto: explicitly forbid
C     sparse RHS and A-1 computation
      IF (id%KEEP(252).EQ.1 .AND. id%MYID.EQ.MASTER) THEN
        IF (id%ICNTL(20).EQ.1) THEN ! out-of-range => 0
C         NB: in doc ICNTL(20) only accessed during solve
C         In practice, will have failed earlier if RHS not allocated.
C         Still it looks safer to keep this test.
          id%INFO(1)=-43
          id%INFO(2)=20
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Sparse RHS is incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(30).NE.0) THEN ! out-of-range => 1
          id%INFO(1)=-43
          id%INFO(2)=30
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: A-1 functionality incompatible with forward',
     &       ' performed during factorization (ICNTL(32)=1)'
        ELSE IF (id%ICNTL(9) .NE. 1) THEN
          id%INFO(1)=-43
          id%INFO(2)=9
          IF (LPOK) WRITE(LP,'(A)')
     &       ' ERROR: Transpose system (ICNTL(9).NE.0) not ',
     &       ' compatible with forward performed during',
     &       ' factorization (ICNTL(32)=1)'
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C
      IF (id%INFO(1).LT.0) GOTO 530
C
C     The memory allowed is given by ICNTL(23) in Mbytes
C     0 means that nothing is provided.
C     Save memory available, ICNTL(23) in KEEP8(4)
C
      IF ( ICNTL(23) .GT. 0 ) THEN
        ITMP = 1
      ELSE
        ITMP = 0
      ENDIF
      CALL MPI_ALLREDUCE( ITMP, JTMP, 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
      IF ( id%MYID.EQ.MASTER ) THEN
C       Negative values considered 0
        ITMP = max(ICNTL(23),0)
      END IF
      CALL MPI_BCAST( ITMP, 1, MPI_INTEGER,
     &                MASTER, id%COMM, IERR )
C     JTMP: nb of procs with nonzero ICNTL(23)
C     ITMP: value of ICNTL(23) on master
      IF ( ITMP .GT. 0 .AND. JTMP .EQ. 1 ) THEN
C       ICNTL(23)>0 only on master
      ELSE
C       Local values of ICNTL(23) are used, note that
C       they could all be zeros
        ITMP = ICNTL(23)
      ENDIF
C 
      ITMP8 = int(ITMP, 8)
      id%KEEP8(4) = ITMP8 * 1000000_8   ! convert to nb of bytes
C     Compute \sum of memories allowed
      CALL MPI_REDUCE( id%KEEP8(4), ITMP8, 1, MPI_INTEGER8,
     &                  MPI_SUM, MASTER, id%COMM, IERR )
      ITMP8 = ITMP8 / 1000000_8 ! Use to print \sum_{ICNTL(23)}
      IF ( PROKG ) THEN
          NWORKING = id%NSLAVES
          CALL MUMPS_SETI8TOI4( id%KEEP8(129), TMP_INFOG_4)
          WRITE( MPG, 172 ) 
     &    NWORKING, id%ICNTL(22), KEEP_486_FOR_PRINT,
     &    KEEP(489),
     &    id%ICNTL(49),
     &    KEEP(12), 
     &    id%KEEP8(111), TMP_INFOG_4, KEEP(127), KEEP(28),
     &    id%KEEP8(4)/1000000_8, ITMP8, LWK_USER_SUM8, CNTL1
          IF (KEEP(252).GT.0) 
     &    WRITE(MPG,173) KEEP(253)
          IF (KEEP(269).NE.0) 
     &    WRITE(MPG,174) KEEP(269)
      ENDIF
      IF (KEEP(201).LE.0) THEN
C       In-core version or no factors
        KEEP(IXSZ)=XSIZE_IC
      ELSE IF (KEEP(201).EQ.2) THEN
C       OOC version, no panels
        KEEP(IXSZ)=XSIZE_OOC_NOPANEL
      ELSE IF (KEEP(201).EQ.1) THEN
C     Panel versions:
        IF (KEEP(50).EQ.0) THEN
          KEEP(IXSZ)=XSIZE_OOC_UNSYM
        ELSE
          KEEP(IXSZ)=XSIZE_OOC_SYM
        ENDIF
      ENDIF
      IF ( KEEP(486) .NE. 0 ) THEN !LR is activated
C       Stats initialization for LR
        CALL INIT_STATS_GLOBAL(id)
       END IF
C
*     **********************************
*     Begin intializations regarding the
*     computation of the determinant
*     **********************************
      IF (id%MYID.EQ.MASTER) KEEP(258)=ICNTL(33)
      CALL MPI_BCAST(KEEP(258), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF (KEEP(258) .NE. 0) THEN
        KEEP(259) = 0      ! Initial exponent of the local determinant
        KEEP(260) = 1      ! Number of permutations
        id%DKEEP(6)  = 1.0E0  ! real part of the local determinant
      ENDIF
*     ********************************
*     End intializations regarding the
*     computation of the determinant
*     ********************************
C
         CALL MUMPS_STOP_ON_USER_REQUEST(id%KEEP,id%KEEP8, id%ICNTL,
     &   id%INFO, id%MYID)
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0) GOTO 530
*     **********************
*     Begin of Scaling phase
*     **********************
C
C     SCALING MANAGEMENT
C     * Options 1, 3, 4 centralized only
C  
C     * Options 7, 8  : also works for distributed matrix
C
C     At this point, we have the scaling arrays allocated
C     on the master. They have been allocated on the master
C     inside the main MUMPS driver.
C
      CALL MPI_BCAST(KEEP(52), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      LSCAL = ((KEEP(52) .GT. 0) .AND. (KEEP(52) .LE. 8))
      IF (LSCAL) THEN
C
        IF ( id%MYID.EQ.MASTER ) THEN
          CALL MUMPS_SECDEB(TIMEET)
        ENDIF
C       -----------------------
C       Retrieve parameters for
C       simultaneous scaling
C       -----------------------
        IF (KEEP(52) .EQ. 7) THEN 
C       -- Cheap setting of SIMSCALING (it is the default in 4.8.4)
           K231= KEEP(231)
           K232= KEEP(232)
           K233= KEEP(233)
        ELSEIF (KEEP(52) .EQ. 8) THEN
C       -- More expensive setting of SIMSCALING (it was the default in 4.8.1,2,3)
           K231= KEEP(239)
           K232= KEEP(240)
           K233= KEEP(241)
        ENDIF
        CALL MPI_BCAST(id%DKEEP(3),1,MPI_REAL,MASTER,
     &       id%COMM,IERR)
C
        IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND. 
     &       KEEP(54).NE.0 ) THEN
C{        ------------------------------
C         Scaling for distributed matrix
C         We need to allocate scaling
C         arrays on all processors, not
C         only the master.
C         ------------------------------
           IF ( id%MYID .NE. MASTER ) THEN
              IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
              IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
            ALLOCATE( id%ROWSCA(N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=N
            ENDIF
         ENDIF
         M = N
         BUMAXMN=M
         IF(N > BUMAXMN) BUMAXMN = N
         LIWK = 4_8*int(BUMAXMN,8)
         ALLOCATE (IWK(LIWK), BURP(M), BUCP(N),
     &            BURS(2* id%NPROCS), BUCS(2* id%NPROCS),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(LIWK+int(M,8)+int(N,8)+
     &      4_8*int(id%NPROCS,8) , id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 1
C        -- LWK not used
         LWK_REAL = 1_8
         LWRKR_TH = 1_8
         LWRKC_TH = 1_8
         NOMP_MAX   = 1
         ALLOCATE(WK_REAL(LWK_REAL),
     &            WRKR_TH(LWRKR_TH,NOMP_MAX),
     &            WRKC_TH(LWRKC_TH,NOMP_MAX),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(
     &             LWK_REAL+ 
     &             LWRKR_TH*int(NOMP_MAX,8)+
     &             LWRKC_TH*int(NOMP_MAX,8),
     &             id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
         CALL SMUMPS_SIMSCALEABS(
     &        IRN_loc_PTR(1), JCN_loc_PTR(1), A_loc_PTR(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL, 
     &        WRKR_TH, LWRKR_TH, 
     &        WRKC_TH, LWRKC_TH,
     &        NOMP_MAX, id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         IF(LIWK < BUINTSZ) THEN
            DEALLOCATE(IWK)
            LIWK = BUINTSZ
            ALLOCATE(IWK(LIWK), stat=allocok)
            IF (allocok > 0) THEN
               id%INFO(1)=-13
               CALL MUMPS_SET_IERROR(LIWK,id%INFO(2))
            ENDIF
         ENDIF
         DEALLOCATE(WK_REAL, WRKR_TH, WRKC_TH)
         LWK_REAL     = BURESZ
C
C        -- Set NOMP_MAX from KEEP(281)
         CALL SMUMPS_SET_NOMP_MAX(id%KEEP(281), id%KEEP(361), 
     &                            N, NOMP_MAX)
C
         IF (NOMP_MAX.LE.1) THEN
C          temp array per thread not used
           LWRKR_TH = 1
           LWRKC_TH = 1
         ELSE
           LWRKR_TH = N
           IF (id%KEEP(50).NE.0) THEN
C            WRKC_TH not used on symmetric matrices
             LWRKC_TH = 1
           ELSE
             LWRKC_TH = N
           ENDIF
         ENDIF
         ALLOCATE(WK_REAL(LWK_REAL),
     &            WRKR_TH(LWRKR_TH,max(NOMP_MAX,1)),
     &            WRKC_TH(LWRKC_TH,max(NOMP_MAX,1)),
     &            stat=allocok)
         IF (allocok > 0) THEN
            id%INFO(1)=-13
            CALL MUMPS_SET_IERROR(
     &             LWK_REAL+ 
     &             LWRKR_TH*int(max(NOMP_MAX,1),8)+
     &             LWRKC_TH*int(max(NOMP_MAX,1),8),
     &             id%INFO(2))
         ENDIF
C        --- Propagate enventual error
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &        id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 517
C        -- estimation of memory and construction of partvecs
         BUJOB = 2
         CALL SMUMPS_SIMSCALEABS(
     &        IRN_loc_PTR(1), JCN_loc_PTR(1), A_loc_PTR(1),
     &        id%KEEP8(29),
     &        M, N,  id%NPROCS, id%MYID, id%COMM,
     &        BURP, BUCP,
     &        BURS, BUCS, BUREGISTRE,
     &        IWK, LIWK,
     &        BUINTSZ, BURESZ, BUJOB,
     &        id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &        WRKR_TH, LWRKR_TH, 
     &        WRKC_TH, LWRKC_TH,
     &        NOMP_MAX, id%KEEP(50),
     &        K231, K232, K233, 
     &        id%DKEEP(3),
     &        SCONEERR, SCINFERR)
         id%DKEEP(4) = SCONEERR
         id%DKEEP(5) = SCINFERR
CXXXX 
         DEALLOCATE(IWK, BURP,BUCP,BURS, BUCS)
         DEALLOCATE(WK_REAL, WRKR_TH, WRKC_TH)
C}
        ELSE IF ( KEEP(54) .EQ. 0 ) THEN
C{        ------------------
C         Centralized matrix
C         ------------------
          IF ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8))  THEN
C             -------------------------------
C             Create a communicator of size 1
C             -------------------------------
              IF (id%MYID.EQ.MASTER) THEN
                COLOUR = 0
              ELSE
                COLOUR = MPI_UNDEFINED
              ENDIF
              CALL MPI_COMM_SPLIT( id%COMM, COLOUR, 0,
     &             COMM_FOR_SCALING, IERR )
              IF (id%MYID.EQ.MASTER) THEN
                 M = N
                 BUMAXMN=N
                 IF(N > BUMAXMN) BUMAXMN = N
                 LIWK = 1_8
                 ALLOCATE(IWK(LIWK), BURP(1), BUCP(1),
     &                BURS(1), BUCS(1),
     &                stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    CALL MUMPS_SET_IERROR( LIWK+4_8, id%INFO(2) )
                    GOTO 400
                 ENDIF
                 LWK_REAL = int(M,8) + int(N,8)  
C
C                -- Set NOMP_MAX from KEEP(281)
                 CALL SMUMPS_SET_NOMP_MAX(id%KEEP(281), id%KEEP(361), 
     &                                    N, NOMP_MAX)
C
                 IF (NOMP_MAX.LE.1) THEN
C                  temp array per thread not used
                   LWRKR_TH = 1
                   LWRKC_TH = 1
                 ELSE
                   LWRKR_TH = N
                   IF (id%KEEP(50).NE.0) THEN
C                    WRKC_TH not used on symmetric matrices
                     LWRKC_TH = 1
                   ELSE
                     LWRKC_TH = N
                   ENDIF
                 ENDIF
                 ALLOCATE(WK_REAL(LWK_REAL),
     &                    WRKR_TH(LWRKR_TH,max(NOMP_MAX,1)),
     &                    WRKC_TH(LWRKC_TH,max(NOMP_MAX,1)),
     &                    stat=allocok)
                 IF (allocok > 0) THEN
                    id%INFO(1)=-13
                    CALL MUMPS_SET_IERROR(
     &                     LWK_REAL+ 
     &                     LWRKR_TH*int(max(NOMP_MAX,1),8)+
     &                     LWRKC_TH*int(max(NOMP_MAX,1),8),
     &                     id%INFO(2))
                 ENDIF
                 CALL MPI_COMM_RANK(COMM_FOR_SCALING, SCMYID, IERR)
                 CALL MPI_COMM_SIZE(COMM_FOR_SCALING, SCNPROCS, IERR)
                 BUJOB = 1
                 CALL SMUMPS_SIMSCALEABS(
     &                id%IRN(1), id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                WRKR_TH, LWRKR_TH, 
     &                WRKC_TH, LWRKC_TH,
     &                NOMP_MAX, id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 IF(LWK_REAL < BURESZ) THEN
                    id%INFO(1) = -136
                    GOTO 400
                 ENDIF
                 BUJOB = 2
                 CALL SMUMPS_SIMSCALEABS(id%IRN(1),
     &                id%JCN(1), id%A(1),
     &                id%KEEP8(28),
     &                M, N,  SCNPROCS, SCMYID, COMM_FOR_SCALING,
     &                BURP, BUCP,
     &                BURS, BUCS, BUREGISTRE,
     &                IWK, LIWK,
     &                BUINTSZ, BURESZ, BUJOB,
     &                id%ROWSCA(1), id%COLSCA(1), WK_REAL, LWK_REAL,
     &                WRKR_TH, LWRKR_TH, 
     &                WRKC_TH, LWRKC_TH,
     &                NOMP_MAX, id%KEEP(50),
     &                K231, K232, K233, 
     &                id%DKEEP(3),
     &                SCONEERR, SCINFERR)
                 id%DKEEP(4) = SCONEERR
                 id%DKEEP(5) = SCINFERR
 400             CONTINUE
                 IF (allocated(WK_REAL)) DEALLOCATE(WK_REAL)
                 IF (allocated(WRKR_TH)) DEALLOCATE(WRKR_TH)
                 IF (allocated(WRKC_TH)) DEALLOCATE(WRKC_TH)
                 IF (allocated(IWK))     DEALLOCATE(IWK)
                 IF (allocated(BURP))    DEALLOCATE(BURP)
                 IF (allocated(BUCP))    DEALLOCATE(BUCP)
                 IF (allocated(BURS))    DEALLOCATE(BURS)
                 IF (allocated(BUCS))    DEALLOCATE(BUCS)
              ENDIF
C             Centralized matrix: make DKEEP(4:5) available to all processors
              CALL MPI_BCAST( id%DKEEP(4),2,MPI_REAL,
     &                        MASTER, id%COMM, IERR )
              IF (id%MYID.EQ.MASTER) THEN
C               Communicator should only be
C               freed on the master process
                CALL MPI_COMM_FREE(COMM_FOR_SCALING, IERR)
              ENDIF
              CALL MUMPS_PROPINFO(ICNTL(1), id%INFO(1),
     &             id%COMM, id%MYID)
              IF (id%INFO(1).LT.0) GOTO 517
          ELSE IF (id%MYID.EQ.MASTER) THEN
C           -------------------
C           Centralized scaling
C           -------------------
            IF (KEEP(52).GT.0 .AND. KEEP(52).LE.6) THEN
C             ---------------------
C             Allocate temporary
C             workspace for scaling
C             ---------------------
              IF (KEEP(52) .eq. 1 ) THEN
C               No workspace indeed needed
                LWK = 1_8
                LWK_REAL = 1_8
              ELSE IF ( KEEP(52) .eq. 3 ) THEN 
                LWK = 1_8
                LWK_REAL = int(N,8)
              ELSE IF ( KEEP(52) .eq. 4 ) THEN
C               Options 3 or 4
                LWK = 1_8
                LWK_REAL = 2_8*int(N,8)
              END IF
C             Real workarray
              ALLOCATE( WK_REAL( LWK_REAL ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                CALL MUMPS_SET_IERROR( LWK_REAL, id%INFO(2) )
                GOTO 137
              END IF
C             Real/complex workarray
              ALLOCATE( WK( LWK ), stat = IERR )
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
                GOTO 137
              END IF
              CALL SMUMPS_FAC_A(N, id%KEEP8(28), KEEP(52), id%A(1),
     &             id%IRN(1), id%JCN(1),
     &             id%COLSCA(1), id%ROWSCA(1),
     &             WK, LWK, WK_REAL, LWK_REAL, ICNTL(1), id%INFO(1) )
              DEALLOCATE( WK_REAL )
              DEALLOCATE( WK )
            ENDIF
          ENDIF
C}
        ENDIF ! Scaling distributed matrices or centralized
        IF (KEEP(125).NE.0) THEN
C         ------------------------
C           If we enable the scaling of the |A11 A12| block
C           we set to 1 the scaling corresponding to the Schur
C           complement matrix A22
C         ------------------------
            IF ((KEEP(60).GT.0) .and. (KEEP(116).GT.0)) THEN
C            Schur is active, reset Schur entries to ONE
             IF ( ((KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)) .AND.
     &          KEEP(54).NE.0 ) THEN
C             Scaling available on all procs
              DO I=1, N
               IF (id%SYM_PERM(I).GT.id%N-KEEP(116)) THEN
                 id%COLSCA(I) = ONE
                 id%ROWSCA(I) = ONE
               ENDIF
              ENDDO
             ELSE IF ( id%MYID .EQ. MASTER) THEN
C             Scaling available on master
              DO I=1, N
               IF (id%SYM_PERM(I).GT.id%N-KEEP(116)) THEN
                 id%COLSCA(I) = ONE
                 id%ROWSCA(I) = ONE
               ENDIF
              ENDDO
             ENDIF
            ENDIF
        ENDIF 
        IF (id%MYID.EQ.MASTER) THEN
            CALL MUMPS_SECFIN(TIMEET)
            id%DKEEP(92)=real(TIMEET)
          IF (PROKG) WRITE( MPG, 140 ) TIMEET
C         Print inf-norm after last KEEP(233) iterations of
C         scaling option KEEP(52)=7 or 8 (SimScale)
C
          IF (PROKG.AND.(KEEP(52).EQ.7.OR.KEEP(52).EQ.8) 
     &             .AND. (K233+K231+K232).GT.0) THEN
           IF (K232.GT.0) WRITE(MPG, 166) id%DKEEP(4)
          ENDIF
        ENDIF
      ENDIF ! LSCAL
C
C       scaling might also be provided by the user
        LSCAL = (LSCAL .OR. (KEEP(52) .EQ. -1) .OR. KEEP(52) .EQ. -2)
        IF (LSCAL .AND. KEEP(258).NE.0 .AND. id%MYID .EQ. MASTER) THEN
          DO I = 1, id%N
            CALL SMUMPS_UPDATEDETER_SCALING(id%ROWSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
          ENDDO
          IF (KEEP(50) .EQ. 0) THEN ! unsymmetric
            DO I = 1, id%N
              CALL SMUMPS_UPDATEDETER_SCALING(id%COLSCA(I),
     &           id%DKEEP(6),    ! determinant
     &           KEEP(259))   ! exponent of the determinant
            ENDDO
          ELSE
C           -----------------------------------------
C           In this case COLSCA = ROWSCA
C           Since determinant was initialized to 1,
C           compute square of the current determinant
C           rather than going through COLSCA.
C           -----------------------------------------
            CALL SMUMPS_DETER_SQUARE(id%DKEEP(6), KEEP(259))
          ENDIF
C         Now we should have taken the
C         inverse of the scaling vectors
          CALL SMUMPS_DETER_SCALING_INVERSE(id%DKEEP(6), KEEP(259))
        ENDIF
C
C       ********************
C       End of Scaling phase
C       At this point: either (matrix is distributed and KEEP(52)=7 or 8)
C       in which case scaling arrays are allocated on all processors,
C       or scaling arrays are only on the host processor.
C       In case of distributed matrix input, we will free the scaling
C       arrays on procs with MYID .NE. 0 after the all-to-all distribution
C       of the original matrix.
C       ********************
C
 137  CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
      IF (id%INFO(1).LT.0) GOTO 517
C     Fwd in facto: in case of repeated factorizations
C     with different Schur options we prefer to free
C     systematically this array now than waiting for
C     the root node. We rely on the fact that it is
C     allocated or not during the solve phase so if
C     it was allocated in a 1st call to facto and not
C     in a second, we don't want the solve to think
C     it was allocated in the second call.
      IF (associated(id%root%RHS_CNTR_MASTER_ROOT)) THEN
        DEALLOCATE (id%root%RHS_CNTR_MASTER_ROOT)
        NULLIFY (id%root%RHS_CNTR_MASTER_ROOT)
      ENDIF
C     Fwd in facto: check that id%NRHS has not changed
      IF ( id%MYID.EQ.MASTER.AND. KEEP(252).EQ.1 .AND.
     &      id%NRHS .NE. id%KEEP(253) ) THEN
C         Error: NRHS should not have
C         changed since the analysis
          id%INFO(1)=-42
          id%INFO(2)=id%KEEP(253)
      ENDIF
C     Fwd in facto: allocate and broadcast RHS_MUMPS
C     to make it available on all processors.
      IF (id%KEEP(252) .EQ. 1) THEN
          IF ( id%MYID.NE.MASTER ) THEN
            id%KEEP(254) = N              ! Leading dimension
            id%KEEP8(85) = int(N,8)*int(id%KEEP(253),8) ! Tot size
            ALLOCATE(RHS_MUMPS(id%KEEP8(85)),stat=IERR)
            IF (IERR > 0) THEN
               id%INFO(1)=-13
               CALL MUMPS_SET_IERROR(id%KEEP8(85), id%INFO(2))
               IF (LPOK)
     &         WRITE(LP,*) 'ERROR while allocating RHS on a slave'
               NULLIFY(RHS_MUMPS)
            ENDIF
            RHS_MUMPS_ALLOCATED = .TRUE.
          ELSE 
C           Case of non working master
            id%KEEP(254)=id%LRHS              ! Leading dimension
            id%KEEP8(85)=int(id%LRHS,8)*int(id%KEEP(253)-1,8)+
     &                   int(id%N,8) ! Tot size
            RHS_MUMPS=>id%RHS
            RHS_MUMPS_ALLOCATED = .FALSE.
            IF (LSCAL) THEN
C             Scale before broadcast: apply row
C             scaling (remark that we assume no
C             transpose).
              DO K=1, id%KEEP(253)
                DO I=1, N
                  RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          = RHS_MUMPS( id%KEEP(254) * (K-1) + I )
     &          * id%ROWSCA(I)
                ENDDO
              ENDDO
            ENDIF
          ENDIF
      ELSE
          id%KEEP8(85)=1_8
          ALLOCATE(RHS_MUMPS(1),stat=IERR)
          IF (IERR > 0) THEN
             id%INFO(1)=-13
             id%INFO(2)=1
             IF (LPOK)
     &            WRITE(LP,*) 'ERREUR while allocating RHS on a slave'
             NULLIFY(RHS_MUMPS)
          ENDIF
          RHS_MUMPS_ALLOCATED = .TRUE.
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 516
      IF (KEEP(252) .EQ. 1) THEN
C
C         Broadcast the columns of the right-hand side
C         one by one. Leading dimension is keep(254)=N
C         on procs with MYID > 0 but may be larger on
C         the master processor.
          DO I= 1, id%KEEP(253)
            CALL MPI_BCAST(RHS_MUMPS((I-1)*id%KEEP(254)+1), N,
     &           MPI_REAL, MASTER,id%COMM,IERR)
          END DO
      ENDIF
      IF (id%MYID.EQ. MASTER) THEN
C       Copy the value of ICNTL(24) and make it
C       available on all working processors.
        KEEP(110)=id%ICNTL(24)
C       KEEP(110) defaults to 0 for out of range values
        IF (KEEP(110).NE.1) KEEP(110)=0
      ENDIF
      CALL MPI_BCAST(KEEP(110), 1, MPI_INTEGER,
     &               MASTER, id%COMM, IERR)
      IF ((KEEP(50).EQ.2).AND.(KEEP(219).NE.0)) THEN
       CALL SMUMPS_BUF_MAX_ARRAY_MINSIZE(max(KEEP(108),1),IERR)
       IF (IERR .NE. 0) THEN
C      ------------------------
C      Error allocating SMUMPS_BUF
C      ------------------------
          id%INFO(1) = -13
          id%INFO(2) = max(KEEP(108),1)
       END IF
      ENDIF
C     -----------------------------------------------
C     Depending on the option used for 
C       -detecting null pivots (ICNTL(24)/KEEP(110))
C         CNTL(3) is used to set DKEEP(1)
C               ( A row is considered as null if ||row|| < DKEEP(1) )
C         CNTL(5) is then used to define if a large 
C                 value is set on the diagonal or if a 1 is set
C                 and other values in the row are reset to zeros.
C       -rank revealing on the Schur (ICNTL(56)/KEEP(19))
C         SEUIL* corresponds to the minimum required 
C                absolute value of pivot.
C         SEUIL_LDLT_NIV2 is used only in the 
C                case of SYM=2 within a niv2 node for which 
C                we have only a partial view of the fully summed rows.
      IF (id%MYID .EQ. MASTER) CNTL3 = id%CNTL(3)
      CALL MPI_BCAST(CNTL3, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL5 = id%CNTL(5)
      CALL MPI_BCAST(CNTL5, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) CNTL6 = id%CNTL(6)
      CALL MPI_BCAST(CNTL6, 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      IF (id%MYID .EQ. MASTER) id%DKEEP(8) = id%CNTL(7)
      CALL MPI_BCAST(id%DKEEP(8), 1, MPI_REAL,
     &               MASTER, id%COMM, IERR)
      id%DKEEP(11) = id%DKEEP(8)/id%KEEP(461)
      id%DKEEP(12) = id%DKEEP(8)/id%KEEP(462)
      IF (KEEP(486).EQ.0) id%DKEEP(8) = ZERO
      COMPUTE_ANORMINF = .FALSE.
      IF ( (KEEP(486) .NE. 0).AND. (id%DKEEP(8).LT.ZERO)) THEN
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(19).NE.0) THEN
C       Rank revealing factorisation
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (KEEP(110).NE.0) THEN
C       Null pivot detection       
        COMPUTE_ANORMINF = .TRUE.
      ENDIF
      IF (id%DKEEP(8).LT.ZERO) THEN
C       Experimental setting of CNTL(7)
        IF (COMPUTE_ANORMINF) THEN
         EFF_SIZE_SCHUR = 0
         CALL SMUMPS_ANORMINF(  id , ANORMINF, LSCAL, EFF_SIZE_SCHUR )
C        If no schur ANORMINF fine for other cases
        ELSE
         ANORMINF = ZERO
        ENDIF
        id%DKEEP(8) = abs(id%DKEEP(8))*ANORMINF
C       ANORMINF need be recomputed in case of schur
        IF ((KEEP(60).GT.0).AND.KEEP(116).GT.0) ANORMINF=ZERO
      ENDIF
C     -------------------------------------------------------
C        We compute ANORMINF, when needed, based on
C        the infinite norm of Rowsca *A*Colsca
C        and make it available on all working processes.
      IF (COMPUTE_ANORMINF) THEN
         EFF_SIZE_SCHUR = 0
         IF (KEEP(60).GT.0) EFF_SIZE_SCHUR = KEEP(116)
         CALL SMUMPS_ANORMINF(  id , ANORMINF, LSCAL, EFF_SIZE_SCHUR )
      ELSE
         ANORMINF = ZERO
      ENDIF
C
       IF ((KEEP(19).NE.0).OR.(KEEP(110).NE.0)) THEN 
         IF (PROKG) THEN 
           IF (KEEP(19).NE.0) THEN
           WRITE(MPG,'(A,1PD16.4)')
     &    ' CNTL(3) for null pivot rows/singularities  =',CNTL3      
           ELSE
           WRITE(MPG,'(A,1PD16.4)')
     &    ' CNTL(3) for null pivot row detection       =',CNTL3      
           ENDIF
         ENDIF
       ENDIF
       IF (KEEP(19).EQ.0) THEN 
C        -- RR is off
         SEUIL = ZERO
         id%DKEEP(9) = ZERO
       ELSE
C        -- RR is on
C
C      CNTL(3) is the threshold used in the following to compute  
C      DKEEP(9) the threshold under which the sing val. are considered
C      as null and from which we start to look for a gap between two
C      sing val.
         IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(9) = abs(CNTL(3))
         ELSE IF  (CNTL3 .GT. ZERO) THEN
           id%DKEEP(9) = CNTL3*ANORMINF
         ELSE  !  (CNTL(3) .EQ. ZERO) THEN
           CALL MUMPS_NPIV_CRITICAL_PATH(
     &     N, KEEP(28), id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &     id%NA(1), id%LNA, id%NE_STEPS(1), NPIV_CRITICAL_PATH )
           id%DKEEP(9) = sqrt(REAL(NPIV_CRITICAL_PATH))*EPS*ANORMINF  
         ENDIF
         IF (PROKG) THEN 
           WRITE(MPG, '(A,I16)')
     &    ' ICNTL(56) rank revealing effective value   =',KEEP(19)   
           WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Threshold for singularities on the root =',id%DKEEP(9)
         ENDIF
C        RR postponing considers that pivot rows with norm smaller 
C        than SEUIL should be postponed. 
C        SEUIL should be bigger than DKEEP(9), this means that 
C        DKEEP(13) should be bigger than 1. 
         Thresh_Seuil = id%DKEEP(13)
         IF (id%DKEEP(13).LT.1)  Thresh_Seuil = 10
         SEUIL = id%DKEEP(9)*Thresh_Seuil
         IF (PROKG) WRITE(MPG,'(A,1PD16.4)')
     &   ' ...Threshold for postponing                =',SEUIL
       ENDIF !end KEEP(19)
       SEUIL_LDLT_NIV2 = SEUIL
C     -------------------------------
C     -- Null pivot row detection
C     -------------------------------
       IF (KEEP(110).EQ.0) THEN 
C        -- Null pivot is off
C        Initialize DKEEP(1) to a negative value
C        in order to avoid detection of null pivots
C        (test max(AMAX,RMAX,abs(PIVOT)).LE.PIVNUL
C        in SMUMPS_FAC_I, where PIVNUL=DKEEP(1))
         id%DKEEP(1) = -1.0E0
         id%DKEEP(2) = ZERO
       ELSE
C        -- Null pivot is on
        IF (KEEP(19).NE.0) THEN
C        -- RR is on
C      RR postponing considers that pivot rows of norm smaller that SEUIL 
C      should be postponed, but pivot rows smaller than DKEEP(1) are 
C      directly added to null space and thus considered as null pivot rows. 
         IF ((id%DKEEP(10).LE.0).OR.(id%DKEEP(10).GT.1)) THEN
C           DKEEP(10) is out of range, set to the default value 10-1
           id%DKEEP(1) = id%DKEEP(9)*1E-1
         ELSE
           id%DKEEP(1) = id%DKEEP(9)*id%DKEEP(10) 
         ENDIF    
        ELSE
C        -- RR is off
C        -- only Null pivot detection
C         We keep strategy currently used in MUMPS 4.10.0
          IF (CNTL3 .LT. ZERO) THEN
           id%DKEEP(1)  = abs(CNTL(3))
          ELSE IF  (CNTL3 .GT. ZERO) THEN
            id%DKEEP(1)  = CNTL3*ANORMINF
          ELSE !  (CNTL(3) .EQ. ZERO) THEN
c          id%DKEEP(1)  = NPIV_CRITICAL_PATH*EPS*ANORMINF
           CALL MUMPS_NPIV_CRITICAL_PATH(
     &         N, KEEP(28), id%STEP(1), id%FRERE_STEPS(1), id%FILS(1),
     &         id%NA(1), id%LNA, id%NE_STEPS(1), NPIV_CRITICAL_PATH )
           id%DKEEP(1)  = sqrt(REAL(NPIV_CRITICAL_PATH))*EPS*ANORMINF
          ENDIF
        ENDIF ! fin rank revealing        
        IF ((KEEP(110).NE.0).AND.(PROKG)) THEN
           WRITE(MPG, '(A,I16)')
     &    ' ICNTL(24) null pivot rows detection        =',KEEP(110)   
           WRITE(MPG,'(A,1PD16.4)')
     &    ' ...Zero pivot detection threshold          =',id%DKEEP(1)
        ENDIF 
        IF (CNTL5.GT.ZERO) THEN
            id%DKEEP(2) = CNTL5 * ANORMINF
            IF (PROKG) WRITE(MPG,'(A,1PD10.3)')
     &    ' ...Fixation for null pivots               =',id%DKEEP(2)
         ELSE
            IF (PROKG) WRITE(MPG,*) '...Infinite fixation '
            IF (id%KEEP(50).EQ.0) THEN
C             Unsym
            ! the user let us choose a fixation. set in NEGATIVE
            ! to detect during facto when to set row to zero !
             id%DKEEP(2) = -max(1.0E10*ANORMINF, 
     &                sqrt(huge(ANORMINF))/1.0E8)
            ELSE
C             Sym
            id%DKEEP(2) = ZERO
            ENDIF
         ENDIF
       ENDIF ! fin null pivot detection.
C     Find id of root node if RR is on 
      IF (KEEP(53).NE.0) THEN
        ID_ROOT =MUMPS_PROCNODE(id%PROCNODE_STEPS(id%STEP(KEEP(20))),
     &                          id%KEEP(199))
        IF ( KEEP( 46 )  .NE. 1 ) THEN
          ID_ROOT = ID_ROOT + 1
        END IF
      ENDIF
C Second pass:  set parameters for null pivot detection
C Allocate PIVNUL_LIST_STRUCT in case of null pivot detection
C and in case of rank revealing
      KEEP(109) = 0
      LPN_LIST  = 0
      IF(KEEP(110) .EQ. 1) THEN
         LPN_LIST = 100
      ENDIF
      IF (KEEP(19).NE.0 .AND.
     &   (ID_ROOT.EQ.id%MYID .OR. id%MYID.EQ.MASTER)) THEN
         LPN_LIST = 100
      ENDIF
      IF (LPN_LIST.GT.0) THEN
        PIVNUL_LIST_STRUCT%SIZE_PIVNUL_LIST = LPN_LIST
        ALLOCATE( PIVNUL_LIST_STRUCT%PIVNUL_LIST(LPN_LIST),stat = IERR )
        IF ( IERR .GT. 0 ) THEN
           id%INFO(1)=-13
           id%INFO(2)=LPN_LIST
        END IF
        PIVNUL_LIST_STRUCT%PIVNUL_LIST(1:LPN_LIST) = 0
      ENDIF
C end set parameter for null pivot detection
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 516
C   --------------------------------------------------------------
C   STATIC PIVOTING
C     -- Static pivoting only when RR and Null pivot detection OFF
C   --------------------------------------------------------------
      KEEP(97) = 0
      IF ((KEEP(19).EQ.0).AND.(KEEP(110).EQ.0)) THEN
        IF (id%MYID .EQ. MASTER) CNTL4 = id%CNTL(4)
        CALL MPI_BCAST( CNTL4, 1, MPI_REAL,
     &                MASTER, id%COMM, IERR )
C 
        IF ( CNTL4 .GE. ZERO ) THEN
         KEEP(97) = 1
         IF ( CNTL4 .EQ. ZERO ) THEN
C           -- set seuil to sqrt(eps)*||A||
            IF(ANORMINF .EQ. ZERO) THEN
               EFF_SIZE_SCHUR = 0
               IF (KEEP(60).GT.0) EFF_SIZE_SCHUR = KEEP(116)
               CALL SMUMPS_ANORMINF( id , ANORMINF, LSCAL, 
     &              EFF_SIZE_SCHUR )
            ENDIF
            SEUIL = sqrt(EPS) * ANORMINF
         ELSE
            SEUIL = CNTL4
         ENDIF
         SEUIL_LDLT_NIV2 = SEUIL
        ELSE 
         SEUIL = ZERO
        ENDIF
      ENDIF
C     set number of tiny pivots / 2x2 pivots in types 1 /
C     2x2 pivots in types 2, to zero. This is because the
C     user can call the factorization step several times.
      KEEP(98)  = 0
      KEEP(103) = 0
      KEEP(105) = 0
      MAXS      = 1_8
*
*     Start allocations
*     *****************
*
C
C  The slaves can now perform the factorization
C
C
C  Allocate id%S on all nodes
C  or point to user provided data WK_USER when LWK_USER>0
C  =======================
C
C     Compute BLR_STRAT and a first estimation 
C     of MAXS, the size of id%S
      CALL  SMUMPS_SET_BLRSTRAT_AND_MAXS_K8 (
     &           MAXS_BASE8, MAXS_BASE_RELAXED8,
     &           BLR_STRAT,
     &           id%KEEP(1), id%KEEP8(1))
C   
      MAXS = MAXS_BASE_RELAXED8
      IF (WK_USER_PROVIDED) THEN
C       -- Set MAXS to size of WK_USER_
        MAXS = id%KEEP8(24)
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 516
      ENDIF
C
      id%KEEP8(75) = huge(id%KEEP8(75))
      id%KEEP8(76) = huge(id%KEEP8(76))
      IF (I_AM_SLAVE) THEN
C
        IF (id%KEEP8(4) .NE. 0_8) THEN
C
          IF ( .NOT. WK_USER_PROVIDED ) THEN
C           Set MAXS given BLR_STRAT, KEEP(201) and MAXS_BASE_RELAXED8
            CALL SMUMPS_MEM_ALLOWED_SET_MAXS (
     &          MAXS,
     &          BLR_STRAT, id%KEEP(201), MAXS_BASE_RELAXED8,
     &          id%KEEP(1), id%KEEP8(1), id%MYID, id%N, id%NELT, 
     &          id%NA(1), id%LNA, id%NSLAVES, 
     &          KEEP464COPY,  KEEP465COPY,
     &          id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
C           Given  MAXS and  max memory allowed KEEP8(4)
C           compute in KEEP8(75) the number of real/complex 
C           available for dynamic allocations 
            CALL  SMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           MAXS, id%MYID, 
     &           .FALSE., ! UNDER_L0_OMP
     &           N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &           BLR_STRAT, id%KEEP(201), 
     &           id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
          ELSE
C            KEEP8(75) dow not include MAXS, since WK_USER is provided
             CALL  SMUMPS_MEM_ALLOWED_SET_K75 ( 
     &           0_8, id%MYID, 
     &           .FALSE., ! UNDER_L0_OMP
     &           N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &           BLR_STRAT, id%KEEP(201), 
     &           id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &           , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &           size(id%I8_L0_OMP,2)
     &           )
          ENDIF
          IF (KEEP(400) .GT.0) THEN
C            ------------------------------
C            compute KEEP8(75) under L0_OMP
C            ------------------------------
C            Save KEEP8(75) above L0_OMP to reset KEEP8(75) 
C            when starting FAC_PAR_M
             id%KEEP8(76) = id%KEEP8(75)
             CALL  SMUMPS_MEM_ALLOWED_SET_K75 ( 
     &             0_8, ! MAXS=0_8 
     &             id%MYID, 
     &             .TRUE., ! UNDER_L0_OMP
     &             id%N, id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &             BLR_STRAT, id%KEEP(201), 
     &             id%KEEP, id%KEEP8, id%INFO(1), id%INFO(2)
     &             , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &             size(id%I8_L0_OMP,2)
     &             )
C            KEEP8(75) holds the number of entries that 
C            can be allocated underL0.
C            It will be used during SMUMPS_FAC_L0_OMP to adjust the 
C            the size of MUMPS_TPS_ARR(ITH)%LA 
          ENDIF
        ENDIF ! MEM_ALLOWED
C
      ENDIF ! I_AM_SLAVE THEN
C
      IF (I_AM_SLAVE) THEN
        IF ( (KEEP(400).GT.0) .AND. (KEEP(406).EQ.2) ) THEN
C         Compute KEEP8(77) the peak authorized used by 
C         SMUMPS_PERFORM_COPIES
          CALL SMUMPS_L0_COMPUTE_PEAK_ALLOWED(
     &          id%MYID, id%N,
     &          id%NELT, id%NA(1), id%LNA, id%NSLAVES,
     &          BLR_STRAT, id%KEEP(201),
     &          id%KEEP(1), id%KEEP8(1), id%INFO(1), id%INFO(2)
     &          , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &          size(id%I8_L0_OMP,2)
     &           )
        ENDIF
      ENDIF ! I_AM_SLAVE)
C
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) .LT. 0) THEN
        GOTO 516
      ENDIF
      CALL MUMPS_SETI8TOI4(MAXS, id%INFO(39))
      CALL SMUMPS_AVGMAX_STAT8(PROKG, MPG, MAXS, id%NSLAVES,
     & PRINT_MAXAVG,
     & id%COMM, " Effective size of S     (based on INFO(39))=   ")
C
      IF ( I_AM_SLAVE ) THEN
C       ------------------
C       Dynamic scheduling
C       ------------------
        CALL SMUMPS_LOAD_SET_INICOST( dble(id%COST_SUBTREES),
     &        KEEP(64), id%DKEEP(15), KEEP(375), MAXS )
        MEMORY_MD_ARG = min(int(PERLU,8) * ( MAXS_BASE8 / 100_8 + 1_8 ),
C       Restrict freedom from dynamic scheduler when
C       MEM_ALLOWED=ICNTL(23) is small (case where KEEP8(4)-MAXS_BASE8
C       is negative after call to SMUMPS_MAX_MEM)
     &                      max(0_8, MAXS-MAXS_BASE8))
        CALL SMUMPS_LOAD_INIT( id, MEMORY_MD_ARG, MAXS )
C
C       Out-Of-Core (OOC) issues. Case where we ran one factorization OOC
C       and the second one is in-core: we try to free OOC
C       related data from previous factorization.
C
        CALL SMUMPS_CLEAN_OOC_DATA(id, IERR)
        IF (IERR < 0) THEN
          id%INFO(1) = -90
          id%INFO(2) = 0
          GOTO 112
        ENDIF
        IF (KEEP(201) .GT. 0) THEN
C          -------------------
C          OOC initializations
C          -------------------
           IF (KEEP(201).EQ.1 !PANEL Version
     &         .AND.KEEP(50).EQ.0 ! Unsymmetric
     &         .AND.KEEP(251).NE.2 ! Store L to disk
     &         ) THEN
              id%OOC_NB_FILE_TYPE=2 ! declared in MUMPS_OOC_COMMON
           ELSE
              id%OOC_NB_FILE_TYPE=1 ! declared in MUMPS_OOC_COMMON
           ENDIF
C          ------------------------------
C          Dimension IO buffer, KEEP(100)
C          ------------------------------
           IF (KEEP(205) .GT. 0) THEN
             KEEP(100) = KEEP(205)
           ELSE
             IF (KEEP(201).EQ.1) THEN ! PANEL version
               I8TMP = int(id%OOC_NB_FILE_TYPE,8) *
     &               2_8 * int(KEEP(226),8)
             ELSE
               I8TMP = 2_8 * id%KEEP8(119)
             ENDIF
             I8TMP = I8TMP +  int(max(KEEP(12),0),8) *
     &               (I8TMP/100_8+1_8)
C            we want to avoid too large IO buffers.
C            12M corresponds to 100Mbytes given to buffers.
             I8TMP = min(I8TMP, 12000000_8)
             KEEP(100)=int(I8TMP)
           ENDIF
           IF (KEEP(201).EQ.1) THEN
C            Panel version. Force the use of a buffer. 
             IF ( KEEP(99) < 3 ) THEN
               KEEP(99) = KEEP(99) + 3
             ENDIF
           ENDIF
C          --------------------------
C          Reset KEEP(100) to 0 if no
C          buffer is used for OOC.
C          --------------------------
           IF (KEEP(99) .LT.3) KEEP(100)=0
           IF((dble(KEEP(100))*dble(KEEP(35))/dble(2)).GT.
     &          (dble(1999999999)))THEN
             IF (PROKG) THEN
               WRITE(MPG,*)id%MYID,': Warning: DIM_BUF_IO might be
     &  too big for Filesystem'
             ENDIF
           ENDIF
           ALLOCATE (id%OOC_INODE_SEQUENCE(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_INODE_SEQUENCE)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_TOTAL_NB_NODES(id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE
              NULLIFY(id%OOC_TOTAL_NB_NODES)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_SIZE_OF_BLOCK(KEEP(28),
     &          id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_SIZE_OF_BLOCK)
              GOTO 112
           ENDIF
           ALLOCATE (id%OOC_VADDR(KEEP(28),id%OOC_NB_FILE_TYPE),
     &          stat=IERR)
           IF ( IERR .GT. 0 ) THEN
              id%INFO(1) = -13
              id%INFO(2) = id%OOC_NB_FILE_TYPE*KEEP(28)
              NULLIFY(id%OOC_VADDR)
              GOTO 112
           ENDIF
        ENDIF
      ENDIF
 112  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF (id%INFO(1) < 0) THEN
C       LOAD_END must be done but not OOC_END_FACTO
        GOTO 513
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (KEEP(201) .GT. 0) THEN
           IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
             CALL SMUMPS_OOC_INIT_FACTO(id,MAXS)
          ELSE
             WRITE(*,*) "Internal error in SMUMPS_FAC_DRIVER"
             CALL MUMPS_ABORT()
           ENDIF
           IF(id%INFO(1).LT.0)THEN
              GOTO 111
           ENDIF
        ENDIF
C       First increment corresponds to the number of
C       floating-point operations for subtrees allocated
C       to the local processor.
        CALL SMUMPS_LOAD_UPDATE(0,.FALSE.,dble(id%COST_SUBTREES),
     &          id%KEEP(1),id%KEEP8(1))
        IF (id%INFO(1).LT.0) GOTO 111
      END IF
C     -----------------------
C     Manage main workarray S
C     -----------------------
      EARLYT3ROOTINS = KEEP(200) .EQ.0
     &        .OR.  ( KEEP(200) .LT. 0 .AND. KEEP(400) .EQ. 0 )
#if defined (LARGEMATRICES)
      IF ( id%MYID .ne. MASTER ) THEN
#endif
      IF (.NOT.WK_USER_PROVIDED) THEN
        IF ( EARLYT3ROOTINS ) THEN
C         Standard allocation strategy
          ALLOCATE (id%S(MAXS),stat=IERR)
          id%KEEP8(23) = MAXS
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(MAXS, id%INFO(2))
C           On some platforms (IBM for example), an
C           allocation failure returns a non-null pointer.
C           Therefore we nullify S
            NULLIFY(id%S)
            id%KEEP8(23)=0_8
          ENDIF
        ENDIF
      ELSE
       id%S => id%WK_USER(1:id%KEEP8(24))
       id%KEEP8(23) = 0_8
      ENDIF
#if defined (LARGEMATRICES)
      END IF
#endif
C
 111  CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 514
C     --------------------------
C     Initialization of modules
C     related to data management
C     --------------------------
      NB_ACTIVE_FRONTS_ESTIM = 3
      NB_THREADS = 1
!$    NB_THREADS = OMP_GET_MAX_THREADS()
C     
      NB_ACTIVE_FRONTS_ESTIM = 3*NB_THREADS
      IF (I_AM_SLAVE) THEN
C
        CALL MUMPS_FDM_INIT('A',NB_ACTIVE_FRONTS_ESTIM, id%INFO)
C
        IF ( (KEEP(486).EQ.2) 
     &     .OR. ((KEEP(489).NE.0).AND.(KEEP(400).GT.1))
     &    ) THEN
C         In case of LRSOLVE or CompressCB, 
C         initialize nb of handlers to nb of BLR
C         nodes estimated at analysis
          NB_FRONTS_F_ESTIM = KEEP(470)
        ELSE
          IF (KEEP(489).NE.0) THEN
C          Compress CB and no L0 OMP (or 1 thread under L0): 
C          NB_ACTIVE_FRONTS_ESTIM is too small, 
C          to limit nb of reallocations make it twice larger
           NB_FRONTS_F_ESTIM = 2*NB_ACTIVE_FRONTS_ESTIM
          ELSE
           NB_FRONTS_F_ESTIM = NB_ACTIVE_FRONTS_ESTIM
          ENDIF
        ENDIF
        CALL MUMPS_FDM_INIT('F',NB_FRONTS_F_ESTIM, id%INFO )
        IF (id%INFO(1) .LT. 0 ) GOTO 114
#if ! defined(NO_FDM_DESCBAND)
C         Storage of DESCBAND information
          CALL MUMPS_FDBD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
#if ! defined(NO_FDM_MAPROW)
C         Storage of MAPROW and ROOT2SON information
          CALL MUMPS_FMRD_INIT( NB_ACTIVE_FRONTS_ESTIM, id%INFO )
#endif
        CALL SMUMPS_BLR_INIT_MODULE( NB_FRONTS_F_ESTIM, id%INFO )
 114    CONTINUE
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                       id%COMM, id%MYID )
C     GOTO 500: one of the above module initializations failed
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C
C  Allocate space for matrix in arrowhead form
C  ===========================================
C
C  CASE 1 : Matrix is assembled
C  CASE 2 : Matrix is elemental
C
      IF ( KEEP(55) .eq. 0 ) THEN
C       ------------------------------------
C       Space has been allocated already for
C       the integer part during analysis
C       Only slaves need the arrowheads.
C       ------------------------------------
        IF ( I_AM_SLAVE .and. id%KEEP8(26) .ne. 0_8 ) THEN
          ALLOCATE( DBLARR( id%KEEP8(26) ),
     &              INTARR( id%KEEP8(27) ), stat = IERR )
        ELSE
          ALLOCATE( DBLARR( 1 ),
     &              INTARR( 1 ), stat =IERR )
        END IF
        IF ( IERR .NE. 0 ) THEN
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for DBLARR(',id%KEEP8(26),')+INTARR(',
     &      id%KEEP8(27),')'
          ENDIF
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR( max(id%KEEP8(26),1_8)+
     &                           max(id%KEEP8(27),1_8),
     &                            id%INFO(2) )
          GOTO 100
        END IF
      ELSE
C        --------------------------------
C        Allocate element variables lists
C        --------------------------------
         IF ( I_AM_SLAVE .and. id%KEEP8(27) .ne. 0_8 ) THEN
           ALLOCATE( INTARR( id%KEEP8(27) ), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             CALL MUMPS_SET_IERROR(id%KEEP8(27), id%INFO(2))
             GOTO 100
           END IF
         ELSE
C          INTARR also allocated of size 1 on non-working master
           ALLOCATE( INTARR(1),stat=allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             GOTO 100
           END IF
         END IF
C        ---------------------------------------
C        Allocate DBLARR to hold possibly scaled
C        copies of elemental matrices
C        On a working master (hybrid host) and
C        no scaling, avoid the copy and point
C        directly to user data instead.
C        ---------------------------------------
         IF ( I_AM_SLAVE ) THEN
           IF (      id%MYID_NODES .eq. MASTER
     &       .AND.   KEEP(46)   .eq. 1
     &       .AND.   KEEP(52)   .eq. 0 ) THEN
C            -------------------
C            Pointer association
C            -------------------
             DBLARR => id%A_ELT
           ELSE
C            ----------
C            Allocation
C            ----------
             ALLOCATE( DBLARR( max(id%KEEP8(26),1_8) ), stat=allocok )
             IF ( allocok .GT. 0 ) THEN
               id%INFO(1) = -13
               CALL MUMPS_SET_IERROR(max(id%KEEP8(26),1_8), id%INFO(2))
               NULLIFY(DBLARR)
               GOTO 100
             END IF
           END IF
         ELSE
           ALLOCATE( DBLARR(1), stat = allocok )
           IF ( allocok .GT. 0 ) THEN
             id%INFO(1) = -13
             id%INFO(2) = 1
             NULLIFY(DBLARR)
             GOTO 100
           END IF
         END IF
      END IF
 100  CONTINUE
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                     id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C     ------------------------------------------
C     Prepare stuff for the 2D block-cyclic root
C     ------------------------------------------
      IF ( KEEP(38).NE.0 ) THEN
          ALLOCATE(id%root%RG2L(id%N), stat=allocok)
          IF ( allocok .NE. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%N
          END IF
        IF ( id%INFO(1) .GE.0 ) THEN
          CALL SMUMPS_INIT_ROOT_FAC( id%N, id%MYID,
     &    id%root, id%FILS(1), id%KEEP(1) )
        ENDIF
      ENDIF
C     ----------------
C     Check for errors
C     ----------------
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                     id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C       -----------------------------------
C
C       DISTRIBUTION OF THE ORIGINAL MATRIX
C
C       -----------------------------------
C
C     TIMINGS: computed (and printed) on the host
C     Next line: global time for distrib(arrowheads,elts)
C     on the host. Synchronization has been performed.
      IF (id%MYID.EQ.MASTER) CALL MUMPS_SECDEB(TIME)
C     -------------------------------------------
C     S_PTR_ARG / MAXS_ARG will be used for id%S
C     argument to arrowhead/element distribution
C     routines: if id%S is not allocated, we pass
C     S_DUMMY_ARG instead, which is not accessed.
C     -------------------------------------------
      IF (EARLYT3ROOTINS) THEN
            S_PTR_ARG => id%S
            MAXS_ARG = MAXS
      ELSE
            S_PTR_ARG => S_DUMMY_ARG
            MAXS_ARG = 1
      ENDIF
C
      IF ( KEEP( 55 ) .eq. 0 ) THEN
C     {
C       ----------------------------
C       Original matrix is assembled
C       Arrowhead format to be used.
C       ----------------------------
C       KEEP8(26) and KEEP8(27) hold the number of entries for real/integer
C       for the matrix in arrowhead format. They have been set by the
C       analysis phase (SMUMPS_ANA_F and SMUMPS_ANA_G)
C
C       ------------------------------------------------------------------
C       Blocking is used for sending arrowhead records (I,J,VAL)
C              buffer(1) is used to store number of bytes already packed
C              buffer(2) number of records already packed
C       KEEP(39) : Number of records (blocking factor)
C       ------------------------------------------------------------------
C
C     ---------------------------------------------
C     In case of parallel root compute minimum
C     size of workspace to receive arrowheads
C     of root node. Will be used to check that
C     MAXS is large enough for arrowheads (case
C     of EARLYT3ROOTINS (KEEP(200)=0); if .NOT.
C     EARLYT3ROOTINS (KEEP(200)=1), root will
C     be assembled into id%S later and size of
C     id%S will be checked later)
C     ---------------------------------------------
      IF (EARLYT3ROOTINS .AND. KEEP(38).NE.0 .AND.
     &    KEEP(60) .EQ.0 .AND. I_AM_SLAVE) THEN
        LWK = int(MUMPS_NUMROC( id%root%ROOT_SIZE, id%root%MBLOCK,
     &             id%root%MYROW, 0, id%root%NPROW ),8)
        LWK = max( 1_8, LWK )
        LWK = LWK*
     &        int(MUMPS_NUMROC( id%root%ROOT_SIZE, id%root%NBLOCK,
     &        id%root%MYCOL, 0, id%root%NPCOL ),8)
        LWK = max( 1_8, LWK )
      ELSE
        LWK = 1_8
      ENDIF
C     MAXS must be at least 1, and in case of
C     parallel root, large enough to receive
C     arrowheads of root.
      IF (MAXS .LT. int(LWK,8)) THEN
           id%INFO(1) = -9
           CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
      IF ( id%INFO(1).LT.0 ) GOTO 500
C
C     =====================================================
      IF (KEEP(399).GE.2) THEN
C{    Multihtreaded algorithm taking into account all cases
C     =====================================================
C
C        Timing on master.
         IF (id%MYID.EQ.MASTER) THEN
           CALL MUMPS_SECDEB(TIME)
         END IF
C
         NZ_loc8   = 0_8
         NBRECORDS = KEEP(39)
         SIZESCAL = id%N
C        Set  NZ_loc8, A_loc_PTR, IRN_loc_PTR, JCN_loc_PTR
C        and update NBRECORDS
         IF (KEEP(54).EQ.0) THEN
C         centralized matrix
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
              NBRECORDS = int(id%KEEP8(28))
          ENDIF
          IF (id%MYID.EQ.MASTER) THEN
           NZ_loc8 = id%KEEP8(28)
           A_loc_PTR   => id%A
           IRN_loc_PTR => id%IRN
           JCN_loc_PTR => id%JCN
           IF (LSCAL) THEN
             COLSCA_PTR => id%COLSCA
             ROWSCA_PTR => id%ROWSCA
           ELSE
             COLSCA_PTR  => DUMMYSCA
             ROWSCA_PTR  => DUMMYSCA
             SIZESCAL  = 1
           ENDIF
          ELSE
           A_loc_PTR   => DUMMYA_loc
           IRN_loc_PTR => DUMMYIRN_loc
           JCN_loc_PTR => DUMMYJCN_loc
           COLSCA_PTR  => DUMMYSCA
           ROWSCA_PTR  => DUMMYSCA
           SIZESCAL  = 1
          ENDIF
         ELSE
C         distributed matrix
C         ---------------------------------------------------
C         In order to have possibly IRN_loc/JCN_loc/A_loc
C         of size 0, avoid to pass them inside REDISTRIBUTION
C         and pass id instead
C         NZ_locMAX8 gives as a maximum buffer size (send/recv) used
C          an upper bound to limit buffers on small matrices
C         ---------------------------------------------------
          CALL MPI_ALLREDUCE(id%KEEP8(29), NZ_locMAX8, 1, 
     &              MPI_INTEGER8, MPI_MAX, id%COMM, IERR)
          IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
              NBRECORDS = int(NZ_locMAX8)
          ENDIF
          NZ_loc8 = id%KEEP8(29)
          LSCAL   = (KEEP(52).EQ.7).OR.(KEEP(52).EQ.8)
C         available on all MPI
          IF (LSCAL) THEN
            COLSCA_PTR => id%COLSCA
            ROWSCA_PTR => id%ROWSCA
          ELSE
           COLSCA_PTR  => DUMMYSCA
           ROWSCA_PTR  => DUMMYSCA
           SIZESCAL  = 1
          ENDIF
         ENDIF
#if ! defined(PCPRET)
          IF (id%KEEP(72).EQ.1) THEN
           NBRECORDS = max(3,NBRECORDS/10)
          ENDIF
         CALL SMUMPS_FAC_DIST_ARROWHEADS_OMP ( id%N,
     &     NZ_loc8, 
C          replace id by:
     &        A_loc_PTR(1), IRN_loc_PTR(1), JCN_loc_PTR(1),
     &        SIZESCAL, LSCAL, COLSCA_PTR(1), ROWSCA_PTR(1), 
     &     DBLARR(1), id%KEEP8(26), INTARR(1),
     &     id%KEEP8(27),
     &     id%PTR8ARR(1), id%NINCOLARR(1),
     &     id%NINROWARR(1), id%PTRDEBARR(1),
     &     id%FILS(1),
     &     KEEP(1), id%KEEP8(1), id%MYID, id%COMM, NBRECORDS,
     &     S_PTR_ARG(1), MAXS_ARG, id%root, id%PROCNODE_STEPS(1),
     &     id%NPROCS, id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &     id%ICNTL(1), id%INFO(1), NSEND8, NLOCAL8,
     &     id%ISTEP_TO_INIV2(1),
     &     id%CANDIDATES(1,1) 
     &     )
         CALL MPI_BARRIER(id%COMM, IERR)
#else 
#endif
C     
         IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C          ------------------------------
C          The host is not working -> had
C          no data from initial matrix
C          ------------------------------
           NSEND8  = 0_8
           NLOCAL8 = 0_8
         END IF
C        -------------------------
C        Check for possible errors
C        -------------------------
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
         IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
         CALL MPI_REDUCE( NSEND8, NSEND_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
         CALL MPI_REDUCE( NLOCAL8, NLOCAL_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
         IF ( PROKG ) THEN
           WRITE(MPG,125) NLOCAL_TOT8, NSEND_TOT8
         END IF
C}        
      ELSE
C{   =======================================================
      IF ( KEEP(54) .eq. 0 ) THEN
C     {
C       ================================================
C       FIRST CASE : MATRIX IS NOT INITIALLY DISTRIBUTED
C       ================================================
C       A small integer workspace is needed to
C       send the arrowheads.
        IF ( id%MYID .eq. MASTER ) THEN
#if defined(LARGEMATRICES)
          ALLOCATE (WK(LWK),stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1) = -13
            CALL MUMPS_SET_IERROR(LWK, id%INFO(2))
            write(6,*) ' PB1 ALLOC LARGEMAT'
          ENDIF
#endif
        ENDIF
        CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
        IF ( id%INFO(1).LT.0 ) GOTO 500
        IF ( id%MYID .eq. MASTER ) THEN
C
C         --------------------------------
C         MASTER sends arowheads using the
C         global communicator with ranks
C         also in global communicator
C         IWK is used as temporary
C         workspace of size N.
C         --------------------------------
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
#if defined(LARGEMATRICES)
          CALL SMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28), id%A(1),
     &      id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &      LSCAL, id%COLSCA(1), id%ROWSCA(1),
     &      id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &      NBRECORDS,
     &      id%COMM, id%root, KEEP,id%KEEP8,
     &      id%FILS(1),
     &
     &      INTARR(1), id%KEEP8(27), DBLARR(1), id%KEEP8(26),
     &      id%PTR8ARR(1), id%NINCOLARR(1),
     &      id%NINROWARR(1), id%PTRDEBARR(1),
     &      id%FRERE_STEPS(1), id%STEP(1), WK(1), LWK,
     &      id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &      id%CANDIDATES(1,1), id%ICNTL(1), id%INFO(1) )
C         write(6,*) '!!! A,IRN,JCN are freed during factorization '
          DEALLOCATE (id%A)
          NULLIFY(id%A)
          DEALLOCATE (id%IRN)
          NULLIFY (id%IRN)
          DEALLOCATE (id%JCN)
          NULLIFY (id%JCN)
          IF (.NOT.WK_USER_PROVIDED) THEN
            IF (EARLYT3ROOTINS) THEN
              ALLOCATE (id%S(MAXS),stat=IERR)
              id%KEEP8(23) = MAXS
              IF ( IERR .GT. 0 ) THEN
                id%INFO(1) = -13
                id%INFO(2) = MAXS
                NULLIFY(id%S)
                id%KEEP8(23)=0_8
                write(6,*) ' PB2 ALLOC LARGEMAT',MAXS
                CALL MUMPS_ABORT()
               ENDIF
             ENDIF
            ENDIF
          ELSE
            id%S => id%WK_USER(1:id%KEEP8(24))
          ENDIF
          IF (EARLYT3ROOTINS) THEN
            id%S(MAXS-LWK+1_8:MAXS) = WK(1_8:LWK)
          ENDIF
          DEALLOCATE (WK)
#else
          CALL SMUMPS_FACTO_SEND_ARROWHEADS(id%N, id%KEEP8(28),
     &    id%A(1), id%IRN(1), id%JCN(1), id%SYM_PERM(1),
     &    LSCAL, id%COLSCA(1), id%ROWSCA(1),   
     &    id%MYID, id%NSLAVES, id%PROCNODE_STEPS(1),
     &    NBRECORDS,
     &    id%COMM, id%root, KEEP(1),id%KEEP8(1),
     &    id%FILS(1),
     &
     &    INTARR(1), id%KEEP8(27), DBLARR(1), id%KEEP8(26),
     &    id%PTR8ARR(1), id%NINCOLARR(1),
     &    id%NINROWARR(1), id%PTRDEBARR(1),
     &    id%FRERE_STEPS(1), id%STEP(1), S_PTR_ARG(1), MAXS_ARG,
     &    id%ISTEP_TO_INIV2(1), id%I_AM_CAND(1),
     &    id%CANDIDATES(1,1), id%ICNTL(1), id%INFO(1) )
#endif
        ELSE
          NBRECORDS = KEEP(39)
          IF (id%KEEP8(28) .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(id%KEEP8(28))
          ENDIF
          CALL SMUMPS_FACTO_RECV_ARROWHD2( id%N,
     &       DBLARR(1), id%KEEP8(26),
     &       INTARR(1), id%KEEP8(27),
     &       id%PTR8ARR(1), id%NINCOLARR(1),
     &       id%NINROWARR(1), id%PTRDEBARR(1),
     &       KEEP( 1 ), id%KEEP8(1), id%FILS(1), id%MYID, id%COMM,
     &       NBRECORDS,
     &
     &       S_PTR_ARG(1), MAXS_ARG,
     &       id%root,
     &       id%PROCNODE_STEPS(1), id%NSLAVES,
     &       id%SYM_PERM(1), id%FRERE_STEPS(1), id%STEP(1),
     &       id%ICNTL(1), id%INFO(1) )
        ENDIF
C     }
      ELSE
C     {
C     =============================================
C     SECOND CASE : MATRIX IS INITIALLY DISTRIBUTED
C     =============================================
C     Timing on master.
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECDEB(TIME)
      END IF
      IF ( I_AM_SLAVE ) THEN
C     {
C       ---------------------------------------------------
C       In order to have possibly IRN_loc/JCN_loc/A_loc
C       of size 0, avoid to pass them inside REDISTRIBUTION
C       and pass id instead
C       NZ_locMAX8 gives as a maximum buffer size (send/recv) used
C        an upper bound to limit buffers on small matrices
C       ---------------------------------------------------
        CALL MPI_ALLREDUCE(id%KEEP8(29), NZ_locMAX8, 1, MPI_INTEGER8,
     &                     MPI_MAX, id%COMM_NODES, IERR)
        NBRECORDS = KEEP(39)
        IF (NZ_locMAX8 .LT. int(NBRECORDS,8)) THEN
            NBRECORDS = int(NZ_locMAX8)
        ENDIF
        CALL SMUMPS_REDISTRIBUTION( id%N,
     &  id%KEEP8(29),
     &  id,
     &  DBLARR(1), id%KEEP8(26), INTARR(1),
     &  id%KEEP8(27),
     &  id%PTR8ARR(1), id%NINCOLARR(1),
     &  id%NINROWARR(1), id%PTRDEBARR(1),
     &  id%FILS(1),
     &  KEEP(1), id%KEEP8(1), id%MYID_NODES,
     &  id%COMM_NODES, NBRECORDS,
     &  S_PTR_ARG(1), MAXS_ARG, id%root, id%PROCNODE_STEPS(1),
     &  id%NSLAVES, id%SYM_PERM(1), id%STEP(1),
     &  id%ICNTL(1), id%INFO(1), NSEND8, NLOCAL8,
     &  id%ISTEP_TO_INIV2(1),
     &  id%CANDIDATES(1,1) )
        IF ( ( KEEP(52).EQ.7 ).OR. (KEEP(52).EQ.8) ) THEN
C         -------------------------------------------------
C         In that case, scaling arrays have been allocated
C         on all processors. They were useful for matrix
C         distribution. But we now really only need them
C         on the host. In case of distributed solution, we
C         will have to broadcast either ROWSCA or COLSCA
C         (depending on MTYPE) but this is done later.
C
C         In other words, on exit from the factorization,
C         we want to have scaling arrays available only
C         on the host.
C         -------------------------------------------------
          IF ( id%MYID > 0 ) THEN
            IF (associated(id%ROWSCA)) THEN
              DEALLOCATE(id%ROWSCA)
              NULLIFY(id%ROWSCA)
            ENDIF
            IF (associated(id%COLSCA)) THEN
              DEALLOCATE(id%COLSCA)
              NULLIFY(id%COLSCA)
            ENDIF
          ENDIF
        ENDIF
#if defined(LARGEMATRICES)
C      deallocate id%IRN_loc, id%JCN(loc) to free extra space
C      Note that in this case IRN_loc cannot be used
C      anymore during the solve phase for IR and Error analysis.
        IF (associated(id%IRN_loc)) THEN
           DEALLOCATE(id%IRN_loc)
           NULLIFY(id%IRN_loc)
        ENDIF
        IF (associated(id%JCN_loc)) THEN
           DEALLOCATE(id%JCN_loc)
           NULLIFY(id%JCN_loc)
        ENDIF
        IF (associated(id%A_loc)) THEN
           DEALLOCATE(id%A_loc)
           NULLIFY(id%A_loc)
        ENDIF
        write(6,*) ' Warning :', 
     &        ' id%A_loc, IRN_loc, JCN_loc deallocated !!! '
#endif
        IF (PROK) THEN
          WRITE(MP,120) NLOCAL8, NSEND8
        END IF
C     }
      END IF ! I_AM_SLAVE
C     -------------------------
C     Check for possible errors
C     -------------------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
      IF ( KEEP(46) .eq. 0 .AND. id%MYID.eq.MASTER ) THEN
C       ------------------------------
C       The host is not working -> had
C       no data from initial matrix
C       ------------------------------
        NSEND8  = 0_8
        NLOCAL8 = 0_8
      END IF
C     --------------------------
C     Put into some info/infog ?
C     --------------------------
      CALL MPI_REDUCE( NSEND8, NSEND_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( NLOCAL8, NLOCAL_TOT8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
        WRITE(MPG,125) NLOCAL_TOT8, NSEND_TOT8
      END IF
C     }
      ENDIF ! distributed matrix
C     }
      ENDIF !  "old" multithreaded algorithm 
C     }
      ELSE
C     {
C       -------------------
C       Matrix is elemental,
C       provided on the
C       master only
C       -------------------
        IF ( id%MYID.eq.MASTER)
     &   CALL SMUMPS_MAXELT_SIZE( id%ELTPTR(1),
     &                        id%NELT,
     &                        MAXELT_SIZE )
C
C         Perform the distribution of the elements.
C         A this point,
C           PTRAIW/PTRARW have been computed.
C           INTARR/DBLARR have been allocated
C           ELTPROC gives the mapping of elements
C
        CALL SMUMPS_ELT_DISTRIB( id%N, id%NELT, id%KEEP8(30),
     &     id%COMM, id%MYID,
     &     id%NSLAVES, id%PTRAR(1),
     &     id%PTRAR(id%NELT+2),
     &     INTARR(1), DBLARR(1), id%KEEP8(27), id%KEEP8(26),
     &     id%KEEP(1), id%KEEP8(1), MAXELT_SIZE,
     &     id%FRTPTR(1), id%FRTELT(1),
     &     S_PTR_ARG(1), MAXS_ARG, id%FILS(1),
     &     id, id%root )
C       ----------------
C       Broadcast errors
C       ----------------
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C     }
      END IF ! Element entry
C     ------------------------
C     Time the redistribution:
C     ------------------------
      IF ( id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(93) = real(TIME)
        IF (PROKG) WRITE(MPG,160) id%DKEEP(93)
      END IF
C     -------------------------------------
C     Small memory optimizaiton: we can now
C     free RG2L on the non working host,
C     -------------------------------------
      IF (id%KEEP(38) .NE. 0 .AND. .NOT. I_AM_SLAVE) THEN
          IF (associated(id%root%RG2L)) THEN
            DEALLOCATE(id%root%RG2L)
            NULLIFY   (id%root%RG2L)
          ENDIF
      ENDIF
      IF ( KEEP(400) .GT. 0 ) THEN
C       L0-OMP was active at analysis and
C       thus will be active at factorization
C       We check the number of threads.
        NOMP=1
!$      NOMP = omp_get_max_threads()
          IF ( NOMP .NE. KEEP(400) ) THEN
            id%INFO(1)=-58
            id%INFO(2)=KEEP(400)
            IF (LPOK) WRITE(LP,'(A,A,I5,A,I5)')
     &" FAILURE DETECTED IN FACTORIZATION: #threads for multithreaded",
     &" tree parallelism changed from",KEEP(400)," at analysis to", 
     &  NOMP
          ENDIF
C       error check
        CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
      ENDIF
C
C     TIMINGS:
C     Next line: elapsed time for factorization
      IF (id%MYID.EQ.MASTER)  CALL MUMPS_SECDEB(TIME)
C
C  Allocate buffers on the workers
C  ===============================
C
      IF ( I_AM_SLAVE )  THEN
        CALL SMUMPS_BUF_INI_MYID(id%MYID_NODES)
C
C  Some buffers are required to pack/unpack data and for
C  receiving MPI messages.
C  For packing/unpacking : the buffer must be large
C  enough to send several messages while receives might not
C  be posted yet.
C  It is assumed that the size of an integer is held in KEEP(34)
C  while the size of a complex is held in KEEP(35).
C  BUFR and LBUFR are declared of type integer, since byte is not
C  a standard datatype.
C  We now use KEEP(43) or KEEP(379) and KEEP(44) or KEEP(380)
C  as estimated at analysis to allocate appropriate buffer sizes
C
C  Reception buffer
C  ----------------
        IF (KEEP(486).NE.0) THEN
          SMUMPS_LBUFR_BYTES8 = int(KEEP( 380 ),8) * int(KEEP( 35 ),8)
        ELSE
          SMUMPS_LBUFR_BYTES8 = int(KEEP( 44 ),8) * int(KEEP( 35 ),8)
        ENDIF
C       ---------------------------------------
C       Ensure a reasonable minimal buffer size
C       ---------------------------------------
        SMUMPS_LBUFR_BYTES8 = max( SMUMPS_LBUFR_BYTES8,
     &                      100000_8 )
        IF ((KEEP(50).NE.0).AND.(KEEP(489).GT.0).AND.
     &       (id%NSLAVES.GE.2)) THEN
C        ----------------------------------------------------------
C        Ensure large enough receive buffer in case of BLR with
C        CompressCB for symmetric matrices. 
C        -----------------------------------------------------------
          RATIOK465 = real(KEEP465COPY)/real(1000)
          SMUMPS_LBUFR_BYTES8 = max (SMUMPS_LBUFR_BYTES8, 
     &     int( 
     &       RATIOK465*
     &       real(
     &        int(KEEP(2)+1,8)*int(KEEP(142),8)*int(KEEP(35),8) 
     &           ) 
     &       ,8)
     &                              )
        ENDIF
C
C  If there is pivoting, size of the message might still increase.
C  We use a relaxation (so called PERLU) to increase the estimate.
C
C  Note: PERLU is a global estimate for pivoting. 
C  It may happen that one large contribution block size is increased 
C  by more than that.
C  This is why we use an extra factor 2 relaxation coefficient for 
C  the relaxation of
C  the reception buffer in the case where pivoting is allowed.
C  A more dynamic strategy could be applied: if message to
C  be received is larger than expected, reallocate a larger
C  buffer. (But this won't work with IRECV.)
C  Finally, one may want (as we are currently doing it for 
C  most messages)
C  to cut large messages into a series of smaller ones.
C
        IF (KEEP(48).EQ.5) THEN
          MIN_PERLU = 2
        ELSE
          MIN_PERLU = 0
        ENDIF
C
        SMUMPS_LBUFR_BYTES8 = SMUMPS_LBUFR_BYTES8
     &        + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &        real(SMUMPS_LBUFR_BYTES8)/100E0, 8)
        SMUMPS_LBUFR_BYTES8 = min(SMUMPS_LBUFR_BYTES8,
     &                            int(huge(I4)-100,8))
        SMUMPS_LBUFR_BYTES  = int( SMUMPS_LBUFR_BYTES8 )
C       SMUMPS_LBUFR is the size of the buffer as a number of integers,
C       we round SMUMPS_LBUFR (size in #integers) above to have at least
C       SMUMPS_LBUFR_BYTES available in the buffer.
        SMUMPS_LBUFR = (SMUMPS_LBUFR_BYTES+KEEP(34)-1)/KEEP(34)
C       Finally, make SMUMPS_LBUFR_BYTES a multiple of KEEP(34) by setting
C       SMUMPS_LBUFR_BYTES to the number of bytes that will be allocated
        SMUMPS_LBUFR_BYTES = SMUMPS_LBUFR*KEEP(34)
        IF (KEEP(48)==5) THEN
C          Since the buffer is going to be allocated, use
C          it as the constraint for memory/granularity
C          in hybrid scheduler
C
           id%KEEP8(21) = id%KEEP8(22) +
     &        int( real(max(PERLU,MIN_PERLU))*
     &        real(id%KEEP8(22))/100E0,8)
        ENDIF
C
C  Now estimate the size for the buffer for asynchronous
C  sends of contribution blocks (so called CB). We want to be able to send at
C  least KEEP(213)/100 (two in general) messages at the
C  same time.
C
C   Send buffer
C   -----------
        IF (KEEP(486).NE.0) THEN
         SMUMPS_LBUF8 = int( real(KEEP(213)) / 100.0E0 *
     &                      real(KEEP(379)) * real(KEEP(35)), 8  )
        ELSE
         SMUMPS_LBUF8 = int( real(KEEP(213)) / 100.0E0 *
     &                      real(KEEP(43)) * real(KEEP(35)), 8  )
        ENDIF
        SMUMPS_LBUF8 = max( SMUMPS_LBUF8, 100000_8 )
        SMUMPS_LBUF8 = SMUMPS_LBUF8
     &                 + int( 2.0E0 * real(max(PERLU,MIN_PERLU))*
     &                   real(SMUMPS_LBUF8)/100E0, 8)
C       Make SMUMPS_LBUF8 small enough to be stored in a standard integer
        SMUMPS_LBUF8 = min(SMUMPS_LBUF8, int(huge(I4)-100,8))
C
C       No reason to have send buffer smaller than receive buffer.
C       This should never occur with the formulas above but just
C       in case:
        SMUMPS_LBUF8 = max(SMUMPS_LBUF8, SMUMPS_LBUFR_BYTES8+3*KEEP(34))
        SMUMPS_LBUF  = int(SMUMPS_LBUF8)
        IF(id%KEEP(48).EQ.4)THEN
           SMUMPS_LBUFR_BYTES=SMUMPS_LBUFR_BYTES*5
           SMUMPS_LBUF=SMUMPS_LBUF*5
        ENDIF
C
C  Estimate size of buffer for small messages 
C  Each node can send ( NSLAVES - 1 ) messages to (NSLAVES-1) nodes
C
C  KEEP(56) is the number of nodes of level II.
C  Messages will be sent for the symmetric case
C  for synchronisation issues.
C
C  We take an upperbound
C
        SMUMPS_LBUF_INT = ( KEEP(56) + id%NSLAVES * id%NSLAVES ) * 5
     &               * KEEP(34)
        IF ( KEEP( 38 ) .NE. 0 ) THEN
C
C
          KKKK = MUMPS_PROCNODE( id%PROCNODE_STEPS(id%STEP(KEEP(38))),
     &                           id%KEEP(199) )
          IF ( KKKK .EQ. id%MYID_NODES ) THEN
             SMUMPS_LBUF_INT = SMUMPS_LBUF_INT + 4 * KEEP(34) *
     &         ( id%NSLAVES + id%NE_STEPS(id%STEP(KEEP(38)))
     &      + min(KEEP(56), id%NE_STEPS(id%STEP(KEEP(38)))) * id%NSLAVES
     &         )
          END IF
        END IF
C       At this point, SMUMPS_LBUFR_BYTES, SMUMPS_LBUF
C       and SMUMPS_LBUF_INT have been computed (all
C       are in numbers of bytes).
        IF ( PROK ) THEN
          WRITE( MP, 9999 ) SMUMPS_LBUFR_BYTES,
     &                      SMUMPS_LBUF, SMUMPS_LBUF_INT
        ELSE 
         IF (PROKG) THEN
          WRITE( MPG, 9999 ) SMUMPS_LBUFR_BYTES,
     &                      SMUMPS_LBUF, SMUMPS_LBUF_INT
         ENDIF
        END IF
 9999   FORMAT( /,' Allocated buffers',/,' ------------------',/,
     &  ' Size of reception buffer in bytes ...... = ', I12,
     &  /,
     &  ' Size of async. emission buffer (bytes).. = ', I12,/,
     &  ' Small emission buffer (bytes) .......... = ', I12)
C       --------------------------
C       Allocate small send buffer
C       required for SMUMPS_FAC_B
C       --------------------------
        CALL SMUMPS_BUF_ALLOC_SMALL_BUF( SMUMPS_LBUF_INT, IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)= -13
C         convert to size in integer  id%INFO(2)= SMUMPS_LBUF_INT
          id%INFO(2)= (SMUMPS_LBUF_INT+KEEP(34)-1)/KEEP(34)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &     ':Allocation error in SMUMPS_BUF_ALLOC_SMALL_BUF'
     &     ,id%INFO(2)
          ENDIF
          GO TO 110
        END IF
C
C       --------------------------------------
C       Allocate reception buffer on all procs
C       This is done now.
C       --------------------------------------
        ALLOCATE( BUFR( SMUMPS_LBUFR ),stat=IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1) = -13
          id%INFO(2) = SMUMPS_LBUFR
          IF (LPOK) THEN
            WRITE(LP,*)
     &     ': Allocation error for BUFR(', SMUMPS_LBUFR,
     &     ') on MPI process',id%MYID
          ENDIF
          GO TO 110
        END IF
C       -----------------------------------------
C       Estimate MAXIS. IS will be allocated in
C       SMUMPS_FAC_B. It will contain factors and
C       contribution blocks integer information
C       -----------------------------------------
C       Relax integer workspace based on PERLU
        PERLU          = KEEP( 12 )
        IF (KEEP(201).GT.0) THEN
C         OOC panel or non panel (note that
C         KEEP(15)=KEEP(225) if non panel)
          MAXIS_ESTIM   = KEEP(225)
        ELSE
C         In-core or reals for factors not stored
          MAXIS_ESTIM   = KEEP(15)
        ENDIF
        MAXIS = max( 1, int( min(  int(huge(MAXIS),8),
     &          int(MAXIS_ESTIM,8) + 3_8 * max(int(PERLU,8),10_8) *
     &          ( int(MAXIS_ESTIM,8) / 100_8 + 1_8 )
     &                          )   ! min
     &                     )  ! int
     &             ) !max
C       ----------------------------
C       Allocate PTLUST_S and PTRFAC
C       They will be used to access
C       factors in the solve phase.
C       They are also needed for
C       SMUMPS_FAC_L0_OMP.
C       ----------------------------
        ALLOCATE( id%PTLUST_S( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTLUST_S(', id%KEEP(28),')'
          ENDIF
          NULLIFY(id%PTLUST_S)
          GOTO 110
        END IF
        ALLOCATE( id%PTRFAC( id%KEEP(28) ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=id%KEEP(28)
          NULLIFY(id%PTRFAC)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for id%PTRFAC(', id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C       -----------------------------
C       Reserve temporary workspace :
C       IPOOL, PTRWB, ITLOC, PTRIST
C       PTRWB will be subdivided again
C       in routine SMUMPS_FAC_B
C       -----------------------------
        PTRIST = 1
        PTRWB  = PTRIST + id%KEEP(28)
        ITLOC  = PTRWB  + 2 * id%KEEP(28)
C Fwd in facto: ITLOC of size id%N + id%KEEP(253)
        IPOOL  = ITLOC  + id%N + id%KEEP(253)
C
C       --------------------------------
C       NA(1) is an upperbound for LPOOL
C       --------------------------------
C       Structure of the pool:
C     ____________________________________________________
C    | Subtrees   |         | Top nodes           | 1 2 3 |
C     ----------------------------------------------------
        LPOOL = MUMPS_GET_POOL_LENGTH(id%NA(1), id%KEEP(1),id%KEEP8(1))
        ALLOCATE( IWK(  IPOOL + LPOOL - 1 ), stat = IERR )
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=IPOOL + LPOOL - 1
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWK(',IPOOL+LPOOL-1,')'
          ENDIF
          GOTO 110
        END IF
        ALLOCATE(IWK8( 2 * id%KEEP(28)), stat = IERR)
        IF ( IERR .NE. 0 ) THEN
          id%INFO(1)=-13
          id%INFO(2)=2 * id%KEEP(28)
          IF (LPOK) THEN
            WRITE(LP,*) id%MYID,
     &      ': Allocation error for IWKB(', 2*id%KEEP(28),')'
          ENDIF
          GOTO 110
        END IF
C
C  Return to SPMD
C
      ENDIF
C
 110  CONTINUE
C     ----------------
C     Broadcast errors
C     ----------------
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C
      IF ( I_AM_SLAVE )  THEN
C       Store size of receive buffers in SMUMPS_LBUF module
        CALL SMUMPS_BUF_DIST_IRECV_SIZE( SMUMPS_LBUFR_BYTES )
        IF (PROK) THEN
          WRITE( MP, 170 ) MAXS, MAXIS, MAXS_BASE8, KEEP(15),
     &    id%KEEP8(26), id%KEEP8(27), id%KEEP8(11), KEEP(26), KEEP(27)
        ENDIF
      END IF
C     ===============================================================
C     Before calling the main driver, SMUMPS_FAC_B,
C     some statistics should be initialized to 0,
C     even on the host node because they will be
C     used in REDUCE operations afterwards.
C     --------------------------------------------
C     Size of factors written. It will be set to POSFAC in
C     IC, otherwise we accumulate written factors in it.
      id%KEEP8(31)= 0_8
C     Size of factors under L0 will be returned
C     in id%KEEP8(64), not included in KEEP8(31))
C     Number of entries in factors
      id%KEEP8(10) = 0_8
C     KEEP8(8) will hold the volume of extra copies due to
C              in-place stacking in fac_mem_stack.F
      id%KEEP8(8)=0_8
      id%INFO(9:14)=0
      RINFO(2:3)=ZERO
      IF ( I_AM_SLAVE ) THEN
C       ------------------------------------
C       Call effective factorization routine
C       ------------------------------------
        IF ( KEEP(55) .eq. 0 ) THEN
          LDPTRAR = 1 ! PTRAR no longer used (of size 2)
        ELSE
          LDPTRAR = id%NELT + 1
        END IF
        IF ( id%KEEP(55) .NE. 0 ) THEN
          NELT_arg = id%NELT
        ELSE
C         ------------------------------
C         Use size 1 to avoid complaints
C         when using check bound options
C         ------------------------------
          NELT_arg = 1
        END IF
      ENDIF
      IF (I_AM_SLAVE) THEN
        IF (associated(id%L0_OMP_MAPPING))
     &    DEALLOCATE(id%L0_OMP_MAPPING)
        IF (KEEP(400) .GT. 0) THEN
          id%LL0_OMP_MAPPING = KEEP(28)
        ELSE
          id%LL0_OMP_MAPPING = 1
        ENDIF
        ALLOCATE(id%L0_OMP_MAPPING(id%LL0_OMP_MAPPING), stat=allocok)
        IF ( allocok > 0) THEN
            write(*,*) "Problem allocating L0_OMP_MAPPING",
     &              IERR, KEEP(28)
            GOTO 115
        ENDIF
        IF (associated(id%L0_OMP_FACTORS)) THEN
          CALL SMUMPS_FREE_L0_OMP_FACTORS(id%L0_OMP_FACTORS)
        ENDIF
        IF (KEEP(400) .GT. 0) THEN
          id%LL0_OMP_FACTORS = KEEP(400)
        ELSE
          id%LL0_OMP_FACTORS = 1
        ENDIF
        ALLOCATE(id%L0_OMP_FACTORS(id%LL0_OMP_FACTORS),stat = allocok)
        IF (allocok > 0) THEN
           id%INFO(1)=-7
           id%INFO(2)=NB_THREADS
           GOTO 111
        ENDIF
        CALL SMUMPS_INIT_L0_OMP_FACTORS(id%L0_OMP_FACTORS)
      ENDIF
 115  CONTINUE
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO( 1 ) .LT. 0 ) GOTO 500
C     Compute DKEEP(17) 
      AVG_FLOPS    =  RINFOG(1)/(real(id%NSLAVES))
      id%DKEEP(17) = max ( id%DKEEP(18), AVG_FLOPS/real(50) )
     &               
      IF (PROK.AND.id%MYID.EQ.MASTER) THEN
        IF (id%NSLAVES.LE.1) THEN
         WRITE(MP,'(/A,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops (RINFOG(1))                         = ', 
     &   RINFOG(1)
        ELSE
         WRITE(MP,'(/A,A,1PD10.3,A,1PD10.3)')
     &' Start factorization with total', 
     &' estimated flops RINFOG(1) / Average per MPI proc    = ', 
     &      RINFOG(1), ' / ', AVG_FLOPS
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
C       IS/S pointers passed to SMUMPS_FAC_B with
C       implicit interface through intermediate
C       structure S_IS_POINTERS. IS will be allocated
C       during SMUMPS_FAC_B.
C       In case of L0OMP, id%IS and id%S are allocated during
C       SMUMPS_FAC_B, and only after L0OMP nodes are processed,
C       in order to limit the global memory peak. 
        S_IS_POINTERS%IW => id%IS; NULLIFY(id%IS)
        S_IS_POINTERS%A  => id%S ; NULLIFY(id%S)
        CALL SMUMPS_FAC_B(id%N,S_IS_POINTERS,MAXS,MAXIS,id%SYM_PERM(1),
     &  id%NA(1),id%LNA,id%NE_STEPS(1),id%ND_STEPS(1), id%FILS(1),
     &  id%STEP(1),id%FRERE_STEPS(1),id%DAD_STEPS(1),id%CANDIDATES(1,1),
     &  id%ISTEP_TO_INIV2(1),id%TAB_POS_IN_PERE(1,1), id%PTRAR(1),
     &LDPTRAR,id%PTR8ARR(1),id%NINCOLARR(1),id%NINROWARR(1),id%PTRDEBARR
     &  (1), IWK(PTRIST),id%PTLUST_S(1),id%PTRFAC(1),IWK(PTRWB),IWK8,
     &  IWK(ITLOC),RHS_MUMPS(1),IWK(IPOOL),LPOOL,CNTL1,ICNTL(1),
     &  id%INFO(1), RINFO(1),KEEP(1),id%KEEP8(1),id%PROCNODE_STEPS(1),
     &  id%NSLAVES,id%COMM_NODES,id%MYID,id%MYID_NODES,BUFR,SMUMPS_LBUFR
     &  , SMUMPS_LBUFR_BYTES, SMUMPS_LBUF, INTARR(1), DBLARR(1),
     &  id%root, NELT_arg, id%FRTPTR(1), id%FRTELT(1),id%COMM_LOAD,
     &  id%ASS_IRECV, SEUIL, SEUIL_LDLT_NIV2, id%MEM_DIST(0),
     &  id%DKEEP(1), PIVNUL_LIST_STRUCT, id%LRGROUPS(1)
     &  ,id%IPOOL_B_L0_OMP(1),id%LPOOL_B_L0_OMP,
     &  id%IPOOL_A_L0_OMP(1),id%LPOOL_A_L0_OMP,id%L_VIRT_L0_OMP,
     &  id%VIRT_L0_OMP(1), id%VIRT_L0_OMP_MAPPING(1),id%L_PHYS_L0_OMP,
     &  id%PHYS_L0_OMP(1), id%PERM_L0_OMP(1), id%PTR_LEAFS_L0_OMP(1),
     &  id%L0_OMP_MAPPING(1),id%LL0_OMP_MAPPING,
     &  id%THREAD_LA, id%L0_OMP_FACTORS(1), id%LL0_OMP_FACTORS,
     &  id%I4_L0_OMP(1,1), size(id%I4_L0_OMP,1), size(id%I4_L0_OMP,2),
     &  id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1), size(id%I8_L0_OMP,2)
     &     )
        id%IS => S_IS_POINTERS%IW; NULLIFY(S_IS_POINTERS%IW)
        id%S  => S_IS_POINTERS%A ; NULLIFY(S_IS_POINTERS%A)
C
C       ------------------------------
C       Deallocate temporary workspace
C       ------------------------------
        DEALLOCATE( IWK  )
        DEALLOCATE( IWK8 )
      ENDIF
C     ---------------------------------
C     Free some workspace corresponding
C     to the original matrix in
C     arrowhead or elemental format.
C                  -----
C     Note : DBLARR may be a pointer
C     in case of element-entry.
C     ---------------------------------
      IF (allocated( INTARR )) DEALLOCATE( INTARR )
      IF ( KEEP(55) .eq. 0 ) THEN
C       ----------------
C       Assembled matrix
C       ----------------
        IF (associated( DBLARR )) THEN
          DEALLOCATE( DBLARR )
          NULLIFY   ( DBLARR )
        ENDIF
C
      ELSE
C       ----------------
C       Elemental matrix
C       ----------------
C       For the master from an hybrid host
C       execution without scaling, then scalar
C       values were not copied !
        IF (      id%MYID_NODES .eq. MASTER
     &      .AND.   KEEP(46)   .eq. 1
     &      .AND.   KEEP(52)   .eq. 0 ) THEN
            NULLIFY( DBLARR )
        ELSE
            IF (associated( DBLARR )) THEN
              DEALLOCATE( DBLARR )
              NULLIFY   ( DBLARR )
            ENDIF
        END IF
      END IF
C     We also free RG2L now
      IF ( KEEP(38) .NE. 0) THEN
        IF (associated(id%root%RG2L)) THEN
          DEALLOCATE(id%root%RG2L)
          NULLIFY(id%root%RG2L)
        ENDIF
      ENDIF
C
C     Memory statistics
C     -----------------------------------
C     If QR (Keep(19)) is not zero, and if
C     the host does not have the information
C     (ie is not slave), send information
C     computed on the slaves during facto
C     to the host.
C     -----------------------------------
C     Note the KEEP(17), KEEP(143) have been bcasted during fac_par_m
      IF ( KEEP(19) .NE. 0 ) THEN
        IF ( KEEP(46) .NE. 1 ) THEN
C         Host was not working during facto_root
C         Send him the information
          IF ( id%MYID .eq. MASTER ) THEN
            CALL MPI_RECV( KEEP(17), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
            CALL MPI_RECV( KEEP(143), 1, MPI_INTEGER, 1, DEFIC_TAG,
     &                   id%COMM, STATUS, IERR )
          ELSE IF ( id%MYID .EQ. 1 ) THEN
            CALL MPI_SEND( KEEP(17), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
            CALL MPI_SEND( KEEP(143), 1, MPI_INTEGER, 0, DEFIC_TAG,
     &                   id%COMM, IERR )
          END IF
        END IF
      END IF
C     --------------------------------
C     Deallocate communication buffers
C     They will be reallocated
C     in the solve.
C     --------------------------------
      IF (allocated(BUFR)) DEALLOCATE(BUFR)
      CALL SMUMPS_BUF_DEALL_SMALL_BUF( IERR )
C//PIV
      IF ((KEEP(50).EQ.2).AND.(KEEP(219).NE.0)) THEN
      CALL SMUMPS_BUF_DEALL_MAX_ARRAY()
      ENDIF
C
C     Check for errors.
C     After SMUMPS_FAC_B every slave is aware of an error. 
C     If master is included in computations, the call below should
C     not be necessary.
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
C
      CALL SMUMPS_EXTRACT_SCHUR_REDRHS(id)
C     return to user singular values
      IF (id%KEEP(19) .NE.0) THEN
       CALL SMUMPS_EXTRACT_SINGULAR_VALUES(id)
      ENDIF
      IF (KEEP(201) .GT. 0) THEN
         IF ((KEEP(201).EQ.1) .OR. (KEEP(201).EQ.2)) THEN
            IF ( I_AM_SLAVE ) THEN
               CALL SMUMPS_OOC_CLEAN_PENDING(IERR)
               IF(IERR.LT.0)THEN
                  id%INFO(1)=IERR
                  id%INFO(2)=0
               ENDIF
            ENDIF
            CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &           id%COMM, id%MYID )
C             We want to collect statistics even in case of 
C             error to understand if it is due to numerical
C             issues
CC            IF ( id%INFO(1) < 0 ) GOTO 500
         END IF
      END IF
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SECFIN(TIME)
        id%DKEEP(94)=real(TIME)
        IF (KEEP(400).GT.0) THEN
C         Facto time above L0_OMP = total time - facto time under L0_OMP
          id%DKEEP(96)=id%DKEEP(94)-id%DKEEP(95)
        ENDIF
      ENDIF
C     Time to process root node:
      CALL MPI_REDUCE( id%DKEEP(99), TMPTIME, 1,
     &                 MPI_REAL,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
      id%DKEEP(99)=TMPTIME
C     =====================================================================
C     COMPUTE MEMORY ALLOCATED BY MUMPS, INFO(16)
C     ---------------------------------------------
      MEM_EFF_ALLOCATED = .TRUE.
      CALL SMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, .TRUE., TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &    , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &     size(id%I8_L0_OMP,2)
     &     )
      IF (KEEP(400) .GT. 0 ) THEN ! L0 activated
        CALL SMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .TRUE., id%KEEP(201),
     &     BLR_STRAT, .TRUE., TOTAL_BYTES_UNDER_L0,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .TRUE.   ! UNDER_L0_OMP
     &     , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &       size(id%I8_L0_OMP,2)
     &        )
        TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
        TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
      ENDIF
       IF (id%KEEP8(24).NE.0) THEN
C         WK_USER is not part of memory allocated by MUMPS
C         and is not counted, id%KEEP8(23) should be zero
          id%INFO(16) = TOTAL_MBYTES
       ELSE
C         Note that even for the case of ICNTL(23)>0 
C         we report here the memory effectively allocated
C         that can be smaller than ICNTL(23) !
          id%INFO(16) = TOTAL_MBYTES
       ENDIF
C       ----------------------------------------------------
C       Centralize memory statistics on the host
C       id%INFOG(18) = size of mem in Mbytes for facto,
C                   for the processor using largest memory
C       id%INFOG(19) = size of mem in Mbytes for facto,
C                   sum over all processors
C       ----------------------------------------------------
        CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                           id%INFO(16), id%INFOG(18), IRANK )
        CALL SMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, id%INFO(16), id%INFOG(18), id%INFOG(19),
     &       id%NSLAVES, IRANK,
     &       id%KEEP(1) )
C     If WK_USER is provided, this excludes WK_USER
      IF (PROK ) THEN
          WRITE(MP,'(A,I12) ')
     & ' ** Eff. min. Space MBYTES for facto             (INFO(16)):',
     &                TOTAL_MBYTES
      ENDIF
C     ========================(INFO(16) RELATED)======================
C     ---------------------------------------
C     COMPUTE EFFECTIVE MEMORY USED INFO(22)
C     ---------------------------------------
      PERLU_ON = .TRUE.
      MEM_EFF_ALLOCATED = .FALSE.
      CALL SMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &     id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &     id%KEEP8(30),
     &     id%NSLAVES, TOTAL_MBYTES, .TRUE., id%KEEP(201),
     &     BLR_STRAT, PERLU_ON, TOTAL_BYTES,
     &     IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &     , .FALSE.   ! UNDER_L0_OMP
     &     , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &       size(id%I8_L0_OMP,2)
     &     )
      IF (KEEP(400) .GT. 0 ) THEN ! L0 activated
        CALL SMUMPS_MAX_MEM( id%KEEP(1),id%KEEP8(1),
     &    id%MYID, N, id%NELT, id%NA(1), id%LNA, id%KEEP8(28),
     &    id%KEEP8(30),
     &    id%NSLAVES, TOTAL_MBYTES_UNDER_L0, .TRUE., id%KEEP(201),
     &    BLR_STRAT, PERLU_ON, TOTAL_BYTES_UNDER_L0,
     &    IDUMMY, BDUMMY, MEM_EFF_ALLOCATED
     &    , .TRUE. ! UNDER_L0_OMP
     &    , id%I8_L0_OMP(1,1), size(id%I8_L0_OMP,1),
     &      size(id%I8_L0_OMP,2)
     &  )
        TOTAL_MBYTES = max (TOTAL_MBYTES,TOTAL_MBYTES_UNDER_L0) 
        TOTAL_BYTES  = max (TOTAL_BYTES, TOTAL_BYTES_UNDER_L0)
      ENDIF
C     -- TOTAL_BYTES and TOTAL_MBYTES includes both static 
C     -- (MAXS) and BLR structures computed as the SUM of the PEAKS 
C     -- (KEEP8(67) + KEEP8(70))
      id%KEEP8(7) = TOTAL_BYTES
C     -- INFO(22) holds the effective space (in Mbytes) used by MUMPS
C     -- (it includes part of WK_USER used if provided by user)
      id%INFO(22) = TOTAL_MBYTES
C     ----------------------------------------------------
C     Centralize memory statistics on the host
C       INFOG(21) = size of effective mem (Mbytes) for facto,
C                   for the processor using largest memory
C       INFOG(22) = size of effective mem (Mbytes) for facto,
C                   sum over all processors
C     ----------------------------------------------------
      CALL MUMPS_MEM_CENTRALIZE( id%MYID, id%COMM,
     &                    id%INFO(22), id%INFOG(21), IRANK )
      IF ( PROKG ) THEN
       IF (PRINT_MAXAVG) THEN
        WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, max in  Mbytes     (INFOG(21)):',
     &  id%INFOG(21)
       ENDIF
       WRITE( MPG,'(A,I12) ')
     &  ' ** Memory effectively used, total in Mbytes    (INFOG(22)):',
     &  id%INFOG(22)
      END IF
      SUM_INFO22_THIS_NODE=0
      CALL MPI_REDUCE( id%INFO(22), SUM_INFO22_THIS_NODE, 1,
     &                 MPI_INTEGER,
     &                 MPI_SUM, 0, id%KEEP(411), IERR )
      CALL MPI_REDUCE( SUM_INFO22_THIS_NODE, MAX_SUM_INFO22_THIS_NODE,
     &                 1, MPI_INTEGER, MPI_MAX, 0, id%COMM, IERR )
      IF (PROKG .AND. PRINT_NODEINFO) THEN
        WRITE(MPG,'(A,I12)')
     &  ' ** Max. effective space per compute node, in MBytes       :',
     &  MAX_SUM_INFO22_THIS_NODE
      ENDIF
C
      IF (I_AM_SLAVE) THEN
       K67 = id%KEEP8(67)
       K68 = id%KEEP8(68)
       K70 = id%KEEP8(70)
       K74 = id%KEEP8(74)
       K75 = id%KEEP8(75)
      ELSE
       K67 = 0_8
       K68 = 0_8
       K70 = 0_8
       K74 = 0_8
       K75 = 0_8
      ENDIF
C     -- Save the number of entries effectively used
C        in main working array S
      CALL MUMPS_SETI8TOI4(K67,id%INFO(21))
C
C
      IF (KEEP(400) .GT.0 ) THEN
        IF (.NOT. I_AM_SLAVE) THEN
          id%DKEEP(95) = 0.0E0
          id%DKEEP(16) = 0.0E0
        ENDIF
        IF (id%NPROCS .GT. 1) THEN
C         Compute average and max (across MPI's)
          CALL MPI_REDUCE(id%DKEEP(95), TMPTIME, 1,
     &    MPI_REAL, MPI_SUM, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) TIMEAVG=dble(TMPTIME)
          CALL MPI_REDUCE(id%DKEEP(16), TMPFLOP, 1,
     &    MPI_REAL, MPI_SUM, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) FLOPAVG=dble(TMPFLOP)
          IF (id%MYID.EQ.MASTER) THEN
            TIMEAVG = TIMEAVG / id%NSLAVES
            FLOPAVG = FLOPAVG / id%NSLAVES
          ENDIF
          CALL MPI_REDUCE(id%DKEEP(95), TMPTIME, 1, 
     &    MPI_REAL, MPI_MAX, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) TIMEMAX=dble(TMPTIME)
          CALL MPI_REDUCE(id%DKEEP(16), TMPFLOP, 1, 
     &    MPI_REAL, MPI_MAX, MASTER, id%COMM, IERR)
          IF (id%MYID.EQ.MASTER) FLOPMAX=dble(TMPFLOP)
C         (PROKG may only be true on master)
          IF ( PROKG ) THEN
            WRITE(MPG,190) FLOPAVG, FLOPMAX
            WRITE(MPG,188) TIMEAVG, TIMEMAX
          ENDIF
        ELSE
C          Print DKEEP(95) directly without reduction
           IF ( PROKG ) THEN
             WRITE(MPG,189) id%DKEEP(16)
             WRITE(MPG,187) id%DKEEP(95)
           ENDIF
        ENDIF
      ENDIF
      IF ( PROKG ) THEN
          IF ( ( KEEP(38).NE.0 .OR. KEEP(20).NE.0 ) .AND.
     &        KEEP(60) .EQ. 0 ) THEN
            WRITE(MPG,186) id%DKEEP(99)
          ENDIF
C         Elapsed time for factorization:
          IF (id%INFO(1) .GE.0) THEN
            WRITE(MPG,180) id%DKEEP(94)
          ELSE
            WRITE(MPG,185) id%DKEEP(94)
          ENDIF
      ENDIF
C
C  Sum RINFO(2) : total number of flops for assemblies
C  Sum RINFO(3) : total number of flops for eliminations
C     Initialize RINFO(4) in case BLR was not activated
      RINFO(4) = RINFO(3)
C 
C  Should work even if the master does some work
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
C     Reduce needed to dimension small working array
C     on all procs during SMUMPS_GATHER_SOLUTION
      KEEP(247) = 0
      CALL MPI_REDUCE( KEEP(246), KEEP(247), 1, MPI_INTEGER, 
     &                 MPI_MAX, MASTER, id%COMM, IERR)
C
C     Reduce compression times: get max compression times
      CALL MPI_REDUCE( id%DKEEP(97), id%DKEEP(98), 1,
     &     MPI_REAL,
     &     MPI_MAX, MASTER, id%COMM, IERR)
C
      CALL MPI_REDUCE( RINFO(2), RINFOG(2), 2,
     &                 MPI_REAL,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MUMPS_REDUCEI8( id%KEEP8(31)+id%KEEP8(64),id%KEEP8(6),
     &                     MPI_SUM, MASTER, id%COMM )
C
      IF (id%MYID.EQ.0) THEN
C      In MegaBytes
       RINFOG(16) = real(id%KEEP8(6)*int(KEEP(35),8))/real(1E6)
       IF (KEEP(201).LE.0) THEN
        RINFOG(16) = ZERO
       ENDIF
      ENDIF
      CALL MUMPS_REDUCEI8( id%KEEP8(48),id%KEEP8(148), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(148), INFOG(9))
C
      CALL MPI_REDUCE( int(id%INFO(10),8), id%KEEP8(128),
     &                 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      IF (id%MYID.EQ.MASTER) THEN
        CALL MUMPS_SETI8TOI4(id%KEEP8(128), id%INFOG(10))
      ENDIF
C     Use MPI_MAX for this one to get largest front size
      CALL MPI_ALLREDUCE( id%INFO(11), INFOG(11), 1, MPI_INTEGER,
     &                 MPI_MAX, id%COMM, IERR)
C     make maximum effective frontal size available on all procs
C     for solve phase
C     (Note that INFO(11) includes root size on root master)
      KEEP(133) = INFOG(11)
      CALL MPI_REDUCE( id%INFO(12), INFOG(12), 3, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( id%INFO(40), INFOG(50), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      CALL MPI_REDUCE( KEEP(103), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(229) = INFOG(25)
      CALL MPI_REDUCE( KEEP(105), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, MASTER, id%COMM, IERR)
      KEEP(230) = INFOG(25)
C
      id%INFO(25) = KEEP(98)
      CALL MPI_ALLREDUCE( id%INFO(25), INFOG(25), 1, MPI_INTEGER,
     &                 MPI_SUM, id%COMM, IERR)
C     Extra copies due to in-place stacking
      CALL MUMPS_REDUCEI8( id%KEEP8(8), id%KEEP8(108), MPI_SUM,
     &                     MASTER, id%COMM )
C     Entries in factors
      CALL MUMPS_SETI8TOI4(id%KEEP8(10), id%INFO(27))
      CALL MUMPS_REDUCEI8( id%KEEP8(10),id%KEEP8(110), MPI_SUM,
     &                     MASTER, id%COMM )
      CALL MUMPS_SETI8TOI4(id%KEEP8(110), INFOG(29))
C     Initialize INFO(28)/INFOG(35) in case BLR not activated
      id%INFO(28)  = id%INFO(27)
      INFOG(35)    = INFOG(29)
C     ==============================
C     LOW-RANK
C     ==============================
      IF ( KEEP(486) .NE. 0 ) THEN  !LR is activated
C Compute and Save local amount of flops in case of BLR
            RINFO(4) = real(FLOP_FRFRONTS + FLOP_FACTO_FR - FLOP_LRGAIN
     &                 + FLOP_COMPRESS  + FLOP_FRFRONTS)
C
C Compute and Save local number of entries in compressed factors
C 
            ITMP8 =  id%KEEP8(10) - int(MRY_LU_LRGAIN,8)
            CALL MUMPS_SETI8TOI4( ITMP8, id%INFO(28))
C
            CALL MPI_REDUCE( MRY_LU_LRGAIN, TMP_MRY_LU_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_LU_FR, TMP_MRY_LU_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_CB_FR, TMP_MRY_CB_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( MRY_CB_LRGAIN, TMP_MRY_CB_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_LRGAIN, TMP_FLOP_LRGAIN
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM_FR, TMP_FLOP_TRSM_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM_LR, TMP_FLOP_TRSM_LR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_FR, TMP_FLOP_UPDATE_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_LR, TMP_FLOP_UPDATE_LR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FRSWAP_COMPRESS,
     &                       TMP_FLOP_FRSWAP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_MIDBLK_COMPRESS,
     &                       TMP_FLOP_MIDBLK_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_UPDATE_LRLR3, TMP_FLOP_UPDATE_LRLR3
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(FLOP_ACCUM_COMPRESS, TMP_FLOP_ACCUM_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_TRSM, TMP_FLOP_TRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_PANEL, TMP_FLOP_PANEL
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FRFRONTS, TMP_FLOP_FRFRONTS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_COMPRESS, TMP_FLOP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_DECOMPRESS, TMP_FLOP_DECOMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_CB_COMPRESS, TMP_FLOP_CB_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_CB_DECOMPRESS,TMP_FLOP_CB_DECOMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( FLOP_FACTO_FR, TMP_FLOP_FACTO_FR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( CNT_NODES,TMP_CNT_NODES
     &                      , 1, MPI_INTEGER,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            IF (id%NPROCS.GT.1) THEN
              FLOP_FACTO_LR = FLOP_FACTO_FR - FLOP_LRGAIN
     &                          + FLOP_COMPRESS + FLOP_FRFRONTS
              CALL MPI_REDUCE( FLOP_FACTO_LR, AVG_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_SUM, MASTER, id%COMM, IERR)
              IF (id%MYID.EQ.MASTER) THEN
                AVG_FLOP_FACTO_LR = AVG_FLOP_FACTO_LR/id%NPROCS
              ENDIF
              CALL MPI_REDUCE( FLOP_FACTO_LR, MIN_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_MIN, MASTER, id%COMM, IERR)
              CALL MPI_REDUCE( FLOP_FACTO_LR, MAX_FLOP_FACTO_LR
     &                        , 1, MPI_DOUBLE_PRECISION,
     &                         MPI_MAX, MASTER, id%COMM, IERR)
            ENDIF ! NPROCS > 1
            CALL MPI_REDUCE( TIME_UPDATE, TMP_TIME_UPDATE
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR1, TMP_TIME_UPDATE_LRLR1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR2, TMP_TIME_UPDATE_LRLR2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_LRLR3, TMP_TIME_UPDATE_LRLR3
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_FRLR, TMP_TIME_UPDATE_FRLR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_UPDATE_FRFR, TMP_TIME_UPDATE_FRFR
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DIAGCOPY, TMP_TIME_DIAGCOPY
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_COMPRESS,TMP_TIME_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_MIDBLK_COMPRESS, 
     &                       TMP_TIME_MIDBLK_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRSWAP_COMPRESS,
     &                       TMP_TIME_FRSWAP_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_CB_COMPRESS, TMP_TIME_CB_COMPRESS
     &                      , 1, MPI_DOUBLE_PRECISION, 
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP, TMP_TIME_DECOMP
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_DECOMP_UCFS, TMP_TIME_DECOMP_UCFS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_LRASM_NIV1, TMP_TIME_LRASM_NIV1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_LRASM_LOCASM2, TMP_TIME_LRASM_LOCASM2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_LRASM_MAPLIG1, TMP_TIME_LRASM_MAPLIG1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_LRASM_CONTRIB2, TMP_TIME_LRASM_CONTRIB2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_FRASM_LOCASM2, TMP_TIME_FRASM_LOCASM2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_FRASM_MAPLIG1, TMP_TIME_FRASM_MAPLIG1
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE(TIME_FRASM_CONTRIB2, TMP_TIME_FRASM_CONTRIB2
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_PANEL, TMP_TIME_PANEL
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_I, TMP_TIME_FAC_I
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_MQ, TMP_TIME_FAC_MQ
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FAC_SQ, TMP_TIME_FAC_SQ
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_LRTRSM, TMP_TIME_LRTRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRTRSM, TMP_TIME_FRTRSM
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_FRFRONTS, TMP_TIME_FRFRONTS
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
            CALL MPI_REDUCE( TIME_LR_MODULE, TMP_TIME_LR_MODULE
     &                      , 1, MPI_DOUBLE_PRECISION,
     &                       MPI_SUM, MASTER, id%COMM, IERR)
        IF (id%MYID.EQ.MASTER) THEN
         IF (id%NPROCS.GT.1) THEN
C rename the stat variable so that COMPUTE_GLOBAL_GAINS can work for any
C number of procs
            MRY_LU_FR            = TMP_MRY_LU_FR
            MRY_LU_LRGAIN        = TMP_MRY_LU_LRGAIN
            MRY_CB_FR            = TMP_MRY_CB_FR
            MRY_CB_LRGAIN        = TMP_MRY_CB_LRGAIN
            FLOP_LRGAIN          = TMP_FLOP_LRGAIN
            FLOP_PANEL           = TMP_FLOP_PANEL
            FLOP_TRSM            = TMP_FLOP_TRSM
            FLOP_TRSM_FR         = TMP_FLOP_TRSM_FR
            FLOP_TRSM_LR         = TMP_FLOP_TRSM_LR
            FLOP_UPDATE_FR       = TMP_FLOP_UPDATE_FR
            FLOP_UPDATE_LR       = TMP_FLOP_UPDATE_LR
            FLOP_UPDATE_LRLR3    = TMP_FLOP_UPDATE_LRLR3
            FLOP_COMPRESS        = TMP_FLOP_COMPRESS
            FLOP_MIDBLK_COMPRESS = TMP_FLOP_MIDBLK_COMPRESS
            FLOP_FRSWAP_COMPRESS = TMP_FLOP_FRSWAP_COMPRESS
            FLOP_ACCUM_COMPRESS  = TMP_FLOP_ACCUM_COMPRESS
            FLOP_CB_COMPRESS     = TMP_FLOP_CB_COMPRESS
            FLOP_DECOMPRESS      = TMP_FLOP_DECOMPRESS
            FLOP_CB_DECOMPRESS   = TMP_FLOP_CB_DECOMPRESS
            FLOP_FRFRONTS        = TMP_FLOP_FRFRONTS
            FLOP_FACTO_FR        = TMP_FLOP_FACTO_FR
            CNT_NODES            = TMP_CNT_NODES
            TIME_UPDATE          = TMP_TIME_UPDATE         /id%NPROCS
            TIME_UPDATE_LRLR1    = TMP_TIME_UPDATE_LRLR1   /id%NPROCS
            TIME_UPDATE_LRLR2    = TMP_TIME_UPDATE_LRLR2   /id%NPROCS
            TIME_UPDATE_LRLR3    = TMP_TIME_UPDATE_LRLR3   /id%NPROCS
            TIME_UPDATE_FRLR     = TMP_TIME_UPDATE_FRLR    /id%NPROCS
            TIME_UPDATE_FRFR     = TMP_TIME_UPDATE_FRFR    /id%NPROCS
            TIME_COMPRESS        = TMP_TIME_COMPRESS       /id%NPROCS
            TIME_MIDBLK_COMPRESS = TMP_TIME_MIDBLK_COMPRESS/id%NPROCS
            TIME_FRSWAP_COMPRESS = TMP_TIME_FRSWAP_COMPRESS/id%NPROCS
            TIME_DIAGCOPY        = TMP_TIME_DIAGCOPY       /id%NPROCS
            TIME_CB_COMPRESS     = TMP_TIME_CB_COMPRESS    /id%NPROCS
            TIME_PANEL           = TMP_TIME_PANEL          /id%NPROCS
            TIME_FAC_I           = TMP_TIME_FAC_I          /id%NPROCS
            TIME_FAC_MQ          = TMP_TIME_FAC_MQ         /id%NPROCS
            TIME_FAC_SQ          = TMP_TIME_FAC_SQ         /id%NPROCS
            TIME_LRTRSM          = TMP_TIME_LRTRSM         /id%NPROCS
            TIME_FRTRSM          = TMP_TIME_FRTRSM         /id%NPROCS
            TIME_FRFRONTS        = TMP_TIME_FRFRONTS       /id%NPROCS
            TIME_LR_MODULE       = TMP_TIME_LR_MODULE      /id%NPROCS
            TIME_DECOMP          = TMP_TIME_DECOMP         /id%NPROCS
            TIME_DECOMP_UCFS     = TMP_TIME_DECOMP_UCFS    /id%NPROCS
            TIME_LRASM_NIV1      = TMP_TIME_LRASM_NIV1     /id%NPROCS
            TIME_LRASM_LOCASM2   = TMP_TIME_LRASM_LOCASM2  /id%NPROCS
            TIME_LRASM_MAPLIG1   = TMP_TIME_LRASM_MAPLIG1  /id%NPROCS
            TIME_LRASM_CONTRIB2  = TMP_TIME_LRASM_CONTRIB2 /id%NPROCS
            TIME_FRASM_LOCASM2   = TMP_TIME_FRASM_LOCASM2  /id%NPROCS
            TIME_FRASM_MAPLIG1   = TMP_TIME_FRASM_MAPLIG1  /id%NPROCS
            TIME_FRASM_CONTRIB2  = TMP_TIME_FRASM_CONTRIB2 /id%NPROCS
         ENDIF
         CALL COMPUTE_GLOBAL_GAINS(id%KEEP8(110),id%RINFOG(3),
     &        id%KEEP8(49), PROKG, MPG)
C         Number of entries in factor  INFOG(35) in 
C         compressed form is updated as long as 
C         BLR is activated,  this independently of the 
C         fact that factors are saved in LR.
          CALL MUMPS_SETI8TOI4(id%KEEP8(49),  id%INFOG(35))
          FRONTWISE = 0
C         WRITE gains also compute stats stored in DKEEP array
          IF (LPOK) THEN
            IF (CNTL(7) < 0.0E0) THEN
C           Warning : using negative values is an experimental and 
C            non recommended setting.
             WRITE(LP,'(/A/,A/,A/,A,A)') 
     &  ' WARNING in BLR input setting',
     &  '          CNTL(7) < 0 is experimental: ',
     &  '          RRQR precision = |CNTL(7| x ||A_pre||, ',
     &  '          where A_pre is the preprocessed matrix as defined',
     &  ' in the Users guide '
            ENDIF
          ENDIF
          CALL SAVEandWRITE_GAINS(FRONTWISE,
     &                KEEP(489), id%DKEEP, N,  id%ICNTL(36),
     &                KEEP(487), KEEP(488), KEEP(490),
     &                KEEP(491), KEEP(50), KEEP(486), 
     &                KEEP(249)*max(KEEP(381), 1),
     &                KEEP(472), KEEP(475), KEEP(478), 
     &                KEEP(480), KEEP(481),
     &                KEEP(483), KEEP(484), 
     &                id%KEEP8(110), id%KEEP8(49),
     &                KEEP(28), id%NPROCS, MPG, PROKG)
C           flops when BLR activated
            RINFOG(14) = id%DKEEP(56)
          ELSE
            RINFOG(14) = 0.0E00
          ENDIF
        IF (id%MYID .eq. MASTER) THEN
          KEEP(399) = KEEP399_SAVE
        ENDIF
      ENDIF
C     ==============================
C     NULL PIVOTS AND RANK-REVEALING
C     ==============================
      IF(KEEP(110) .EQ. 1) THEN
C        -- make available to users the local number of null pivots detected 
C        -- with ICNTL(24) = 1.
         id%INFO(18) = KEEP(109)
         CALL MPI_ALLREDUCE( KEEP(109), KEEP(112), 1, MPI_INTEGER,
     &        MPI_SUM, id%COMM, IERR)
      ELSE
         id%INFO(18)  = 0
         KEEP(109) = 0
         KEEP(112) = 0
      ENDIF
      IF ( associated( id%PIVNUL_LIST) ) THEN
        DEALLOCATE(id%PIVNUL_LIST)
        NULLIFY(id%PIVNUL_LIST)
      ENDIF
C     set INFOG(28) even in case of error 
      IF (id%MYID.EQ.MASTER) THEN
C       INFOG(28) deficiency resulting from ICNTL(24) and ICNTL(56)
        INFOG(28)=KEEP(112)
        IF (KEEP(17).GT.0) THEN
          INFOG(28)=KEEP(112)+KEEP(17)
        ENDIF
      ENDIF
C
      IF (id%INFO(1).GE.0) THEN
C{     PIVNUL_LIST not meaningful in case of error
C      (do not allocate)
       IF (id%MYID.EQ.MASTER) THEN
        IF ( INFOG(28) .GT. 0 ) THEN
          ALLOCATE(id%PIVNUL_LIST(INFOG(28)), stat=IERR)
          IF ( IERR .GT. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=INFOG(28)
          END IF
        ENDIF
       ELSE
C       id%PIVNUL_LIST(1:KEEP(109)) used during sol_driver on slaves 
C       to initialize id%RHSINTR
        IF (KEEP(109).GT.0) THEN
          ALLOCATE(id%PIVNUL_LIST(KEEP(109)), stat=IERR)
          IF ( IERR .GT. 0 ) THEN
           id%INFO(1)=-13
           id%INFO(2)=INFOG(28)
          END IF
        ENDIF
       ENDIF
C}
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1).lt.0 ) GOTO 490
      IF ( (KEEP(19).NE.0) .AND. (KEEP(143) .NE. KEEP(17)) ) THEN
C
C         Raise a warning (on all MPI processes) since determinant or
C         inertia or null pivot list is not consistent with deficiency
C         computed with ICNTL(46)=1.
C
        id%INFO(1) = id%INFO(1)+16
        id%INFO(2) = KEEP(112)+KEEP(143)
        IF (KEEP(118) .GE. 40) THEN
          IF ( PROKG ) THEN
            WRITE(MPG,'(/A,A/,A,A,I8/,A,A,I8/)')
     &      " WARNING: in the context of rank-revealing,",
     &      " the inertia, determinant and pivnul list",  
     &      " are computed with RR (rank-revealing)-LU,", 
     &      " but the deficiency found by RR-LU: ", 
     &      id%INFO(2),
     &      " is different from the deficiency computed",
     &      " with ICNTL(56)>0: ",  KEEP(112)+KEEP(17)
          ENDIF
        ELSE
          IF ( LP .GT. 0 ) THEN
            WRITE(LP,'(/A,A/,A/)')
     &      " ERROR : in the context of rank-revealing,",
     &      " the inertia, determinant and pivnul list",  
     &      " are not correct because RR LU not called "
          ENDIF
        ENDIF
      ENDIF
C     ========================================
C     We now provide to the host the part of
C     PIVNUL_LIST resulting from the processing
C     of the root node and we update id%INFO(18)
C     on the processor holding the root to
C     include null pivots relative to the root
C     ========================================
      IF ( KEEP(109).GT.0 ) THEN
          DO I=1, KEEP(109) 
            id%PIVNUL_LIST(I)= 
     &             PIVNUL_LIST_STRUCT%PIVNUL_LIST(I)
          ENDDO
      ENDIF
      IF (KEEP(17) .NE. 0) THEN
        IF (id%MYID .EQ. ID_ROOT) THEN
C         Include in id%INFO(18) null pivots resulting
C         from deficiency on the root. In this way,
C         the sum of all id%INFO(18) is equal to INFOG(28).
          id%INFO(18)=id%INFO(18)+KEEP(17)
        ENDIF
        IF (ID_ROOT .EQ. MASTER) THEN
          IF (id%MYID.EQ.MASTER) THEN
C           --------------------------------------------------
C           Null pivots of root have been stored in
C           PIVNUL_LIST_STRUCT%PIVNUL_LIST(
C                         KEEP(109)+1:KEEP(109)+KEEP(17) )
C           Shift them at the end of the list because:
C           * this is what we need to build the null space
C           * we would otherwise overwrite them on the host
C             when gathering null pivots from other processors
C           --------------------------------------------------
            DO I= KEEP(17), 1, -1
              id%PIVNUL_LIST(KEEP(112)+I)= 
     &               PIVNUL_LIST_STRUCT%PIVNUL_LIST(KEEP(109)+I)
            ENDDO
          ENDIF
        ELSE
C         ---------------------------------
C         Null pivots of root must be sent
C         from the processor responsible of
C         the root to the host (or MASTER).
C         ---------------------------------
          IF (id%MYID .EQ. ID_ROOT) THEN
            CALL MPI_SEND(
     &               PIVNUL_LIST_STRUCT%PIVNUL_LIST(KEEP(109)+1), 
     &               KEEP(17), MPI_INTEGER, MASTER, ZERO_PIV,
     &               id%COMM, IERR)
          ELSE IF (id%MYID .EQ. MASTER) THEN
            CALL MPI_RECV(id%PIVNUL_LIST(KEEP(112)+1), KEEP(17),
     &               MPI_INTEGER, ID_ROOT, ZERO_PIV,
     &               id%COMM, STATUS, IERR )
          ENDIF
        ENDIF
      ENDIF
C     ===========================
C     gather zero pivots indices
C     on the host node
C     ===========================
C     In case of non working host, the following code also
C     works considering that KEEP(109) is equal to 0 on
C     the non-working host
      IF(KEEP(110) .EQ. 1) THEN
         ALLOCATE(ITMP2(id%NPROCS),stat = IERR )  ! deallocated in 490
         IF ( IERR .GT. 0 ) THEN
            id%INFO(1)=-13
            id%INFO(2)=id%NPROCS
         END IF
         CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
         IF (id%INFO(1).LT.0) GOTO 490
         CALL MPI_GATHER ( KEEP(109),1, MPI_INTEGER, 
     &        ITMP2(1), 1, MPI_INTEGER, 
     &        MASTER, id%COMM, IERR)
         IF(id%MYID .EQ. MASTER) THEN
            POSBUF = ITMP2(1)+1
C           First null pivot of master is in
C           position 1 of global list
            KEEP(220)=1
            DO I = 1,id%NPROCS-1
              IF (ITMP2(I+1).GT.0) THEN
                CALL MPI_RECV(id%PIVNUL_LIST(POSBUF), ITMP2(I+1), 
     &              MPI_INTEGER,I, 
     &              ZERO_PIV, id%COMM, STATUS, IERR)
C              Send position POSBUF of first null pivot of proc I
C              in global list. Will allow to quickly identify during
C              the solve step if one is concerned by a global position
C              K, 0 <= K <= INFOG(28).
               CALL MPI_SEND(POSBUF, 1, MPI_INTEGER, I, ZERO_PIV,
     &              id%COMM, IERR)
               POSBUF = POSBUF + ITMP2(I+1)
              ENDIF
            ENDDO
         ELSE
           IF (KEEP(109).GT.0) THEN
            CALL MPI_SEND( 
     &           PIVNUL_LIST_STRUCT%PIVNUL_LIST(1), KEEP(109), 
     &           MPI_INTEGER,
     &           MASTER,ZERO_PIV, id%COMM, IERR)
            CALL MPI_RECV( KEEP(220), 1, MPI_INTEGER, MASTER, ZERO_PIV,
     &           id%COMM, STATUS, IERR )
           ENDIF
         ENDIF
      ENDIF
      IF (associated( PIVNUL_LIST_STRUCT%PIVNUL_LIST)) THEN
        DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
        NULLIFY(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
      ENDIF
C     =====================================
C     Statistics concerning the determinant
C     =====================================
C
C     1/ on the host better take into account null pivots if scaling:
C
C     Since null pivots are excluded from the computation
C     of the determinant, we also exclude the corresponding
C     scaling entries. Since those entries have already been
C     taken into account before the factorization, we multiply
C     the determinant on the host by the scaling values corresponding
C     to pivots in PIVNUL_LIST.
      IF (id%MYID.EQ.MASTER .AND. LSCAL. AND. KEEP(258).NE.0) THEN
        K = min(KEEP(143), KEEP(17))
        K = max(K, 0)
        DO I = 1, KEEP(112)+ K
c       DO I = 1, id%INFOG(28)  ! all null pivots + singular values
          CALL SMUMPS_UPDATEDETER_SCALING(
     &                            id%ROWSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
          CALL SMUMPS_UPDATEDETER_SCALING(
     &                            id%COLSCA(id%PIVNUL_LIST(I)),
     &                            id%DKEEP(6), KEEP(259))
        ENDDO
      ENDIF
C
C     2/ Swap signs depending on pivoting on each proc
C
      IF (KEEP(258).NE.0) THEN
C       Return the determinant in INFOG(34) and RINFOG(12/13)
C       In case of real arithmetic, initialize
C       RINFOG(13) to 0 (no imaginary part and
C       not touched by SMUMPS_DETER_REDUCTION)
        RINFOG(13)=0.0E0
        IF (KEEP(260).EQ.-1) THEN ! Local to each processor
          id%DKEEP(6)=-id%DKEEP(6)
        ENDIF
C
C       3/ Perform a reduction
C
        CALL SMUMPS_DETER_REDUCTION(
     &           id%COMM, id%DKEEP(6), KEEP(259),
     &           RINFOG(12), INFOG(34), id%NPROCS)
C
C       4/ Swap sign if needed
C
        IF (id%KEEP(50).EQ.0 .AND. id%MYID.EQ. MASTER) THEN
C         Modify sign of determinant according
C         to unsymmetric permutation (max-trans
C         of max-weighted matching)
          IF (id%KEEP(23).NE.0) THEN
            CALL SMUMPS_DETER_SIGN_PERM(
     &           RINFOG(12), id%N,
     &           id%UNS_PERM(1) )
C           Remark that RINFOG(12/13) are modified only
C           on the host but will be broadcast on exit
C           from MUMPS (see SMUMPS_DRIVER)
          ENDIF
        ENDIF
      ENDIF
 490  IF (allocated(ITMP2)) DEALLOCATE(ITMP2)
C
C     =====================================
C     Statistics relative to min/max pivots
C     =====================================
      CALL MPI_REDUCE( id%DKEEP(19), RINFOG(19), 1, 
     &                 MPI_REAL,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(20), RINFOG(20), 1, 
     &                 MPI_REAL,
     &                 MPI_MIN, MASTER, id%COMM, IERR )
      CALL MPI_REDUCE( id%DKEEP(21), RINFOG(21), 1, 
     &                 MPI_REAL,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
C     =========================================
C     Centralized number of swaps for pivoting
C     =========================================
      CALL MPI_REDUCE( id%KEEP8(80), ITEMP8, 1, MPI_INTEGER8,
     &                 MPI_SUM, MASTER, id%COMM, IERR )
      IF (id%MYID .EQ. MASTER) THEN
        CALL MUMPS_SETI8TOI4(ITEMP8,id%INFOG(48))
      ENDIF
C     ==========================================
C     Centralized largest increase of panel size
C     ==========================================
      CALL MPI_REDUCE( id%KEEP(425), id%INFOG(49), 1, MPI_INTEGER,
     &                 MPI_MAX, MASTER, id%COMM, IERR )
      IF ( PROKG ) THEN
C{    -----------------------------
C     PRINT STATISTICS  (on master)
C     -----------------------------
          WRITE(MPG,99984) RINFOG(2),RINFOG(3),KEEP(52), 
     &                    id%KEEP8(148),
     &                    id%KEEP8(128), INFOG(11), id%KEEP8(110)
          IF (id%KEEP(50) == 1 .OR. id%KEEP(50) == 2) THEN
            ! negative pivots
            WRITE(MPG, 99987) INFOG(12)
          END IF
          IF (id%KEEP(50) == 0) THEN
            ! off diag pivots
            WRITE(MPG, 99985) INFOG(12)
          END IF
          IF (id%KEEP(50) .NE. 1) THEN
            ! delayed pivots
            WRITE(MPG, 99982) INFOG(13)
          END IF
          IF (KEEP(97) .NE. 0) THEN
            ! tiny pivots
            WRITE(MPG, '(A,D16.4)') 
     &     ' Effective static pivoting thresh., CNTL(4) =', SEUIL
            WRITE(MPG, 99986) INFOG(25)
          ENDIF
          IF (id%KEEP(50) == 2) THEN
            !number of 2x2 pivots in type 1 nodes
             WRITE(MPG, 99988) KEEP(229)
            !number of 2x2 pivots in type 2 nodes
             WRITE(MPG, 99989) KEEP(230)
          ENDIF
          !number of zero pivots
          IF (KEEP(110) .NE.0) THEN
              WRITE(MPG, 99991) KEEP(112)
          ENDIF
          !Deficiency on root
          IF ( KEEP(19) .ne. 0 )
c         IF ( KEEP(17) .ne. 0 )
     &    WRITE(MPG, 99983) KEEP(17)
          !Total deficiency
          IF (KEEP(110).NE.0.OR.KEEP(19).NE.0)
     &    WRITE(MPG, 99992) KEEP(17)+KEEP(112) ! = INFOG(28)
          IF (id%KEEP(50) .EQ. 1 .OR. id%KEEP(50) .EQ. 2) THEN
            IF (KEEP(110) .NE. 0 .OR .KEEP(19).NE.0) THEN
              WRITE(MPG, 99997) INFOG(50)
            ENDIF
          ENDIF
          !Smallest pivot with also null pivots in abs value
          WRITE(MPG, 99995) RINFOG(19)
          !Smallest pivot in abs value
          WRITE(MPG, 99993) RINFOG(20)
          !Largest pivot in abs value
          WRITE(MPG, 99994) RINFOG(21)
          !value of ICNTL(12) that was effectively used.
          WRITE(MPG, 99996) INFOG(24)
          ! Memory compress
          WRITE(MPG, 99981) INFOG(14)
          ! Extra copies due to ip stack in unsym case
          ! in core case (or OLD_OOC_PANEL)
          IF (id%KEEP8(108) .GT. 0_8) THEN
            WRITE(MPG, 99980) id%KEEP8(108)
          ENDIF
          IF  ((KEEP(60).NE.0) .AND. INFOG(25).GT.0) THEN
          !  Schur on and tiny pivots set in last level 
          ! before the Schur if KEEP(114)=0
           WRITE(MPG, '(A)') 
     & " ** Warning Static pivoting was necessary"
           WRITE(MPG, '(A)') 
     & " ** to factor interior variables with Schur ON"
          ENDIF
          IF (KEEP(258).NE.0) THEN
            WRITE(MPG,99978) RINFOG(12)
            WRITE(MPG,99977) INFOG(34)
          ENDIF
C}
      END IF
* ==========================================
*
*  End of Factorization Phase
*
* ==========================================
C
C  Goto 500 is done when
C  LOAD_INIT
C  OOC_INIT_FACTO
C  MUMPS_FDM_INIT
#if ! defined(NO_FDM_DESCBAND)
C  MUMPS_FDBD_INIT
#endif
#if ! defined(NO_FDM_MAPROW)
C  MUMPS_FMRD_INIT
#endif
C  are all called.
C
 500  CONTINUE
C     Redo free INTARR and DBLARR in case an error occurred
C     after allocating them and before freeing them.
      IF (associated(DBLARR)) THEN
        DEALLOCATE(DBLARR)
        NULLIFY(DBLARR)
      ENDIF
      IF (allocated(INTARR)) THEN
        DEALLOCATE(INTARR)
      ENDIF
      IF ( KEEP(38) .NE. 0) THEN
        IF (associated(id%root%RG2L)) THEN
          DEALLOCATE(id%root%RG2L)
          NULLIFY(id%root%RG2L)
        ENDIF
      ENDIF
#if ! defined(NO_FDM_DESCBAND)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDBD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
#if ! defined(NO_FDM_MAPROW)
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FMRD_END(id%INFO(1))  ! INFO(1): input only
      ENDIF
#endif
      IF (I_AM_SLAVE) THEN
C       Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
C         Store pointer to BLR_ARRAY in MUMPS structure
C         (requires successful factorization otherwise module is freed)
          CALL SMUMPS_BLR_MOD_TO_STRUC(id%BLRARRAY_ENCODING)
        ELSE
C         INFO(1) positive or negative
          CALL SMUMPS_BLR_END_MODULE(id%INFO(1), id%KEEP8, id%KEEP(34))
        ENDIF
      ENDIF
      IF (I_AM_SLAVE) THEN
        CALL MUMPS_FDM_END('A')
C       Terminate BLR module except if it is still needed for solve
        IF ( 
     &       (
     &         (KEEP(486).EQ.2) 
     &       )
     &      .AND. id%INFO(1).GE.0
     &     ) THEN
           CALL MUMPS_FDM_MOD_TO_STRUC('F', id%FDM_F_ENCODING,
     &        id%INFO(1))
           IF (.NOT. associated(id%FDM_F_ENCODING)) THEN
             WRITE(*,*) "Internal error 2 in SMUMPS_FAC_DRIVER"
           ENDIF
        ELSE
           CALL MUMPS_FDM_END('F')
        ENDIF
      ENDIF
C
C  Goto 514 is done when an
C  error occurred in MUMPS_FDM_INIT
C  or (after FDM_INIT but before
C  OOC_INIT)
C
 514  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         IF ((KEEP(201).EQ.1).OR.(KEEP(201).EQ.2)) THEN
            CALL SMUMPS_OOC_END_FACTO(id,IERR)
            IF (id%ASSOCIATED_OOC_FILES) THEN
              id%ASSOCIATED_OOC_FILES = .FALSE.
            ENDIF
            IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
         ENDIF
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE IF (KEEP(201).NE.0) THEN
C           ----------------------------------------
C           In OOC or if KEEP(201).EQ.-1 we always
C           free S at end of factorization. As id%S
C           may be unassociated in case of error
C           during or before the allocation of id%S,
C           we only free S when it was associated.
C           ----------------------------------------
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            id%KEEP8(23)=0_8
         ELSE
C
C         Fwd in facto: free RHS_MUMPS in case
C         it was allocated.
          IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
          NULLIFY(RHS_MUMPS)
          RHS_MUMPS_ALLOCATED = .FALSE.
          IF (.NOT.WK_USER_PROVIDED) THEN
C{
            IF (id%KEEP8(31).EQ.0) THEN
C            No factors stored in id%S
             IF (associated(id%S)) THEN
               DEALLOCATE (id%S)
               NULLIFY(id%S)
               id%KEEP8(23) = 0
             ENDIF
            ELSE
C{            Factors stored in id%S
               TMPpeak = id%KEEP8(73) + id%KEEP8(31) -
     &           (SMUMPS_LBUFR_BYTES+SMUMPS_LBUF)/int(KEEP(35),8)
               Compact_S_Authorized = .FALSE.
C              Set Compact_S_Authorized 
               IF (id%KEEP8(4).GT.0_8) THEN
                 IF (TMPpeak.LT.id%KEEP8(75)) 
     &               Compact_S_Authorized=.TRUE.
               ELSE
                 Compact_S_Authorized = .TRUE.
               ENDIF
              IF (ICNTL49_LOC.EQ.1.AND..NOT.Compact_S_Authorized) THEN
                  IF (PROK) THEN
                   WRITE(MP,'(A,I4,A,I3,A,/A,I8)') 
     &              " ** WARNING ** on MPI proc= ", id%MYID, 
     &              " ICNTL(49)= ", ICNTL49_LOC, 
     &              ", but not enough memory to compact id%S due to ",
     &              "    memory limitation given by ICNTL(23)=",
     &              ICNTL(23)
                  ENDIF
              ELSE IF (
     &           (ICNTL49_LOC.EQ.1.AND.Compact_S_Authorized) 
     &            .OR.
     &           (ICNTL49_LOC.EQ.2) 
     &           ) THEN
C{             Try to compact id%S of size MAXS
                 ALLOCATE(S_IS_POINTERS%A(id%KEEP8(31)), stat=IERR)
                 IF (IERR .GT. 0 ) THEN
                  IF (PROK) THEN
                   WRITE(MP,'(A,I4,A,I3,A)') 
     &              " ** WARNING ** on MPI proc= ", id%MYID, 
     &              " ICNTL(49)= ", ICNTL49_LOC, 
     &              ", but not enough memory to compact id%S "
                  ENDIF
                   id%INFO(1) = id%INFO(1) + 4
                   GOTO 513
                 ENDIF
C
!$      NOMP = OMP_GET_MAX_THREADS()
!$      CHUNK8=int(KEEP(361),8)
!$OMP   PARALLEL DO PRIVATE(I8) SCHEDULE(STATIC, CHUNK8)
!$OMP&  IF ( id%KEEP8(31) > int(KEEP(361),8) .AND. NOMP .GT. 1)
#if defined(__ve__)
!NEC$ IVDEP
#endif
                 DO I8=1_8, id%KEEP8(31)
                   S_IS_POINTERS%A(I8) = id%S(I8)
                 ENDDO
!$OMP   END PARALLEL DO
                 DEALLOCATE (id%S)
                 id%S  => S_IS_POINTERS%A;  NULLIFY(S_IS_POINTERS%A)
                 id%KEEP8(23) = id%KEEP8(31)
C}
              ENDIF
C}
            ENDIF
C}
          ENDIF
         ENDIF
      ELSE  ! host not working
         IF (WK_USER_PROVIDED) THEN
C     at the end of a phase S is always freed when WK_USER provided
            NULLIFY(id%S)
         ELSE
            IF (associated(id%S))  DEALLOCATE(id%S)
            NULLIFY(id%S)   ! in all cases
            id%KEEP8(23)=0_8
         END IF
      END IF
C
C     Goto 513 is done in case of error where LOAD_INIT was
C     called but not the scaling nor OOC_INIT_FACTO.
 513  CONTINUE
      IF ( I_AM_SLAVE ) THEN
         CALL SMUMPS_LOAD_END( id%INFO(1), id%NSLAVES, IERR )
         IF (IERR.LT.0 .AND. id%INFO(1) .GE. 0) id%INFO(1) = IERR
      ENDIF
      CALL MUMPS_PROPINFO( ICNTL(1), id%INFO(1),
     &     id%COMM, id%MYID )
C
C     Goto 516 is done in case of error when GPU initialiwqtion
C     has been performed and scaling was optionally computed but
C     not LOAD_INIT nor OOC_INIT_FACTO. We can then extract
C     scaling arrays in case of error.
 516  CONTINUE
C     --------------------------------------------
C     We now build id%ROWSCA_loc and id%COLSCA_loc
C     in case of successful factorization, in the
C     numbering associated to the fully summed
C     variables of the frontal matrices.
C     This requires the factorization to be
C     successful because otherwise we do not have
C     the final lists of pivots associated to
C     the fronts, including delayed pivots and
C     symmetric/unsymmetric permutations done
C     during the factorization process.
C     --------------------------------------------
      IF (LSCAL .AND. id%INFO(1).GE.0) THEN
        CALL SMUMPS_EXTRACT_SCALING(id)
C       occurs during scaling extraction, keep the error.
        IF ( id%INFO(1) .LT. 0 ) THEN
          IF (LPOK) THEN
            WRITE(LP,*) "Facto OK but error during EXTRACT_SCALING",
     &      id%INFO(1:2)
          ENDIF
        ENDIF
      ENDIF
C
C     Goto 517 is done when an error occurs when GPU initialization
C     has been performed but not LOAD_INIT or OOC_INIT_FACTO, e.g.
C     when an error occurred during the scaling.
 517  CONTINUE
      IF (associated( PIVNUL_LIST_STRUCT%PIVNUL_LIST)) THEN
        DEALLOCATE(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
        NULLIFY(PIVNUL_LIST_STRUCT%PIVNUL_LIST)
      ENDIF
C
C     Goto 530 is done when an error occurs before
C     the calls to GPU_INIT, LOAD_INIT and OOC_INIT_FACTO
 530  CONTINUE
C  Fwd in facto: free RHS_MUMPS in case
C  it was allocated.
      IF (RHS_MUMPS_ALLOCATED) DEALLOCATE(RHS_MUMPS)
      NULLIFY(RHS_MUMPS)
      RHS_MUMPS_ALLOCATED = .FALSE.
C
      id%KEEP8(26) = KEEP826_SAVE
      RETURN
 120  FORMAT(/' Local redistrib: data local/sent           =',I16,I16)
 125  FORMAT(/' Redistrib: total data local/sent           =',I16,I16)
 130  FORMAT(//'****** FACTORIZATION STEP ********'/)
 140  FORMAT(/' Statistics on the scaling phase'
     & /' Elapsed time for scaling                   =',F12.4)
 160  FORMAT(
     & ' Elapsed time to reformat/distribute matrix =',F12.4/)
 166  FORMAT(' Max difference from 1 after scaling the entries',
     &       ' for ONE-NORM (option 7/8)   =',D9.2/)
 170  FORMAT(' STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Size of internal working array S           =',I16/
     &        ' Size of internal working array IS          =',I16/
     &        ' Minimum (ICNTL(14)=0) size of S            =',I16/
     &        ' Minimum (ICNTL(14)=0) size of IS           =',I16/
     &        ' Real space for original matrix             =',I16/
     &        ' Integer space for original matrix          =',I16/
     &        ' INFO(3) Real space for factors (estimated) =',I16/ 
     &        ' INFO(4) Integer space for factors (estim.) =',I16/
     &        ' Maximum frontal size (estimated)           =',I16)
 172  FORMAT(' GLOBAL STATISTICS PRIOR NUMERICAL FACTORIZATION ...'/
     &        ' Number of working processes                =',I16/
     &        ' ICNTL(22) Out-of-core option               =',I16/
     &        ' ICNTL(35) BLR activation (eff. choice)     =',I16/
     &        ' ICNTL(37) BLR CB compression (eff. choice) =',I16/
     &        ' ICNTL(49) Compact workarray S (end facto.) =',I16/
     &        ' ICNTL(14) Memory relaxation                =',I16/
     &        ' INFOG(3) Real space for factors (estimated)=',I16/
     &        ' INFOG(4) Integer space for factors (estim.)=',I16/
     &        ' Maximum frontal size (estimated)           =',I16/
     &        ' Number of nodes in the tree                =',I16/
     &        ' ICNTL(23) Memory allowed (value on host)   =',I16/
     &        '           Sum over all procs               =',I16/
     &        ' Memory provided by user, sum of LWK_USER   =',I16/
     &        ' Effective threshold for pivoting, CNTL(1)  =',D16.4)
 173  FORMAT( ' Perform forward during facto, NRHS         =',I16)
 174  FORMAT( ' KEEP(268) Relaxed pivoting effective value =',I16)
 180  FORMAT(/' Elapsed time for factorization                     =',
     &          F12.4)
 185  FORMAT(/' Elapsed time for (failed) factorization            =',
     & F12.4)
 186  FORMAT(/' Elapsed time to process root node                  =',
     & F12.4)
 187  FORMAT( ' Elapsed time under L0                      =',F12.4)
 188  FORMAT( ' Elapsed time under L0 (avg/max across MPI) =',
     &          F12.4,F12.4)
 189  FORMAT(/' Flops under L0 layer                       =',1PD12.3)
 190  FORMAT(/' Flops under L0 Layer  (avg/max across MPI) =',
     &          1PD12.3,1PD12.3)
99977 FORMAT( ' INFOG(34)  Determinant (base 2 exponent)           =',
     &          I16)
99978 FORMAT( ' RINFOG(12) Determinant (real part)                 =',
     &          F16.8)
99980 FORMAT( ' Extra copies due to In-Place stacking              =',
     &          I16)
99981 FORMAT( ' INFOG (14) Number of memory compress               =',
     &          I16)
99982 FORMAT( ' INFOG (13) Number of delayed pivots                =',
     &      I16)
99983 FORMAT( ' Nb of singularities detected by ICNTL(56)          =',
     &      I16)
99991 FORMAT( ' Nb of null pivots detected by ICNTL(24)            =',
     &      I16)
99992 FORMAT( ' INFOG (28) Estimated deficiency                    =',
     &      I16)
99997 FORMAT( ' INFOG (50) Number of negative AND null pivots      =',
     &      I16)
99995 FORMAT( ' RINFOG(19) Smallest pivot WITH perturbed pivots    =',
     &     1PD10.3)
99993 FORMAT( ' RINFOG(20) Smallest pivot WITHOUT perturbed pivots =',
     &     1PD10.3)
99994 FORMAT( ' RINFOG(21) Largest pivot in absolute value         =',
     &     1PD10.3)
99996 FORMAT( ' INFOG (24) Effective value of ICNTL(12)            =',
     &     I16)
99984 FORMAT(/'Leaving factorization with ...'/
     &        ' RINFOG (2) Operations in node assembly             =',
     &          1PD10.3/                                   
     &        ' ------ (3) Operations in node elimination          =',
     &          1PD10.3/                                       
     &        ' ICNTL  (8) Scaling effectively used                =',
     &          I16/                                           
     &        ' INFOG  (9) Real space for factors                  =',
     &          I16/                                           
     &        ' INFOG (10) Integer space for factors               =',
     &          I16/                                           
     &        ' INFOG (11) Maximum front size                      =',
     &          I16/                                           
     &        ' INFOG (29) Number of entries in factors            =',
     &          I16)                                           
99985 FORMAT( ' INFOG (12) Number of off diagonal pivots           =',
     &          I16)                                           
99986 FORMAT( ' INFOG (25) Number of tiny pivots(static)           =',
     &          I16)                                           
99987 FORMAT( ' INFOG (12) Number of negative pivots               =',
     &          I16)
99988 FORMAT( ' Number of 2x2 pivots in type 1 nodes               =',
     &          I16)
99989 FORMAT( ' Number of 2x2 pivots in type 2 nodes               =',
     &          I16)
      END SUBROUTINE SMUMPS_FAC_DRIVER
C
      SUBROUTINE SMUMPS_PRINT_ALLOCATED_MEM( PROK, PROKG, PRINT_MAXAVG,
     &       MP, MPG, INFO16, INFOG18, INFOG19, NSLAVES, IRANK, KEEP )
      IMPLICIT NONE
C
C  Purpose:
C  =======
C     Print memory allocated during factorization
C     - called at beginning of factorization in full-rank
C     - called at end of factorization in low-rank (because
C       of dynamic allocations)
C
      LOGICAL, INTENT(IN) :: PROK, PROKG, PRINT_MAXAVG
      INTEGER, INTENT(IN) :: MP, MPG, INFO16, INFOG18, INFOG19
      INTEGER, INTENT(IN) :: IRANK, NSLAVES
      INTEGER, INTENT(IN) :: KEEP(500)
C
      IF ( PROKG ) THEN
        IF (PRINT_MAXAVG) THEN
         WRITE( MPG,'(A,I12) ')
     &   ' ** Memory allocated, max in Mbytes             (INFOG(18)):',
     &   INFOG18
        ENDIF
        WRITE( MPG,'(/A,I12) ')
     &   ' ** Memory allocated, total in Mbytes           (INFOG(19)):',
     &    INFOG19
      END IF
      RETURN
      END SUBROUTINE SMUMPS_PRINT_ALLOCATED_MEM
      SUBROUTINE SMUMPS_AVGMAX_STAT8(PROKG, MPG, VAL, NSLAVES,
     &     PRINT_MAXAVG, COMM, MSG)
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      LOGICAL, intent(in) :: PROKG
      INTEGER, intent(in) :: MPG
      INTEGER(8), intent(in) :: VAL
      INTEGER, intent(in) :: NSLAVES
      LOGICAL, intent(in) :: PRINT_MAXAVG
      INTEGER, intent(in) :: COMM
      CHARACTER*48 MSG 
C  Local
      INTEGER(8) MAX_VAL
      INTEGER IERR, MASTER
      REAL LOC_VAL, AVG_VAL
      PARAMETER(MASTER=0)
C
      CALL MUMPS_REDUCEI8( VAL, MAX_VAL, MPI_MAX, MASTER, COMM)
      LOC_VAL = real(VAL)/real(NSLAVES)
      CALL MPI_REDUCE( LOC_VAL, AVG_VAL, 1, MPI_REAL,
     &                 MPI_SUM, MASTER, COMM, IERR )
      IF (PROKG) THEN
        IF (PRINT_MAXAVG) THEN
          WRITE(MPG,100) " Average", MSG, int(AVG_VAL,8)
        ELSE
          WRITE(MPG,110)  MSG, MAX_VAL
        ENDIF
      ENDIF
      RETURN
 100  FORMAT(A8,A48,I18)
 110  FORMAT(A48,I18)
      END SUBROUTINE SMUMPS_AVGMAX_STAT8
C
C
C     ==================================================================
C
      SUBROUTINE SMUMPS_EXTRACT_SCALING(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Extract distributed scaling arrays from SMUMPS_EXTRACT_SCALING
C     In case of unsymmetric permutation, ROWSCA and COLSCA correspond
C     to Dr and Dc, in the expression Dr A Q Dc. In other terms, Dc
C     is compatbile with the front column indices, it does not
C     correspond to the column indices of A, meaning that Q is not
C     needed to just extract the scaling values.
C
C
      TYPE(SMUMPS_STRUC) :: id
      INTEGER, EXTERNAL  :: MUMPS_PROCNODE
C
C     MPI
C     ===
C
      INCLUDE 'mpif.h'
C
C     Local declarations
C     ==================
C
      REAL, DIMENSION(:), POINTER :: COLSCA
      REAL, DIMENSION(:), POINTER :: ROWSCA
      INTEGER, PARAMETER :: MASTER = 0
C
      INTEGER    :: ISTEP, NPIV, LIELL
      INTEGER    :: IERR_MPI, allocok
      INTEGER    :: ISCA
      INTEGER    :: JROW, JCOL, IPOS, JJ ! access to IS
      INTEGER    :: LIW_PASSED
      INTEGER(8) :: LALLOC
C
C     Free and reallocate distributed scaling arrays :
C     - in symmetric, COLSCA_loc points on ROWSCA_loc.
C     - not allocated if KEEP(89)=0
C
      NULLIFY(ROWSCA)
      NULLIFY(COLSCA)
      IF (associated(id%ROWSCA_loc)) THEN
        DEALLOCATE(id%ROWSCA_loc)
        NULLIFY(id%ROWSCA_loc)
      ENDIF
      IF (associated(id%COLSCA_loc)) THEN
        IF (id%KEEP(50) .EQ. 0) THEN
          DEALLOCATE(id%COLSCA_loc)
        ENDIF
        NULLIFY(id%COLSCA_loc)
      ENDIF
C     id%ROWSCA and id%COLSCA are available on master,
C     allocate ROWSCA and COLSCA of order N on other procs
      IF ( id%MYID .EQ. MASTER ) THEN
        ROWSCA => id%ROWSCA
        COLSCA => id%COLSCA
        IF (.NOT. associated(ROWSCA)) THEN
          WRITE(*,*) "Internal error 1 in SMUMPS_EXTRACT_SCALING"
          CALL MUMPS_ABORT()
        ENDIF
        IF (.NOT. associated(COLSCA)) THEN
          WRITE(*,*) "Internal error 2 in SMUMPS_EXTRACT_SCALING"
          CALL MUMPS_ABORT()
        ENDIF
      ELSE
        IF (id%KEEP(50).EQ.0) THEN
          ALLOCATE(ROWSCA(id%N),COLSCA(id%N),stat=allocok)
          LALLOC = int(id%N+id%N,8)
        ELSE
          ALLOCATE(ROWSCA(id%N),stat=allocok)
          COLSCA => ROWSCA
          LALLOC = int(id%N,8)
        ENDIF
        IF (allocok .GT. 0) THEN
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR(LALLOC,id%INFO(2))
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C     Jump to 110 in case of error on ROWSCA or COLSCA
C     on one of the MPI processes.
      IF (id%INFO(1) .LT. 0) GOTO 110
C
      IF ( id%KEEP(89) .GT. 0) THEN
        IF (id%KEEP(50).EQ.0) THEN
          ALLOCATE(id%ROWSCA_loc(id%KEEP(89)),
     &             id%COLSCA_loc(id%KEEP(89)),stat=allocok)
          LALLOC = int(id%KEEP(89),8)*2_8
        ELSE
          ALLOCATE(id%ROWSCA_loc(id%KEEP(89)),stat=allocok)
          id%COLSCA_loc => id%ROWSCA_loc
          LALLOC = int(id%KEEP(89),8)
        ENDIF
        IF (allocok .GT. 0) THEN
          id%INFO(1)=-13
          CALL MUMPS_SET_IERROR(LALLOC,id%INFO(2))
        ENDIF
      ENDIF
      CALL MUMPS_PROPINFO( id%ICNTL(1), id%INFO(1),
     &                        id%COMM, id%MYID )
C     Jump to 100 in case of error (we free everything)
      IF (id%INFO(1) .LT. 0) GOTO 100
      CALL MPI_BCAST(ROWSCA(1), id%N, MPI_REAL,
     &               MASTER, id%COMM, IERR_MPI)
      IF (id%KEEP(50) .EQ. 0) THEN
        CALL MPI_BCAST(COLSCA(1), id%N, MPI_REAL,
     &                 MASTER, id%COMM, IERR_MPI)
      ENDIF
      LIW_PASSED = max(id%KEEP(32),1)
      ISCA = 1
      IF ( id%MYID .ne. MASTER .OR.
     &     id%KEEP(46) .eq. 1 ) THEN ! I_AM_SLAVE
        DO ISTEP = 1, id%KEEP(28)
          IF ( id%MYID_NODES.EQ. MUMPS_PROCNODE(
     &         id%PROCNODE_STEPS(ISTEP),
     &         id%KEEP(199) ) ) THEN
              CALL MUMPS_SOL_GET_NPIV_LIELL_IPOS( ISTEP, id%KEEP,
     &        NPIV, LIELL, IPOS,
     &        id%IS(1), LIW_PASSED, id%PTLUST_S(1), id%STEP(1), id%N)
            IF ( id%KEEP(50) .EQ. 0 ) THEN
C             Row indices:
              JROW = IPOS + 1
C             Column indices:
              JCOL = IPOS + 1 + LIELL
            ELSE
C             Use row indices because column indices may have
C             been set to negative to flag 2x2 pivots
              JROW = IPOS + 1
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              DO JJ = 1, NPIV
                id%ROWSCA_loc(ISCA+JJ-1) = ROWSCA(id%IS(JROW+JJ-1))
                id%COLSCA_loc(ISCA+JJ-1) = COLSCA(id%IS(JCOL+JJ-1))
              ENDDO
            ELSE
              DO JJ = 1, NPIV
                id%ROWSCA_loc(ISCA+JJ-1) = ROWSCA(id%IS(JROW+JJ-1))
              ENDDO
            ENDIF
            ISCA = ISCA + NPIV
          ENDIF
        ENDDO
      ENDIF
C     End of EXTRACT_SCALING, we keep id%ROWSCA_loc and id%COLSCA_loc
C     but free ROWSCA and COLSCA
      GOTO 110
      RETURN
 100  CONTINUE
C     Exit with error, free what was allocated
      IF (associated(id%ROWSCA_loc)) THEN
        DEALLOCATE(id%ROWSCA_loc)
        NULLIFY(id%ROWSCA_loc)
      ENDIF
      IF (associated(id%COLSCA_loc)) THEN
        IF (id%KEEP(50) .EQ. 0) THEN
          DEALLOCATE(id%COLSCA_loc)
        ENDIF
        NULLIFY(id%COLSCA_loc)
      ENDIF
 110  CONTINUE
C     Free local ROWSCA and COLSCA arrays
      IF ( id%MYID .NE. 0) THEN
        IF (associated(ROWSCA)) DEALLOCATE(ROWSCA)
        IF ( id%KEEP(50) .EQ. 0 ) THEN
          IF (associated(COLSCA)) DEALLOCATE(COLSCA)
        ENDIF
      ENDIF
      NULLIFY(ROWSCA)
      NULLIFY(COLSCA)
      RETURN
      END SUBROUTINE SMUMPS_EXTRACT_SCALING
C
C     ==================================================================
C
      SUBROUTINE SMUMPS_EXTRACT_SCHUR_REDRHS(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     Extract the Schur and possibly also the reduced right-hand side
C     (if Fwd in facto) from the processor working on Schur and copy
C     it into the user datastructures id%SCHUR and id%REDRHS on the host.
C     This routine assumes that the integer list of the Schur has not
C     been permuted and still corresponds to LISTVAR_SCHUR.
C
C     If the Schur is centralized, the master of the Schur holds the
C     Schur and possibly also the reduced right-hand side.
C     If the Schur is distribued (already built in user's datastructure),
C     then the master of the Schur may hold the reduced right-hand side,
C     in which case it is available in root%RHS_CNTR_MASTER_ROOT.
C     
      TYPE(SMUMPS_STRUC) :: id
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_SCHUR, SIZE_SCHUR, LD_SCHUR, IB, BL4
      INTEGER(4) :: I4 ! 32-bit even in 64-bit version
      INTEGER :: ROW_LENGTH, I
      INTEGER(8) :: SURFSCHUR8, BL8, SHIFT8
      INTEGER(8) :: ISCHUR_SRC, ISCHUR_DEST, ISCHUR_SYM, ISCHUR_UNS
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Schur option off
      IF (id%KEEP(60) .EQ. 0) RETURN
C     Get Schur id
      ID_SCHUR =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(max(id%KEEP(20),id%KEEP(38)))),
     &    id%KEEP(199))
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_SCHUR = ID_SCHUR + 1
      END IF
C     Get size of Schur
      IF (id%MYID.EQ.ID_SCHUR) THEN
        IF (id%KEEP(60).EQ.1) THEN
C         Sequential Schur
          LD_SCHUR =
     &    id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))+2+id%KEEP(IXSZ))
          SIZE_SCHUR = LD_SCHUR - id%KEEP(253)
        ELSE
C         Parallel Schur
          LD_SCHUR   = -999999 ! not used
          SIZE_SCHUR = id%root%TOT_ROOT_SIZE
        ENDIF
      ELSE IF (id%MYID .EQ. MASTER) THEN
        SIZE_SCHUR = id%KEEP(116)
        LD_SCHUR = -44444 ! Not used
      ELSE
C       Proc is not concerned with Schur, return
        RETURN
      ENDIF
      SURFSCHUR8 = int(SIZE_SCHUR,8)*int(SIZE_SCHUR,8)
C     =================================
C     Case of parallel Schur: if REDRHS
C     was requested, obtain it directly
C     from id%root%RHS_CNTR_MASTER_ROOT
C     =================================
      IF (id%KEEP(60) .GT. 1) THEN
        IF (id%KEEP(221).EQ.1 .AND. id%KEEP(252).GT.0) THEN
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
              CALL scopy(SIZE_SCHUR,
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1), 1,
     &             id%REDRHS((I-1)*id%LREDRHS+1), 1)
            ELSE
              IF (id%MYID.EQ.ID_SCHUR) THEN
C               Send
                CALL MPI_SEND(
     &             id%root%RHS_CNTR_MASTER_ROOT((I-1)*SIZE_SCHUR+1),
     &             SIZE_SCHUR,
     &             MPI_REAL,
     &             MASTER, TAG_SCHUR,
     &             id%COMM, IERR )
              ELSE ! MYID.EQ.MASTER
C               Receive
                CALL MPI_RECV( id%REDRHS((I-1)*id%LREDRHS+1),
     &             SIZE_SCHUR,
     &             MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &             id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
          ENDDO
C         ------------------------------
C         In case of parallel Schur, we
C         free root%RHS_CNTR_MASTER_ROOT
C         ------------------------------
          IF (id%MYID.EQ.ID_SCHUR) THEN
            DEALLOCATE(id%root%RHS_CNTR_MASTER_ROOT)
            NULLIFY   (id%root%RHS_CNTR_MASTER_ROOT)
          ENDIF
        ENDIF
C       return because this is all we need to do
C       in case of parallel Schur complement
        RETURN
      ENDIF
C     ============================
C     Centralized Schur complement
C     ============================
C     PTRAST has been freed at the moment of calling this
C     routine. Schur is available through
C     PTRFAC(IW( PTLUST_S( STEP(KEEP(20)) ) + 4 +KEEP(IXSZ) ))
      IF (id%KEEP(252).EQ.0) THEN
C       CASE 1 (ORIGINAL CODE):
C       Schur is contiguous on ID_SCHUR
        IF ( ID_SCHUR .EQ. MASTER ) THEN ! Necessarily equals id%MYID
C         ---------------------
C         Copy Schur complement
C         ---------------------
          CALL SMUMPS_COPYI8SIZE( SURFSCHUR8,
     &      id%S(id%PTRFAC(id%STEP(id%KEEP(20)))),
     &      id%SCHUR(1) )
        ELSE
C         -----------------------------------------
C         The processor responsible of the Schur
C         complement sends it to the host processor
C         Use blocks to avoid too large messages.
C         -----------------------------------------
          BL8=int(huge(I4)/id%KEEP(35)/10,8)
          DO IB=1, int((SURFSCHUR8+BL8-1_8) / BL8)
            SHIFT8 = int(IB-1,8) * BL8                ! Where to send
            BL4    = int(min(BL8,SURFSCHUR8-SHIFT8)) ! Size of block
            IF ( id%MYID .eq. ID_SCHUR ) THEN
C             Send Schur complement
              CALL MPI_SEND( id%S( SHIFT8 +
     &          id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ)))),
     &          BL4,
     &          MPI_REAL,
     &          MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
            ELSE IF ( id%MYID .eq. MASTER ) THEN
C             Receive Schur complement
              CALL MPI_RECV( id%SCHUR(1_8 + SHIFT8),
     &                     BL4,
     &                     MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                     id%COMM, STATUS, IERR )
            END IF
          ENDDO
        END IF
      ELSE
C       CASE 2 (Fwd in facto): Schur is not contiguous on ID_SCHUR,
C       process it row by row.
C
C       2.1: We first centralize Schur complement into id%SCHUR
        ISCHUR_SRC = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &               +4+id%KEEP(IXSZ)))
        ISCHUR_DEST= 1_8
        DO I=1, SIZE_SCHUR
          ROW_LENGTH = SIZE_SCHUR
          IF (ID_SCHUR.EQ.MASTER) THEN ! Necessarily = id%MYID
            CALL scopy(ROW_LENGTH, id%S(ISCHUR_SRC), 1,
     &                 id%SCHUR(ISCHUR_DEST),1)
          ELSE
            IF (id%MYID.EQ.ID_SCHUR) THEN
C             Send
              CALL MPI_SEND( id%S(ISCHUR_SRC), ROW_LENGTH,
     &        MPI_REAL,
     &        MASTER, TAG_SCHUR,
     &        id%COMM, IERR )
            ELSE
C             Recv
              CALL MPI_RECV( id%SCHUR(ISCHUR_DEST),
     &                   ROW_LENGTH,
     &                   MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &                   id%COMM, STATUS, IERR )
            ENDIF
          ENDIF
          ISCHUR_SRC = ISCHUR_SRC+int(LD_SCHUR,8)
          ISCHUR_DEST= ISCHUR_DEST+int(SIZE_SCHUR,8)
        ENDDO
C       2.2: Get REDRHS on host
C       2.2.1: Symmetric => REDRHS is available in last KEEP(253)
C              rows of Schur structure on ID_SCHUR
C       2.2.2: Unsymmetric => REDRHS corresponds to last KEEP(253)
C              columns. However it must be transposed.
        IF (id%KEEP(221).EQ.1) THEN ! Implies Fwd in facto
          ISCHUR_SYM = id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8) *
     &                    int(LD_SCHUR,8)
          ISCHUR_UNS =
     &                 id%PTRFAC(id%IS(id%PTLUST_S(id%STEP(id%KEEP(20)))
     &                    +4+id%KEEP(IXSZ))) + int(SIZE_SCHUR,8)
          ISCHUR_DEST = 1_8
          DO I = 1, id%KEEP(253)
            IF (ID_SCHUR .EQ. MASTER) THEN ! necessarily = id%MYID
              IF (id%KEEP(50) .EQ. 0) THEN
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ELSE
                CALL scopy(SIZE_SCHUR, id%S(ISCHUR_SYM), 1,
     &                     id%REDRHS(ISCHUR_DEST), 1)
              ENDIF
            ELSE
              IF (id%MYID .NE. MASTER) THEN
                IF (id%KEEP(50) .EQ. 0) THEN
C                 Use id%S(ISCHUR_SYM) as temporary contig. workspace
C                 of size SIZE_SCHUR. 
                  CALL scopy(SIZE_SCHUR, id%S(ISCHUR_UNS), LD_SCHUR,
     &            id%S(ISCHUR_SYM), 1)
                ENDIF
                CALL MPI_SEND(id%S(ISCHUR_SYM), SIZE_SCHUR,
     &          MPI_REAL, MASTER, TAG_SCHUR,
     &          id%COMM, IERR )
              ELSE
                CALL MPI_RECV(id%REDRHS(ISCHUR_DEST),
     &          SIZE_SCHUR, MPI_REAL, ID_SCHUR, TAG_SCHUR,
     &          id%COMM, STATUS, IERR )
              ENDIF
            ENDIF
            IF (id%KEEP(50).EQ.0) THEN
              ISCHUR_UNS = ISCHUR_UNS + int(LD_SCHUR,8)
            ELSE
              ISCHUR_SYM = ISCHUR_SYM + int(LD_SCHUR,8)
            ENDIF
            ISCHUR_DEST = ISCHUR_DEST + int(id%LREDRHS,8)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_EXTRACT_SCHUR_REDRHS
      SUBROUTINE SMUMPS_EXTRACT_SINGULAR_VALUES(id)
      USE SMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose
C  =======
C
C     
      TYPE(SMUMPS_STRUC) :: id
C
C  Local variables
C  ===============
C
      INCLUDE 'mpif.h'
      INCLUDE 'mumps_tags.h'
      INCLUDE 'mumps_headers.h'
      INTEGER :: STATUS(MPI_STATUS_SIZE)
      INTEGER :: IERR
      INTEGER, PARAMETER :: MASTER = 0
      INTEGER :: ID_ROOT, ALLOCOK
C
C  External functions
C  ==================
C
      INTEGER MUMPS_PROCNODE
      EXTERNAL MUMPS_PROCNODE
C     Quick return in case factorization did not terminate correctly
      IF (id%INFO(1) .LT. 0) RETURN
C     Quick return if Postponing + rank revealing option off
      IF (id%KEEP(19) .EQ. 0) RETURN
C     Get Root id
      ID_ROOT =MUMPS_PROCNODE(
     &    id%PROCNODE_STEPS(id%STEP(id%KEEP(20))),
     &    id%KEEP(199))
      IF ( id%KEEP( 46 )  .NE. 1 ) THEN
        ID_ROOT = ID_ROOT + 1
      END IF  
C     =================================
C     Singular values are stored in 
C     root%SINGULAR_VALUES
C     We copy it to id%SINGULAR_VALUES
C     =================================
      IF ((ID_ROOT.EQ.id%MYID).AND.(id%MYID.EQ.MASTER)) THEN
C        write(6,*) " singular_values already on host"
         IF (associated(id%SINGULAR_VALUES))
     &        DEALLOCATE(id%SINGULAR_VALUES)
         NULLIFY(id%SINGULAR_VALUES)
         id%NB_SINGULAR_VALUES=id%root%NB_SINGULAR_VALUES
         ALLOCATE(id%SINGULAR_VALUES(id%NB_SINGULAR_VALUES)
     &        , stat=ALLOCOK)
         IF (ALLOCOK.GT.0) THEN
            id%INFO(1)=-13
            id%INFO(2)= id%NB_SINGULAR_VALUES
            RETURN
         END IF
         CALL scopy(id%NB_SINGULAR_VALUES,
     &        id%root%SINGULAR_VALUES(1), 1,
     &        id%SINGULAR_VALUES(1), 1)
      ELSE
         IF (id%MYID.EQ.ID_ROOT) THEN
C     Send
C           write(6,*) " id%MYID sends singular_values "
            CALL MPI_SEND(
     &           id%root%NB_SINGULAR_VALUES,
     &           1,
     &           MPI_INTEGER,
     &           MASTER, TAG_ROOT1,
     &           id%COMM, IERR )
            CALL MPI_SEND(
     &           id%root%SINGULAR_VALUES(1),
     &           id%root%NB_SINGULAR_VALUES,
     &           MPI_REAL,
     &           MASTER, TAG_ROOT2,
     &           id%COMM, IERR )
         ELSEIF (id%MYID.EQ.MASTER) THEN
C     Receive
            CALL MPI_RECV( id%NB_SINGULAR_VALUES,
     &           1,
     &           MPI_INTEGER, ID_ROOT, TAG_ROOT1,
     &           id%COMM, STATUS, IERR )
            IF (associated(id%SINGULAR_VALUES))
     &           DEALLOCATE(id%SINGULAR_VALUES)
            NULLIFY(id%SINGULAR_VALUES)
            ALLOCATE(id%SINGULAR_VALUES(id%NB_SINGULAR_VALUES)
     &           , stat=ALLOCOK)
            IF (ALLOCOK.GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)= id%NB_SINGULAR_VALUES
               RETURN
            END IF
            CALL MPI_RECV( id%SINGULAR_VALUES(1),
     &           id%NB_SINGULAR_VALUES,
     &           MPI_REAL, ID_ROOT, TAG_ROOT2,
     &           id%COMM, STATUS, IERR )
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_EXTRACT_SINGULAR_VALUES
      SUBROUTINE SMUMPS_SET_NOMP_MAX(KEEP281, KEEP361,
     &                               N, NOMP_MAX)
!$    USE OMP_LIB
C
C  Purpose
C  =======
C     set NOMP_MAX from KEEP(281)
C     on output NOMP_MAX >=0
C
C  Parameters
C  ==========
C
      INTEGER, INTENT(IN)    :: KEEP281, KEEP361, N
      INTEGER, INTENT(OUT)   :: NOMP_MAX
C     
C     Local variables
C   
      INTEGER :: NOMP
C
C     out-of-range entries treated as -1
      NOMP_MAX= max(-1, KEEP281)
      NOMP    = 1
!$    NOMP = omp_get_max_threads()
      IF (NOMP_MAX.EQ.-1) THEN
C     automatic setting
       IF (N.LE.KEEP361) THEN
         NOMP_MAX = 0
         RETURN
       ENDIF
       IF (NOMP.GT.1) THEN
C       conservative because of memory allocation
        NOMP_MAX = min(NOMP, 10)
       ELSE
C       no multithreading and all parallel do suppressed
        NOMP_MAX = 0
       ENDIF
      ELSE
C     NOMP_MAX >=0
C     use provided value
       NOMP_MAX = min(NOMP_MAX, NOMP)
      ENDIF
C
      RETURN
      END SUBROUTINE SMUMPS_SET_NOMP_MAX
