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
      SUBROUTINE SMUMPS_SIMSCALEABS(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &     WRKR_TH, LWRKR_TH,
     &     WRKC_TH, LWRKC_TH,
     &     NOMP_MAX, SYM, NB1, NB2, NB3, EPS,
     &     ONENORMERR,INFNORMERR)
C----------------------------------------------------------------------
C    IF SYM=0 CALLs unsymmetric variant SMUMPS_SIMSCALEABSUNS.
C    IF SYM=2 CALLS symmetric variant where only one of a_ij and a_ji 
C         is stored. SMUMPS_SIMSCALEABSSYM
C---------------------------------------------------------------------
C    For details, see the two subroutines below
C         SMUMPS_SIMSCALEABSUNS and SMUMPS_SIMSCALEABSSYM
C ---------------------------------------------------------------------
C
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) NZ_loc
      INTEGER(8) :: IWRKSZ
      INTEGER(8) :: ISZWRKRC, LWRKR_TH, LWRKC_TH
      INTEGER    :: NOMP_MAX
      INTEGER M, N, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER(8) :: INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      REAL A_loc(NZ_loc)
      INTEGER RPARTVEC(M)
      INTEGER RSNDRCVSZ(2*NUMPROCS)
      INTEGER CPARTVEC(N) 
      INTEGER CSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER(8) :: REGISTRE(12)
      REAL ROWSCA(M)
      REAL COLSCA(N)
      REAL WRKRC(ISZWRKRC)
      REAL WRKR_TH(LWRKR_TH,max(NOMP_MAX,1))
      REAL WRKC_TH(LWRKC_TH,max(NOMP_MAX,1))
      REAL ONENORMERR,INFNORMERR
C     LOCALS
C     FOR the scaling phase
      INTEGER SYM, NB1, NB2, NB3
      REAL EPS
C     EXTERNALS
      EXTERNAL SMUMPS_SIMSCALEABSUNS,SMUMPS_SIMSCALEABSSYM, 
     &         SMUMPS_INITREAL
C     MUST HAVE IT
      INTEGER I
      INTEGER, PARAMETER :: K361 = 2048
!$    INTEGER    :: NOMP 
!$    INTEGER    :: CHUNK
!$    IF (NOMP_MAX.GT.0) THEN
!$     NOMP  = omp_get_max_threads()
!$     CHUNK= max(K361/2, (N+NOMP-1) / NOMP )
!$    ENDIF
      IF(SYM.EQ.0) THEN
         CALL SMUMPS_SIMSCALEABSUNS(IRN_loc, JCN_loc, A_loc, 
     &        NZ_loc,
     &        M, N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, CPARTVEC,
     &        RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, COLSCA, WRKRC, ISZWRKRC, 
     &        WRKR_TH, LWRKR_TH,
     &        WRKC_TH, LWRKC_TH, NOMP_MAX,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
      ELSE
         CALL SMUMPS_SIMSCALEABSSYM(IRN_loc, JCN_loc, A_loc, 
     &        NZ_loc,
     &        N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, 
     &        RSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, WRKRC, ISZWRKRC,  
     &        WRKR_TH, LWRKR_TH, NOMP_MAX,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
         IF (OP.EQ.2) THEN
          IF (NOMP_MAX.LE.0) THEN
            DO I=1,N
             COLSCA(I) = ROWSCA(I)
            ENDDO
          ELSE 
!$OMP      PARALLEL DO PRIVATE(I)
!$OMP&     SCHEDULE(STATIC,CHUNK)
!$OMP&     IF ( N > K361 .AND. NOMP .GT. 1)
            DO I=1,N
             COLSCA(I) = ROWSCA(I)
            ENDDO
          ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SIMSCALEABS
      SUBROUTINE SMUMPS_SIMSCALEABSUNS(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC, 
     &     WRKR_TH, LWRKR_TH,
     &     WRKC_TH, LWRKC_TH, NOMP_MAX,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
C----------------------------------------------------------------------
C Input parameters:
C      M, N: size of matrix (in general M=N, but the algorithm 
C            works for rectangular matrices as well (norms other than
C            inf-norm are not possible mathematically in this case).
C      NUMPROCS, MYID, COMM: guess what are those
C      RPARTVEC:  row partvec to be filled when OP=1
C      CPARTVEC:  col partvec to be filled when OP=1
C      RSNDRCVSZ: send recv sizes for row operations. 
C                 to be filled when OP=1
C      CSNDRCVSZ: send recv sizes for col operations. 
C                 to be filled when OP=1
C      REGISTRE:  to store some pointers (size etc)
C      IWRK: working space. when OP=1 IWRKSZ.GE.4*MAXMN
C            when OP=2 INTSZ portion is used. Thus, IWRKSZ>INTSZ 
C            when OP=2
C      IWRKSZ: size
C      INTSZ: to be computed when OP=1, necessary integer space to run 
C             scaling algo when OP=2
C      RESZ:  to be computed when OP=1, necessary real space to run 
C             scaling algo when OP=2
C      OP: 
C          =1 estimation of memory and construction of partvecs
C           writes into RPARTVEC,CPARTVEC,RSNDRCVSZ,CSNDRCVSZ,REGISTRE
C           does not access WRKRC, uses IWRK as workspace
C           computes INTSZ and RESZ.
C          =2 Compute scalings 
C           restores pointers from REGISTRE, 
C           stores communication structure in IWRK (from the start). 
C
C      ROWSCA: space for row scaling factor; has size M
C      COLSCA: space for col scaling factor; has size N
C      WRKRC: real working space. when OP=1, is not accessed. Thus, it
C             can be declared to be of size 1 at OP=1 call.
C      ISZWRKRC: size
C      SYM: is matrix symmetric
C      NB1, NB2, NB3: algo runs 
C                     NB1 iters of inf-norm (default  1/1), 
C                     NB2 iters of 1-norm   (default  3/10),
C                     NB3 iters of inf-norm (default  3/10).
C          in succession.
C      EPS: tolerance for concergence. 
C           IF EPS < 0.R0 then does not test convergence.
C           If convergence occured during the first set of inf-norm
C           iterations, we start performing one-norm iterations.
C           If convergence occured during the one-norm iterations,
C           we start performing the second set of inf-norm iterations.
C           If convergence occured during the second set of inf-norm,
C           we prepare to return.
C     ONENORMERR : error in one norm scaling (associated with the scaling 
C                  arrays of the previous iterations), 
C     INFNORMERR : error in inf norm scaling (associated with the scaling 
C                  arrays of the previous iterations).
C---------------------------------------------------------------------
C On input:
C      OP=1==>Requirements
C             IWRKSZ.GE.4*MAXMN
C             RPARTVEC  of size M
C             CPARTVEC  of size N
C             RSNDRCVSZ of size 2*NUMPROCS
C             CSNDRCVSZ of size 2*NUMPROCS
C             REGISTRE  of size 12
C             
C      OP=2==>Requirements
C             INTSZ .GE. REGISTRE(11)
C             RESZ  .GE. REGISTRE(12) 
C---------------------------------------------------------------------
C On output:
C     ROWSCA and COLSCA 
C            at processor 0 of COMM: complete factors.
C            at other processors   : only the ROWSCA(i) or COLSCA(j) 
C            for which there is a nonzero a_i* or a_*j are useful.
C     ONENORMERR : error in one norm scaling 
C                = -1.0 if iter2=0.
C     INFNORMERR : error in inf norm scaling 
C                = inf norm error at iter3 if iter3 > 0
C                = inf norm error at iter1 if iter1 > 0, iter3=0
C                = -1.0 if iter1=iter3=0
C ---------------------------------------------------------------------
C References:
C     The scaling algorithms are based on those discussed in
C     [1] D. Ruiz, "A scaling algorithm to equilibrate both rows and 
C         columns norms in matrices", Tech. Rep. Rutherford 
C         Appleton Laboratory, Oxon, UK and ENSEEIHT-IRIT, 
C         Toulouse, France, RAL-TR-2001-034 and RT/APO/01/4, 2001.
C     [2] D. Ruiz and B. Ucar, "A symmetry preserving algorithm for
C         matrix scaling", in preparation as of Jan'08.
C
C     The parallelization approach is discussed in
C     [3] P. R. Amestoy, I. S. Duff, D. Ruiz, and B. Ucar,
C         "A parallel matrix scaling algorithm".
C         In proceedings of VECPAR'08-International Meeting-High 
C         Performance Computing for Computational Science, Jan'08.
C     and was supported by ANR-SOLSTICE project (ANR-06-CIS6-010)
C ---------------------------------------------------------------------
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) :: NZ_loc
      INTEGER(8) :: IWRKSZ, INTSZ
      INTEGER(8) :: ISZWRKRC, LWRKR_TH, LWRKC_TH
      INTEGER    :: M, N, OP
      INTEGER    :: NUMPROCS, MYID, COMM, NOMP_MAX
      INTEGER(8) :: RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      REAL A_loc(NZ_loc)
      INTEGER RPARTVEC(M) 
      INTEGER CPARTVEC(N) 
      INTEGER RSNDRCVSZ(2*NUMPROCS)
      INTEGER CSNDRCVSZ(2*NUMPROCS)
      INTEGER(8) :: REGISTRE(12)
      INTEGER IWRK(IWRKSZ)
      REAL ROWSCA(M)
      REAL COLSCA(N)
      REAL WRKRC(ISZWRKRC)
      REAL WRKR_TH(LWRKR_TH,max(NOMP_MAX,1))
      REAL WRKC_TH(LWRKC_TH,max(NOMP_MAX,1))
      REAL ONENORMERR,INFNORMERR
C     LOCALS
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER ICSNDRCVNUM, OCSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER ICSNDRCVVOL, OCSNDRCVVOL
      INTEGER  INUMMYR, INUMMYC
C IMPORTANT POINTERS
      INTEGER(8) :: IMYRPTR,IMYCPTR 
      INTEGER(8) :: IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER(8) :: ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER(8) :: ICNGHBPRCS, ICSNDRCVIA,ICSNDRCVJA
      INTEGER(8) :: OCNGHBPRCS, OCSNDRCVIA,OCSNDRCVJA
      INTEGER(8) :: ISTATUS, REQUESTS, TMPWORK
      INTEGER(8) :: ITDRPTR, ITDCPTR, ISRRPTR
      INTEGER(8) :: OSRRPTR, ISRCPTR, OSRCPTR
C     FOR the scaling phase
      INTEGER NB1, NB2, NB3
      REAL EPS
C     Iteration vars 
      INTEGER ITER, IR, IC
      INTEGER(8) :: NZIND
      REAL ELM
C     COMM TAGS....
      INTEGER TAG_COMM_COL
      PARAMETER(TAG_COMM_COL=100)
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
C     FUNCTIONS
      EXTERNAL SMUMPS_CREATEPARTVEC,
     &     SMUMPS_NUMVOLSNDRCV, 
     &     SMUMPS_SETUPCOMMS,
     &     SMUMPS_FILLMYROWCOLINDICES,
     &     SMUMPS_INITREAL,
     &     SMUMPS_INITREALLST,
     &     SMUMPS_DOCOMMINF,
     &     SMUMPS_DOCOMM1N
      REAL SMUMPS_ERRSCALOC
      REAL SMUMPS_ERRSCA1
      INTRINSIC abs
      REAL RONE, RZERO
      PARAMETER(RONE=1.0E0,RZERO=0.0E0)
C     TMP VARS
      INTEGER(8) :: RESZR, RESZC
      INTEGER(8) :: INTSZR, INTSZC
      INTEGER MAXMN
      INTEGER I, IERROR, IOMP
      REAL ONEERRROW, ONEERRCOL, ONEERRL, ONEERRG
      REAL INFERRROW, INFERRCOL, INFERRL, INFERRG
      LOGICAL OORANGEIND
      INTEGER, PARAMETER :: K361 = 2048
!$    INTEGER    :: NOMP
!$    INTEGER    :: CHUNK, CHUNK_NZ
!$    ! Too large => pb with cache L3 ?
!$    ! INTEGER(8) :: CHUNK8 
!$    IF (NOMP_MAX.GT.0) THEN
!$     NOMP     = omp_get_max_threads()
!$     CHUNK    = max(K361/2, (N+NOMP-1) / NOMP )
!$    ! CHUNK8= (NZ_loc+int(NOMP-1,8) / int(NOMP,8) ) )
!$     CHUNK_NZ = max(K361/2, (N+NOMP_MAX-1) / NOMP_MAX )
!$    ENDIF
C
      OORANGEIND = .FALSE.
      INFERRG = -RONE
      ONEERRG = -RONE
      MAXMN = M
      IF(MAXMN < N) MAXMN = N
C     Create row partvec and col partvec
      IF(OP == 1) THEN
         IF(NUMPROCS > 1) THEN
C     Check done outside
C     IF(IWRKSZ.LT.4*MAXMN) THEN   ERROR.... 
            CALL SMUMPS_CREATEPARTVEC(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           RPARTVEC, M, N,
     &           IWRK, IWRKSZ, INUMMYR, NOMP_MAX)
            CALL SMUMPS_CREATEPARTVEC(MYID, NUMPROCS, COMM,
     &           JCN_loc, IRN_loc,  NZ_loc,
     &           CPARTVEC, N, M,
     &           IWRK, IWRKSZ, INUMMYC, NOMP_MAX)
C     Compute sndrcv sizes, store them for later use           
            CALL SMUMPS_NUMVOLSNDRCV(MYID, NUMPROCS, M, RPARTVEC,
     &           NZ_loc, IRN_loc, N, JCN_loc,
     &           IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM,ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
            CALL SMUMPS_NUMVOLSNDRCV(MYID, NUMPROCS, N, CPARTVEC,
     &           NZ_loc, JCN_loc, M, IRN_loc, 
     &           ICSNDRCVNUM,ICSNDRCVVOL,
     &           OCSNDRCVNUM,OCSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS), COMM)         
            INTSZR =  int(IRSNDRCVNUM,8) + int(ORSNDRCVNUM,8) + 
     &           int(IRSNDRCVVOL,8) + int(ORSNDRCVVOL,8) +
     &           2_8*int(NUMPROCS+1,8) + int(INUMMYR,8)
            INTSZC = int(ICSNDRCVNUM,8) + int(OCSNDRCVNUM,8) + 
     &           int(ICSNDRCVVOL,8) + int(OCSNDRCVVOL,8) +
     &           2_8*int(NUMPROCS+1,8) + int(INUMMYC,8)
            INTSZ = INTSZR + INTSZC + int(MAXMN,8) + 
     &           int(MPI_STATUS_SIZE +1,8) * int(NUMPROCS,8)
         ELSE
C     NUMPROCS IS 1
             IRSNDRCVNUM = 0
             ORSNDRCVNUM = 0
             IRSNDRCVVOL = 0
             ORSNDRCVVOL = 0
             INUMMYR = 0
             ICSNDRCVNUM  = 0
             OCSNDRCVNUM  = 0
             ICSNDRCVVOL = 0
             OCSNDRCVVOL  = 0
             INUMMYC = 0
             INTSZ = 0_8
          ENDIF
C     CALCULATE NECESSARY REAL SPACE
          RESZR = int(M,8) + int(IRSNDRCVVOL,8) + int(ORSNDRCVVOL,8)
          RESZC = int(N,8) + int(ICSNDRCVVOL,8) + int(OCSNDRCVVOL,8)
          RESZ = RESZR  + RESZC  
C     CALCULATE NECESSARY INT SPACE
C     The last maxmn is tmpwork for setup comm and fillmyrowcol
          REGISTRE(1)  = int(IRSNDRCVNUM,8)
          REGISTRE(2)  = int(ORSNDRCVNUM,8)
          REGISTRE(3)  = int(IRSNDRCVVOL,8)
          REGISTRE(4)  = int(ORSNDRCVVOL,8)
          REGISTRE(5)  = int(ICSNDRCVNUM,8)
          REGISTRE(6)  = int(OCSNDRCVNUM,8)
          REGISTRE(7)  = int(ICSNDRCVVOL,8)
          REGISTRE(8)  = int(OCSNDRCVVOL,8)
          REGISTRE(9)  = int(INUMMYR,8)
          REGISTRE(10) = int(INUMMYC,8)
          REGISTRE(11) = INTSZ
          REGISTRE(12) = RESZ
       ELSE
C     else of op=1. That is op=2 now.
C     restore the numbers
          IRSNDRCVNUM = int(REGISTRE(1))
          ORSNDRCVNUM = int(REGISTRE(2))
          IRSNDRCVVOL = int(REGISTRE(3))
          ORSNDRCVVOL = int(REGISTRE(4))
          ICSNDRCVNUM = int(REGISTRE(5))
          OCSNDRCVNUM = int(REGISTRE(6))
          ICSNDRCVVOL = int(REGISTRE(7))
          OCSNDRCVVOL = int(REGISTRE(8))
          INUMMYR     = int(REGISTRE(9))
          INUMMYC     = int(REGISTRE(10))
          IF(NUMPROCS > 1) THEN
C     Check done outsize
C             IF(INTSZ < REGISTRE(11)) THEN ERROR
C             IF(RESZ < REGISTRE(12)) THEN ERROR
C     Fill up myrows and my colsX
             CALL SMUMPS_FILLMYROWCOLINDICES(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            RPARTVEC, CPARTVEC, M, N,
     &            IWRK(1_8), INUMMYR,
     &            IWRK(1_8+int(INUMMYR,8)), INUMMYC,     
     &            IWRK(1_8+int(INUMMYR,8)+int(INUMMYC,8)),
     &            IWRKSZ-int(INUMMYR,8)-int(INUMMYC,8), NOMP_MAX )
             IMYRPTR = 1_8
             IMYCPTR = IMYRPTR + int(INUMMYR,8)
C     Set up comm and run.
C     set pointers in iwrk (4 parts)
C     
C     ROWS    [---------------------------------------------]
             IRNGHBPRCS = IMYCPTR    + int(INUMMYC    ,8)
             IRSNDRCVIA = IRNGHBPRCS + int(IRSNDRCVNUM,8)
             IRSNDRCVJA = IRSNDRCVIA + int(NUMPROCS+1 ,8)
             ORNGHBPRCS = IRSNDRCVJA + int(IRSNDRCVVOL,8)
             ORSNDRCVIA = ORNGHBPRCS + int(ORSNDRCVNUM,8)
             ORSNDRCVJA = ORSNDRCVIA + int(NUMPROCS+1 ,8)
C     COLS    [---------------------------------------------]
             ICNGHBPRCS = ORSNDRCVJA + int(ORSNDRCVVOL,8)
             ICSNDRCVIA = ICNGHBPRCS + int(ICSNDRCVNUM,8)
             ICSNDRCVJA = ICSNDRCVIA + int(NUMPROCS+1 ,8)
             OCNGHBPRCS = ICSNDRCVJA + int(ICSNDRCVVOL,8)
             OCSNDRCVIA = OCNGHBPRCS + int(OCSNDRCVNUM,8)
             OCSNDRCVJA = OCSNDRCVIA + int(NUMPROCS+1 ,8)
C     
C     MPI     [-----------------]
             REQUESTS = OCSNDRCVJA + int(OCSNDRCVVOL,8)
             ISTATUS = REQUESTS + int(NUMPROCS,8)
C     
C     TMPWRK  [-----------------]
             TMPWORK = ISTATUS + int(MPI_STATUS_SIZE,8)*int(NUMPROCS,8)
             CALL SMUMPS_SETUPCOMMS(MYID, NUMPROCS, M, RPARTVEC,
     &            NZ_loc, IRN_loc,N, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL SMUMPS_SETUPCOMMS(MYID, NUMPROCS, N, CPARTVEC,
     &            NZ_loc, JCN_loc, M, IRN_loc,
     &            ICSNDRCVNUM, ICSNDRCVVOL, 
     &            IWRK(ICNGHBPRCS),
     &            IWRK(ICSNDRCVIA),
     &            IWRK(ICSNDRCVJA),
     &            OCSNDRCVNUM, OCSNDRCVVOL, 
     &            IWRK(OCNGHBPRCS),IWRK(OCSNDRCVIA),IWRK(OCSNDRCVJA),
     &            CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS),  IWRK(REQUESTS),
     &            TAG_COMM_COL, COMM)
             CALL SMUMPS_INITREAL(ROWSCA, M, RZERO, NOMP_MAX)
             CALL SMUMPS_INITREAL(COLSCA, N, RZERO, NOMP_MAX)
             CALL SMUMPS_INITREALLST(ROWSCA, M, 
     &            IWRK(IMYRPTR),INUMMYR, RONE, NOMP_MAX)
             CALL SMUMPS_INITREALLST(COLSCA, N, 
     &            IWRK(IMYCPTR),INUMMYC, RONE, NOMP_MAX)   
          ELSE
             CALL SMUMPS_INITREAL(ROWSCA, M, RONE, NOMP_MAX)
             CALL SMUMPS_INITREAL(COLSCA, N, RONE, NOMP_MAX)            
          ENDIF
          ITDRPTR = 1_8
          ITDCPTR = ITDRPTR + int(M,8)
C     
          ISRRPTR = ITDCPTR + int(N,8)
          OSRRPTR = ISRRPTR + int(IRSNDRCVVOL,8)
C     
          ISRCPTR = OSRRPTR + int(ORSNDRCVVOL,8)
          OSRCPTR = ISRCPTR + int(ICSNDRCVVOL,8)
C     To avoid bound check errors...
          IF(NUMPROCS == 1)THEN
             OSRCPTR = OSRCPTR - 1_8
             ISRCPTR = ISRCPTR - 1_8
             OSRRPTR = OSRRPTR - 1_8
             ISRRPTR = ISRRPTR - 1_8
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1_8
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1_8
             IF(ICSNDRCVVOL == 0) ISRCPTR = ISRCPTR - 1_8
             IF(OCSNDRCVVOL == 0) OSRCPTR = OSRCPTR - 1_8
          ENDIF
          ITER = 1
          DO WHILE (ITER.LE.NB1+NB2+NB3)
C{
C           -------------------------
C           CLEAR temporary Dr and Dc
C           -------------------------
            IF (NOMP_MAX.GT.1 .AND. 
     &        (ITER.GT.NB1 .AND.ITER.LE.NB1+NB2)
     &      ) THEN
C{
!$OMP         PARALLEL
!$OMP&         PRIVATE(IOMP)
!$OMP&         NUM_THREADS(NOMP_MAX)
               IOMP  = 1
!$             IOMP = OMP_GET_THREAD_NUM() + 1
               IF(NUMPROCS > 1) THEN
                 CALL SMUMPS_ZEROOUT(WRKR_TH(1,IOMP),N,
     &                 IWRK(IMYRPTR),INUMMYR, 0)
                 CALL SMUMPS_ZEROOUT(WRKC_TH(1,IOMP),N,
     &                 IWRK(IMYCPTR),INUMMYC, 0)
               ELSE
                 CALL SMUMPS_INITREAL(WRKR_TH(1,IOMP),N, RZERO, 
     &                                0)
                 CALL SMUMPS_INITREAL(WRKC_TH(1,IOMP),N, RZERO, 
     &                                0)
               ENDIF
!$OMP         END PARALLEL
C}
            ELSE
C{
               IF(NUMPROCS > 1) THEN
                  CALL SMUMPS_ZEROOUT(WRKRC(ITDRPTR),M,
     &                 IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
                  CALL SMUMPS_ZEROOUT(WRKRC(ITDCPTR),N,
     &                 IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
               ELSE
                  CALL SMUMPS_INITREAL(WRKRC(ITDRPTR),M, RZERO, 
     &                                NOMP_MAX)
                  CALL SMUMPS_INITREAL(WRKRC(ITDCPTR),N, RZERO, 
     &                                NOMP_MAX)
               ENDIF
C}
             ENDIF
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
C     ------------------
C     INF-NORM ITERATION
C     ------------------
               IF (NOMP_MAX.LE.0) THEN
                IF((ITER.EQ.1).OR.(OORANGEIND)) THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                         IF(WRKRC(IR)<ELM) WRKRC(IR)= ELM
                         IF(WRKRC(ITDCPTR-1_8+int(IC,8))<ELM) 
     &                      WRKRC(ITDCPTR-1_8+int(IC,8))= ELM
                      ELSE
                         OORANGEIND = .TRUE.
                      ENDIF
                   ENDDO
                ELSEIF(.NOT.OORANGEIND) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                      IF(WRKRC(IR)<ELM) WRKRC(IR)= ELM
                      IF(WRKRC(ITDCPTR-1_8+int(IC,8))<ELM) 
     &                   WRKRC(ITDCPTR-1_8+int(IC,8))= ELM
                   ENDDO
                ENDIF
               ELSE
                IF((ITER.EQ.1).OR.(OORANGEIND)) THEN
!$OMP   PARALLEL DO PRIVATE(NZIND,IR,IC,ELM) 
!$OMP&  SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&  IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
!$OMP&  REDUCTION(.OR.:OORANGEIND)
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
!$OMP                    ATOMIC UPDATE
                         WRKRC(IR)= max (ELM, WRKRC(IR))
!$OMP                    END ATOMIC
C                         IF(WRKRC(ITDRPTR-1_8+int(IR,8))<ELM) THEN
C                            WRKRC(ITDRPTR-1_8+int(IR,8))= ELM
C                         ENDIF
C                        IF(WRKRC(ITDCPTR-1_8+int(IC,8))<ELM) THEN
C                           WRKRC(ITDCPTR-1_8+int(IC,8))= ELM
C                         ENDIF
C                        ITDCPTR+IC might be > int4
!$OMP                    ATOMIC UPDATE
                         WRKRC(ITDCPTR-1_8+int(IC,8)) = 
     &                   max (ELM,WRKRC(ITDCPTR-1_8+int(IC,8)))
!$OMP                    END ATOMIC
                      ELSE
                         OORANGEIND = .TRUE.
                      ENDIF
                   ENDDO
!$OMP   END PARALLEL DO
                ELSEIF(.NOT.OORANGEIND) THEN
!$OMP   PARALLEL DO PRIVATE(NZIND,IR,IC,ELM) 
!$OMP&  SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&  IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
C                      IF(WRKRC(ITDRPTR-1_8+int(IR,8))<ELM) THEN
C                         WRKRC(ITDRPTR-1_8+int(IR,8))= ELM
C                      ENDIF
!$OMP                 ATOMIC UPDATE
                      WRKRC(IR)= max (ELM, WRKRC(IR))
!$OMP                 END ATOMIC
C                      IF(WRKRC(ITDCPTR-1_8+int(IC,8))<ELM) THEN
C                         WRKRC(ITDCPTR-1_8+int(IC,8))= ELM
C                      ENDIF
!$OMP                 ATOMIC UPDATE
                      WRKRC(ITDCPTR-1_8+int(IC,8)) = 
     &                   max (ELM,WRKRC(ITDCPTR-1_8+int(IC,8)))
!$OMP                 END ATOMIC
                   ENDDO
!$OMP   END PARALLEL DO
                ENDIF
C}
               ENDIF
                IF(NUMPROCS > 1) THEN
                   CALL SMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                  ICSNDRCVNUM,IWRK(ICNGHBPRCS),
     &                  ICSNDRCVVOL,IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                  WRKRC(ISRCPTR),
     &                  OCSNDRCVNUM,IWRK(OCNGHBPRCS),
     &                  OCSNDRCVVOL,IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                  WRKRC( OSRCPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
C     
                  CALL SMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = SMUMPS_ERRSCALOC(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
C     find error for the cols
                     INFERRCOL = SMUMPS_ERRSCALOC(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
C     get max of those two errors
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL ) THEN
                        INFERRL = INFERRROW                   
                     ENDIF
C     
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_REAL,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL SMUMPS_UPDATESCALE(COLSCA,  
     &                       WRKRC(ITDCPTR),N,
     &                       IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
                        CALL SMUMPS_UPDATESCALE(ROWSCA,  
     &                       WRKRC(ITDRPTR),M,
     &                       IWRK(IMYRPTR),INUMMYR, NOMP_MAX)         
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF                  
               ELSE
C     SINGLE PROCESSOR CASE: INF-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = SMUMPS_ERRSCA1(ROWSCA, 
     &                    WRKRC(ITDRPTR), M, NOMP_MAX)
C     find error for the cols
                     INFERRCOL = SMUMPS_ERRSCA1(COLSCA,  
     &                    WRKRC(ITDCPTR), N, NOMP_MAX)
C     get max of those two errors
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL) THEN
                        INFERRL = INFERRROW                    
                     ENDIF                     
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL SMUMPS_UPSCALE1(COLSCA, WRKRC(ITDCPTR), N, 
     &                                       NOMP_MAX)
                        CALL SMUMPS_UPSCALE1(ROWSCA, WRKRC(ITDRPTR), M, 
     &                                       NOMP_MAX)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
               ENDIF
            ELSE
C     ----------------------------------------
C     WE HAVE ITER.GT.NB1 AND ITER.LE.NB1+NB2. 
C     ONE-NORM ITERATION
C     ----------------------------------------
             IF (NOMP_MAX.LE.1) THEN
               IF((ITER .EQ.1).OR.(OORANGEIND))THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
C                       WRKRC(ITDRPTR-1_8+int(IR,8)) =
C     &                 WRKRC(ITDRPTR-1_8+int(IR,8)) + ELM
                        WRKRC(IR) =  WRKRC(IR) + ELM
                        WRKRC(ITDCPTR-1_8+int(IC,8)) =
     &                  WRKRC(ITDCPTR-1_8+int(IC,8)) + ELM
                     ELSE
                        OORANGEIND = .TRUE.
                     ENDIF
                  ENDDO
               ELSEIF(.NOT.OORANGEIND) THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
C                    WRKRC(ITDRPTR-1_8+int(IR,8)) =
C     &              WRKRC(ITDRPTR-1_8+int(IR,8)) + ELM
                     WRKRC(IR) =  WRKRC(IR) + ELM
                     WRKRC(ITDCPTR-1_8+int(IC,8)) =
     &               WRKRC(ITDCPTR-1_8+int(IC,8)) + ELM
                  ENDDO
               ENDIF
C}
             ELSE  ! NOMP_MAX>1
               IF((ITER .EQ.1).OR.(OORANGEIND))THEN
!$OMP   PARALLEL PRIVATE(IOMP)
!$OMP&  NUM_THREADS(NOMP_MAX)
!$OMP&    IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                    IOMP  = 1
!$                  IOMP = OMP_GET_THREAD_NUM() + 1
!$OMP     DO PRIVATE(NZIND,IR,IC,ELM)
!$OMP&    SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&    REDUCTION(.OR.:OORANGEIND)
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                        WRKR_TH(IR,IOMP) = WRKR_TH(IR,IOMP) + ELM
                        IF (IR.NE.IC)
     &                  WRKC_TH(IC,IOMP) = WRKC_TH(IC,IOMP) + ELM
                     ELSE
                        OORANGEIND = .TRUE.
                     ENDIF
                  ENDDO
!$OMP     END DO
!$OMP   END PARALLEL
               ELSEIF(.NOT.OORANGEIND) THEN
!$OMP   PARALLEL PRIVATE(IOMP)
!$OMP&  NUM_THREADS(NOMP_MAX)
!$OMP&    IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                    IOMP  = 1
!$                  IOMP = OMP_GET_THREAD_NUM() + 1
!$OMP     DO PRIVATE(NZIND,IR,IC,ELM)
!$OMP&    SCHEDULE(STATIC,CHUNK_NZ)
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                     WRKR_TH(IR,IOMP) = WRKR_TH(IR,IOMP) + ELM
                     IF (IR.NE.IC)
     &               WRKC_TH(IC,IOMP) = WRKC_TH(IC,IOMP) + ELM
                  ENDDO
!$OMP     END DO
!$OMP   END PARALLEL
               ENDIF
C
C              For all i on MYID:
C              Build WRKRC(i) = Sum       (WRKR_TH(i,IOMP)
C                              IOMP \in [1:NOMP_MAX]
               IF(NUMPROCS > 1) THEN
                 CALL SMUMPS_REDUCE_WRK_MPI (WRKRC, N, WRKR_TH, 
     &                 NOMP_MAX, 
     &                 IWRK(IMYRPTR),INUMMYR)
                 CALL SMUMPS_REDUCE_WRK_MPI (WRKRC(ITDCPTR), 
     &                 N, WRKC_TH, NOMP_MAX, 
     &                 IWRK(IMYCPTR),INUMMYC)
               ELSE
                 CALL SMUMPS_REDUCE_WRK (WRKRC, N, WRKR_TH, NOMP_MAX)
                 CALL SMUMPS_REDUCE_WRK (WRKRC(ITDCPTR), 
     &                 N, WRKC_TH, NOMP_MAX)
               ENDIF
C}
             ENDIF
               IF(NUMPROCS > 1) THEN                 
                  CALL SMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                 ICSNDRCVNUM, IWRK(ICNGHBPRCS),
     &                 ICSNDRCVVOL, IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                 WRKRC(ISRCPTR),
     &                 OCSNDRCVNUM, IWRK(OCNGHBPRCS),
     &                 OCSNDRCVVOL, IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                 WRKRC( OSRCPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
C     
                  CALL SMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = SMUMPS_ERRSCALOC(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
C     find error for the cols
                     ONEERRCOL = SMUMPS_ERRSCALOC(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
C     get max of those two errors
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL ) THEN
                        ONEERRL = ONEERRROW                   
                     ENDIF
C     
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_REAL,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(ONEERRG.LE.EPS) THEN
                        CALL SMUMPS_UPDATESCALE(COLSCA,
     &                       WRKRC(ITDCPTR),N,
     &                       IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
                        CALL SMUMPS_UPDATESCALE(ROWSCA,
     &                       WRKRC(ITDRPTR),M,
     &                       IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF                            
               ELSE
C     SINGLE-PROCESSOR CASE: ONE-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = SMUMPS_ERRSCA1(ROWSCA, 
     &                    WRKRC(ITDRPTR), M, NOMP_MAX)
C     find error for the cols
                     ONEERRCOL = SMUMPS_ERRSCA1(COLSCA,  
     &                    WRKRC(ITDCPTR), N, NOMP_MAX)
C     get max of those two errors
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL) THEN
                        ONEERRL = ONEERRROW                    
                     ENDIF                     
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                       CALL SMUMPS_UPSCALE1(COLSCA, WRKRC(ITDCPTR), N, 
     &                                      NOMP_MAX)
                       CALL SMUMPS_UPSCALE1(ROWSCA, WRKRC(ITDRPTR), M, 
     &                                      NOMP_MAX)
                       ITER = NB1+NB2+1                        
                       CYCLE
                     ENDIF
                  ENDIF                  
               ENDIF 
            ENDIF
            IF(NUMPROCS > 1) THEN               
               CALL SMUMPS_UPDATESCALE(COLSCA,  WRKRC(ITDCPTR), N,
     &              IWRK(IMYCPTR),INUMMYC, NOMP_MAX)
               CALL SMUMPS_UPDATESCALE(ROWSCA,  WRKRC(ITDRPTR), M,
     &              IWRK(IMYRPTR),INUMMYR, NOMP_MAX)               
C   
            ELSE
C     SINGLE PROCESSOR CASE: Conv check and update of sca arrays
               CALL SMUMPS_UPSCALE1(COLSCA, WRKRC(ITDCPTR), N, 
     &                              NOMP_MAX)
               CALL SMUMPS_UPSCALE1(ROWSCA, WRKRC(ITDRPTR), M, 
     &                              NOMP_MAX)
            ENDIF
            ITER = ITER + 1
C}
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
C{
            CALL MPI_REDUCE(ROWSCA, WRKRC(1), M, MPI_REAL,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
               DO I=1, M
                  ROWSCA(I) = WRKRC(I)
               ENDDO
            ENDIF
C     Scaling factors are printed
C     WRITE (6,*) MYID, 'ROWSCA=',ROWSCA
C     WRITE (6,*) MYID, 'COLSCA=',COLSCA
C     CALL FLUSH(6)
c     REduce the whole scaling factors to processor 0 of COMM
            CALL MPI_REDUCE(COLSCA, WRKRC(1_8+int(M,8)), N, 
     &           MPI_REAL,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
C{
              IF (NOMP_MAX.LE.0) THEN
               DO I=1, N
                  COLSCA(I) = WRKRC(int(I,8)+int(M,8))
               ENDDO
              ELSE
!$OMP       PARALLEL DO PRIVATE(I)
!$OMP&      SCHEDULE(STATIC,CHUNK)
!$OMP&      IF ( N > K361 .AND. NOMP .GT. 1)
               DO I=1, N
                  COLSCA(I) = WRKRC(int(I,8)+int(M,8))
               ENDDO
!$OMP       END PARALLEL DO
              ENDIF
C}
            ENDIF         
C}
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SIMSCALEABSUNS
C
C 
C     SEPARATOR: Another function begins
C
C 
      SUBROUTINE SMUMPS_SIMSCALEABSSYM(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     N, NUMPROCS, MYID, COMM,
     &     PARTVEC, 
     &     RSNDRCVSZ, 
     &     REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     SCA, WRKRC, ISZWRKRC,
     &     WRKR_TH, LWRKR_TH, NOMP_MAX,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
C----------------------------------------------------------------------
C Input parameters:
C     N: size of matrix (sym matrix, square).
C     NUMPROCS, MYID, COMM: guess what are those
C     PARTVEC:  row/col partvec to be filled when OP=1
C     RSNDRCVSZ:send recv sizes for row/col operations. 
C               to be filled when OP=1
C     REGISTRE: to store some pointers (size etc). Its size is 12,
C               but we do not use all in this routine.
C     IWRK: working space. when OP=1 IWRKSZ.GE.2*MAXMN
C           when OP=2 INTSZ portion is used. Donc, IWRKSZ>INTSZ 
C           when OP=2
C      IWRKSZ: size
C      INTSZ: to be computed when OP=1, necessary integer space to run 
C             scaling algo when OP=2
C      RESZ:  to be computed when OP=1, necessary real space to run 
C             scaling algo when OP=2
C      OP: 
C          =1 estimation of memory and construction of partvecs
C           writes into PARTVEC,RSNDRCVSZ,REGISTRE
C           does not access WRKRC, uses IWRK as workspace
C           computes INTSZ and RESZ.
C          =2 Compute scalings 
C           restores pointers from REGISTRE, 
C           stores communication structure in IWRK (from the start). 
C
C      SCA: space for row/col scaling factor; has size M
C      WRKRC: real working space. when OP=1, is not accessed. Donc, it
C             can be declared to be of size 1 at OP=1 call.
C      ISZWRKRC: size
C      SYM: is matrix symmetric
C      NB1, NB2, NB3: algo runs 
C                     NB1 iters of inf-norm (default  1/1), 
C                     NB2 iters of 1-norm   (default  3/10),
C                     NB3 iters of inf-norm (default  3/10).
C          in succession.
C      EPS: tolerance for concergence. 
C           IF EPS < 0.R0 then does not test convergence.
C           See comments for the uns case above.      
C     ONENORMERR : error in one norm scaling (see comments for the 
C                  uns case above), 
C     INFNORMERR : error in inf norm scaling (see comments for the 
C                  uns case above).
C---------------------------------------------------------------------
C On input:
C      OP=1==>Requirements
C             IWRKSZ.GE.2*MAXMN   XXXX compare with uns variant.
C             PARTVEC  of size N
C             SNDRCVSZ of size 2*NUMPROCS
C             REGISTRE  of size 12
C             
C      OP=2==>Requirements
C             INTSZ .GE. REGISTRE(11)
C             RESZ  .GE. REGISTRE(12) 
C---------------------------------------------------------------------
C On output:
C     SCA
C            at processor 0 of COMM: complete factors.
C            at other processors   : only the SCA(i) and SCA(j)
C            for which there is a nonzero a_ij.
C     ONENORMERR : error in one norm scaling 
C                = -1.0 if iter2=0.
C     INFNORMERR : error in inf norm scaling 
C                = inf norm error at iter3 if iter3 > 0
C                = inf norm error at iter1 if iter1 > 0, iter3=0
C                = -1.0 if iter1=iter3=0
C ---------------------------------------------------------------------
C NOTE: some variables are named in such a way that they correspond
C       to the row variables in unsym case. They are used for both 
C       row and col communications.
C ---------------------------------------------------------------------
C References:
C     The scaling algorithms are based on those discussed in
C     [1] D. Ruiz, "A scaling algorithm to equilibrate both rows and 
C         columns norms in matrices", Tech. Rep. Rutherford 
C         Appleton Laboratory, Oxon, UK and ENSEEIHT-IRIT, 
C         Toulouse, France, RAL-TR-2001-034 and RT/APO/01/4, 2001.
C     [2] D. Ruiz and B. Ucar, "A symmetry preserving algorithm for
C         matrix scaling", in preparation as of Jan'08.
C
C     The parallelization approach is based on discussion in
C     [3] P. R. Amestoy, I. S. Duff, D. Ruiz, and B. Ucar, "A parallel
C         matrix scaling algorithm", accepted for publication, 
C         In proceedings of VECPAR'08-International Meeting-High 
C         Performance Computing for Computational Science, Jan'08.
C     and was supported by ANR-SOLSTICE project (ANR-06-CIS6-010)
C ---------------------------------------------------------------------
!$    USE OMP_LIB
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) :: NZ_loc 
      INTEGER N, OP
      INTEGER(8) :: IWRKSZ, LWRKR_TH
      INTEGER NUMPROCS, MYID, COMM, NOMP_MAX
      INTEGER(8) :: INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      REAL A_loc(NZ_loc)
      INTEGER PARTVEC(N), RSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER(8) :: REGISTRE(12)
      REAL SCA(N)
      INTEGER(8) :: ISZWRKRC
      REAL WRKRC(ISZWRKRC), 
     &     WRKR_TH(LWRKR_TH, max(NOMP_MAX,1))
C     LOCALS
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER  INUMMYR
C IMPORTANT POINTERS
      INTEGER(8) :: IMYRPTR,IMYCPTR 
      INTEGER(8) :: IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER(8) :: ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER(8) :: ISTATUS, REQUESTS, TMPWORK
      INTEGER(8) :: ITDRPTR, ISRRPTR, OSRRPTR
      REAL ONENORMERR,INFNORMERR
C     FOR the scaling phase  
      INTEGER NB1, NB2, NB3
      REAL EPS
C     Iteration vars 
      INTEGER ITER, IR, IC
      INTEGER(8) :: NZIND
      REAL ELM
C     COMM TAGS....
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
C     FUNCTIONS
      EXTERNAL SMUMPS_CREATEPARTVECSYM,
     &     SMUMPS_NUMVOLSNDRCVSYM, 
     &     SMUMPS_SETUPCOMMSSYM,
     &     SMUMPS_FILLMYROWCOLINDICESSYM,
     &     SMUMPS_DOCOMMINF,
     &     SMUMPS_DOCOMM1N,
     &     SMUMPS_INITREAL,
     &     SMUMPS_INITREALLST
      REAL SMUMPS_ERRSCALOC
      REAL SMUMPS_ERRSCA1
      INTRINSIC abs
      REAL RONE, RZERO
      PARAMETER(RONE=1.0E0,RZERO=0.0E0)
C     TMP VARS
      INTEGER(8) :: INTSZR
      INTEGER MAXMN
      INTEGER I, IERROR
      REAL ONEERRL, ONEERRG
      REAL INFERRL, INFERRG
      LOGICAL OORANGEIND
      INTEGER, PARAMETER :: K361 = 2048
      INTEGER            :: IOMP
!$    INTEGER    :: NOMP
!$    INTEGER    :: CHUNK, CHUNK_NZ
!$    ! Too large => pb with cache L3 ?
!$    ! INTEGER(8) :: CHUNK8 
!$    ! CHUNK8= max(int(K361/2,8), 
!$    !  &   (NZ_loc+int(NOMP-1,8) / int(NOMP,8) ) )
!$    ! CHUNK8   = min(CHUNK8, huge(CHUNK)-1_8)
!$    NOMP  = omp_get_max_threads()
!$    CHUNK= max(K361/2, (N+NOMP-1) / NOMP )
!$    IF (NOMP_MAX.GT.0) THEN
!$      CHUNK_NZ = max(K361/2, (N+NOMP_MAX-1) / NOMP_MAX )
!$    ENDIF
C
      OORANGEIND = .FALSE.
      INFERRG = -RONE
      ONEERRG = -RONE
      MAXMN = N
      IF(OP == 1) THEN
C{
         IF(NUMPROCS > 1) THEN
C     Check done outside
C     IF(IWRKSZ.LT.2*MAXMN) THEN   ERROR.... 
            CALL SMUMPS_CREATEPARTVECSYM(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           PARTVEC, N,
     &           IWRK, IWRKSZ, INUMMYR )
C
C     Check done outside
            CALL SMUMPS_NUMVOLSNDRCVSYM(MYID, NUMPROCS, N, PARTVEC,
     &           NZ_loc, IRN_loc, JCN_loc, IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM, ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
C     
C     
            INTSZR =  int(IRSNDRCVNUM,8) + int(ORSNDRCVNUM,8) + 
     &           int(IRSNDRCVVOL,8) + int(ORSNDRCVVOL,8) +
     &           2_8*int(NUMPROCS+1,8) + int(INUMMYR,8)
            INTSZ = INTSZR + int(N,8) + 
     &           int(MPI_STATUS_SIZE +1,8) * int(NUMPROCS,8)
         ELSE
C     NUMPROCS IS 1
            IRSNDRCVNUM = 0
            ORSNDRCVNUM = 0
            IRSNDRCVVOL = 0 
            ORSNDRCVVOL = 0
            INUMMYR = 0
            INTSZ = 0_8
         ENDIF
C     CALCULATE NECESSARY REAL SPACE
         RESZ = int(N,8) + int(IRSNDRCVVOL,8) + int(ORSNDRCVVOL,8)
         REGISTRE(1)  = int(IRSNDRCVNUM,8)
         REGISTRE(2)  = int(ORSNDRCVNUM,8)
         REGISTRE(3)  = int(IRSNDRCVVOL,8)
         REGISTRE(4)  = int(ORSNDRCVVOL,8)
         REGISTRE(9)  = int(INUMMYR,8)
         REGISTRE(11) = INTSZ
         REGISTRE(12) = RESZ
C}
      ELSE
C{
C     else of op=1. That is op=2 now.
C     restore the numbers
         IRSNDRCVNUM = int(REGISTRE(1))
         ORSNDRCVNUM = int(REGISTRE(2))
         IRSNDRCVVOL = int(REGISTRE(3))
         ORSNDRCVVOL = int(REGISTRE(4))
         INUMMYR     = int(REGISTRE(9))
          IF(NUMPROCS > 1) THEN
C     Check done outsize
C             IF(INTSZ < REGISTRE(11)) THEN ERROR
C             IF(RESZ < REGISTRE(12)) THEN ERROR
C     Fill up myrows and my colsX
             CALL SMUMPS_FILLMYROWCOLINDICESSYM(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            PARTVEC, N,
     &            IWRK(1), INUMMYR,
     &            IWRK(1+INUMMYR), IWRKSZ-int(INUMMYR,8), NOMP_MAX)
             IMYRPTR = 1_8
             IMYCPTR = IMYRPTR + int(INUMMYR,8)
C     Set up comm and run.
C     set pointers in iwrk (3 parts)
C     
C     ROWS    [---------------------------------------------]
             IRNGHBPRCS = IMYCPTR 
             IRSNDRCVIA = IRNGHBPRCS + int(IRSNDRCVNUM,8)
             IRSNDRCVJA = IRSNDRCVIA + int(NUMPROCS+1,8)
             ORNGHBPRCS = IRSNDRCVJA + int(IRSNDRCVVOL,8)
             ORSNDRCVIA = ORNGHBPRCS + int(ORSNDRCVNUM,8)
             ORSNDRCVJA = ORSNDRCVIA + int(NUMPROCS + 1,8)
C     MPI     [-----------------]
             REQUESTS = ORSNDRCVJA + int(ORSNDRCVVOL,8)
             ISTATUS = REQUESTS + int(NUMPROCS,8)
C     TMPWRK  [-----------------]
             TMPWORK = ISTATUS + int(MPI_STATUS_SIZE,8)*int(NUMPROCS,8)
             CALL SMUMPS_SETUPCOMMSSYM(MYID, NUMPROCS, N, PARTVEC,
     &            NZ_loc, IRN_loc, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL SMUMPS_INITREAL(SCA, N, RZERO, NOMP_MAX)
             CALL SMUMPS_INITREALLST(SCA, N, 
     &            IWRK(IMYRPTR),INUMMYR, RONE, NOMP_MAX)
          ELSE
             CALL SMUMPS_INITREAL(SCA, N, RONE, NOMP_MAX)
          ENDIF
          ITDRPTR = 1_8
          ISRRPTR = ITDRPTR + int(N,8)
          OSRRPTR = ISRRPTR + int(IRSNDRCVVOL,8)
C     
C     To avoid bound check errors...
          IF(NUMPROCS == 1)THEN
             OSRRPTR = OSRRPTR - 1
             ISRRPTR = ISRRPTR - 1
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1
          ENDIF
C     computation starts
          ITER = 1
          DO WHILE(ITER.LE.NB1+NB2+NB3)
C{
C            -------------------------
C            CLEAR temporary Dr and Dc 
C            -------------------------
             IF (NOMP_MAX.GT.1 .AND. 
     &         (ITER.GT.NB1 .AND.ITER.LE.NB1+NB2)
     &       ) THEN
C             if one norm iteration and multithreading activated
C             WRKR_TH need be initialized and 
C             WRKRC will be set by reduction of WRKR_TH
!$OMP         PARALLEL
!$OMP&         PRIVATE(IOMP)
!$OMP&         NUM_THREADS(NOMP_MAX)
               IOMP  = 1
!$             IOMP = OMP_GET_THREAD_NUM() + 1
               IF(NUMPROCS > 1) THEN
                 CALL SMUMPS_ZEROOUT(WRKR_TH(1,IOMP),N,
     &                 IWRK(IMYRPTR),INUMMYR, 0)
               ELSE
                 CALL SMUMPS_INITREAL(WRKR_TH(1,IOMP),N, RZERO, 
     &                                0)
               ENDIF
!$OMP         END PARALLEL
             ELSE
              IF(NUMPROCS > 1) THEN
CFIXME         Size N should be adjusted to effective size
                 CALL SMUMPS_ZEROOUT(WRKRC(ITDRPTR),N,
     &                IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
              ELSE
                 CALL SMUMPS_INITREAL(WRKRC(ITDRPTR),N, RZERO,
     &                                NOMP_MAX)
              ENDIF
             ENDIF
C
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
C            ------------------
C{           INF-NORM ITERATION
C            ------------------
               IF (NOMP_MAX.LE.0) THEN
                 IF((ITER .EQ.1).OR.(OORANGEIND))THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                         IF (WRKRC(IR)<ELM) WRKRC(IR)= ELM
                         IF (WRKRC(IC)<ELM) WRKRC(IC)= ELM
                      ELSE
                         OORANGEIND = .TRUE.
                      ENDIF
                   ENDDO
                 ELSEIF(.NOT.OORANGEIND) THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                      IF (WRKRC(IR)<ELM) WRKRC(IR)= ELM
                      IF (WRKRC(IC)<ELM) WRKRC(IC)= ELM
                   ENDDO
                 ENDIF                      
               ELSE
                IF((ITER.EQ.1).OR.(OORANGEIND)) THEN
!$OMP   PARALLEL DO PRIVATE(NZIND,IR,IC,ELM) 
!$OMP&  SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&  IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
!$OMP&  REDUCTION(.OR.:OORANGEIND)
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
!$OMP                    ATOMIC UPDATE
                         WRKRC(IR)= max (ELM, WRKRC(IR))
!$OMP                    END ATOMIC
C                         IF(WRKRC(ITDRPTR-1_8+int(IR,8))<ELM) THEN
C                            WRKRC(ITDRPTR-1_8+int(IR,8))= ELM
C                         ENDIF
!$OMP                    ATOMIC UPDATE
                         WRKRC(IC) = max (ELM,WRKRC(IC))
!$OMP                    END ATOMIC
C                         IF(WRKRC(ITDRPTR-1_8+int(IC,8))<ELM) THEN
C                            WRKRC(ITDRPTR-1_8+int(IC,8))= ELM
C                         ENDIF
                      ELSE
                         OORANGEIND = .TRUE.
                      ENDIF
                   ENDDO
!$OMP   END PARALLEL DO
                ELSEIF(.NOT.OORANGEIND) THEN
!$OMP   PARALLEL DO PRIVATE(NZIND,IR,IC,ELM) 
!$OMP&  SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&  IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
!$OMP                 ATOMIC UPDATE
                      WRKRC(IR)= max (ELM, WRKRC(IR))
!$OMP                 END ATOMIC
C                      IF(WRKRC(ITDRPTR-1_8+int(IR,8))<ELM) THEN
C                         WRKRC(ITDRPTR-1_8+int(IR,8))= ELM
C                      ENDIF
!$OMP                 ATOMIC UPDATE
                      WRKRC(IC) = max (ELM, WRKRC(IC))
!$OMP                 END ATOMIC
C                      IF(WRKRC(ITDRPTR-1_8+int(IC,8))<ELM) THEN
C                         WRKRC(ITDRPTR-1_8+int(IC,8))= ELM
C                      ENDIF
                   ENDDO
!$OMP   END PARALLEL DO
                ENDIF                      
C}
               ENDIF
               IF(NUMPROCS > 1) THEN
C{
                  CALL SMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
CCCC FIXME #if defined(dev_version)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = SMUMPS_ERRSCALOC(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR, NOMP_MAX) 
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_REAL,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL SMUMPS_UPDATESCALE(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
CCC #endif
C}
               ELSE
C{             SINGLE PROCESSOR CASE: INF-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = SMUMPS_ERRSCA1(SCA, 
     &                    WRKRC(ITDRPTR), N, NOMP_MAX)
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL SMUMPS_UPSCALE1(SCA, WRKRC(ITDRPTR), N,
     &                                       NOMP_MAX)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
C}
               ENDIF
C}
            ELSE
C           ----------------------------------------
C{          WE HAVE ITER.GT.NB1 AND ITER.LE.NB1+NB2. 
C           ONE-NORM ITERATION
C           ----------------------------------------
             IF (NOMP_MAX.LE.1) THEN
               IF((ITER.EQ.1).OR.(OORANGEIND))THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                        WRKRC(IR) = WRKRC(IR) + ELM
                        IF(IR.NE.IC) THEN
                           WRKRC(IC) = WRKRC(IC) + ELM
                        ENDIF
                     ELSE
                        OORANGEIND = .TRUE.
                     ENDIF
                  ENDDO
               ELSEIF(.NOT.OORANGEIND)THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                     WRKRC(IR) = WRKRC(IR) + ELM
                     IF(IR.NE.IC) THEN
                        WRKRC(IC) = WRKRC(IC) + ELM
                     ENDIF
                  ENDDO
               ENDIF
             ELSE ! NOMP_MAX>1
               IF((ITER.EQ.1).OR.(OORANGEIND))THEN
!$OMP   PARALLEL PRIVATE(IOMP)
!$OMP&  NUM_THREADS(NOMP_MAX)
!$OMP&    IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                    IOMP  = 1
!$                  IOMP = OMP_GET_THREAD_NUM() + 1
!$OMP     DO PRIVATE(NZIND,IR,IC,ELM)
!$OMP&    SCHEDULE(STATIC,CHUNK_NZ)
!$OMP&    REDUCTION(.OR.:OORANGEIND)
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                        WRKR_TH(IR,IOMP) = WRKR_TH(IR,IOMP) + ELM
                        IF(IR.NE.IC) THEN
                          WRKR_TH(IC,IOMP) = WRKR_TH(IC,IOMP) + ELM
                        ENDIF
                     ELSE
                        OORANGEIND = .TRUE.
                     ENDIF
                  ENDDO
!$OMP     END DO
!$OMP   END PARALLEL
               ELSEIF(.NOT.OORANGEIND)THEN
!$OMP   PARALLEL PRIVATE(IOMP)
!$OMP&  NUM_THREADS(NOMP_MAX)
!$OMP&    IF ( NZ_loc > int(K361,8) .AND. NOMP .GT. 1)
                    IOMP  = 1
!$                  IOMP = OMP_GET_THREAD_NUM() + 1
!$OMP     DO PRIVATE(NZIND,IR,IC,ELM)
!$OMP&    SCHEDULE(STATIC,CHUNK_NZ)
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                     WRKR_TH(IR,IOMP) = WRKR_TH(IR,IOMP) + ELM
                     IF(IR.NE.IC) THEN
                       WRKR_TH(IC,IOMP) = WRKR_TH(IC,IOMP) + ELM
                     ENDIF
                  ENDDO
!$OMP     END DO
!$OMP   END PARALLEL
C}
               ENDIF
C
C              For all i on MYID:
C              Build WRKRC(i) = Sum       (WRKR_TH(i,IOMP)
C                              IOMP \in [1:NOMP_MAX]
               IF(NUMPROCS > 1) THEN
                 CALL SMUMPS_REDUCE_WRK_MPI (WRKRC, N, WRKR_TH, 
     &                 NOMP_MAX, 
     &                 IWRK(IMYRPTR),INUMMYR)
               ELSE
                 CALL SMUMPS_REDUCE_WRK (WRKRC, N, WRKR_TH, NOMP_MAX)
               ENDIF
             ENDIF
               IF(NUMPROCS > 1) THEN
C{
                  CALL SMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = SMUMPS_ERRSCALOC(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR, NOMP_MAX) 
C     mpi allreduce.
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_REAL,
     &                    MPI_MAX, COMM, IERROR)
                     IF(ONEERRG.LE.EPS) THEN
                        CALL SMUMPS_UPDATESCALE(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
C}
               ELSE
C{             SINGLE-PROCESSOR CASE: ONE-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = SMUMPS_ERRSCA1(SCA, 
     &                    WRKRC(ITDRPTR), N, NOMP_MAX)
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                        CALL SMUMPS_UPSCALE1(SCA, WRKRC(ITDRPTR), N,
     &                                       NOMP_MAX)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
               ENDIF
C}
            ENDIF
            IF(NUMPROCS > 1) THEN
               CALL SMUMPS_UPDATESCALE(SCA, WRKRC(ITDRPTR), N,
     &              IWRK(IMYRPTR),INUMMYR, NOMP_MAX)
            ELSE
               CALL SMUMPS_UPSCALE1(SCA, WRKRC(ITDRPTR), N, 
     &                              NOMP_MAX)
            ENDIF     
            ITER = ITER + 1
C}
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
C{
            CALL MPI_REDUCE(SCA, WRKRC(1), N, MPI_REAL,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
              IF (NOMP_MAX.LE.0) THEN
               DO I=1, N
                  SCA(I) = WRKRC(I)
               ENDDO
              ELSE
!$OMP       PARALLEL DO PRIVATE(I)
!$OMP&      SCHEDULE(STATIC,CHUNK)
!$OMP&      IF ( N > K361 .AND. NOMP .GT. 1)
               DO I=1, N
                  SCA(I) = WRKRC(I)
               ENDDO
!$OMP       END PARALLEL DO
              ENDIF
            ENDIF
C}
         ENDIF
C}
      ENDIF
      RETURN
      END SUBROUTINE SMUMPS_SIMSCALEABSSYM
