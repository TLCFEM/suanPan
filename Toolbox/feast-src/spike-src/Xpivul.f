
      SUBROUTINE T_DEC(GBTRFUL)( M, N, KL, KU, AB, LDAB, IPIV, INFO )
      IMPLICIT NONE 
* Inputs and outputs
      INTEGER M,N,KL,KU,LDAB,INFO
      INTEGER IPIV(N)
      MAIN_TYPE AB(LDAB,N)
* Iterators and local variables
      INTEGER J
      INTEGER KLU,KTOT
* Debug
*      REAL START_TIME_SWAP,END_TIME_SWAP,START_TIME_FAC,END_TIME_FAC
#ifdef EXTRA_DEBUG
      DOUBLE PRECISION START_TIME_SWAP,END_TIME_SWAP,START_TIME_FAC,
     $ END_TIME_FAC
      DOUBLE PRECISION OMP_GET_WTIME
      character(len=100) :: format
#endif
* For now, we'll swap all values in AB, including the don't cares inside
* top-left and bottom right -- slightly wasteful but not too bad. 

*      write(format,'(A1,i2,A6)'),'(', n, 'E12.2)'
*      format = '(8E10.2)'
      KLU=MAX(KL,KU)
      KTOT=KLU+KL+KU+1

*      CALL CPU_TIME(START_TIME_SWAP)
*      START_TIME_SWAP = OMP_GET_WTIME()

      DO 25 J=1,N/2
   25   CALL T_DEC(SWAP)(KL+KU+1,AB(KLU+1,J),1,AB(KLU+1,N+1-J),-1)

*     If the number of columns is odd, we still have to swap the middlemost col. 
      IF(MOD(N, 2) .eq. 1)  CALL T_DEC(SWAP)((KL+KU+1)/2,
     $ AB(KLU+1,N/2+1),1,AB(KTOT-(KL+KU+1)/2+1,N/2+1),-1)

*      CALL CPU_TIME(END_TIME_SWAP)
#ifdef EXTRA_DEBUG
      END_TIME_SWAP = OMP_GET_WTIME()
 
*      CALL CPU_TIME(START_TIME_FAC)
      START_TIME_FAC = OMP_GET_WTIME()
#endif     

      CALL T_DEC(GBTRF)(M,N,KU,KL,AB(KLU-KU+1,1),LDAB,IPIV,INFO)

*      CALL CPU_TIME(END_TIME_FAC)
#ifdef EXTRA_DEBUG
      END_TIME_FAC = OMP_GET_WTIME()
#endif
      END SUBROUTINE T_DEC(GBTRFUL)


*  Definition:
*  ===========
*
*       SUBROUTINE DGBTRSL( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*       ..
*  
*
*  =============
*
*
* DGBTRSL performs the first sweep of the LU solve
* This was created by removing the second sweep from 
* default DGBTRS
*    A * X = B  or  A**T * X = B
* with a general band matrix A using the LU factorization computed
* by DGBTRF.
*
*  Arguments:
*  ==========
*
* \param[in] TRANS
* \verbatim
*          TRANS is CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T* X = B  (Transpose)
*          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
* \endverbatim
*
* \param[in] N
* \verbatim
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
* \endverbatim
*
* \param[in] KL
* \verbatim
*          KL is INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
* \endverbatim
*
* \param[in] KU
* \verbatim
*          KU is INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
* \endverbatim
*
* \param[in] NRHS
* \verbatim
*          NRHS is INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
* \endverbatim
*
* \param[in] AB
* \verbatim
*          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
* \endverbatim
*
* \param[in] LDAB
* \verbatim
*          LDAB is INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
* \endverbatim
*
* \param[in] IPIV
* \verbatim
*          IPIV is INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
* \endverbatim
*
* \param[in,out] B
* \verbatim
*          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
* \endverbatim
*
* \param[in] LDB
* \verbatim
*          LDB is INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
* \endverbatim
*
* \param[out] INFO
* \verbatim
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
* \endverbatim
*
*  Authors:
*  ========
*
* \author Univ. of Tennessee 
* \author Univ. of California Berkeley 
* \author Univ. of Colorado Denver 
* \author NAG Ltd. 
*
* \date November 2011
*
* \ingroup doubleGBcomputational
*
*  =====================================================================
      SUBROUTINE T_DEC(GBTRSL)(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, 
     $                         B, LDB, INFO )
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      MAIN_TYPE AB( LDAB, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            LNOTI, NOTRAN
      INTEGER            J, KD, L, LM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           T_DEC(GEMV), T_DEC(SWAP), T_DEC(TBSV), XERBLA
#ifdef SPIKECOMPLEX
      EXTERNAL           T_DEC(GERU)
#else
      EXTERNAL           T_DEC(GER)
#endif
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( MAX(KL,KU)+KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( '?GBTRSL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      KD = KU + KL + 1
      LNOTI = KL.GT.0
*
      IF( NOTRAN ) THEN
*
*        Solve L*X = B, overwriting B with X.
*
*        L is represented as a product of permutations and unit lower
*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
*        where each transformation L(i) is a rank-one modification of
*        the identity matrix.
*
         IF( LNOTI ) THEN
            DO 10 J = 1, N - 1
               LM = MIN( KL, N-J )
               L = IPIV( J )
               IF( L.NE.J ) THEN
                  CALL T_DEC(SWAP)(NRHS, B( L, 1 ), LDB, B( J, 1 ), LDB)
               END IF 
#ifdef SPIKECOMPLEX
               CALL T_DEC(GERU)(LM, NRHS, -ONE_PREC, AB( KD+1, J ), 1,
     $                   B( J, 1 ), LDB, B( J+1, 1 ), LDB )
#else
               CALL T_DEC(GER)(LM, NRHS, -ONE_PREC, AB( KD+1, J ), 1,
     $                    B( J, 1 ), LDB, B( J+1, 1 ), LDB )
#endif
   10       CONTINUE
         END IF
      ELSE
*
*        Solve A**T*X = B.
*
*         DO 30 I = 1, NRHS
*
*           Solve U**T*X = B, overwriting B with X.
*
*            CALL T_DEC(TBSV)( 'Upper', 'Transpose', 'Non-unit', N, KL+KU, AB,
*     $                  LDAB, B( 1, I ), 1 )
*   30    CONTINUE
      CALL T_DEC(TBSM)('U','T','N',N,NRHS,KL+KU,AB,LDAB,B,LDB)
      END IF
      RETURN
*
*     End of DGBTRSL
*
      END






*
*  =========== DOCUMENTATION ===========
*
*       SUBROUTINE DGBTRSU( TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV, B, LDB,
*                          INFO )
* 
*       .. Scalar Arguments ..
*       CHARACTER          TRANS
*       INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*       ..
*       .. Array Arguments ..
*       INTEGER            IPIV( * )
*       DOUBLE PRECISION   AB( LDAB, * ), B( LDB, * )
*       ..
*  
*
*  =============
*
*
* DGBTRSU performs the second sweep in the LU solve
* This was created by removing the first sweep from 
* default DGBTRS
*    A * X = B  or  A**T * X = B
* with a general band matrix A using the LU factorization computed
* by DGBTRF.
*
* Arguments:
* ==========
*
* \param[in] TRANS
* \verbatim
*          TRANS is CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T* X = B  (Transpose)
*          = 'C':  A**T* X = B  (Conjugate transpose = Transpose)
* \endverbatim
*
* \param[in] N
* \verbatim
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
* \endverbatim
*
* \param[in] KL
* \verbatim
*          KL is INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
* \endverbatim
*
* \param[in] KU
* \verbatim
*          KU is INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
* \endverbatim
*
* \param[in] NRHS
* \verbatim
*          NRHS is INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
* \endverbatim
*
* \param[in] AB
* \verbatim
*          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
* \endverbatim
*
* \param[in] LDAB
* \verbatim
*          LDAB is INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
* \endverbatim
*
* \param[in] IPIV
* \verbatim
*          IPIV is INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
* \endverbatim
*
* \param[in,out] B
* \verbatim
*          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
* \endverbatim
*
* \param[in] LDB
* \verbatim
*          LDB is INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
* \endverbatim
*
* \param[out] INFO
* \verbatim
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
* \endverbatim
*
* Authors:
* ========
*
* \author Univ. of Tennessee 
* \author Univ. of California Berkeley 
* \author Univ. of Colorado Denver 
* \author NAG Ltd. 
*
* \date November 2011
*
* \ingroup doubleGBcomputational
*
*  =====================================================================
      SUBROUTINE T_DEC(GBTRSU)(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV,
     $                         B, LDB, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine (version 3.4.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, KL, KU, LDAB, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      MAIN_TYPE AB( LDAB, * ), B( LDB, * )
*     ..
*
*  =====================================================================
*
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
#ifdef EXTRA_DEBUG
      INTEGER            I, J, KD, L, LM
#endif
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           T_DEC(GEMV), T_DEC(SWAP), T_DEC(TBSV), XERBLA
#ifdef SPIKECOMPLEX
      EXTERNAL           T_DEC(GERU)
#else
      EXTERNAL           T_DEC(GER)
#endif

*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KL.LT.0 ) THEN
         INFO = -3
      ELSE IF( KU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDAB.LT.( MAX(KL,KU)+KL+KU+1 ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( '?GBTRSU', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      CALL T_DEC(GBTRS)(TRANS, N, 0, KU+KL, NRHS, AB, LDAB, IPIV, B,LDB,
     $                   INFO )
      END




      SUBROUTINE T_DEC(GBTRSUL)(TRANS, N, KL, KU, NRHS, AB, LDAB, IPIV,
     $                          B, LDB, INFO)
      IMPLICIT NONE 
* Inputs and outputs
      INTEGER N,NRHS,KL,KU,LDAB,LDB,INFO
      INTEGER IPIV(*)
      MAIN_TYPE AB(LDAB,*), B(LDB,*)
      CHARACTER TRANS
* Iterators and local variables
      INTEGER J
      INTEGER KLU
* Debug
#ifdef EXTRA_DEBUG
      INTEGER I, M
      character(len=100) :: format
#endif
*      write(format,'(A1,i1,A6)'),'(', nrhs, 'E12.2)'

*      print *, ""
*      DO 25 I=1,N
*   25    write (*,format), B(I,1:NRHS) 

      KLU=MAX(KL,KU)
      DO 28 J=1,NRHS
   28   CALL T_DEC(SWAP)(N/2,B(1,J),1,B(N-N/2+1,J),-1)

      CALL T_DEC(GBTRSL)(TRANS,N,KU,KL,NRHS,AB(KLU-KU+1,1),
     $                   LDAB,IPIV,B,LDB,INFO)
*      print *, ""
*      DO 25 I=1,N
*   25    write (*,format), B(I,1:NRHS) 
      
      CALL T_DEC(GBTRSU)(TRANS,N,KU,KL,NRHS,AB(KLU-KU+1,1),
     $                   LDAB,IPIV,B,LDB,INFO)
*      CALL T_DEC(GBTRS)(TRANS,N,KU,KL,NRHS,AB(KLU-KU+1,1),
*     $  LDAB,IPIV,B,LDB,INFO)
      DO 34 J=1,NRHS
   34   CALL T_DEC(SWAP)(N/2,B(1,J),1,B(N-N/2+1,J),-1)
      
      END SUBROUTINE T_DEC(GBTRSUL)

















