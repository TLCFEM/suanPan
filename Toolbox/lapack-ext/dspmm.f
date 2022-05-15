* Purpose:  To compute alpha*op(A,B)+beta*C, op(A,B) could be either
*           A*B,
*           B*A,
*           A*B**T,
*           B**T*A.
*
* \param[in] SIDE
*            SIDE is CHARACTER*1
*            Select which problem to compute:
*            'L'    C=alpha*B*A+beta*C,
*            'R'    C=alpha*A*B+beta*C.
*
* \param[in] UPLO
*            UPLO is CHARACTER*1
*            Select which part of A is stored:
*            'U'    Upper Triangle,
*            'L'    Lower Triangle.
*
* \param[in] TRAN
*            TRAN is CHARACTER*1
*            Select if B is transverse:
*            'N'    No transverse,
*            'T'    Transverse.
*
* \param[in] M
*            M is INTEGER
*            The size of square matrix A.
*
* \param[in] N
*            N is INTEGER
*            Another dimension of matrix B.
*            For SIDE='L', B=>(N,M),
*            For SIDE='R', B=>(M,N).
*
* \param[in] ALPHA
*            ALPHA is DOUBLEPRECISION
*            The factor.
*
* \param[in] A
*            A is DOUBLEPRECISION(*) array of DIMENSION ((M+1)*M/2)
*
* \param[in] B
*            B is DOUBLEPRECISION(*,*) array of DIMENSION (M,N) or (N,M)
*
* \param[in] LDB
*            LDB is INTEGER
*            The leading dimension of matrix B, should be at least max(1,M) or max(1,N).
*
* \param[in] BETA
*            BETA is DOUBLEPRECISION
*            The factor.
*
* \param[in/out] C
*                C is DOUBLEPRECISION(*,*) array of DIMENSION (M,N) or (N,M)
*
* \param[in] LDC
*            LDC is INTEGER
*            The leading dimension of matrix C, should be at least max(1,M) or max(1,N) based on SIDE.
*
      SUBROUTINE DSPMM(SIDE,UPLO,TRAN,M,N,A,ALPHA,B,LDB,BETA,C,LDC)

      !...INPUT ARGUMENTS...
      CHARACTER SIDE,UPLO,TRAN
      INTEGER M,N,LDB,LDC
      DOUBLEPRECISION ALPHA,BETA,A(*),B(LDB,*),C(LDC,*)

      !...TEMP VARIABLES...
      INTEGER I,J,K,X,Y,Z,DIMA,DIMB,PTYPE,INFO
      DOUBLEPRECISION TEMPA
      LOGICAL S,U,T

      !...TWO CONSTANTS...
      DOUBLEPRECISION ZERO,ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)

      !...EXTERNAL SUBROUTINES...
      LOGICAL LSAME
      EXTERNAL LSAME
      EXTERNAL XERBLA

      !...FLAGS...
      S=LSAME(SIDE,'R')
      U=LSAME(UPLO,'U')
      T=LSAME(TRAN,'N')

      !...CHECK IF ACCEPTABLE ARGUMENTS ARE GIVEN...
      INFO=0
      IF((.NOT.S).AND.(.NOT.LSAME(SIDE,'L')))THEN
          INFO=1
      ELSEIF((.NOT.U).AND.(.NOT.LSAME(UPLO,'U')))THEN
          INFO=2
      ELSEIF((.NOT.T).AND.(.NOT.LSAME(TRAN,'T')))THEN
          INFO=3
      ENDIF

      IF(INFO.NE.0)THEN
          CALL XERBLA('DSPMM ',INFO)
          RETURN
      ENDIF

      !...SWITCH TO PROPER DIMENSION...
      !...THE DIMENSION OF C IS AWAYS (DIMA,DIMB)...
      IF(S)THEN
          DIMA=M
          DIMB=N
      ELSE
          DIMA=N
          DIMB=M
      ENDIF

      !...QUICK RETURN...
      IF(ALPHA.EQ.ZERO)THEN
          IF(BETA.EQ.ZERO)THEN
              DO 20 J=1,DIMB
                  DO 10 I=1,DIMA
                      C(I,J)=ZERO
   10             CONTINUE
   20         CONTINUE
          ELSEIF(BETA.NE.ONE)THEN
              DO 40 J=1,DIMB
                  DO 30 I=1,DIMA
                      C(I,J)=BETA*C(I,J)
   30             CONTINUE
   40         CONTINUE
          ENDIF
          RETURN
      ENDIF

      !...ALPHA.NE.ZERO...
      !...CHECK beta*C FIRST...
      IF(BETA.EQ.ZERO)THEN
          DO 60 J=1,DIMB
              DO 50 I=1,DIMA
                  C(I,J)=ZERO
   50         CONTINUE
   60     CONTINUE
      ELSEIF(BETA.NE.ONE)THEN
          DO 80 J=1,DIMB
              DO 70 I=1,DIMA
                  C(I,J)=BETA*C(I,J)
   70         CONTINUE
   80     CONTINUE
      ENDIF

      !...ASSIGN PROBLEM TYPE ACCORDING TO GIVEN FLAGS...
      PTYPE=0000
      IF(S)PTYPE=PTYPE+1000
      IF(U)PTYPE=PTYPE+100
      IF(T)PTYPE=PTYPE+10
      IF(ALPHA.EQ.ONE)PTYPE=PTYPE+1

      X=1

      !...U*B
      IF(PTYPE==1111)THEN
          DO J=1,M
              DO I=1,J-1
                  DO K=1,N
                      C(I,K)=C(I,K)+A(X)*B(J,K)
                      C(J,K)=C(J,K)+A(X)*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
              DO K=1,N
                  C(J,K)=C(J,K)+A(X)*B(J,K)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...U*B**T
      IF(PTYPE==1101)THEN
          DO J=1,M
              DO I=1,J-1
                  DO K=1,N
                      C(I,K)=C(I,K)+A(X)*B(K,J)
                      C(J,K)=C(J,K)+A(X)*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
              DO K=1,N
                  C(J,K)=C(J,K)+A(X)*B(K,J)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...B*U
      IF(PTYPE==0111)THEN
          DO J=1,M
              DO I=1,J-1
                  DO K=1,N
                      C(K,I)=C(K,I)+A(X)*B(K,J)
                      C(K,J)=C(K,J)+A(X)*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
              DO K=1,N
                  C(K,J)=C(K,J)+A(X)*B(K,J)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...B**T*U
      IF(PTYPE==0101)THEN
          DO J=1,M
              DO I=1,J-1
                  DO K=1,N
                      C(K,I)=C(K,I)+A(X)*B(J,K)
                      C(K,J)=C(K,J)+A(X)*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
              DO K=1,N
                  C(K,J)=C(K,J)+A(X)*B(J,K)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...L*B
      IF(PTYPE==1011)THEN
          DO J=1,M
              DO K=1,N
                  C(J,K)=C(J,K)+A(X)*B(J,K)
              ENDDO
              X=X+1
              DO I=J+1,M
                  DO K=1,N
                      C(I,K)=C(I,K)+A(X)*B(J,K)
                      C(J,K)=C(J,K)+A(X)*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...L*B**T
      IF(PTYPE==1001)THEN
          DO J=1,M
              DO K=1,N
                  C(J,K)=C(J,K)+A(X)*B(K,J)
              ENDDO
              X=X+1
              DO I=J+1,M
                  DO K=1,N
                      C(I,K)=C(I,K)+A(X)*B(K,J)
                      C(J,K)=C(J,K)+A(X)*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...B*L
      IF(PTYPE==0011)THEN
          DO J=1,M
              DO K=1,N
                  C(K,J)=C(K,J)+A(X)*B(K,J)
              ENDDO
              X=X+1
              DO I=J+1,M
                  DO K=1,N
                      C(K,I)=C(K,I)+A(X)*B(K,J)
                      C(K,J)=C(K,J)+A(X)*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...B**T*L
      IF(PTYPE==0001)THEN
          DO J=1,M
              DO K=1,N
                  C(K,J)=C(K,J)+A(X)*B(J,K)
              ENDDO
              X=X+1
              DO I=J+1,M
                  DO K=1,N
                      C(K,I)=C(K,I)+A(X)*B(J,K)
                      C(K,J)=C(K,J)+A(X)*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...U*B
      IF(PTYPE==1110)THEN
          DO J=1,M
              DO I=1,J-1
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(I,K)=C(I,K)+TEMPA*B(J,K)
                      C(J,K)=C(J,K)+TEMPA*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(J,K)=C(J,K)+TEMPA*B(J,K)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...U*B**T
      IF(PTYPE==1100)THEN
          DO J=1,M
              DO I=1,J-1
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(I,K)=C(I,K)+TEMPA*B(K,J)
                      C(J,K)=C(J,K)+TEMPA*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(J,K)=C(J,K)+TEMPA*B(K,J)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...B*U
      IF(PTYPE==0110)THEN
          DO J=1,M
              DO I=1,J-1
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(K,I)=C(K,I)+TEMPA*B(K,J)
                      C(K,J)=C(K,J)+TEMPA*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(K,J)=C(K,J)+TEMPA*B(K,J)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...B**T*U
      IF(PTYPE==0100)THEN
          DO J=1,M
              DO I=1,J-1
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(K,I)=C(K,I)+TEMPA*B(J,K)
                      C(K,J)=C(K,J)+TEMPA*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(K,J)=C(K,J)+TEMPA*B(J,K)
              ENDDO
              X=X+1
          ENDDO
          RETURN
      ENDIF

      !...L*B
      IF(PTYPE==1010)THEN
          DO J=1,M
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(J,K)=C(J,K)+TEMPA*B(J,K)
              ENDDO
              X=X+1
              DO I=J+1,M
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(I,K)=C(I,K)+TEMPA*B(J,K)
                      C(J,K)=C(J,K)+TEMPA*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...L*B**T
      IF(PTYPE==1000)THEN
          DO J=1,M
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(J,K)=C(J,K)+TEMPA*B(K,J)
              ENDDO
              X=X+1
              DO I=J+1,M
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(I,K)=C(I,K)+TEMPA*B(K,J)
                      C(J,K)=C(J,K)+TEMPA*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...B*L
      IF(PTYPE==0010)THEN
          DO J=1,M
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(K,J)=C(K,J)+TEMPA*B(K,J)
              ENDDO
              X=X+1
              DO I=J+1,M
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(K,I)=C(K,I)+TEMPA*B(K,J)
                      C(K,J)=C(K,J)+TEMPA*B(K,I)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      !...B**T*L
      IF(PTYPE==0000)THEN
          DO J=1,M
              TEMPA=ALPHA*A(X)
              DO K=1,N
                  C(K,J)=C(K,J)+TEMPA*B(J,K)
              ENDDO
              X=X+1
              DO I=J+1,M
                  TEMPA=ALPHA*A(X)
                  DO K=1,N
                      C(K,I)=C(K,I)+TEMPA*B(J,K)
                      C(K,J)=C(K,J)+TEMPA*B(I,K)
                  ENDDO
                  X=X+1
              ENDDO
          ENDDO
          RETURN
      ENDIF

      END
