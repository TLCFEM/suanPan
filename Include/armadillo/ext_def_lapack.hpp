namespace arma {
#ifdef ARMA_BLAS_CAPITALS

#define arma_sgbmv SGBMV
#define arma_dgbmv DGBMV

#define arma_dsymv DSYMV

#define arma_ssbmv SSBMV
#define arma_dsbmv DSBMV

#define arma_sspmv SSPMV
#define arma_dspmv DSPMV

#define arma_dsysv DSYSV
#define arma_dsygvx DSYGVX

#define arma_spbsv SPBSV
#define arma_dpbsv DPBSV
#define arma_spbtrs SPBTRS
#define arma_dpbtrs DPBTRS

#define arma_sspsv SSPSV
#define arma_dspsv DSPSV
#define arma_ssptrf SSPTRF
#define arma_dsptrf DSPTRF
#define arma_ssptrs SSPTRS
#define arma_dsptrs DSPTRS
#define arma_ssptri SSPTRI
#define arma_dsptri DSPTRI

#define arma_sppsv SPPSV
#define arma_dppsv DPPSV
#define arma_spptrf SPPTRF
#define arma_dpptrf DPPTRF
#define arma_spptrs SPPTRS
#define arma_dpptrs DPPTRS
#define arma_spptri SPPTRI
#define arma_dpptri DPPTRI

#define arma_dsgesv DSGESV
#define arma_dsposv DSPOSV

#else

#define arma_sgbmv sgbmv
#define arma_dgbmv dgbmv

#define arma_dsymv dsymv

#define arma_ssbmv ssbmv
#define arma_dsbmv dsbmv

#define arma_sspmv sspmv
#define arma_dspmv dspmv

#define arma_dsysv dsysv
#define arma_dsygvx dsygvx

#define arma_spbsv spbsv
#define arma_dpbsv dpbsv
#define arma_spbtrs spbtrs
#define arma_dpbtrs dpbtrs

#define arma_sspsv sspsv
#define arma_dspsv dspsv
#define arma_ssptrf ssptrf
#define arma_dsptrf dsptrf
#define arma_ssptrs ssptrs
#define arma_dsptrs dsptrs
#define arma_ssptri ssptri
#define arma_dsptri dsptri

#define arma_sppsv sppsv
#define arma_dppsv dppsv
#define arma_spptrf spptrf
#define arma_dpptrf dpptrf
#define arma_spptrs spptrs
#define arma_dpptrs dpptrs
#define arma_spptri spptri
#define arma_dpptri dpptri

#define arma_dsgesv dsgesv
#define arma_dsposv dsposv

#endif

    extern "C" {
    void arma_fortran(arma_sgbmv)(const char* TRANS, const blas_int* M, const blas_int* N, const blas_int* KL, const blas_int* KU, const float* ALPHA, const float* A, const blas_int* LDA, const float* X, const blas_int* INCX, const float* BETA, float* Y, const blas_int* INCY);

    void arma_fortran(arma_dgbmv)(const char* TRANS, const blas_int* M, const blas_int* N, const blas_int* KL, const blas_int* KU, const double* ALPHA, const double* A, const blas_int* LDA, const double* X, const blas_int* INCX, const double* BETA, double* Y, const blas_int* INCY);

    void arma_fortran(arma_dsymv)(char* UPLO, blas_int* N, double* ALPHA, double* A, blas_int* LDA, double* X, blas_int* INCX, double* BETA, double* Y, blas_int* INCY);

    void arma_fortran(arma_ssbmv)(const char* UPLO, const blas_int* N, const blas_int* K, const float* ALPHA, const float* A, const blas_int* LDA, const float* X, const blas_int* INCX, const float* BETA, float* Y, const blas_int* INCY);

    void arma_fortran(arma_dsbmv)(const char* UPLO, const blas_int* N, const blas_int* K, const double* ALPHA, const double* A, const blas_int* LDA, const double* X, const blas_int* INCX, const double* BETA, double* Y, const blas_int* INCY);

    void arma_fortran(arma_sspmv)(const char* UPLO, const blas_int* N, const float* ALPHA, const float* AP, const float* X, const blas_int* INCX, const float* BETA, float* Y, const blas_int* INCY);

    void arma_fortran(arma_dspmv)(const char* UPLO, const blas_int* N, const double* ALPHA, const double* AP, const double* X, const blas_int* INCX, const double* BETA, double* Y, const blas_int* INCY);

    // symmetric matrix
    void arma_fortran(arma_dsysv)(char* UPLO, blas_int* N, blas_int* NRHS, double* A, blas_int* LDA, blas_int* IPIV, double* B, blas_int* LDB, double* WORK, blas_int* LWORK, blas_int* INFO);

    void arma_fortran(arma_dsygvx)(blas_int* ITYPE, char* JOBZ, char* RANGE, char* UPLO, blas_int* N, double* A, blas_int* LDA, double* B, blas_int* LDB, double* VL, double* VU, blas_int* IL, blas_int* IU, double* ABSTOL, blas_int* M, double* W, double* Z, blas_int* LDZ, double* WORK, blas_int* LWORK, blas_int* IWORK, blas_int* IFAIL, blas_int* INFO);

    // symmetric positive definite band matrix
    void arma_fortran(arma_spbsv)(const char* UPLO, const blas_int* N, const blas_int* KD, const blas_int* NRHS, float* AB, const blas_int* LDAB, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dpbsv)(const char* UPLO, const blas_int* N, const blas_int* KD, const blas_int* NRHS, double* AB, const blas_int* LDAB, double* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_spbtrs)(const char* UPLO, const blas_int* N, const blas_int* KD, const blas_int* NRHS, float* AB, const blas_int* LDAB, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dpbtrs)(const char* UPLO, const blas_int* N, const blas_int* KD, const blas_int* NRHS, double* AB, const blas_int* LDAB, double* B, const blas_int* LDB, blas_int* INFO);

    // symmetric matrix in packed format
    void arma_fortran(arma_sspsv)(const char* UPLO, const blas_int* N, const blas_int* NRHS, float* AP, blas_int* IPIV, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dspsv)(const char* UPLO, const blas_int* N, const blas_int* NRHS, double* AP, blas_int* IPIV, double* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_ssptrf)(const char* UPLO, const blas_int* N, float* AP, blas_int* IPIV, blas_int* INFO);

    void arma_fortran(arma_dsptrf)(const char* UPLO, const blas_int* N, double* AP, blas_int* IPIV, blas_int* INFO);

    void arma_fortran(arma_ssptrs)(const char* UPLO, const blas_int* N, const blas_int* NRHS, const float* AP, const blas_int* IPIV, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dsptrs)(const char* UPLO, const blas_int* N, const blas_int* NRHS, const double* AP, const blas_int* IPIV, double* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_ssptri)(const char* UPLO, const blas_int* N, float* AP, const blas_int* IPIV, float* WORK, blas_int* INFO);

    void arma_fortran(arma_dsptri)(const char* UPLO, const blas_int* N, double* AP, const blas_int* IPIV, double* WORK, blas_int* INFO);

    void arma_fortran(arma_sppsv)(const char* UPLO, const blas_int* N, const blas_int* NRHS, float* AP, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dppsv)(const char* UPLO, const blas_int* N, const blas_int* NRHS, double* AP, double* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_spptrf)(const char* UPLO, const blas_int* N, float* AP, blas_int* INFO);

    void arma_fortran(arma_dpptrf)(const char* UPLO, const blas_int* N, double* AP, blas_int* INFO);

    void arma_fortran(arma_spptrs)(const char* UPLO, const blas_int* N, const blas_int* NRHS, const float* AP, float* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_dpptrs)(const char* UPLO, const blas_int* N, const blas_int* NRHS, const double* AP, double* B, const blas_int* LDB, blas_int* INFO);

    void arma_fortran(arma_spptri)(const char* UPLO, const blas_int* N, float* AP, float* WORK, blas_int* INFO);

    void arma_fortran(arma_dpptri)(const char* UPLO, const blas_int* N, double* AP, double* WORK, blas_int* INFO);

    void arma_fortran(arma_dsgesv)(blas_int*, blas_int*, double*, blas_int*, blas_int*, double*, blas_int*, double*, blas_int*, double*, float*, blas_int*, blas_int*);

    void arma_fortran(arma_dsposv)(const char*, blas_int*, blas_int*, double*, blas_int*, double*, blas_int*, double*, blas_int*, double*, float*, blas_int*, blas_int*);
    }
} // namespace arma
