#ifdef ARMA_BLAS_CAPITALS

#define arma_sgbmv SGBMV
#define arma_dgbmv DGBMV

#define arma_dsymv DSYMV

#define arma_ssbmv SSBMV
#define arma_dsbmv DSBMV

#define arma_sspmv SSPMV
#define arma_dspmv DSPMV
#define arma_sspmm SSPMM
#define arma_dspmm DSPMM

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
#define arma_sspmm sspmm
#define arma_dspmm dspmm

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
void arma_fortran(arma_sgbmv)(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU, const float* ALPHA, const float* A, const int* LDA, const float* X, const int* INCX, const float* BETA, float* Y, const int* INCY);

void arma_fortran(arma_dgbmv)(const char* TRANS, const int* M, const int* N, const int* KL, const int* KU, const double* ALPHA, const double* A, const int* LDA, const double* X, const int* INCX, const double* BETA, double* Y, const int* INCY);

void arma_fortran(arma_dsymv)(char* UPLO, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);

void arma_fortran(arma_ssbmv)(const char* UPLO, const int* N, const int* K, const float* ALPHA, const float* A, const int* LDA, const float* X, const int* INCX, const float* BETA, float* Y, const int* INCY);

void arma_fortran(arma_dsbmv)(const char* UPLO, const int* N, const int* K, const double* ALPHA, const double* A, const int* LDA, const double* X, const int* INCX, const double* BETA, double* Y, const int* INCY);

void arma_fortran(arma_sspmv)(const char* UPLO, const int* N, const float* ALPHA, const float* AP, const float* X, const int* INCX, const float* BETA, float* Y, const int* INCY);

void arma_fortran(arma_dspmv)(const char* UPLO, const int* N, const double* ALPHA, const double* AP, const double* X, const int* INCX, const double* BETA, double* Y, const int* INCY);

void arma_fortran(arma_sspmm)(const char* SIDE, const char* UPLO, const char* TRAN, const int* M, const int* N, const float* A, const float* ALPHA, const float* B, const int* LDB, const float* BETA, float* C, const int* LDC);

void arma_fortran(arma_dspmm)(const char* SIDE, const char* UPLO, const char* TRAN, const int* M, const int* N, const double* A, const double* ALPHA, const double* B, const int* LDB, const double* BETA, double* C, const int* LDC);

// symmetric matrix
void arma_fortran(arma_dsysv)(char* UPLO, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, double* WORK, int* LWORK, int* INFO);

void arma_fortran(arma_dsygvx)(int* ITYPE, char* JOBZ, char* RANGE, char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, double* VL, double* VU, int* IL, int* IU, double* ABSTOL, int* M, double* W, double* Z, int* LDZ, double* WORK, int* LWORK, int* IWORK, int* IFAIL, int* INFO);

// symmetric positive definite band matrix
void arma_fortran(arma_spbsv)(const char* UPLO, const int* N, const int* KD, const int* NRHS, float* AB, const int* LDAB, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dpbsv)(const char* UPLO, const int* N, const int* KD, const int* NRHS, double* AB, const int* LDAB, double* B, const int* LDB, int* INFO);

void arma_fortran(arma_spbtrs)(const char* UPLO, const int* N, const int* KD, const int* NRHS, float* AB, const int* LDAB, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dpbtrs)(const char* UPLO, const int* N, const int* KD, const int* NRHS, double* AB, const int* LDAB, double* B, const int* LDB, int* INFO);

// symmetric matrix in packed format
void arma_fortran(arma_sspsv)(const char* UPLO, const int* N, const int* NRHS, float* AP, int* IPIV, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dspsv)(const char* UPLO, const int* N, const int* NRHS, double* AP, int* IPIV, double* B, const int* LDB, int* INFO);

void arma_fortran(arma_ssptrf)(const char* UPLO, const int* N, float* AP, int* IPIV, int* INFO);

void arma_fortran(arma_dsptrf)(const char* UPLO, const int* N, double* AP, int* IPIV, int* INFO);

void arma_fortran(arma_ssptrs)(const char* UPLO, const int* N, const int* NRHS, const float* AP, const int* IPIV, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dsptrs)(const char* UPLO, const int* N, const int* NRHS, const double* AP, const int* IPIV, double* B, const int* LDB, int* INFO);

void arma_fortran(arma_ssptri)(const char* UPLO, const int* N, float* AP, const int* IPIV, float* WORK, int* INFO);

void arma_fortran(arma_dsptri)(const char* UPLO, const int* N, double* AP, const int* IPIV, double* WORK, int* INFO);

void arma_fortran(arma_sppsv)(const char* UPLO, const int* N, const int* NRHS, float* AP, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dppsv)(const char* UPLO, const int* N, const int* NRHS, double* AP, double* B, const int* LDB, int* INFO);

void arma_fortran(arma_spptrf)(const char* UPLO, const int* N, float* AP, int* INFO);

void arma_fortran(arma_dpptrf)(const char* UPLO, const int* N, double* AP, int* INFO);

void arma_fortran(arma_spptrs)(const char* UPLO, const int* N, const int* NRHS, const float* AP, float* B, const int* LDB, int* INFO);

void arma_fortran(arma_dpptrs)(const char* UPLO, const int* N, const int* NRHS, const double* AP, double* B, const int* LDB, int* INFO);

void arma_fortran(arma_spptri)(const char* UPLO, const int* N, float* AP, float* WORK, int* INFO);

void arma_fortran(arma_dpptri)(const char* UPLO, const int* N, double* AP, double* WORK, int* INFO);

void arma_fortran(arma_dsgesv)(int*, int*, double*, int*, int*, double*, int*, double*, int*, double*, float*, int*, int*);

void arma_fortran(arma_dsposv)(const char*, int*, int*, double*, int*, double*, int*, double*, int*, double*, float*, int*, int*);
}
