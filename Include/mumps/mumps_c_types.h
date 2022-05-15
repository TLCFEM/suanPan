#ifndef MUMPS_C_TYPES_H
#define MUMPS_C_TYPES_H

#include <stdint.h>
#ifdef INTSIZE64
#define MUMPS_INT int64_t
#else
#define MUMPS_INT int
#endif

#define MUMPS_INT8 int64_t

#define SMUMPS_COMPLEX float
#define SMUMPS_REAL float

#define DMUMPS_COMPLEX double
#define DMUMPS_REAL double

typedef struct {
	float r, i;
} mumps_complex;

typedef struct {
	double r, i;
} mumps_double_complex;

#define CMUMPS_COMPLEX mumps_complex
#define CMUMPS_REAL float

#define ZMUMPS_COMPLEX mumps_double_complex
#define ZMUMPS_REAL double

#ifndef mumps_ftnlen
# define mumps_ftnlen MUMPS_INT
#endif

#define MUMPS_ARITH_s 1
#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_c 4
#define MUMPS_ARITH_z 8

#define MUMPS_ARITH_REAL   ( MUMPS_ARITH_s | MUMPS_ARITH_d )
#define MUMPS_ARITH_CMPLX  ( MUMPS_ARITH_c | MUMPS_ARITH_z )
#define MUMPS_ARITH_SINGLE ( MUMPS_ARITH_s | MUMPS_ARITH_c )
#define MUMPS_ARITH_DBL    ( MUMPS_ARITH_d | MUMPS_ARITH_z )

#endif
