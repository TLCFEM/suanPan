#ifndef MPI_H
#define MPI_H

#ifdef INTSIZE64
#include <inttypes.h>
#define LIBSEQ_INT int64_t
#else
#define LIBSEQ_INT int
#endif

#if ! defined(LIBSEQ_CALL)
#if defined(_WIN32) && ! defined(__MINGW32__)
#define LIBSEQ_CALL
#else
#define LIBSEQ_CALL
#endif
#endif

#ifndef MUMPS_MPI_H
#define MUMPS_MPI_H

#ifdef __cplusplus
extern "C" {
#endif

typedef LIBSEQ_INT MPI_Comm;
static MPI_Comm MPI_COMM_WORLD = (MPI_Comm)0;

LIBSEQ_INT LIBSEQ_CALL MPI_Init(LIBSEQ_INT* pargc, char*** pargv);
LIBSEQ_INT LIBSEQ_CALL MPI_Comm_rank(LIBSEQ_INT comm, LIBSEQ_INT* rank);
LIBSEQ_INT LIBSEQ_CALL MPI_Finalize();

#ifdef __cplusplus
}
#endif

#endif

#endif
