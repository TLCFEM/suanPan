#ifndef MUMPS_COMPAT_H
#define MUMPS_COMPAT_H

#if defined(_WIN32) && ! defined(__MINGW32__)
#define MUMPS_WIN32 1
#endif

#ifndef MUMPS_CALL
#ifdef MUMPS_WIN32
#define MUMPS_CALL
#else
#define MUMPS_CALL
#endif
#endif

#if (__STDC_VERSION__ >= 199901L)
#define MUMPS_INLINE static inline
#else
#define MUMPS_INLINE
#endif

#endif
