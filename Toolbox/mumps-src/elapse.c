/*
 *
 *  This file is part of MUMPS 5.5.1, released
 *  on Tue Jul 12 13:17:24 UTC 2022
 *
 *
 *  Copyright 1991-2022 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#if defined(_WIN32)
#include "elapse.h"
#include <sys/timeb.h>
#include <time.h>
void MUMPS_CALL mumps_elapse(double* val) {
    time_t ltime;
    struct _timeb tstruct;

    time(&ltime);
    _ftime(&tstruct);
    *val = (double)ltime + (double)tstruct.millitm * (0.001);
}

#else

#include <sys/time.h>
#include "elapse.h"
void mumps_elapse(double* val) {
    struct timeval time;
    gettimeofday(&time, (struct timezone*)0);
    *val = time.tv_sec + time.tv_usec * 1.e-6;
}
#endif
