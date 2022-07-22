/*
 *
 *  This file is part of MUMPS 5.2.1, released
 *  on Fri Jun 14 14:46:05 UTC 2019
 *
 *
 *  Copyright 1991-2019 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license:
 *  http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.html
 *
 */
#if defined(_WIN32)
#include  "elapse.h"
#include  <time.h>
#include  <sys/timeb.h>

void MUMPS_CALL mumps_elapse(double* val) {
	time_t ltime;
	struct _timeb tstruct;

	time(&ltime);
	_ftime(&tstruct);
	*val = (double)ltime + (double)tstruct.millitm * (0.001);
}

#else

#include "elapse.h"
#include <sys/time.h>
void mumps_elapse(double *val)
  {
    struct timeval time;
    gettimeofday(&time,(struct timezone *)0);
    *val=time.tv_sec+time.tv_usec*1.e-6;
  }
#endif
