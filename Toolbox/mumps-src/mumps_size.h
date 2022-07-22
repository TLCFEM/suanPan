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
#ifndef MUMPS_SIZE_H
#define MUMPS_SIZE_H
#include "mumps_common.h"
#include "mumps_c_types.h"
#define MUMPS_SIZE_C \
        F_SYMBOL(size_c, SIZE_C)
void MUMPS_CALL MUMPS_SIZE_C(char* a, char* b, MUMPS_INT8* diff);
#define MUMPS_ADDR_C \
        F_SYMBOL(addr_c, ADDR_C)
void MUMPS_CALL MUMPS_ADDR_C(char* a, MUMPS_INT8* addr);
#endif
