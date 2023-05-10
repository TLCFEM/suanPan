/*
 *
 *  This file is part of MUMPS 5.6.0, released
 *  on Wed Apr 19 15:50:57 UTC 2023
 *
 *
 *  Copyright 1991-2023 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#ifndef MUMPS_ADDR_H
#define MUMPS_ADDR_H
#include "mumps_common.h"
#include "mumps_c_types.h"
#define MUMPS_SIZE_C \
        F_SYMBOL(size_c, SIZE_C)
void  MUMPS_CALL MUMPS_SIZE_C(char *a, char *b, MUMPS_INT8 *diff);
#define MUMPS_ADDR_C \
        F_SYMBOL(addr_c, ADDR_C)
void  MUMPS_CALL MUMPS_ADDR_C(char *a, MUMPS_INT8 *addr);
#define MUMPS_GETVAL_ADDR_C \
        F_SYMBOL(getval_addr_c, GETVAL_ADDR_C)
void  MUMPS_CALL MUMPS_GETVAL_AT_ADDR_C(volatile MUMPS_INT *val, MUMPS_INT8 *addr);
#define MUMPS_SETRVAL_ADDR_C \
        F_SYMBOL(setrval_addr_c, SETRVAL_ADDR_C)
void MUMPS_CALL MUMPS_SETRVAL_ADDR_C(SMUMPS_REAL *val, MUMPS_INT8 *addr);
#define MUMPS_SETDVAL_ADDR_C \
        F_SYMBOL(setdval_addr_c, SETDVAL_ADDR_C)
void MUMPS_CALL MUMPS_SETDVAL_ADDR_C(DMUMPS_REAL *val, MUMPS_INT8 *addr);
#endif
