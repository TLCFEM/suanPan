/*
 *
 *  This file is part of MUMPS 5.7.3, released
 *  on Mon Jul 15 11:44:21 UTC 2024
 *
 *
 *  Copyright 1991-2024 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
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
#define MUMPS_INT_SIZE_C \
        F_SYMBOL(int_size_c, INT_SIZE_C)
void  MUMPS_CALL MUMPS_INT_SIZE_C(MUMPS_INT8 *i);
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
#define MUMPS_CLANGAOCC_C \
        F_SYMBOL(clangaocc_c, CLANGAOCC_C)
void MUMPS_CALL MUMPS_CLANGAOCC_C(MUMPS_INT8 *i8);
#endif
