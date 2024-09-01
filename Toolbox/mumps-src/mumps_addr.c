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
/* Utility to automatically get the sizes of Fortran types */
#include "mumps_addr.h"
void  MUMPS_CALL MUMPS_INT_SIZE_C( MUMPS_INT8 *i)
{
   *i=sizeof(MUMPS_INT);
}
void  MUMPS_CALL MUMPS_SIZE_C(char *a, char *b, MUMPS_INT8 *diff)
{
    *diff = (MUMPS_INT8) (b - a);
}
void MUMPS_CALL MUMPS_ADDR_C(char *a, MUMPS_INT8 *addr)
{
*addr=*(MUMPS_INT8*)&a;
/*
    With the form "*addr=(MUMPS_INT8)a", "(MUMPS_INT8)a"
    and "a" may have different binary representations
    for large addresses. In the above code, "(MUMPS_INT8*)&a"
    is a pointer to the address "a", considering that "a" is
    a MUMPS_INT8 rather than an address. Then the content of
    that pointer is the exact binary representation of the
    address a, but stored in a MUMPS_INT8 (signed 64-bit
    integer).
*/
}
void MUMPS_CALL MUMPS_GETVAL_ADDR_C(volatile MUMPS_INT *val, MUMPS_INT8 *addr)
{
*val=*(MUMPS_INT*)*addr;
}
void MUMPS_CALL MUMPS_SETRVAL_ADDR_C(SMUMPS_REAL *val, MUMPS_INT8 *addr)
{
*(SMUMPS_REAL*)*addr=*val;
}
void MUMPS_CALL MUMPS_SETDVAL_ADDR_C(DMUMPS_REAL *val, MUMPS_INT8 *addr)
{
*(DMUMPS_REAL*)*addr=*val;
}
void  MUMPS_CALL MUMPS_CLANGAOCC_C( MUMPS_INT8 *i8)
{
#if defined(__aocc__)
   *i8=1;
#else
   *i8=0;
#endif
}
