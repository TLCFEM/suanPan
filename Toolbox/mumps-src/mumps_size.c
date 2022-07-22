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
/* Utility to automatically get the sizes of Fortran types */
#include "mumps_size.h"

void MUMPS_CALL MUMPS_SIZE_C(char* a, char* b, MUMPS_INT8* diff) { *diff = (MUMPS_INT8)(b - a); }

void MUMPS_CALL MUMPS_ADDR_C(char* a, MUMPS_INT8* addr) {
	*addr = *(MUMPS_INT8*)&a;
	/*
	    With the form "*addr=(MUMPS_INT8)a", "(MUMPS_INT8)a"
	    and "a" may have different binary representations
	    for large addresses. In the aboce code, "(MUMPS_INT8*)&a"
	    is a pointer to the address "a", considering that "a" is
	    a MUMPS_INT8 rather than an address. Then the content of
	    that pointer is the exact binary representation of the
	    address a, but stored in a MUMPS_INT8 (signed 64-bit
	    integer).
	*/
}
