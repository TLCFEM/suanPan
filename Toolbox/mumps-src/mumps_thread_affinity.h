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
#ifndef MUMPS_THREAD_AFFINITY_H
#define MUMPS_THREAD_AFFINITY_H
#include "mumps_common.h"
#define MUMPS_THREAD_AFFINITY_RETURN \
  F_SYMBOL(thread_affinity_return,THREAD_AFFINITY_RETURN)
void MUMPS_CALL
MUMPS_THREAD_AFFINITY_RETURN();
#endif /* MUMPS_THREAD_AFFINITY_H */
