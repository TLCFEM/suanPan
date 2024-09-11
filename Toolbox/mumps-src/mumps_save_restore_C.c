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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mumps_save_restore_C.h"
#include "mumps_common.h"
#if ! defined(NO_SAVE_RESTORE)
/* Functions */
void MUMPS_CALL
MUMPS_GET_SAVE_DIR_C(MUMPS_INT *len_save_dir, char* save_dir, mumps_ftnlen l1)
{
  char *tmp_save_dir;
  tmp_save_dir = getenv ("MUMPS_SAVE_DIR");
  if (tmp_save_dir==NULL)
    {
      tmp_save_dir = "NAME_NOT_INITIALIZED";
    }
  *len_save_dir = strlen(tmp_save_dir);
  save_dir = strncpy(save_dir, tmp_save_dir, l1); 
}
void MUMPS_CALL
MUMPS_GET_SAVE_PREFIX_C(MUMPS_INT *len_save_prefix, char* save_prefix, mumps_ftnlen l1)
{
  char *tmp_save_prefix;
  tmp_save_prefix = getenv ("MUMPS_SAVE_PREFIX");
  if (tmp_save_prefix==NULL)
    {
      tmp_save_prefix = "NAME_NOT_INITIALIZED";
    }
  *len_save_prefix = strlen(tmp_save_prefix);
  save_prefix = strncpy(save_prefix, tmp_save_prefix, l1); 
}
#else
void MUMPS_CALL
MUMPS_SAVE_RESTORE_RETURN_C()
{
}
#endif
