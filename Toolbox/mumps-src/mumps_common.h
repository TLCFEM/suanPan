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
#ifndef MUMPS_COMMON_H
#define MUMPS_COMMON_H
#include "mumps_compat.h"
#include "mumps_c_types.h"
/**
 * F_SYMBOL is a macro that converts a couple (lower case symbol, upper
 * case symbol) into the symbol defined by the compiler convention.
 * Example: For MUMPS_XXX, first define
 *   #define MUMPS_XXX F_SYMBOL(xxx,XXX) and then use
 *   MUMPS_XXX in the code to get rid of any symbol convention annoyance.
 *
 * NB: We need to provide both upper and lower case versions because to our
 *     knowledge, there is no way to perform the conversion with CPP
 *     directives only.
 */
#if defined(UPPER)
#define F_SYMBOL(lower_case, upper_case) MUMPS_##upper_case
#elif defined(Add_)
#define F_SYMBOL(lower_case, upper_case) mumps_##lower_case##_
#elif defined(Add__)
#define F_SYMBOL(lower_case, upper_case) mumps_##lower_case##__
#else
#define F_SYMBOL(lower_case, upper_case) mumps_##lower_case
#endif
MUMPS_INT *
mumps_get_mapping();
#define MUMPS_ASSIGN_MAPPING \
    F_SYMBOL(assign_mapping, ASSIGN_MAPPING)
void MUMPS_CALL
MUMPS_ASSIGN_MAPPING(MUMPS_INT *f77mapping);
#define MUMPS_NULLIFY_C_MAPPING F_SYMBOL(nullify_c_mapping, NULLIFY_C_MAPPING)
void MUMPS_CALL
MUMPS_NULLIFY_C_MAPPING();
MUMPS_INT *
mumps_get_pivnul_list();
#define MUMPS_ASSIGN_PIVNUL_LIST \
    F_SYMBOL(assign_pivnul_list, ASSIGN_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_ASSIGN_PIVNUL_LIST(MUMPS_INT *f77pivnul_list);
#define MUMPS_NULLIFY_C_PIVNUL_LIST \
    F_SYMBOL(nullify_c_pivnul_list, NULLIFY_C_PIVNUL_LIST)
void MUMPS_CALL
MUMPS_NULLIFY_C_PIVNUL_LIST();
MUMPS_INT *
mumps_get_uns_perm();
#define MUMPS_ASSIGN_UNS_PERM \
    F_SYMBOL(assign_uns_perm, ASSIGN_UNS_PERM)
void MUMPS_CALL
MUMPS_ASSIGN_UNS_PERM(MUMPS_INT *f77sym_perm);
#define MUMPS_NULLIFY_C_UNS_PERM \
    F_SYMBOL(nullify_c_uns_perm, NULLIFY_C_UNS_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_UNS_PERM();
MUMPS_INT *
mumps_get_sym_perm();
#define MUMPS_ASSIGN_SYM_PERM \
    F_SYMBOL(assign_sym_perm, ASSIGN_SYM_PERM)
void MUMPS_CALL
MUMPS_ASSIGN_SYM_PERM(MUMPS_INT *f77sym_perm);
#define MUMPS_NULLIFY_C_SYM_PERM \
    F_SYMBOL(nullify_c_sym_perm, NULLIFY_C_SYM_PERM)
void MUMPS_CALL
MUMPS_NULLIFY_C_SYM_PERM();
MUMPS_INT *
mumps_get_glob2loc_rhs();
#define MUMPS_ASSIGN_GLOB2LOC_RHS \
    F_SYMBOL(assign_glob2loc_rhs, ASSIGN_GLOB2LOC_RHS)
void MUMPS_CALL
MUMPS_ASSIGN_GLOB2LOC_RHS(MUMPS_INT *f77glob2loc_rhs);
#define MUMPS_NULLIFY_C_GLOB2LOC_RHS \
    F_SYMBOL(nullify_c_glob2loc_rhs, NULLIFY_C_GLOB2LOC_RHS)
void MUMPS_CALL
MUMPS_NULLIFY_C_GLOB2LOC_RHS();
MUMPS_INT *
mumps_get_glob2loc_sol();
#define MUMPS_ASSIGN_GLOB2LOC_SOL \
    F_SYMBOL(assign_glob2loc_sol, ASSIGN_GLOB2LOC_SOL)
void MUMPS_CALL
MUMPS_ASSIGN_GLOB2LOC_SOL(MUMPS_INT *f77glob2loc_sol);
#define MUMPS_NULLIFY_C_GLOB2LOC_SOL \
    F_SYMBOL(nullify_c_glob2loc_sol, NULLIFY_C_GLOB2LOC_SOL)
void MUMPS_CALL
MUMPS_NULLIFY_C_GLOB2LOC_SOL();
#define MUMPS_ICOPY_32TO64_64C_IP_C \
    F_SYMBOL(icopy_32to64_64c_ip_c, ICOPY_32TO64_64C_IP_C)
void MUMPS_CALL
MUMPS_ICOPY_32TO64_64C_IP_C(MUMPS_INT *inouttab, MUMPS_INT8 *sizetab);
#define MUMPS_ICOPY_64TO32_64C_IP_C \
    F_SYMBOL(icopy_64to32_64c_ip_c, ICOPY_64TO32_64C_IP_C)
void MUMPS_CALL
MUMPS_ICOPY_64to32_64C_IP_C(MUMPS_INT8 *inouttab, MUMPS_INT8 *sizetab);
#define MUMPS_MALLOC_C \
    F_SYMBOL(malloc_c, MALLOC_C)
void MUMPS_CALL
MUMPS_MALLOC_C(MUMPS_INT8 *address, MUMPS_INT8 *size);
#define MUMPS_FREE_C \
    F_SYMBOL(free_c, FREE_C)
void MUMPS_CALL
MUMPS_FREE_C(void *address);
#endif /* MUMPS_COMMON_H */
