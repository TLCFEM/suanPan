/*******************************************************************************
 * Copyright (C) 2017-2025 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "FEAST.h"

#include <Domain/DomainBase.h>
#include <Domain/Factory.hpp>
#include <Solver/Integrator/Integrator.h>

extern "C" {
void new_dfeast_srci_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hrci_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* zAq, double* zSq, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_srcix_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hrcix_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* zAq, double* zSq, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_srci_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* qr, la_it* mode, double* resr, la_it* info);
void new_dfeast_grci_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_grci_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_srcix_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* qr, la_it* mode, double* resr, la_it* info, double* Zne, double* Wne);
void new_dfeast_grcix_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_grcix_(char* ijob, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);

void new_zfeast_grcipevx_(char* ijob, la_it* d, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_grcipev_(char* ijob, la_it* d, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_zfeast_srcipevx_(char* ijob, la_it* d, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_srcipev_(char* ijob, la_it* d, la_it* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_dfeast_sbgvx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_sbgv_(const char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_sbevx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_sbev_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hbgvx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hbgv_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hbevx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hbev_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gbgvx_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* klb, la_it* kub, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gbgv_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* klb, la_it* kub, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gbevx_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gbev_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_sbgvx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_sbgv_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* klb, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_sbevx_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_sbev_(char* UPLO, la_it* N, la_it* kla, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gbgvx_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* klb, la_it* kub, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gbgv_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* klb, la_it* kub, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gbevx_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gbev_(la_it* N, la_it* kla, la_it* kua, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_dfeast_sygvx_(char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_sygv_(const char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_syevx_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_syev_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hegvx_(char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hegv_(char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_heevx_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_heev_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gegvx_(la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gegv_(la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_geevx_(la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_geev_(la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_sygvx_(char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_sygv_(char* UPLO, la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_syevx_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_syev_(char* UPLO, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gegvx_(la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gegv_(la_it* N, double* A, la_it* LDA, double* B, la_it* LDB, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_geevx_(la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_geev_(la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_dfeast_sypevx_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_sypev_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hepevx_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hepev_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gepevx_(la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gepev_(la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_sypevx_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_sypev_(char* UPLO, la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gepevx_(la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gepev_(la_it* d, la_it* N, double* A, la_it* LDA, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_dfeast_scsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_scsrgv_(const char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_scsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_scsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hcsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hcsrgv_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hcsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hcsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gcsrgvx_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gcsrgv_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gcsrevx_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gcsrev_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_scsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_scsrgv_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_scsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_scsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gcsrgvx_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gcsrgv_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gcsrevx_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gcsrev_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_difeast_scsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_scsrgv_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_difeast_scsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_scsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_hcsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_hcsrgv_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_hcsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_hcsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emin, double* Emax, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_gcsrgvx_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_gcsrgv_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_gcsrevx_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_gcsrev_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_scsrgvx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_scsrgv_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_scsrevx_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_scsrev_(char* UPLO, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_difeast_gcsrgvx_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_gcsrgv_(la_it* N, double* sa, la_it* isa, la_it* jsa, double* sb, la_it* isb, la_it* jsb, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_difeast_gcsrevx_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_gcsrev_(la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_dfeast_scsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_scsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_hcsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_hcsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_gcsrpevx_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_gcsrpev_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zfeast_scsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zfeast_scsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_dfeast_gcsrpevx_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_dfeast_gcsrpev_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_difeast_scsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_scsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_hcsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_hcsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_gcsrpevx_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_gcsrpev_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_zifeast_scsrpevx_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_zifeast_scsrpev_(char* UPLO, la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);
void new_difeast_gcsrpevx_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info, double* Zne, double* Wne);
void new_difeast_gcsrpev_(la_it* d, la_it* N, double* sa, la_it* isa, la_it* jsa, la_it* fpm, double* epsout, la_it* loop, double* Emid, double* r, la_it* M0, double* lambda, double* q, la_it* mode, double* res, la_it* info);

void new_feastinit_(la_it* fpm);
void new_feastinit_driver_(la_it* fpm, la_it* N);
void new_cfeast_customcontour_(la_it* fpm2, la_it* N, la_it* Nedge, la_it* Tedge, float* Zedge, float* Zne, float* Wne);
void new_zfeast_customcontour_(la_it* fpm2, la_it* N, la_it* Nedge, la_it* Tedge, double* Zedge, double* Zne, double* Wne);
void new_zfeast_contour_(double* Emin, double* Emax, la_it* fpm2, la_it* fpm16, la_it* fpm18, double* Zne, double* Wne);
void new_cfeast_contour_(float* Emin, float* Emax, la_it* fpm2, la_it* fpm16, la_it* fpm18, float* Zne, float* Wne);
void new_zfeast_gcontour_(double* Emid, double* r, la_it* fpm2, la_it* fpm17, la_it* fpm19, double* Zne, double* Wne);
void new_dfeast_rational_(double* Emin, double* Emax, la_it* fpm2, la_it* fpm16, la_it* fpm18, double* Eig, la_it* M0, double* f);
void new_sfeast_rational_(float* Emin, float* Emax, la_it* fpm2, la_it* fpm16, la_it* fpm18, float* Eig, la_it* M0, float* f);
void new_sfeast_rationalx_(float* Zne, float* Wne, la_it* fpm2, float* Eig, la_it* M0, float* f);
void new_zfeast_grational_(double* Emid, double* r, la_it* fpm2, la_it* fpm17, la_it* fpm19, double* Eig, la_it* M0, double* f);
void new_cfeast_grational_(float* Emid, float* r, la_it* fpm2, la_it* fpm17, la_it* fpm19, float* Eig, la_it* M0, float* f);
void new_zfeast_grationalx_(double* Zne, double* Wne, la_it* fpm2, double* Eig, la_it* M0, double* f);
void new_cfeast_grationalx_(float* Zne, float* Wne, la_it* fpm2, float* Eig, la_it* M0, float* f);
void new_dfeast_rationalx_(double* Zne, double* Wne, la_it* fpm2, double* Eig, la_it* M0, double* f);
void new_cfeast_gcontour_(float* Emid, float* r, la_it* fpm2, la_it* fpm17, la_it* fpm19, float* Zne, float* Wne);
}

int FEAST::linear_solve(const shared_ptr<LongFactory>& W) const {
    auto& mass = W->get_mass();
    auto& stiffness = W->get_stiffness();

    std::vector fpm(64, la_it{0});

    new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
    fpm[0] = 1;
#endif

    std::vector output(4, la_it{0});
    std::vector input(4, 0.);
    input[1] = centre - radius; // centre
    input[2] = centre + radius; // radius

    output[1] = static_cast<la_it>(eigen_num);

    auto N = static_cast<la_it>(W->get_size());

    auto M = static_cast<la_it>(eigen_num);
    std::vector R(M, 0.);
    std::vector E(M, 0.);
    M *= N;
    std::vector X(M, 0.);

    if(const auto scheme = W->get_storage_scheme(); StorageScheme::FULL == scheme) new_dfeast_sygv_(&UPLO, &N, stiffness->memptr(), &N, mass->memptr(), &N, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    else if(StorageScheme::SPARSE == scheme || StorageScheme::SPARSESYMM == scheme) {
        auto fs = std::async([&] { return csr_form<double, la_it>(to_triplet_form<double, la_it>(stiffness), SparseBase::ONE); });
        auto fm = std::async([&] { return csr_form<double, la_it>(to_triplet_form<double, la_it>(mass), SparseBase::ONE); });

        auto t_stiff = fs.get();
        auto t_mass = fm.get();

        if(0 == t_mass.n_elem) {
            suanpan_error("The global mass matrix is empty, check if the element type supports mass and density settings.\n");
            return SUANPAN_FAIL;
        }

        new_dfeast_scsrgv_(&UPLO, &N, t_stiff.val_mem(), t_stiff.row_mem(), t_stiff.col_mem(), t_mass.val_mem(), t_mass.row_mem(), t_mass.col_mem(), fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    }
    else if(StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) {
        fpm[41] = 0;

        const auto [l, u] = W->get_bandwidth();
        auto KL = static_cast<la_it>(l);
        const auto KU = static_cast<la_it>(u);
        auto LD = KL + KU + 1;

        new_dfeast_sbgv_(&UPLO, &N, &KL, stiffness->memptr(), &LD, &KL, mass->memptr(), &LD, fpm.data(), &input[3], &output[0], &input[1], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);
    }
    else {
        suanpan_error("The current matrix storage scheme is not supported by the FEAST solver.\n");
        return SUANPAN_FAIL;
    }

    if(0 != output[3]) {
        suanpan_error("Error code {} received.\n", output[3]);
        return SUANPAN_FAIL;
    }

    W->modify_eigenvalue() = vec(E.data(), output[2]);
    W->modify_eigenvector() = mat(X.data(), N, output[2]);

    return SUANPAN_SUCCESS;
}

int FEAST::quadratic_solve(const shared_ptr<LongFactory>& W) const {
    std::vector fpm(64, la_it{0});

    new_feastinit_(fpm.data());

#ifdef SUANPAN_DEBUG
    fpm[0] = 1;
#endif

    auto N = static_cast<la_it>(W->get_size());

    std::vector output(4, la_it{0});
    std::vector input(4, 0.);
    input[0] = centre; // centre
    input[1] = 0.;     // centre
    input[2] = radius; // radius

    output[1] = static_cast<la_it>(eigen_num);

    auto M = 2 * static_cast<la_it>(eigen_num);
    std::vector R(M, 0.);
    std::vector E(M, 0.);
    M *= 2 * N;
    std::vector X(M, 0.);

    la_it P = 2;

    auto fk = std::async([&] { return csr_form<double, la_it>(to_triplet_form<double, la_it>(W->get_stiffness()), SparseBase::ONE); });
    auto fd = std::async([&] { return csr_form<double, la_it>(to_triplet_form<double, la_it>(W->get_damping()), SparseBase::ONE); });
    auto fm = std::async([&] { return csr_form<double, la_it>(to_triplet_form<double, la_it>(W->get_mass()), SparseBase::ONE); });

    const auto t_stiff = fk.get();
    const auto t_damping = fd.get();
    const auto t_mass = fm.get();

    size_t n_elem1 = std::max(std::max(t_stiff.n_elem, t_damping.n_elem), t_mass.n_elem);
    size_t n_elem2 = n_elem1 + n_elem1;
    size_t n_elem3 = n_elem1 + n_elem2;
    size_t n_size1 = N;
    size_t n_size2 = n_size1 + n_size1;
    size_t n_size3 = n_size2 + n_size1;

    std::vector A(n_elem3, 0.);
    std::vector JA(n_elem3, la_it{0});
    std::vector IA(n_size3, la_it{0});

    suanpan::for_each(t_stiff.n_elem, [&](const auto I) {
        A[I] = t_stiff.val_mem()[I];
        JA[I] = t_stiff.col_mem()[I];
    });
    suanpan::for_each(t_damping.n_elem, [&](const auto I) {
        A[I + n_elem1] = t_damping.val_mem()[I];
        JA[I + n_elem1] = t_damping.col_mem()[I];
    });
    suanpan::for_each(t_mass.n_elem, [&](const auto I) {
        A[I + n_elem2] = t_mass.val_mem()[I];
        JA[I + n_elem2] = t_mass.col_mem()[I];
    });

    suanpan::for_each(N, [&](const auto I) {
        JA[I] = t_stiff.row_mem()[I];
        JA[I + n_size1] = t_damping.row_mem()[I];
        JA[I + n_size2] = t_mass.row_mem()[I];
    });

    new_dfeast_gcsrpev_(&P, &N, A.data(), IA.data(), JA.data(), fpm.data(), &input[3], &output[0], &input[0], &input[2], &output[1], E.data(), X.data(), &output[2], R.data(), &output[3]);

    if(0 != output[3]) {
        suanpan_error("Error code {} received.\n", output[3]);
        return SUANPAN_FAIL;
    }

    auto& eigval = W->modify_eigenvalue();
    eigval.set_size(output[2]);

    for(uword I = 0; I < eigval.n_elem; ++I) eigval(I) = E[2 * I];

    auto& eigvec = W->modify_eigenvector();
    eigvec.resize(N, output[2]);

    for(uword I = 0; I < eigvec.n_elem; ++I) eigvec(I) = X[2 * I];

    return SUANPAN_SUCCESS;
}

FEAST::FEAST(const unsigned T, const unsigned N, const double C, const double R, const bool Q)
    : Solver(T)
    , quadratic(Q)
    , eigen_num(N)
    , centre(C)
    , radius(R) {}

int FEAST::initialize() {
    auto& G = get_integrator();

    if(nullptr == G) {
        suanpan_error("A valid integrator is required.\n");
        return SUANPAN_FAIL;
    }

    auto& W = G->get_domain()->get_factory();

    if(const auto scheme = W->get_storage_scheme(); StorageScheme::SYMMPACK == scheme) {
        suanpan_error("The symmetric pack storage is not supported.\n");

        return SUANPAN_FAIL;
    }
    else if((StorageScheme::BAND == scheme || StorageScheme::BANDSYMM == scheme) && !W->contain_solver_type(SolverType::SPIKE)) {
        suanpan_error("The SPIKE system solver is required for banded storage.\n");

        return SUANPAN_FAIL;
    }

    return SUANPAN_SUCCESS;
}

int FEAST::analyze() {
    auto& G = get_integrator();
    const auto D = G->get_domain();
    auto& W = D->get_factory();

    if(SUANPAN_SUCCESS != G->process_modifier()) return SUANPAN_FAIL;

    D->assemble_trial_mass();
    D->assemble_trial_stiffness();
    if(quadratic) D->assemble_trial_damping();

    // if(SUANPAN_SUCCESS != G->process_load()) return SUANPAN_FAIL;
    if(SUANPAN_SUCCESS != G->process_constraint()) return SUANPAN_FAIL;

    return quadratic ? quadratic_solve(W) : linear_solve(W);
}

void FEAST::print() {
    suanpan_info("An eigen solver using FEAST solver.\n");
}
