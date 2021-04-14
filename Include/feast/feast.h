/*
!=========================================================================================
!Copyright (c) 2009-2019, The Regents of the University of Massachusetts, Amherst.
!Developed by E. Polizzi
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, 
!are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!   and the following disclaimer.
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!   and the following disclaimer in the documentation and/or other materials provided with the distribution.
!3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!    products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
!BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
!ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
!EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
!LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
!IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!==========================================================================================
*/

#ifndef FEAST_C_H
#define FEAST_C_H

#ifdef __cplusplus
extern "C" {
#endif
void dfeast_srci_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hrci_(char* ijob, int* N, double* Ze, double* work, double* workc, double* zAq, double* zSq, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_srcix_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hrcix_(char* ijob, int* N, double* Ze, double* work, double* workc, double* zAq, double* zSq, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_srci_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* qr, int* mode, double* resr, int* info);
void dfeast_grci_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_grci_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_srcix_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* qr, int* mode, double* resr, int* info, double* Zne, double* Wne);
void dfeast_grcix_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_grcix_(char* ijob, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);

void zfeast_grcipevx_(char* ijob, int* d, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_grcipev_(char* ijob, int* d, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void zfeast_srcipevx_(char* ijob, int* d, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_srcipev_(char* ijob, int* d, int* N, double* Ze, double* work, double* workc, double* Aq, double* Sq, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void dfeast_sbgvx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_sbgv_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_sbevx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_sbev_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hbgvx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hbgv_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hbevx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hbev_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gbgvx_(int* N, int* kla, int* kua, double* A, int* LDA, int* klb, int* kub, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gbgv_(int* N, int* kla, int* kua, double* A, int* LDA, int* klb, int* kub, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gbevx_(int* N, int* kla, int* kua, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gbev_(int* N, int* kla, int* kua, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_sbgvx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_sbgv_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* klb, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_sbevx_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_sbev_(char* UPLO, int* N, int* kla, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gbgvx_(int* N, int* kla, int* kua, double* A, int* LDA, int* klb, int* kub, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gbgv_(int* N, int* kla, int* kua, double* A, int* LDA, int* klb, int* kub, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gbevx_(int* N, int* kla, int* kua, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gbev_(int* N, int* kla, int* kua, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void dfeast_sygvx_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_sygv_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_syevx_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_syev_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hegvx_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hegv_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_heevx_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_heev_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gegvx_(int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gegv_(int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_geevx_(int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_geev_(int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_sygvx_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_sygv_(char* UPLO, int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_syevx_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_syev_(char* UPLO, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gegvx_(int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gegv_(int* N, double* A, int* LDA, double* B, int* LDB, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_geevx_(int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_geev_(int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void dfeast_sypevx_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_sypev_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hepevx_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hepev_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gepevx_(int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gepev_(int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_sypevx_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_sypev_(char* UPLO, int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gepevx_(int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gepev_(int* d, int* N, double* A, int* LDA, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void dfeast_scsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_scsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_scsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_scsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hcsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hcsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hcsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hcsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gcsrgvx_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gcsrgv_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gcsrevx_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gcsrev_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_scsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_scsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_scsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_scsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gcsrgvx_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gcsrgv_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gcsrevx_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gcsrev_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void difeast_scsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_scsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void difeast_scsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_scsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_hcsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_hcsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_hcsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_hcsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emin, double* Emax, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_gcsrgvx_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_gcsrgv_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_gcsrevx_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_gcsrev_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_scsrgvx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_scsrgv_(char* UPLO, int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_scsrevx_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_scsrev_(char* UPLO, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void difeast_gcsrgvx_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_gcsrgv_(int* N, double* sa, int* isa, int* jsa, double* sb, int* isb, int* jsb, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void difeast_gcsrevx_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_gcsrev_(int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

// non-linear polynomial
void dfeast_scsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_scsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_hcsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_hcsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_gcsrpevx_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_gcsrpev_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zfeast_scsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zfeast_scsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void dfeast_gcsrpevx_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void dfeast_gcsrpev_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void difeast_scsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_scsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_hcsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_hcsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_gcsrpevx_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_gcsrpev_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void zifeast_scsrpevx_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void zifeast_scsrpev_(char* UPLO, int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);
void difeast_gcsrpevx_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info, double* Zne, double* Wne);
void difeast_gcsrpev_(int* d, int* N, double* sa, int* isa, int* jsa, int* fpm, double* epsout, int* loop, double* Emid, double* r, int* M0, double* lambda, double* q, int* mode, double* res, int* info);

void feastinit_(int* fpm);
void feastinit_driver_(int* fpm, int* N);
void cfeast_customcontour_(int* fpm2, int* N, int* Nedge, int* Tedge, float* Zedge, float* Zne, float* Wne);
void zfeast_customcontour_(int* fpm2, int* N, int* Nedge, int* Tedge, double* Zedge, double* Zne, double* Wne);
void zfeast_contour_(double* Emin, double* Emax, int* fpm2, int* fpm16, int* fpm18, double* Zne, double* Wne);
void cfeast_contour_(float* Emin, float* Emax, int* fpm2, int* fpm16, int* fpm18, float* Zne, float* Wne);
void zfeast_gcontour_(double* Emid, double* r, int* fpm2, int* fpm17, int* fpm19, double* Zne, double* Wne);
void dfeast_rational_(double* Emin, double* Emax, int* fpm2, int* fpm16, int* fpm18, double* Eig, int* M0, double* f);
void sfeast_rational_(float* Emin, float* Emax, int* fpm2, int* fpm16, int* fpm18, float* Eig, int* M0, float* f);
void sfeast_rationalx_(float* Zne, float* Wne, int* fpm2, float* Eig, int* M0, float* f);
void zfeast_grational_(double* Emid, double* r, int* fpm2, int* fpm17, int* fpm19, double* Eig, int* M0, double* f);
void cfeast_grational_(float* Emid, float* r, int* fpm2, int* fpm17, int* fpm19, float* Eig, int* M0, float* f);
void zfeast_grationalx_(double* Zne, double* Wne, int* fpm2, double* Eig, int* M0, double* f);
void cfeast_grationalx_(float* Zne, float* Wne, int* fpm2, float* Eig, int* M0, float* f);
void dfeast_rationalx_(double* Zne, double* Wne, int* fpm2, double* Eig, int* M0, double* f);
void cfeast_gcontour_(float* Emid, float* r, int* fpm2, int* fpm17, int* fpm19, float* Zne, float* Wne);
#ifdef __cplusplus
}
#endif

#endif // FEAST_C_H
