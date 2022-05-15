!==========================================================================================
! Copyright (c) 2009-2018, The Regents of the University of Massachusetts, Amherst.
! E. Polizzi research lab
! All rights reserved.
! 
! Redistribution and use in source and binary forms, with or without modification, 
! are permitted provided that the following conditions are met:
! 
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions 
!    and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions 
!    and the following disclaimer in the documentation and/or other materials provided with the distribution.
! 3. Neither the name of the University nor the names of its contributors may be used to endorse or promote
!     products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
! BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING 
! IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!===========================================================================================

! Denotes the type for lapack-style calls
#define T d
! Types for declarations.
#define MAIN_TYPE double precision
#define REAL_TYPE double precision

! One, in the working precision 
#define ONE_PREC 1.0d0
#define ZERO_PREC 0.0d0
! The real part of one. 
#define ONE_REAL_PREC 1.0d0
#define ZERO_REAL_PREC 0.0d0
#define NZERO 1.0d-13
#define SPM_PREC_ENTRY 12

#undef SPIKECOMPLEX
! To decorate lapack/blas style function calls. 
! 'norm dec' is necessary because complex norm functions return reals, 
! and so they have names like dznrm2. 
#define T_DEC(A) glue_helper(d,A)
#define T_REAL_DEC(A) glue_helper(d,A)
#define T_NORM_DEC(A) glue_helper(d,A)
#define T_INORM_DEC(A) glue_helper(id,A)
#define T_INORM_REAL_DEC(A) glue_helper(id,A)

#define T_ALLOC(A) glue_helper(A,d)
#define T_REAL_ALLOC(A) glue_helper(A,d)
#define T_REAL_PRINT(A) glue_helper(A,d)
