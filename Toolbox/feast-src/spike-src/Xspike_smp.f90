
!=========================================================================================
!Copyright (c) 2009-2018, The Regents of the University of Massachusetts, Amherst.
!E. Polizzi research lab
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! SPIKE-SMP - REAL DOUBLE-PRECISION                !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! List of subroutines:

!!!!!! Documented !!!!!!!!

! Xspike_gbtrf
! Xspike_gbtrfp
! Xspike_gbtrs
! Xspike_gbtrsp
! Xspike_gbtrsi
! Xspike_gbsv
! Xspike_tune

!!!! Auxiliary !!!!!!!!!

! Xspike_gbtrf2
! Xspike_gbtrs2
! Xspike_gbtrst2
! Xspike_itrefinement
! Xspike_matmul
! Xspike_multi
! Xspike_invred
! Xspike_prep_recn
! Xspike_solve_recn
! Xspike_solve_recn_transpose
! Xspike_multi_transpose

! Xspike_GBTRFk2
! Xspike_GBTRSk2
! Xspike_simple_solve
! Xspike_simple_solve_transpose
! Xspike_vector_norm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine T_DEC(spike_gbtrf)(spm,n,kl,ku,A,lda,work,info)
! Purpose
! -------
! This subroutine performs the SPIKE DS factorization
! NOTE : This subroutine naturally destroys A, so if you 
! might want to perform iterative refinement later, or
! use A for something else, you should save it before
! calling this subroutine.
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! work (in/out)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
! 
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
  !=====================================================================
  !  Braegan Spring - Eric Polizzi - 2018
  !=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer, dimension(64) :: spm_default
integer :: n,kl,ku,lda,klu,info,i
MAIN_TYPE, dimension(*) :: work
MAIN_TYPE, dimension(*) :: A
integer, dimension(:), pointer :: Ajmin
integer :: infoloc
integer(8),parameter :: fout=6
integer :: ipiv_dummy

klu=max(kl,ku)

if (spm(1)==1) then
  call wwrite_n(fout)
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '*********** SPIKE-SMP -BEGIN ******************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout)
!!!!!!!!!!!! Print the Spike parameters which has been changed from default
  call wwrite_s(fout, 'List of input parameters spm(1:64)-- if different from default')
  call wwrite_n(fout)

  call spikeinit_default(spm_default)
  do i=1,19
      if ((i/=10).and.(spm(i)/=spm_default(i))) then
      call wwrite_s(fout, '   spm(')
      call wwrite_i(fout, i)
      call wwrite_s(fout, ')=')
      call wwrite_i(fout, spm(i))
      call wwrite_n(fout)
    endif
  enddo
  call wwrite_s(fout, 'Size system')  
  call wwrite_t(fout) 
  call wwrite_i(fout,N)
  call wwrite_n(fout)
  call wwrite_s(fout, 'kl,ku      ')  
  call wwrite_t(fout) 
  call wwrite_i(fout,kl)
  call wwrite_t(fout) 
  call wwrite_i(fout,ku)
  call wwrite_n(fout)
  call wwrite_s(fout, '#Threads (Total) available           ')  
  call wwrite_i(fout,spm(25))
  call wwrite_n(fout)
  if (spm(25)/=spm(22)) then
    call wwrite_s(fout, '**ATTENTION** --- SPIKE cannot use all the threads for this problem') 
    call wwrite_n(fout)
    call wwrite_s(fout, '#Threads effectively used by SPIKE   ') 
    call wwrite_i(fout,spm(22))
    call wwrite_n(fout)
  end if
  call wwrite_i(fout, spm(20) )
  call wwrite_s(fout, " Partitions and " )
  call wwrite_i(fout, spm(22) )
  call wwrite_s(fout, " Threads" )
  call wwrite_n(fout)
  call wwrite_s(fout, "Partition ratios: " )
  call wwrite_d(fout, spm(4) * .1d0)
  call wwrite_t(fout)
  call wwrite_d(fout, spm(5) * .1d0)
  call wwrite_n(fout)
  call wwrite_s(fout, "#Partitions with 1 and 2 threads" )
  call wwrite_n(fout)
  call wwrite_i(fout, spm(24) )
  call wwrite_t(fout)
  call wwrite_i(fout, spm(23) )
  call wwrite_n(fout)

  call wallocate_1i(Ajmin,spm(20)+1,infoloc)
  call spikerl_calc_size_partitions(spm,Ajmin,n)

  call wwrite_s(fout, "Partition sizes: " )
  do i=1,spm(20)
    call wwrite_i(fout, Ajmin(i+1) - Ajmin(i))
    call wwrite_t(fout)
  enddo
  call wwrite_n(fout)

end if

if( .not. ((spm(3) .eq. 0) .or. (spm(3) .eq. 2))) then
  spm(3) = 0
endif
call T_DEC(spike_gbtrf2)&
(spm,n,kl,ku,A,lda,work(1),work(2*klu*klu*spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+1),ipiv_dummy,info)

if (spm(1)==1) call wdeallocate_1i(Ajmin)


end subroutine T_DEC(spike_gbtrf)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine T_DEC(spike_gbtrfp)(spm,n,kl,ku,A,lda,work,ipiv,info)
! Purpose
! -------
! This subroutine performs the SPIKE DS factorization with pivoting
! NOTE : This subroutine naturally destroys A, so if you 
! might want to perform iterative refinement later, or
! use A for something else, you should save it before
! calling this subroutine.
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! work (in/out)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
! 
! ipiv (out) 
! The pivoting array, which represents the permutations applied during partial pivoting.
! Size is n.
!
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
  !=====================================================================
  !  Braegan Spring - Eric Polizzi - 2018
  !=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer, dimension(64) :: spm_default
integer :: n,kl,ku,lda,klu,info,i
MAIN_TYPE, dimension(*) :: work
MAIN_TYPE, dimension(*) :: A
integer, dimension(:), pointer :: Ajmin
integer :: infoloc
integer(8),parameter :: fout=6
integer, dimension(*) :: ipiv

klu=max(kl,ku)

if (spm(1)==1) then
  call wwrite_n(fout)
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '********** SPIKE-SMP-Ppivoting Begin **********')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout)
!!!!!!!!!!!! Print the Spike parameters which has been changed from default
  call wwrite_s(fout, 'List of input parameters spm(1:64)-- if different from default')
  call wwrite_n(fout)

  call spikeinit_default(spm_default)
  do i=1,19
      if ((i/=10).and.(spm(i)/=spm_default(i))) then
      call wwrite_s(fout, '   spm(')
      call wwrite_i(fout, i)
      call wwrite_s(fout, ')=')
      call wwrite_i(fout, spm(i))
      call wwrite_n(fout)
    endif
  enddo
  call wwrite_s(fout, 'Size system')  
  call wwrite_t(fout) 
  call wwrite_i(fout,N)
  call wwrite_n(fout)
  call wwrite_s(fout, 'kl,ku      ')  
  call wwrite_t(fout) 
  call wwrite_i(fout,kl)
  call wwrite_t(fout) 
  call wwrite_i(fout,ku)
  call wwrite_n(fout)
  call wwrite_s(fout, '#Threads (Total) available           ')  
  call wwrite_i(fout,spm(25))
  call wwrite_n(fout)
  if (spm(25)/=spm(22)) then
    call wwrite_s(fout, '**ATTENTION** --- SPIKE cannot use all the threads for this problem') 
    call wwrite_n(fout)
    call wwrite_s(fout, '#Threads effectively used by SPIKE   ') 
    call wwrite_i(fout,spm(22))
    call wwrite_n(fout)
  end if
  call wwrite_i(fout, spm(20) )
  call wwrite_s(fout, " Partitions and " )
  call wwrite_i(fout, spm(22) )
  call wwrite_s(fout, " Threads" )
  call wwrite_n(fout)
  call wwrite_s(fout, "Partition ratios: " )
  call wwrite_d(fout, spm(4) * .1d0)
  call wwrite_t(fout)
  call wwrite_d(fout, spm(5) * .1d0)
  call wwrite_n(fout)
  call wwrite_s(fout, "#Partitions with 1 and 2 threads" )
  call wwrite_n(fout)
  call wwrite_i(fout, spm(24) )
  call wwrite_t(fout)
  call wwrite_i(fout, spm(23) )
  call wwrite_n(fout)

  call wallocate_1i(Ajmin,spm(20)+1,infoloc)
  call spikerl_calc_size_partitions(spm,Ajmin,n)

  call wwrite_s(fout, "Partition sizes: " )
  do i=1,spm(20)
    call wwrite_i(fout, Ajmin(i+1) - Ajmin(i))
    call wwrite_t(fout)
  enddo
  call wwrite_n(fout)

end if

spm(3) = 1
call T_DEC(spike_gbtrf2)(spm,n,kl,ku,A,lda,work(1),work(2*klu*klu*spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),ipiv,info)

if (spm(1)==1) call wdeallocate_1i(Ajmin)


end subroutine T_DEC(spike_gbtrfp)



subroutine T_DEC(spike_gbtrsp) &
(spm,n,kl,ku,nrhs,A,lda,work,ipiv,B,ldb)
! Purpose
! -------
! This subroutine performs the SPIKE solve
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The number of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! A (in) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
!  
! work (in)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
!
! ipiv (in) 
! The pivoting array, which represents the permutations applied during partial pivoting.
! Size is n.
!
! B (in/out)
! The collection of vectors on which the solve is performed
!
! ldb (in)
! The leading dimension of B. Usually this is n, but if you feel clever you can try other things. 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs
integer,  dimension(*) :: ipiv
MAIN_TYPE, dimension(*) :: work
MAIN_TYPE, dimension(lda,*) :: A
MAIN_TYPE, dimension(ldb,*) :: B
integer(8),parameter :: fout=6

klu=max(kl,ku)

spm(3) = 1
call T_DEC(spike_gbtrs2)&
(spm,n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),ipiv,B,ldb) 


if((spm(1) .eq. 1) .and. (spm(31) .eq. 0)) then
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '*********** SPIKE-SMP Pivoting END ************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_n(fout)
endif

end subroutine T_DEC(spike_gbtrsp)



subroutine T_DEC(spike_gbtrs) &
(spm,trans,n,kl,ku,nrhs,A,lda,work,B,ldb)
! Purpose
! -------
! This subroutine performs the SPIKE solve
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The number of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! A (in) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
!  
! work (in)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
!
! B (in/out)
! The collection of vectors on which the solve is performed
!
! ldb (in)
! The leading dimension of B. Usually this is n, but if you feel clever you can try other things. 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs
MAIN_TYPE, dimension(*) :: work
MAIN_TYPE, dimension(lda,*) :: A
MAIN_TYPE, dimension(ldb,*) :: B
character :: trans
integer(8),parameter :: fout=6
integer :: ipiv

klu=max(kl,ku)

spm(3) = 0
if((trans .eq. 'n') .or. (trans .eq. 'N')) then
  call T_DEC(spike_gbtrs2) &
(spm,n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),ipiv,B,ldb)
else
 ! call T_DEC(spike_gbtrst2)(spm,trans,n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  if((trans .eq. 't') .or. (trans .eq. 'T')) then
    call T_DEC(spike_gbtrst2) &
(spm,'T',n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  else
    call T_DEC(spike_gbtrst2) &
(spm,'C',n,kl,ku,lda,nrhs,A,work(1),work(2*klu*klu* spm(21)*spm(20) + 1),work(2*2*klu*klu*spm(21)*spm(20)+ 1),B,ldb)
  endif
endif


if((spm(1) .eq. 1) .and. (spm(31) .eq. 0)) then
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '************SPIKE-SMP END *********************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_n(fout)
endif

end subroutine T_DEC(spike_gbtrs)



subroutine T_DEC(spike_gbtrsi)(spm,trans,n,kl,ku,nrhs,C,ldc,A,lda,work,B,ldb)
! Purpose
! -------
! This subroutine performs the SPIKE solve with iterative refinement
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The number of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! C (in)
! A copy of the origional version of A, used for iterative refinement
! 
! ldc (in) 
! The leading dimensions of C. ldc=kl+ku+1
!
! A (in) 
! The banded matrix on which the SPIKE factorization was performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
!  
! work (in)
! This is a large work array, which contains parts of the reduced system which must be communicated to the solve stage
! The required size is (klu^2)*4*(#levels*#partitions)+(#partitions-1))
! The value 4*(#levels*#partitions)+(#partitions-1)) is placed in spm(10) by spikeinit.
! So, you may instead allocate work to be (klu*klu)*spm(10)
!
! B (in/out)
! The collection of vectors on which the solve is performed
!
! ldb (in)
! The leading dimension of B. Usually this is n, but if you feel clever you can try other things. 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,klu,ldb,nrhs,ldc
MAIN_TYPE, dimension(*) :: work
MAIN_TYPE, dimension(lda,*) :: A
MAIN_TYPE, dimension(ldb,*) :: B
MAIN_TYPE, dimension(ldc,*) :: C
character :: trans
MAIN_TYPE, dimension(:,:), pointer :: oB
integer :: infoloc
integer(8),parameter :: fout=6
!This is just to make the arguments work out for T_DEC(spike_gbtrs)
call T_ALLOC(wallocate_2)(oB,ldb,nrhs,infoloc)
call T_DEC(LACPY)( 'F',n,nrhs,B, N, oB, N )

klu=max(kl,ku)

!Spike solve
spm(31) = 1
call T_DEC(spike_gbtrs)(spm,trans,n,kl,ku,nrhs,A,lda,work,B,ldb)

!Applying iterative refinement
call T_DEC(spike_itrefinement)(spm,trans,n,kl,ku,nrhs,A,work,lda,C,ldc,B,oB,ldb)

call T_ALLOC(wdeallocate_2)(oB)
spm(31) = 0
if((spm(1) .eq. 1) .and. (spm(31) .eq. 0)) then
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_s(fout, '************SPIKE-SMP END *********************')
  call wwrite_n(fout) 
  call wwrite_s(fout, '***********************************************')  
  call wwrite_n(fout) 
  call wwrite_n(fout)
endif


end subroutine T_DEC(spike_gbtrsi)




subroutine T_DEC(spike_gbsv)(spm,N,KL,KU,NRHS,A,LDA,B,LDB,INFO)
! Purpose
! -------
! This is the do-it-all factorize and solve subroutine. 
! A matrix and a collection of vectors are entered, and the 
! problem  A^(-1)B => B is solved. 
! The residual value for B is checked, and simple iterative 
! refinement is performed if it is above the given tolerance.
!
! Overall, this is similar to the lapack Xgbsv, with the exception that
! the matrix A is returned to the initial state upon return.
! 
! Optimization note: 
! If the optimization flag spm(2)=2 is set, this subroutine 
! will attempt to find appropriate ratios for the partition sizes. 
! These ratios are dependent on the characteristics of the matrix used 
! (specificly the relationship between klu and nrhs) and the hardware 
! on which SPIKE is run. The hardware characteristics are embodied
! in the tuning constant K=spm(7).
! This subroutine may attempt to find K, but finding K requires a 
! large single threaded factorization and solve. To avoid this cost,
! it will check if there is a non-zero value in spm(7), which 
! indicates that the constant has already been found.
! Because this constant should only depend on the the hardware, 
! one may manually set the tuning constant by setting spm(7) after calling
! spikeinit. 
! (dspike_tune will help with this process)
! This need not be run every time Xspike_gbsv is called, or even every
! time the calling program is run (although communicating the constant
! between program runs is an exercise left to the reader [I don't want
! to muck up your environment variables]).  
! 
! -------
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! nrhs (in)
! The number of columns in B
! 
! A (in) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! B (in/out) 
! The set of vectors on which the matrix solve will be performed; solution on exit
! 
! ldb (in)
! The leading dimensions of B
! 
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer :: n,kl,ku,klu,nrhs
integer :: ldc,lda,ldb,info,i,p
MAIN_TYPE, dimension(lda,*) :: A
MAIN_TYPE, dimension(ldb,*) :: B
MAIN_TYPE, dimension(:,:),pointer :: C
integer,  dimension(:),pointer :: ajmin,ipiv
integer, dimension(64) :: spm
character :: trans
MAIN_TYPE,  dimension(:),pointer :: work
double precision :: K,klud,nrhsd
integer :: infoloc
integer(8),parameter :: fout=6

spm(31) = 1 !Indicate that factorize and solve are done in one call.
klu = max(kl,ku)

trans='n'
klud = 1.0d0*klu 
nrhsd = 1.0d0*nrhs
if((spm(2) .eq. 2)) then
  K = 0.1d0*spm(7)
  spm(4) = int(10*(klud/(klud+K*nrhsd)*.5d0 + (      klud +       nrhsd)/(klud/K + nrhsd)))
  spm(5) = int(10*(klud/(klud+K*nrhsd)      + (1.5d0*klud + 2.0d0*nrhsd)/(klud/K + nrhsd)))
endif
call T_ALLOC(wallocate_1)(work,klu*klu*spm(10),infoloc)
call wallocate_1i(ipiv,n,infoloc)

if(spm(3) .eq. 1) then
  ldc=klu+kl+ku+1
else
  ldc=kl+ku+1
endif

!Give the system a hint that A and  C should be 'owned' by the various processors by copying them in a parallel section
call T_ALLOC(wallocate_2)(C,ldc,n,infoloc)
call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call spikerl_calc_size_partitions(spm,Ajmin,n)

!$OMP PARALLEL DO default(shared) private(p)
do i=1,spm(20)
  p=i
  call T_DEC(LACPY)('A',lda,Ajmin(p+1)-Ajmin(p),A(1,Ajmin(p)),lda,C(ldc-(kl+ku),Ajmin(p)),ldc)
enddo

if(spm(3) .eq. 1) then
  call T_DEC(spike_gbtrfp)(spm,n,kl,ku,C,ldc,work,ipiv,info)
  call T_DEC(spike_gbtrsp)(spm,n,kl,ku,nrhs,C,ldc,work,ipiv,B,ldb)
else
  call T_DEC(spike_gbtrf)(spm,n,kl,ku,C,ldc,work,info)
  if(spm(3) .eq. 0) then
    call T_DEC(spike_gbtrs)(spm,trans,n,kl,ku,nrhs,C,ldc,work,B,ldb)
  else
    call T_DEC(spike_gbtrsi)(spm,trans,n,kl,ku,nrhs,A,lda,C,lda,work,B,ldb)
  endif
endif

call T_ALLOC(wdeallocate_2)(C)
call wdeallocate_1i(Ajmin)

call wdeallocate_1i(ipiv)
call T_ALLOC(wdeallocate_1)(work)

end subroutine T_DEC(spike_gbsv) 





subroutine T_DEC(spike_tune)(spm)
! Purpose
! -------
! Auto-tuning- Calculate the partition size ratios, for COMPLE DOUBLE PRECISION values
! These ratios are described by a tuning constant, which is dependent on the Big-O behavior of the system banded factorization and solve.
! Discovering this constant requires a large, single threaded fatorize and solve.
! As a result, this function is fairly slow. 
! -------
!
! Arguments  
! ---------
! spm (out)
! An array of 64 integers, which control various features of the spike factorization and solve.
! These parameters are described in the description of spikeinit
! For this function we only modify the following values:
!
! spm(4)  = out. This describes the ratio between the size of the large first and last partitions, and the medium 2-thread spike partitions. It is given as 10X the ratio
! spm(5)  = out. This describes the ratio between the size of the large first and last partitions, and the small single-thread LU partitions. It is given as 10X the ratio
! spm(7)  = out. This describes the tuning constant from which these ratios are derived. It is given as 10X the tuning constant. 
! The derived relationship between these values should be:
! spm(4) = spm(7) + 5
! spm(5) = 2*spm(4)
! ---------
!=====================================================================
!  Braegan Spring -  2018
!=====================================================================
use omp_lib
implicit none 
include 'f90_noruntime_interface.fi'
integer n,klu,lda,info
double precision :: t1,t2,t0
integer, dimension(64) :: spm
MAIN_TYPE, dimension(:,:),pointer :: A
MAIN_TYPE, dimension(:,:),pointer :: B
REAL_TYPE :: nzero, norma
integer, dimension(:),    pointer                :: ipiv
double precision :: K
integer :: infoloc
integer(8),parameter :: fout=6
!double precision, external :: OMP_GET_WTIME
logical :: pivot

pivot = spm(3) .eq. 1

n=10000
klu=200
nzero=NZERO

norma=4*klu*ONE_REAL_PREC
if(pivot) then
  lda = 3*klu+1 
else
  lda = 2*klu+1 
endif

call T_ALLOC(wallocate_2)(A,lda,n,infoloc)
call T_ALLOC(wallocate_2)(B,n,klu,infoloc)
call wallocate_1i(ipiv,n,infoloc)

B=ONE_PREC
A=ONE_PREC
A(klu+1,:) = 3*klu*ONE_PREC

if(pivot) then
  t0 = OMP_GET_WTIME()
  call T_DEC(GBTRFUL)(n,n,klu,klu,A,lda,ipiv,info)
  t1 = OMP_GET_WTIME()
!  call DGBTRSL('N',n,klu,klu,klu,A,LDA,ipiv,B,n,info)
  call T_DEC(GBTRSL)('N',n,klu,klu,klu,A,LDA,ipiv,B,n,info)
  call T_DEC(GBTRSU)('N',n,klu,klu,klu,A,LDA,ipiv,B,n,info)
!  call DGBTRS('N',n,klu,klu,klu,A,LDA,ipiv,B,n,info)
  t2 = OMP_GET_WTIME()
  K = 1.0d0*(t2-t1)/(t1-t0)
else
  t0 = OMP_GET_WTIME()
  call T_DEC(GBALU)(n,klu,klu,A,lda,nzero,norma,info)
  t1 = OMP_GET_WTIME()
  call T_DEC(TBSM)('L','N','U',n,klu,klu,A,lda,B,n)
  call T_DEC(TBSM)('U','N','N',n,klu,klu,A,lda,B,n)
  t2 = OMP_GET_WTIME()
  K = 1.0d0*(t2-t1)/(t1-t0)
endif

spm(7) = int(K*10)
!spm(4) = spm(7) +5
spm(5) = spm(7) +spm(7)/2 +10
spm(4) = spm(5)/2

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'Large-bandwith partition ratios found:')
  call wwrite_d(fout, spm(4)*.1d0)
  call wwrite_t(fout)
  call wwrite_d(fout, spm(5)*.1d0)
  call wwrite_n(fout)
endif

call T_ALLOC(wdeallocate_2)(A)
call T_ALLOC(wdeallocate_2)(B)

end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NOTE: Most of the functions below here aren't really intended to be called by an end user.  !
!                          They are not documented                                            !
!                                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





subroutine T_DEC(spike_gbtrf2)(spm,n,kl,ku,A,lda,rV,rW,red,ipiv,info)
! Purpose
! -------
! This is the code that ultimately does the work of the
! SPIKE factorization. However, there is a more user friendly 
! interface. T_DEC(spike_gbtrf) should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (out) 
! Work space, to contain various parts of the reduced system
! Must be communicated to the solve subroutine
!
! info (out)
! Information parameter
! info = 0 -> The banded LU factorization on each partition did not require boosting
! info = 1 -> The banded LU factorization on some partition required boosting. The answer will be approximate. Iterative refinement may be necessary.
! info = 2 -> There was an illegal value in your A matrix.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm 
logical ::invred
integer :: tim
integer :: n,kl,ku,lda,nbpart,nblevel,nthread
integer :: i,j,klu,p,p_adjust,info,max_info
integer :: counter
!integer,external :: omp_get_thread_num
integer, dimension(:), pointer :: Ajmin
!integer, dimension(:,:), pointer :: keys
integer, dimension(4,spm(23)) :: keys
MAIN_TYPE, dimension(lda,*):: A
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku),*):: red
MAIN_TYPE, dimension(2*max(kl,ku),max(kl,ku),spm(21),*):: rV,rW
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*(kl+ku)) :: spike2_grj
MAIN_TYPE, pointer, dimension(:,:) :: spike2_red
MAIN_TYPE, pointer, dimension(:,:) :: spike2_grj
REAL_TYPE :: nzero, norma
MAIN_TYPE,dimension(:,:),pointer::vw
MAIN_TYPE,dimension(:,:),pointer::vw2
double precision, dimension(:),pointer :: timing_array
integer :: infoloc
integer(8),parameter :: fout=6
!double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time, parallel_end_time
character(len=1) :: norm
REAL_TYPE, dimension(:), pointer :: normwork
REAL_TYPE, external :: T_DEC(langb)
logical :: pivot
integer, dimension(*) :: ipiv
integer, dimension(:), pointer :: ipiv_aux
integer :: kd, ibstart, iastartl,iastartu,icstart,mysize
double precision :: start_time_fac, end_time_fac
!character(len=100) :: format
!write(format,'(A1,i2,A6)') '(', 1, 'F6.2)'

pivot = spm(3) .eq. 1

start_time = OMP_GET_WTIME()

if(spm(14) .eq. 0) then
  !Norm infinity
  norm='i'
endif

if(spm(14) .eq. 1) then
  !Norm one
  norm='o'
endif

if(spm(14) .eq. 2) then
  !Norm Frobenius
  norm='f'
endif

call wallocate_1d(timing_array,spm(22),infoloc)

nzero=NZERO
!If we've only been given one thread, do the banded primitive and get out of here.
if(spm(20) == 1) then
  parallel_start_time = OMP_GET_WTIME()
  timing_array(1) = OMP_GET_WTIME()

  if(pivot) then

      klu=max(kl,ku)
      iastartu = klu-kl+1
      call T_DEC(GBTRF)(n,n,kl,ku,A(iastartu,1),lda,ipiv,info)
      !    call T_DEC(GBTRF)(n,n,kl,ku,A,lda,ipiv,info)
  else 
    ! Lapack norm compute requires a work array for infinorm
    if(spm(14) .eq. 0) then
      call T_REAL_ALLOC(wallocate_1)(normwork,n,infoloc)
    else
      !not actually used but passing an unallocated array seems to be a problem 
      call T_REAL_ALLOC(wallocate_1)(normwork,1,infoloc)
    endif
    norma = T_DEC(langb)( norm, n, kl, ku, A, lda, normwork )
    call T_REAL_ALLOC(wdeallocate_1)(normwork)  
    call T_DEC(GBALU)(n,kl,ku,A,lda,nzero,norma,info)
 
  endif

  timing_array(1) = OMP_GET_WTIME() - timing_array(1)
  parallel_end_time = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time   = 0.0d0
else

call T_ALLOC(wallocate_2)(spike2_red, 2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call T_ALLOC(wallocate_2)(spike2_grj, 2*max(kl,ku),spm(23)*(kl+ku),infoloc)

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call spikerl_calc_size_partitions(spm,Ajmin,n)

!Ajmin holds the positions of the sub-arrays. 
nthread = spm(22)
nbpart = spm(20)
nblevel = spm(21)
invred=.false.
klu=max(kl,ku)

!call wallocate_2i(keys,4,spm(23),infoloc)
keys(1,1:spm(23)) = 0
!keys(2:4,1:spm(23)) = -1

! the keys is a shared section of memory that the simplified two-threaded spike partitions will use to communicate 
! when they have finished various steps of their work. 
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

if(spm(23) > 0) then
  call T_ALLOC(wallocate_2)(vw2,Ajmin(3)    -Ajmin(2)       ,(ku+kl)*spm(23),infoloc) 
endif 
if(spm(20)-spm(23) > 0) then
 ! call T_ALLOC(wallocate_2)(vw,Ajmin(nbpart)-Ajmin(nbpart-1),(ku+kl)*spm(20),infoloc) 
  call T_ALLOC(wallocate_2)(vw,Ajmin(nbpart)-Ajmin(nbpart-1),(ku+kl)*(spm(24)-2),infoloc) 
endif

info = 0
max_info=0
parallel_start_time=OMP_GET_WTIME()

! For the pivoting code, We need to pull out the B and C matrices before performing the factorization
! because they will be moved aruond by the swapping code for the UL factorization otherwise.
if(pivot) then
  !Copy the b matrices over
  !$omp parallel do ordered default(none) private(ibstart, icstart, p, mysize) shared(Ajmin,rV,rW,A,keys,nbpart,nthread,spm,kl,ku,klu,lda)
  do i=1,nthread,1

    ibstart = klu+1
    icstart = klu+1+kl+ku
  
    call threadmap(i,p,nbpart,spm(23))
  
  !  rV(:,:,1,p) = ZERO_PREC
    !$omp ordered
    !If we're outside the range of threads that use the two-thread spike, 
    if((p==1) .or. (p > 1 + spm(23) .and. p < nbpart)) then
      call T_DEC(LACPY)( 'L',ku , ku, A(ibstart,Ajmin(p+1)) , lda-1,A(ibstart,Ajmin(p)), lda-1)
    endif
  
    if(p > 1 .and. p <= 1 + spm(23)) then
      !or if we are the appropriate thread in the two-thread range...
      if(keys(3,p-1) == omp_get_thread_num()) then
        !Outside B
        call T_DEC(LACPY)( 'L',ku , ku, A(ibstart,Ajmin(p+1)) , lda-1,A(ibstart,Ajmin(p)), lda-1)
        !Inside  B
        call T_DEC(LACPY)( 'L',ku , ku, A(ibstart,(Ajmin(p)+Ajmin(p+1))/2) , lda-1,A(1,Ajmin(p)), lda-1)
      endif
    endif
    !$omp end ordered
  
  enddo
  !Copy the C matrices over
  !TODO Unfortunately, this loop is almost exactly backwards from how we'd like it to be -- this will probably hurt performance, I'll have to fix it later
  !$omp parallel do ordered default(none) private(ibstart, icstart, p, mysize) shared(Ajmin,rV,rW,A,keys,nbpart,nthread,spm,kl,ku,klu,lda)
  do i=nthread,1,-1
  
    ibstart = klu+1
    icstart = klu+1+kl+ku
  
    call threadmap(i,p,nbpart,spm(23))
  
  !  rW(:,:,1,p) = ZERO_PREC
    !$omp ordered
    !If we're outside the range of threads that use the two-thread spike, 
    if((p > 1 + spm(23) .and. p < nbpart) .or. p==nbpart) then 
      call T_DEC(LACPY)( 'U', kl, kl, A(icstart,Ajmin(p)-kl), lda-1,A(icstart,Ajmin(p+1)-kl), lda-1)
    endif
  
    if(p > 1 .and. p <= 1 + spm(23)) then
      !or if we are the appropriate thread in that range...
  
      if(keys(2,p-1) == i-1) then
        !Inside C
        call T_DEC(LACPY)( 'U', kl, kl, A(icstart,(Ajmin(p)+Ajmin(p+1))/2-kl), lda-1,A(icstart,Ajmin(p+1)-kl), lda-1)
   
        !Outside C
        call T_DEC(LACPY)( 'U', kl, kl, A(icstart,Ajmin(p)-kl), lda-1,A(icstart,(Ajmin(p)+Ajmin(p+1))/2-kl), lda-1)
      endif
    endif
    !$omp end ordered
  
  enddo
endif

!$omp parallel do default(firstprivate) private(kd,ipiv_aux, iastartl, iastartu, ibstart, icstart,i,j,p,p_adjust,tim,counter,info,norma,infoloc,normwork,mysize,start_time_fac,end_time_fac) shared(pivot,nblevel,vw,vw2,spike2_red,spike2_grj,keys,A,rV,rW,spm,ajmin,nthread,max_info,timing_array,lda, ipiv, nbpart) firstprivate(nzero,norm,klu,kl,ku)
do i=1,nthread,1
  timing_array(i) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
  info = 0 
! <<<<<<<<<<< LU Factorize A(1) - A(nbpart-1) and Solving Ai.Vi=bi, A2..nbpart-1.W2..nbpart-1=c2..nbpart-1 >>>>>>>>>>> !

  if(pivot) then
    call wallocate_1i(ipiv_aux,2*klu,infoloc)
    ibstart  = klu+1
    icstart  = klu+1+kl+ku
    iastartl = klu-kl+1
  else   
  ! Lapack norm compute requires a work array for infinorm
    if(spm(14) .eq. 0) then
      call T_REAL_ALLOC(wallocate_1)(normwork,Ajmin(p+1)-Ajmin(p),infoloc)
    else
      !not actually used but passing an unallocated array seems to be a problem 
      call T_REAL_ALLOC(wallocate_1)(normwork,1,infoloc)
    endif
    norma = T_DEC(langb)( norm, Ajmin(p+1)-Ajmin(p), kl, ku, A(1,Ajmin(p)), lda, normwork )
    call T_REAL_ALLOC(wdeallocate_1)(normwork)  
    ibstart  = 1
    icstart  = 1+kl+ku
    iastartl = ku+1
    iastartu = 1
    kd = ku
  endif

  mysize = Ajmin(p+1)-Ajmin(p)

  if(p == 1) then
!      call sleep(1)
    rV(:,:,1,p) = ZERO_PREC

    if(pivot) then
      iastartu = klu-kl+1
      kd = kl+ku
      START_TIME_FAC = OMP_GET_WTIME()
      call T_DEC(GBTRF)(mysize,mysize,kl,ku,A(iastartl,Ajmin(p)),lda,ipiv(Ajmin(p)),info)
      END_TIME_FAC = OMP_GET_WTIME()

      ipiv_aux(1:2*klu) = ipiv(Ajmin(p+1)-2*klu : Ajmin(p+1)-1) - (mysize - 2*klu)
      call T_DEC(LACPY)( 'L',ku , ku, A(ibstart,Ajmin(p)), lda-1, rV(2*klu-ku+1,1,1,p), 2*klu)
      call T_DEC(GBTRSL)('N',2*KLU,KL,KU,KLU,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux,rV(1,1,1,p),2*klu,INFO)
    else
      call T_DEC(LACPY)( 'L',ku , ku, A(ibstart,Ajmin(p+1)), lda-1,rV(2*klu-ku+1,1,1,p), 2*klu)
      call T_DEC(GBALU)(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,nzero,norma,info)
      call T_DEC(TBSM)('L','N','U', klu, ku, kl, A(iastartl,Ajmin(p+1)-klu),lda, rV(klu+1,1,1,p), 2*klu)
    endif

! I need the L^-1b for the transpose case (so save it in an unused portion of rV) 
! This section is unused because there is no top of the V spike for the first partiton
    call T_DEC(LACPY)( 'L',klu , klu, rV(2*klu-ku+1,1,1,p), 2*klu, rV(1,1,1,p), 2*klu)
    call T_DEC(TBSM)('U','N','N', klu, ku, kd, A(iastartu,  Ajmin(p+1)-klu),lda, rV(klu+1,1,1,p), 2*klu)
  endif  

  if(p==nbpart) then
    iastartl = klu-ku+1
    rW(:,:,1,p) = ZERO_PREC

    if(pivot) then
      call T_DEC(GBTRFUL)(mysize,mysize,kl,ku,A(1,Ajmin(p)),lda,ipiv(Ajmin(p)),info)
      ipiv_aux(1:2*klu) = ipiv(Ajmin(p+1)-2*klu : Ajmin(p+1)-1) - (mysize - 2*klu)
      call T_DEC(LACPY)( 'L',kl , kl, A(ibstart,Ajmin(p)), lda-1,rW(2*klu-kl+1,1,1,p), 2*klu)
      call T_DEC(GBTRSL)('N',2*KLU,KU,KL,KLU,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux,rW(1,1,1,p),2*klu,INFO)
      call T_DEC(LACPY)( 'L', kl, kl, rW(2*klu+1-kl,1,1,p), 2*klu, rW(klu+1-kl,1,1,p), 2*klu)
      call T_DEC(GBTRSU)('N',KLU,KU,KL,KLU,A(iastartl,Ajmin(p+1)-KLU),LDA,ipiv_aux,rW(klu+1,1,1,p),2*klu,INFO)

      do j=1,klu/2
        call T_DEC(SWAP)(2*KLU,rW(1,j,1,p),1,rW(1,klu+1-j,1,p),-1)
      enddo
      if(mod(klu,2) .eq. 1) then
        call T_DEC(SWAP)(KLU,rW(1,klu/2+1,1,p),1,rW(klu+1,klu/2+1,1,p),-1)
      endif

    else
      call T_DEC(LACPY)( 'U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, rW(1,klu-kl+1,1,p), 2*klu)

      call T_DEC(GBAUL)(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,NZERO,norma,info)
      call T_DEC(TBSM)('U','N','U', klu, kl, ku,A(1,Ajmin(p)),lda,rW(1,klu-kl+1,1,p),2*klu)
!I need the U^-1c for the transpose case. 
      call T_DEC(LACPY)( 'U', klu, klu, rW(1,klu-kl+1,1,p), 2*klu, rW(klu+1,klu-kl+1,1,p), 2*klu)
      call T_DEC(TBSM)('L','N','N', klu, kl, kl,A(ku+1,Ajmin(p)),lda,rW(1,klu-kl+1,1,p),2*klu)
    endif
  endif

! Two thread partitions 
! Here we will call the simplified two-thread spike to work on the partition.
  if(p > 1 .and. p <= 1 + spm(23)) then 
!copy b and c in to vw

!Copying in to the area of vw worked on by the first thread.
! This is submatrix c
    if(keys(2,p-1) == omp_get_thread_num()) then
      vw2(1:(Ajmin(p+1)-Ajmin(p))/2,(p-2)*(kl+ku)+1:(p-1)*(kl+ku))=ZERO_PREC
      if(pivot) then 
        call T_DEC(LACPY)('U', kl, kl, A(klu+kl+ku+1,(Ajmin(p)+Ajmin(p+1))/2-kl), lda-1, vw2(1,(p-2)*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p))
      else
        call T_DEC(LACPY)('U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, vw2(1,(p-2)*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p))
      endif
    endif

!Copying in to the area of vw worked on by the second thread.
! This is submatrix b
    if(keys(3,p-1) == omp_get_thread_num()) then
      vw2((Ajmin(p+1)-Ajmin(p))/2+1:Ajmin(p+1)-Ajmin(p),(p-2)*(kl+ku)+1:(p-1)*(kl+ku))=ZERO_PREC
      if(pivot) then
        do j=1,ku
          call T_DEC(COPY)(j,A(klu+1,Ajmin(p)+ku-j),1,vw2((Ajmin(p+1)-Ajmin(p))/2+1,(p-1)*(kl+ku)-ku+j),-1)
        enddo
      else
        call T_DEC(LACPY)('L', ku, ku, A(1,Ajmin(p+1))       , lda-1, vw2(Ajmin(p+1)-Ajmin(p)-ku+1,(p-2)*(kl+ku)+1), Ajmin(p+1)-Ajmin(p))
      endif
    endif
    call T_DEC(spike_GBTRFk2)(pivot,mysize,kl,ku,A(1,Ajmin(p)),lda,info,ipiv(Ajmin(p)),keys(1,p-1))
    call T_DEC(spike_GBTRSk2)(pivot,'N',Ajmin(p+1)-Ajmin(p),kl+ku,kl,ku,A(1,Ajmin(p)),lda,vw2(1,(p-2)*(kl+ku)+1),mysize,ipiv(Ajmin(p)),keys(1,p-1), spike2_red(1,(p-2)*2*klu+1), spike2_grj(1,(p-2)*(kl+ku)+1),.true.)
!copy vw out to rV
    if(keys(2,p-1) == omp_get_thread_num()) then
      rV(1:klu,:,1,p) = ZERO_PREC
      rW(1:klu,:,1,p) = ZERO_PREC
!Initializing the space into which the v tips are copied to zero. 

      !Copying in to the area of vw worked on by the first thread... 
      call T_DEC(LACPY)('F', klu, ku, vw2(1,(p-2)*(kl+ku)+1)                           , Ajmin(p+1)-Ajmin(p), rV(1,1,1,p)            ,2*klu)
      call T_DEC(LACPY)('F', klu, kl, vw2(1,(p-2)*(kl+ku)+ku+1)                        , Ajmin(p+1)-Ajmin(p), rW(1,klu-kl+1,1,p)     ,2*klu)
    endif

 
    if(keys(3,p-1) == omp_get_thread_num()) then
      rV(klu+1:2*klu,:,1,p) = ZERO_PREC
      rW(klu+1:2*klu,:,1,p) = ZERO_PREC
      !Copying in to the area of vw worked on by the second thread... 
      if(pivot) then
        do j=1,ku
          call T_DEC(COPY)(klu,vw2((Ajmin(p+1)-Ajmin(p))/2+1,(p-1)*(kl+ku)+1-j),   -1, rV(klu+1,j,1,p)       , 1)
        enddo
        do j=1,kl
          call T_DEC(COPY)(klu,vw2((Ajmin(p+1)-Ajmin(p))/2+1,(p-2)*(kl+ku)+kl+1-j),-1, rW(klu+1,klu-kl+j,1,p), 1)
        enddo
      else
        call T_DEC(LACPY)('F', klu, ku, vw2(Ajmin(p+1)-Ajmin(p)-klu+1,(p-2)*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p), rV(klu+1,1,1,p)        ,2*klu)
        call T_DEC(LACPY)('F', klu, kl, vw2(Ajmin(p+1)-Ajmin(p)-klu+1,(p-2)*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p), rW(klu+1,klu-kl+1,1,p) ,2*klu)
      endif

    endif

  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    p_adjust = p - (spm(23))-2

    if(pivot) then
      icstart = klu+1+kl+ku
      call T_DEC(GBTRF)(mysize,mysize,kl,ku,A(iastartl,Ajmin(p)),lda,ipiv(Ajmin(p)),info)
      vw(:, p_adjust*(kl+ku)+1:(p_adjust+1)*(kl+ku)) = ZERO_PREC
      call T_DEC(LACPY)( 'L', ku, ku, A(ibstart,Ajmin(p)  )   , lda-1, vw(mysize-ku+1, p_adjust*(kl+ku)+1)   , mysize)
      call T_DEC(LACPY)( 'U', kl, kl, A(icstart,Ajmin(p+1)-kl), lda-1, vw(1          , p_adjust*(kl+ku)+ku+1), mysize)

      ipiv_aux(1:2*klu) = ipiv(Ajmin(p+1)-2*klu : Ajmin(p+1)-1) - (mysize - 2*klu)
      call T_DEC(GBTRSL)('N',2*KLU ,KL,KU,KU   ,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux      ,vw(mysize-2*klu+1, p_adjust*(kl+ku)+1)   ,mysize,INFO)
      call T_DEC(GBTRSL)('N',mysize,KL,KU,KL   ,A(iastartl,Ajmin(p)        ),LDA,ipiv(Ajmin(p)),vw(1             , p_adjust*(kl+ku)+ku+1),mysize,INFO)
      call T_DEC(GBTRSU)('N',mysize,KL,KU,KL+KU,A(KLU-KL+1,Ajmin(p))        ,LDA,IPIV(Ajmin(p)),vw(1             , p_adjust*(kl+ku)+1)   ,mysize,INFO)
    else
      call T_DEC(GBALU)(Ajmin(p+1)-Ajmin(p),kl,ku,A(1,Ajmin(p)),lda,nzero,norma,info)
      vw(:, p_adjust*(kl+ku)+1:(p_adjust+1)*(kl+ku)) = ZERO_PREC
      call T_DEC(LACPY)( 'L', ku, ku, A(1,Ajmin(p+1)),        lda-1, vw(Ajmin(p+1)-Ajmin(p)-ku+1,p_adjust*(kl+ku)+1), Ajmin(p+1)-Ajmin(p))
      call T_DEC(LACPY)( 'U', kl, kl, A(kl+ku+1,Ajmin(p)-kl), lda-1, vw(1,p_adjust*(kl+ku)+ku+1),                     Ajmin(p+1)-Ajmin(p))
 
      call T_DEC(TBSM)('L','N','U',klu   , ku   , kl, A(ku+1, Ajmin(p+1)-klu), lda, vw(mysize-klu+1, p_adjust*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p))
      call T_DEC(TBSM)('L','N','U',mysize, kl   , kl, A(ku+1, Ajmin(p)      ), lda, vw(1                        , p_adjust*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p))
      call T_DEC(TBSM)('U','N','N',mysize, ku+kl, ku, A(1   , Ajmin(p)      ), lda, vw(1                        , p_adjust*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p))
    endif

    rV(:,:,1,p) = ZERO_PREC
    rW(:,:,1,p) = ZERO_PREC

    call T_DEC(LACPY)('F', klu, ku, vw(1                        ,p_adjust*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p), rV(1,1,1,p),2*klu)
    call T_DEC(LACPY)('F', klu, kl, vw(1                        ,p_adjust*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p), rW(1,klu-kl+1,1,p),2*klu)
    call T_DEC(LACPY)('F', klu, ku, vw(Ajmin(p+1)-Ajmin(p)-klu+1,p_adjust*(kl+ku)+1)   , Ajmin(p+1)-Ajmin(p), rV(klu+1,1,1,p),2*klu)
    call T_DEC(LACPY)('F', klu, kl, vw(Ajmin(p+1)-Ajmin(p)-klu+1,p_adjust*(kl+ku)+ku+1), Ajmin(p+1)-Ajmin(p), rW(klu+1,klu-kl+1,1,p),2*klu)
  endif

! Handle info
! We want errors to take presedence (as they are most critical -- unrecoverable)

  if(info .eq. 0) then
    !DGBALU was happy -- don't have to do anything
  else
    if(info .gt. 0) then
    !DGBALU had to boost 
      info = 1
    else
    !DGBALU found an error/illegal value 
      info = 2
    endif
  endif
  if(pivot) then
    call wdeallocate_1i(ipiv_aux)
  endif
  !$omp atomic
  max_info=max(info,max_info)
  timing_array(i) = OMP_GET_WTIME() - timing_array(i) 
enddo 

parallel_end_time=OMP_GET_WTIME()

info=max_info

if(spm(23) > 0) then
  call T_ALLOC(wdeallocate_2)(vw2) 
endif 
if(spm(20)-spm(23) > 0) then
  call T_ALLOC(wdeallocate_2)(vw)
endif
!<<<<<<<<<<<<<< Reduced systems initialization >>>>>>>>>>>>>>>>>>!

redsys_start_time = OMP_GET_WTIME()

! Reduced systems initialization - first level !
!$omp parallel do default(shared) private(p)
do p=1,nbpart-1,1
  call T_DEC(spike_invred)(klu,rV(1,1,1,p),rW(1,1,1,p+1),red(1,1,p))
enddo

call T_DEC(spike_prep_recn)(kl,ku,rV(1,1,1,1),rW(1,1,1,1),red(1,1,1),nbpart,nblevel)

redsys_end_time = OMP_GET_WTIME()

!call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
call T_ALLOC(wdeallocate_2)(spike2_red)
call T_ALLOC(wdeallocate_2)(spike2_grj)

endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'---------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Factorization Time|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|   LU or SPIKE on blocks    |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread |   Partition   |          Time         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |       ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'       | ')
    call wwrite_d(fout,timing_array(i))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Factorize |       ')
  call wwrite_d(fout,parallel_end_time-parallel_start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System   |       ')
  call wwrite_d(fout,redsys_end_time-redsys_start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Factorize|       ')
  call wwrite_d(fout,end_time-start_time)
  call wwrite_s(fout,'  |')
  call wwrite_n(fout)
  call wwrite_s(fout,'--------------------------------------------------')
  call wwrite_n(fout)
  endif

call wdeallocate_1d(timing_array)

end subroutine T_DEC(spike_gbtrf2)


subroutine T_DEC(spike_gbtrs2) &
(spm,n,kl,ku,lda,nrhs,A,rV,rW,red,ipiv,f,ldf)
! Purpose
! -------
! This is the code that ultimately does the work of the
! SPIKE non-transpose solve. However, there is a more user friendly 
! interface. T_DEC(spike_gbtrs) should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (in) 
! Work space, to contain various parts of the reduced system
! Came from the factorization 
! 
! f (in/out)
! The collection of vectors B.
!
! ldf (in)
! The leading dimension of f.
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,nrhs, nbpart,nblevel,nthread,ldf
integer :: mysize,mystart
integer, dimension(:), pointer :: Ajmin
integer, dimension(*) :: ipiv
MAIN_TYPE, dimension(lda,*):: A
MAIN_TYPE, dimension(ldf,*):: f
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku),*):: red
MAIN_TYPE, dimension(2*max(kl,ku),max(kl,ku),spm(21),*):: rV,rW
integer, dimension(:,:), pointer :: keys
integer ::i,j,klu,p,info
!integer,external :: omp_get_thread_num
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*nrhs) :: spike2_grj
MAIN_TYPE, pointer, dimension(:,:) :: spike2_red
MAIN_TYPE, pointer, dimension(:,:) :: spike2_grj
MAIN_TYPE, dimension(:,:), pointer:: g
MAIN_TYPE, dimension(:,:), pointer :: g2
MAIN_TYPE, dimension(:,:,:), pointer :: grj
MAIN_TYPE, dimension(:,:), pointer :: aux
logical :: invred
double precision, dimension(:,:),pointer :: timing_array
!character(len=100) :: timing_print_format
integer :: infoloc
integer(8),parameter :: fout=6
!double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: t1,t2,t3,t4
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time_1, parallel_end_time_1
double precision :: parallel_start_time_2, parallel_end_time_2
logical :: pivot
integer :: kd,iastartl,iastartu,ibstart,jbstart,icstart,jcstart,inc,istartf,istartgrj
integer, dimension(:), pointer :: ipiv_aux

pivot = spm(3) .eq. 1

!COPY OPTION
!MAIN_TYPE, pointer, dimension(:,:) :: f_aux1
!MAIN_TYPE, pointer, dimension(:,:) :: f_aux2
call T_ALLOC(wallocate_2)(spike2_red,2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call T_ALLOC(wallocate_2)(spike2_grj,2*max(kl,ku),spm(23)*nrhs,infoloc)

start_time = OMP_GET_WTIME()
call wallocate_2d(timing_array,spm(22),2,infoloc)

if(spm(20) == 1) then
  parallel_start_time_1 = OMP_GET_WTIME()
  timing_array(1,1) = OMP_GET_WTIME()
  if(pivot) then
    klu=max(kl,ku)
    iastartu = klu-kl+1
    call T_DEC(GBTRS)('N',n,kl,ku,nrhs,A(iastartu,1),lda,ipiv,f,ldf,INFO)
  else
    call T_DEC(TBSM)('L','N','U', n, nrhs, kl, A(ku+1,1),lda,f,ldf)
  endif
  timing_array(1,1) = OMP_GET_WTIME() - timing_array(1,1)
  parallel_end_time_1 = OMP_GET_WTIME()

  parallel_start_time_2 = OMP_GET_WTIME()
  timing_array(1,2) = OMP_GET_WTIME()
  if(pivot) then
  else
    call T_DEC(TBSM)('U','N','N', n, nrhs, ku, A        ,lda,f,ldf)
  endif
  timing_array(1,2) = OMP_GET_WTIME() - timing_array(1,2)
  parallel_end_time_2 = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time = 0.0d0
else

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call spikerl_calc_size_partitions(spm,Ajmin,n)
nthread = spm(22)
nbpart = spm(20) 
nblevel = spm(21)
invred=.false.

klu=max(kl,ku)
call T_ALLOC(wallocate_3)(grj,2*klu,nrhs,nbpart,infoloc) !! modified rhs
grj=ZERO_PREC

call wallocate_2i(keys,4,nbpart,infoloc)
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

if(nbpart-(spm(23)+2) > 0) then
  call T_ALLOC(wallocate_2)(g,Ajmin(nbpart)-Ajmin(nbpart-1),nrhs*(nbpart-(spm(23)+2)),infoloc)
endif

if(spm(23) > 0) then
  call T_ALLOC(wallocate_2)(g2,Ajmin(3)-Ajmin(2),nrhs*spm(23),infoloc)
endif

!COPY OPTION
!if(pivot) then
!  call T_ALLOC(wallocate_2)(f_aux1,Ajmin(2)-Ajmin(1),nrhs,infoloc)
!  call T_ALLOC(wallocate_2)(f_aux2,Ajmin(nbpart+1)-Ajmin(nbpart),nrhs,infoloc)
!endif

parallel_start_time_1 = OMP_GET_WTIME()

!COPY OPTION
!NOT $omp parallel do private(t1,t2,t3,t4,kd,mysize,iastartl,iastartu,info,j) firstprivate(p,mystart,nrhs,kl,ku,lda,ldf,klu,nbpart,nthread), default(none), shared(f_aux1,f_aux2,pivot,timing_array,spm,Ajmin,A,f,g2,g,spike2_red,spike2_grj,grj,keys)

!$omp parallel do private(t1,t2,t3,t4,kd,mysize,iastartl,iastartu,info,j) firstprivate(p,mystart,nrhs,kl,ku,lda,ldf,klu,nbpart,nthread), default(none), shared(pivot,timing_array,spm,Ajmin,A,f,g2,g,spike2_red,spike2_grj,grj,keys)
do i=1,nthread,1

  if(pivot) then
    iastartu = klu-kl+1
    kd = kl+ku
  else
    iastartl = ku+1
    iastartu = 1
    kd = ku
  endif

  timing_array(i,1) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
  mysize = Ajmin(p+1)-Ajmin(p)

!! Solving Ajgj=fj
  if (p<=nbpart-1) then
    if (p==1) then
      if(pivot) then
        iastartl = klu-kl+1
!COPY OPTION
!        call T_DEC(LACPY)('F',mysize,nrhs,f(Ajmin(p),1),ldf,f_aux1(1,1),mysize)
!        call T_DEC(GBTRSL)('N',mysize,KL,KU,nrhs,A(iastartl,Ajmin(p)),LDA,ipiv(Ajmin(p)),f_aux1(1,1),mysize,INFO)
!        call T_DEC(LACPY)('F', klu, nrhs, f_aux1(mysize+1-klu,1),mysize, grj(klu+1,1,p), 2*klu )
        call T_DEC(GBTRSL)('N',Ajmin(p+1)-Ajmin(p),KL,KU,nrhs,A(iastartl,Ajmin(p)),LDA,ipiv(Ajmin(p)),f(Ajmin(p),1),ldf,INFO)
        call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1),ldf, grj(klu+1,1,p), 2*klu )
      else
        call T_DEC(TBSM)('L','N','U', Ajmin(p+1)-Ajmin(p), nrhs,kl,A(iastartl,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
        call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
      endif
      call T_DEC(TBSM)('U','N','N', klu, nrhs, kd, A(iastartu,Ajmin(p+1)-klu), lda, grj(klu+1,1,p),2*klu)
      grj(1:klu,1:nrhs,p)=ZERO_PREC

    else
      if(p <= spm(23)+1) then
    ! Spike version case
        mysize = Ajmin(3) - Ajmin(2)
        mystart = nrhs*(p-2)+1

        if(keys(2,p-1) == omp_get_thread_num()) then
          call T_DEC(LACPY)('F',mysize/2,nrhs,f(Ajmin(p),1),ldf,g2(1,mystart),mysize)
        endif

        if(keys(3,p-1) == omp_get_thread_num()) then
          if(pivot) then
!            do j=1,nrhs
!              call T_DEC(COPY)(mysize/2,f(Ajmin(p)+mysize/2,j),-1,g2(1+mysize/2,mystart+j-1),1)
!            enddo

            do j=1,nrhs/2
              call T_DEC(SWAP)(mysize/2,f(ajmin(p)+mysize/2,j),1,f(ajmin(p)+mysize/2,nrhs+1-j),-1)
            enddo
            if(mod(nrhs,2) .eq. 1) then
              call T_DEC(SWAP)(mysize/4,f(ajmin(p)+mysize/2,nrhs/2+1),1,f(ajmin(p)+mysize/2+(mysize/2+1)/2,nrhs/2+1),-1)
            endif
          endif
!          else
            call T_DEC(LACPY)('F',mysize/2,nrhs,f(Ajmin(p)+mysize/2,1),ldf,g2(1+mysize/2,mystart),mysize)
!          endif
        endif

        call T_DEC(spike_GBTRSk2)(pivot,'N',mysize,nrhs,kl,ku,A(1,Ajmin(p)),lda,g2(1,mystart),mysize,ipiv(Ajmin(p)),keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*(nrhs)+1),.false.)

        if(keys(2,p-1) == omp_get_thread_num()) then
          call T_DEC(LACPY)('F',klu,nrhs,g2(1,mystart),mysize,grj(1,1,p),2*klu)
        endif
        if(keys(3,p-1) == omp_get_thread_num()) then
          if(pivot) then
            do j=1,nrhs
              call T_DEC(COPY)(klu,g2(1+mysize/2,mystart+j-1),-1,grj(klu+1,nrhs-j+1,p),1)
            enddo
          else
            call T_DEC(LACPY)('F',klu,nrhs,g2(mysize-klu+1,mystart),mysize,grj(klu+1,1,p),2*klu)
          endif

        endif
      else 
! Single threaded partition 
        mysize = Ajmin(nbpart)-Ajmin(nbpart-1)

        mystart = nrhs*(p-(2+spm(23)))+1
        t1 = OMP_GET_WTIME()
        call T_DEC(LACPY)('F',mysize,nrhs,f(Ajmin(p),1),ldf,g(1, mystart),mysize)
        t2 = OMP_GET_WTIME()
        if(pivot) then
          iastartl = klu-kl+1
!          call T_DEC(GBTRSL)('N',mysize,KL,KU,nrhs, A(iastartl,Ajmin(p)), LDA,ipiv(Ajmin(p)),g(1,mystart),mysize,INFO)
          t3 = OMP_GET_WTIME()
!           call T_DEC(TBSM)('U','N','N', mysize, nrhs, kl+ku, A(klu-kl+1, Ajmin(p)), lda, g(1,mystart), mysize)
!          call T_DEC(GBTRSU)('N',mysize,KL,KU,nrhs, A(iastartl,Ajmin(p)), LDA,IPIV(Ajmin(p)),g(1,mystart),mysize,INFO)
          call T_DEC(GBTRS)('N',mysize,KL,KU,nrhs, A(iastartl,Ajmin(p)), LDA,ipiv(Ajmin(p)),g(1,mystart),mysize,INFO)
          t4 = OMP_GET_WTIME()
!print *, t2-t1, t3-t2, t4-t3
        else
          call T_DEC(TBSM)('L','N','U', mysize, nrhs, kl,A(ku+1,Ajmin(p)),lda,g(1, mystart),mysize)
          call T_DEC(TBSM)('U','N','N', mysize, nrhs, ku,A(1,Ajmin(p)),lda,   g(1, mystart),mysize)
        endif
        call T_DEC(LACPY)('F',klu,nrhs,g(1, mystart), mysize,grj(1,1,p),2*klu)
        call T_DEC(LACPY)('F',klu,nrhs,g(mysize-klu+1, mystart),mysize,grj(klu+1,1,p),2*klu)

      endif
    endif
  endif
!! solving Anbpart gnbpart =f nbpart
  if (p==nbpart) then
    if(pivot) then
      iastartu = klu+1
      iastartl = klu-ku+1

!QfQ
!COPY OPTION
!      do j=1,nrhs,1
!        call T_DEC(COPY)(mysize,f(Ajmin(p),j),1,f_aux2(1,nrhs+1-j),-1)
!      enddo

      do j=1,nrhs/2
        call T_DEC(SWAP)(mysize,f(ajmin(p),j),1,f(ajmin(p),nrhs+1-j),-1)
      enddo
      if(mod(nrhs,2) .eq. 1) then
        call T_DEC(SWAP)(mysize/2,f(ajmin(p),nrhs/2+1),1,f(ajmin(p)+(mysize+1)/2,nrhs/2+1),-1)
      endif
!      do j=1,nrhs,1
!        call T_DEC(COPY)(mysize,f(Ajmin(p),j),1,f_aux2(1,nrhs+1-j),-1)
!      enddo


      call T_DEC(GBTRSL)('N',mysize,KU,KL,NRHS,A(iastartl,ajmin(p)),LDA,IPIV(ajmin(p)),f(ajmin(p),1),LDf,INFO)
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1),ldf, grj(1,1,p), 2*klu )
!      call T_DEC(GBTRSL)('N',mysize,KU,KL,NRHS,A(iastartl,ajmin(p)),LDA,IPIV(ajmin(p)),f_aux2(1,1),mysize,INFO)
!      call T_DEC(LACPY)('F', klu, nrhs, f_aux2(mysize+1-klu,1),mysize, grj(1,1,p), 2*klu )

! IPIV doesn't really matter in this case because the U sweep doesn't use pivoting -- that is just in L. 
! So, that is why IPIV is set to 1.
      call T_DEC(GBTRSU)('N',klu,KU,KL,NRHS,A(iastartl,ajmin(p+1)-klu),LDA,IPIV(1),grj(1,1,p),2*klu,INFO)
      do j=1,nrhs/2
        call T_DEC(SWAP)(klu,grj(1,j,p),1,grj(1,nrhs+1-j,p),-1)
      enddo
      if(mod(nrhs,2) .eq. 1) then
        call T_DEC(SWAP)(klu/2,grj(1,nrhs/2+1,p),1,grj((klu+1)/2+1,nrhs/2+1,p),-1)
      endif
    else
      call T_DEC(TBSM)('U','N','U', mysize, nrhs, ku,A(iastartu,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1),ldf, grj(1,1,p), 2*klu )
      call T_DEC(TBSM)('L','N','N', klu, nrhs, kl,A(ku+1,Ajmin(p)),lda,grj(1,1,p),2*klu)
    endif 

    grj(klu+1:,1:nrhs,p)=ZERO_PREC
  endif
  timing_array(i,1) = OMP_GET_WTIME() - timing_array(i,1)
enddo
 
parallel_end_time_1 = OMP_GET_WTIME()

if(spm(23) > 0) then
  call T_ALLOC(wdeallocate_2)(g2)
endif

redsys_start_time = OMP_GET_WTIME()

call T_DEC(spike_solve_recn)(kl,ku,nrhs,rV(1,1,1,1),rW(1,1,1,1),red(1,1,1),nbpart,nblevel,grj(1,1,1))

!<<<<<<<<<<<<<<<<<<<<<<<< Retrieval  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>!
!COPY OPTION
!NOT $omp parallel do default(none) shared(f_aux1,f_aux2,A,grj,lda,ldf,nbpart,nrhs,kl,ku,klu,ajmin,pivot,spm), private(istartgrj,istartf,info,inc,mysize,aux,ipiv_aux,ibstart,jbstart,iastartl,icstart,jcstart,infoloc)

!$omp parallel do default(none) shared(A,grj,lda,ldf,nbpart,nrhs,kl,ku,klu,ajmin,pivot,spm), private(istartgrj,istartf,info,inc,mysize,aux,ipiv_aux,ibstart,jbstart,iastartl,icstart,jcstart,infoloc)
do p=1,nbpart,1
!!! Aixi=fi-Bi*x(i+1)b-Ci*x(i-1)t
  mysize = Ajmin(p+1)-Ajmin(p)
  if(p==1) then
    if(pivot) then
      ibstart = klu+1
      jbstart = Ajmin(p)
      iastartl = klu-kl+1
    else
      ibstart = 1
      jbstart = Ajmin(p+1)
    endif

    call T_DEC(TRMM)('L','L','N','N',ku,nrhs,-ONE_PREC,A(ibstart,jbstart),lda-1,grj(1,1,p+1),2*klu)

    if(pivot) then

      call wallocate_1i(ipiv_aux,2*klu,infoloc)
      ipiv_aux(1:2*klu) = ipiv(Ajmin(p+1)-2*klu : Ajmin(p+1)-1) - (mysize - 2*klu)

      call T_ALLOC(wallocate_2)(aux,2*klu,nrhs,infoloc)
      aux(1:2*klu,1:nrhs) = ZERO_PREC
      
      call T_DEC(LACPY)('A',ku,nrhs,grj(1,1,p+1),2*klu,aux(2*klu-ku+1,1),2*klu)
!      call T_DEC(GBTRSL)('N',2*KLU,KL,KU,NRHS,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux,aux,2*klu,INFO)
      call T_DEC(GBTRSL)('N',2*KLU,KL,KU,NRHS,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux,aux,2*klu,INFO)
! COPY OPTION
      do i=1,nrhs
        call T_DEC(AXPY)(2*klu,ONE_PREC,aux(1,i),1,f(Ajmin(p+1)-2*klu,i),1) 
!        call T_DEC(AXPY)(2*klu,ONE_PREC,aux(1,i),1,f_aux1(mysize+1-2*klu,i),1) 
      enddo

      call wdeallocate_1i(ipiv_aux)      
      call T_ALLOC(wdeallocate_2)(aux)

    else

      call T_DEC(TBSM)('L','N','U',ku, nrhs,kl,A(ku+1,Ajmin(p+1)-ku),lda,grj(1,1,p+1),2*klu) ! finish U

      do i=1,nrhs
        call T_DEC(AXPY)(ku,ONE_PREC,grj(1,i,p+1),1,f(Ajmin(p+1)-ku,i),1) 
      enddo
    endif

  endif
  if(p==nbpart) then
    if(pivot) then
      !icstart = 2*klu+1-ku
      icstart = klu+1
      jcstart = Ajmin(p)
      iastartl = klu-ku+1

      call wallocate_1i(ipiv_aux,2*klu,infoloc)
      ipiv_aux(1:2*klu) = ipiv(Ajmin(p+1)-2*klu : Ajmin(p+1)-1) - (mysize - 2*klu)

      call T_ALLOC(wallocate_2)(aux,2*klu,nrhs,infoloc)
      aux(1:2*klu,1:nrhs) = ZERO_PREC

      do j=1,nrhs/2
        call T_DEC(SWAP)(kl,grj(2*klu-kl+1,j,p-1),1,grj(2*klu-kl+1,nrhs+1-j,p-1),-1)
      enddo
      if(mod(nrhs,2) .eq. 1) then
        call T_DEC(SWAP)(kl/2,grj((2*klu-kl+1),nrhs/2+1,p-1),1,grj((2*klu-kl+1)+(kl+1)/2,nrhs/2+1,p-1),-1)
      endif

      call T_DEC(TRMM)('L','L','N','N',kl,nrhs,-ONE_PREC,A(icstart,jcstart),lda-1,grj(2*klu-kl+1,1,p-1),2*klu)
      call T_DEC(LACPY)('A',kl,nrhs,grj(2*klu-kl+1,1,p-1),2*klu,aux(2*klu-kl+1,1),2*klu)
      call T_DEC(GBTRSL)('N',2*KLU,KU,KL,NRHS,A(iastartl,Ajmin(p+1)-2*KLU),LDA,ipiv_aux,aux,2*klu,INFO)

      do i=1,nrhs
        call T_DEC(AXPY)(2*klu,ONE_PREC,aux(1,i),1,f(Ajmin(p+1)-2*klu,i),1) 
      enddo

      call wdeallocate_1i(ipiv_aux)      
      call T_ALLOC(wdeallocate_2)(aux)

    else
      call T_DEC(TRMM)('L','U','N','N',kl,nrhs,-ONE_PREC,A(kl+ku+1,Ajmin(p)-kl),lda-1,grj(2*klu-kl+1,1,p-1),2*klu)
      call T_DEC(TBSM)('U','N','U', kl, nrhs, ku,A(1,Ajmin(p)),lda,grj(2*klu-kl+1,1,p-1),2*klu) ! finish L
      do i=1,nrhs
        call T_DEC(AXPY)(kl,ONE_PREC,grj(2*klu-kl+1,i,p-1),1,f(Ajmin(p),i),1) 
      enddo


    endif
  endif

  if (p<=nbpart-1 .and. p>=2) then
    if(pivot) then
      ibstart = klu+1
      icstart = klu+ku+1+kl
      jbstart = Ajmin(p)

      if((p > 1 + spm(23) .and. p < nbpart) .or. p==nbpart) then 
        !normal parition
        jcstart = Ajmin(p+1) - kl
        inc = 1
        istartf = Ajmin(p+1)-ku
        istartgrj=0;
      else
        !2x2 spike parititon
        jcstart = (Ajmin(p)+Ajmin(p+1))/2-kl
        istartf = (Ajmin(p)+Ajmin(p+1))/2
        inc = -1
        istartgrj=nrhs+1;
      endif
    else
      inc = 1
      ibstart = 1
      icstart = kl+ku+1
      jbstart = Ajmin(p+1)
      jcstart = Ajmin(p) - kl
      istartf = Ajmin(p+1)-ku
      istartgrj=0;
    endif 

    
    call T_DEC(TRMM)('L','L','N','N',ku,nrhs,-ONE_PREC,A(ibstart, jbstart),lda-1,grj(1         ,1,p+1),2*klu)
    call T_DEC(TRMM)('L','U','N','N',kl,nrhs,-ONE_PREC,A(icstart, jcstart),lda-1,grj(2*klu-kl+1,1,p-1),2*klu)

    do i=1,nrhs
      call T_DEC(AXPY)(kl,ONE_PREC,grj(2*klu-kl+1,i,p-1),1,f(Ajmin(p),i),1) 
!      call T_DEC(AXPY)(ku,ONE_PREC,grj(1,i,p+1),1,f(Ajmin(p+1)-ku,i),1) 
      call T_DEC(AXPY)(ku,ONE_PREC,grj(1,istartgrj + i*inc,p+1),inc,f(istartf,i),1) 
    enddo

  end if

enddo


redsys_end_time = OMP_GET_WTIME() 

!allocate(keys(4,nbpart))
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(3,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

!! Solving Ajgj=fj
parallel_start_time_2 = OMP_GET_WTIME()
!$omp parallel do  private(t1,t2,p,mysize,j,mystart), shared(g,spike2_red,spike2_grj,keys)
do i=1,nthread,1
  call threadmap(i,p,nbpart,spm(23))
  mysize = Ajmin(p+1)-Ajmin(p)
  timing_array(i,2) = OMP_GET_WTIME()
  if(p==1) then
    if(pivot) then
!COPY OPTION
!      call T_DEC(GBTRSU)('N',mysize,KL,KU,NRHS,A(KLU-KL+1,Ajmin(p)),LDA,IPIV(Ajmin(p)),f_aux1(1,1),mysize,INFO)
!      call T_DEC(LACPY)('F',mysize,nrhs,f_aux1(1,1),mysize,f(Ajmin(p),1),ldf)
      call T_DEC(GBTRSU)('N',mysize,KL,KU,NRHS,A(KLU-KL+1,Ajmin(p)),LDA,IPIV(Ajmin(p)),f(Ajmin(p),1),ldf,INFO)
    else
      call T_DEC(TBSM)('U','N','N', mysize, nrhs, ku,A(1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
    endif
  endif
  if(p==nbpart) then
    if(pivot) then
!COPY OPTION
      call T_DEC(GBTRSU)('N',mysize,KU,KL,NRHS,A(KLU-KU+1,Ajmin(p)),LDA,IPIV(Ajmin(p)),f(Ajmin(p),1),ldf,INFO)
      do j=1,nrhs/2
        call T_DEC(SWAP)(mysize,f(ajmin(p),j),1,f(ajmin(p),nrhs+1-j),-1)
      enddo
      if(mod(nrhs,2) .eq. 1) then
        call T_DEC(SWAP)(mysize/2,f(ajmin(p),nrhs/2+1),1,f(ajmin(p)+(mysize+1)/2,nrhs/2+1),-1)
      endif

!      call T_DEC(GBTRSU)('N',mysize,KU,KL,NRHS,A(KLU-KU+1,Ajmin(p)),LDA,IPIV(Ajmin(p)),f_aux2(1,1),mysize,INFO)
!      do j=1,nrhs,1
!        call T_DEC(COPY)(mysize,f_aux2(1,nrhs+1-j),1,f(Ajmin(p),j),-1)
!      enddo

    else 
      call T_DEC(TBSM)('L','N','N', mysize, nrhs, kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)
    endif
  endif

  if(p > 1 .and. p <= 1 + spm(23)) then 

    call T_DEC(spike_GBTRSk2)(pivot,'N',mysize, nrhs, kl, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf, ipiv(ajmin(p)), keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1),.false.)

    if(pivot) then
      if(keys(3,p-1) == omp_get_thread_num()) then
        do j=1,nrhs/2
          call T_DEC(SWAP)(mysize/2,f(ajmin(p)+mysize/2,j),1,f(ajmin(p)+mysize/2,nrhs+1-j),-1)
        enddo
        if(mod(nrhs,2) .eq. 1) then
          call T_DEC(SWAP)(mysize/4,f(ajmin(p)+mysize/2,nrhs/2+1),1,f(ajmin(p)+mysize/2+(mysize/2+1)/2,nrhs/2+1),-1)
        endif
      endif
    endif
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    if(pivot) then
t1 = OMP_GET_WTIME()
      mystart = nrhs*(p-2)+1
!COPY OPTION
!      call T_DEC(LACPY)('F',mysize,nrhs,f(Ajmin(p),1),ldf,g(1, mystart),mysize)
!      call T_DEC(GBTRS)('N',mysize,KL,KU,nrhs, A(klu-kl+1,Ajmin(p)), LDA,ipiv(Ajmin(p)),g(1,mystart),mysize,INFO)
!      call T_DEC(LACPY)('F',mysize,nrhs,g(1, mystart),mysize,f(Ajmin(p),1),ldf)

      call T_DEC(GBTRS)('N',mysize,KL,KU,nrhs, A(klu-kl+1,Ajmin(p)), LDA,ipiv(Ajmin(p)),f(Ajmin(p),1),ldf,INFO)

!print *, t2-t1, t3-t2
    else
      call T_DEC(TBSM)('L','N','U', Ajmin(p+1)-Ajmin(p), nrhs, kl, A(ku+1,Ajmin(p)), lda, f(Ajmin(p),1), ldf)
      call T_DEC(TBSM)('U','N','N', Ajmin(p+1)-Ajmin(p) ,nrhs, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf)
    endif
  endif

!! solving Anbpart gnbpart =f nbpart
  timing_array(i,2) = OMP_GET_WTIME() - timing_array(i,2)
enddo

parallel_end_time_2 = OMP_GET_WTIME()
call T_ALLOC(wdeallocate_3)(grj)
call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
call T_ALLOC(wdeallocate_2)(spike2_red)
call T_ALLOC(wdeallocate_2)(spike2_grj)
if(nbpart-(2+spm(23)) > 0) then
  call T_ALLOC(wdeallocate_2)(g)
endif


endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Solve Time|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|  Solve Sweeps on blocks    |----------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread | Partition |         Time 1        |          Time 2       |')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |     ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'     | ')
    call wwrite_d(fout,timing_array(i,1))
    call wwrite_s(fout,' | ')
    call wwrite_d(fout,timing_array(i,2))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Solve        | ')
  call wwrite_d(fout,(parallel_end_time_1-parallel_start_time_1))
  call wwrite_s(fout,' | ')
  call wwrite_d(fout,(parallel_end_time_2-parallel_start_time_2))
  call wwrite_s(fout,' |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System Solve| ')
  call wwrite_d(fout,(redsys_end_time-redsys_start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Solve       | ')
  call wwrite_d(fout,(end_time-start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
endif

call wdeallocate_2d(timing_array)


end subroutine T_DEC(spike_gbtrs2)





subroutine T_DEC(spike_gbtrst2) &
(spm,trans,n,kl,ku,lda,nrhs,A,rV,rW,red,f,ldf)
! Purpose
! -------
! This is the code that ultimately does the work of the
! SPIKE transpose solve. However, there is a more user friendly 
! interface. T_DEC(spike_gbtrs) should be used instead. 
! 
! Arguments 
! ---------
! spm (in) 
! The spike parameter array, as described in spikeinit
! 
! n (in) 
! The nuber of rows and columns of the matrix A
!
! kl, ku (in)
! The lower and upper bandwidths of A
! 
! A (in/out) 
! The banded matrix on which the SPIKE factorization will be performed
! For best results, A should be diagonally dominant.
! 
! lda (in) 
! The leading dimensions of A. lda=kl+ku+1
! 
! rV, rW, red, (out) 
! Work space, to contain various parts of the reduced system
! Came from the factorization 
! 
! f (in/out)
! The collection of vectors B.
!
! ldf (in)
! The leading dimension of f.
!=====================================================================
!  Braegan Spring - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: n,kl,ku,lda,nrhs, nbpart,nblevel,nthread
integer :: mysize,mystart
integer, dimension(:), pointer :: Ajmin
MAIN_TYPE, dimension(lda,*):: A
MAIN_TYPE, dimension(ldf,*):: f
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku),*):: red
MAIN_TYPE, dimension(2*max(kl,ku),max(kl,ku),spm(21),*):: rV,rW
integer, dimension(:,:), pointer :: keys
integer ::i,j,klu,p,ldf
!integer,external :: omp_get_thread_num
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*2*max(kl,ku)) :: spike2_red
!MAIN_TYPE, dimension(2*max(kl,ku),spm(23)*nrhs) :: spike2_grj
MAIN_TYPE, pointer, dimension(:,:) :: spike2_red
MAIN_TYPE, pointer, dimension(:,:) :: spike2_grj
MAIN_TYPE, dimension(:,:), pointer :: g
MAIN_TYPE, dimension(:,:), pointer:: g2
MAIN_TYPE, dimension(:,:,:), pointer :: grj
MAIN_TYPE, dimension(:,:), pointer :: gaux
integer, dimension(1) :: ipiv_dummy 
logical :: invred
integer :: t1,t2,tim
double precision, dimension(:,:),pointer :: timing_array
integer :: infoloc
integer(8),parameter :: fout=6
!double precision, external :: OMP_GET_WTIME
double precision :: start_time,end_time
double precision :: redsys_start_time, redsys_end_time
double precision :: parallel_start_time_1, parallel_end_time_1
double precision :: parallel_start_time_2, parallel_end_time_2
character :: trans

start_time = OMP_GET_WTIME()

call T_ALLOC(wallocate_2)(spike2_red,2*max(kl,ku),spm(23)*2*max(kl,ku),infoloc)
call T_ALLOC(wallocate_2)(spike2_grj,2*max(kl,ku),spm(23)*nrhs,infoloc)
call wallocate_2d(timing_array,spm(22),2,infoloc)

if(spm(20) == 1) then
  parallel_start_time_1 = OMP_GET_WTIME()
  timing_array(1,1) = OMP_GET_WTIME()
  call T_DEC(TBSM)('U',trans,'N', n, nrhs, ku, A,lda,f,ldf)
  timing_array(1,1) = OMP_GET_WTIME() - timing_array(1,1)
  parallel_end_time_1 = OMP_GET_WTIME()

  parallel_start_time_2 = OMP_GET_WTIME()
  timing_array(1,2) = OMP_GET_WTIME()
  call T_DEC(TBSM)('L',trans,'U', n, nrhs, kl, A(ku+1,1),lda,f,ldf)
  timing_array(1,2) = OMP_GET_WTIME() - timing_array(1,2)
  parallel_end_time_2 = OMP_GET_WTIME()
  redsys_start_time = 0.0d0
  redsys_end_time = 0.0d0
else

klu = max(kl,ku)

nbpart = spm(20)
nthread = spm(22)
nblevel = spm(21)
invred=.false.

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call T_ALLOC(wallocate_3)(grj,2*klu,nrhs,nbpart,infoloc) !! modified rhs
call T_ALLOC(wallocate_2)(gaux,2*klu,nrhs*nbpart,infoloc) !! modified rhs

call spikerl_calc_size_partitions(spm,Ajmin,n)
grj=ZERO_PREC
gaux=ZERO_PREC

!manually allocate 'keys' for the partitions on which we will use spike
call wallocate_2i(keys,4,nbpart,infoloc)
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

parallel_start_time_1 = OMP_GET_WTIME()
!Set up the zero augmented right hand sides for the reduced system
!$omp parallel do default(shared) private(i,p,mysize,tim,t2,t1)
do i=1,nthread,1
  call threadmap(i,p,nbpart,spm(23))
  timing_array(i,1) = OMP_GET_WTIME()
  mysize = Ajmin(p+1)-Ajmin(p)
  if(p==1) then
    ! A1 is LU factorized
    !Augment the bottom of f1 with zero

    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
    f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = ZERO_PREC

    ! Sweeps over the right hand side
    call T_DEC(TBSM)('U',trans,'N',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)

    ! Perform multiplication 
    ! We will need to unmodified values for f at a later point, so we pull out the ones that will need to work with
    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
    call T_DEC(TRMM)('L', 'L', trans,'N',ku,nrhs,ONE_PREC,rV(1,1,1,1),2*klu,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)

  endif

  if(p==nbpart) then
   !Augment the top of the final part of f with zeroes

    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
    f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = ZERO_PREC

    call T_DEC(TBSM)('L',trans,'N',mysize,nrhs,kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)

    !w's
    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
    call T_DEC(TRMM)('L','U',trans,'N',kl,nrhs,ONE_PREC, rW(klu+1,klu-kl+1,1,p),2*klu,gaux(1,nrhs*(p-1)+1),2*klu)
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 

    !augment with 0
    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
    f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = ZERO_PREC


    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
    f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = ZERO_PREC

    ! Sweeps over the right hand side
    call T_DEC(TBSM)('U',trans,'N',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)
    call T_DEC(TBSM)('L',trans,'U',mysize,nrhs,kl,A(ku+1,Ajmin(p)),lda,f(Ajmin(p),1),ldf)

    !v's
    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
    call T_DEC(TRMM)('L', 'L', trans,'N',ku,nrhs,ONE_PREC,A(1,Ajmin(p+1)),lda-1,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)

    !w's    
    call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
    call T_DEC(TRMM)('L','U',trans,'N',kl,nrhs,ONE_PREC,A(lda,Ajmin(p)-kl),lda-1,gaux(1,nrhs*(p-1)+1),2*klu)

  endif

   
!Two thread partitions 
  if(p > 1 .and. p <= 1 + spm(23)) then 
    !augment with 0
    if(keys(3,p-1) == omp_get_thread_num()) then
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, grj(klu+1,1,p), 2*klu)
      f(Ajmin(p+1)-klu:Ajmin(p+1)-1,1:nrhs) = ZERO_PREC
    endif
    if(keys(2,p-1) == omp_get_thread_num()) then
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, grj(1,1,p), 2*klu)
      f(Ajmin(p):Ajmin(p)+klu-1,1:nrhs) = ZERO_PREC
    endif

    call T_DEC(spike_GBTRSk2)(.false.,trans,Ajmin(p+1)-Ajmin(p), nrhs, kl, ku, A(1,Ajmin(p)), lda, f(Ajmin(p),1), ldf, ipiv_dummy, keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1),.false.)

    !v's
    if(keys(3,p-1) == omp_get_thread_num()) then
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p+1)-klu,1), ldf, gaux(klu+1,nrhs*(p-1)+1), 2*klu)
      call T_DEC(TRMM)('L', 'L', trans,'N',ku,nrhs,ONE_PREC,A(1,Ajmin(p+1)),lda-1,gaux(klu+1+klu-ku,nrhs*(p-1)+1),2*klu)
    endif

    !w's    
    if(keys(2,p-1) == omp_get_thread_num()) then
      call T_DEC(LACPY)('F', klu, nrhs, f(Ajmin(p),1), ldf, gaux(1,nrhs*(p-1)+1), 2*klu)
      call T_DEC(TRMM)('L','U',trans,'N',kl,nrhs,ONE_PREC,A(lda,Ajmin(p)-kl),lda-1,gaux(1,nrhs*(p-1)+1),2*klu)
    endif

  endif

  timing_array(i,1) = OMP_GET_WTIME() - timing_array(i,1)
enddo
parallel_end_time_1 = OMP_GET_WTIME()

redsys_start_time = OMP_GET_WTIME()
!$omp parallel do default(shared) private(p,mysize,t1,t2,tim)
do p=1,nbpart,1
!  call threadmap(i,p,nbpart,spm(23))
! I know p-1+1 is silly, but I want to make it clear that this is inter-thread communication -- 
! the p-1 follows with what was used before, then the +1 stands for the next partition over
! So this works with the w spike from the next partition over
! so for grj, on the first index, 
! (p-1)*2*klu represents that we're working in partition p
! + klu indicates that we're writing to the bottom half
! + klu - ku indicates that we have implied 0's in gaux, because the dgemm that filled gaux
! began filling at index 1, rather than index ku
  if(p == 1) then
    do j=1,nrhs,1
      call T_DEC(AXPY)(kl,-ONE_PREC,gaux(1,nrhs*(p-1+1)+j),1,grj(2*klu-kl+1,j,p),1)
    enddo 
  endif
 
  if(p == nbpart) then
! This works with the v spike from the previous partition
    do j=1,nrhs,1
      call T_DEC(AXPY)(ku,-ONE_PREC,gaux(2*klu-ku+1,nrhs*(p-1-1)+j),1,grj(1,j,p),1)
    enddo
  endif 

  if(p > 1 .and. p < nbpart) then
    do j=1,nrhs,1
      call T_DEC(AXPY)(kl,-ONE_PREC,gaux(1,nrhs*(p-1+1)+j),1,grj(2*klu-kl+1,j,p),1)
    enddo 
    do j=1,nrhs,1
      call T_DEC(AXPY)(ku,-ONE_PREC,gaux(2*klu-ku+1,nrhs*(p-1-1)+j),1,grj(1,j,p),1)
    enddo
  endif
enddo

call T_ALLOC(wdeallocate_2)(gaux) !! modified rhs

!Reduced system
call T_DEC(spike_solve_recn_transpose)(trans,kl,ku,nrhs,rV,rW,red,nbpart,nblevel,grj(1,1,1)) 


redsys_end_time = OMP_GET_WTIME()

 
keys(1,1:nbpart) = 0
do i=1,spm(23)
  keys(2,i) = 2*(i)-1
  keys(3,i) = 2*(i)
enddo

! We did the sweeps A over f[middle] in the S stage, so we need only do the A sweeps on [yb;0;yt]
! That is, the vector we are solving for has the top and bottom tips equal to the things we got
! from the reduced system, and the middle all zeroes. 
if(nbpart-(spm(23)+2) > 0) then
  call T_ALLOC(wallocate_2)(g,Ajmin(nbpart)-Ajmin(nbpart-1),nrhs*(nbpart-(spm(23)+2)),infoloc)
endif

!g2 is the temp vector for the 2spike threads. It must be shared, because it is used to communicate between them.
if(spm(23) > 0) then
  call T_ALLOC(wallocate_2)(g2,Ajmin(3)-Ajmin(2),nrhs*spm(23),infoloc)
endif

parallel_start_time_2 = OMP_GET_WTIME()
!$omp parallel do default(shared) shared(g2,g) private(j,mystart,p,mysize,t1,t2,tim)
do i=1,nthread,1
  timing_array(i,2) = OMP_GET_WTIME()
  call threadmap(i,p,nbpart,spm(23))
  mysize = Ajmin(p+1)-Ajmin(p)
  if(p == 1) then
    call T_DEC(TBSM)('U',trans,'N', klu, nrhs, ku, A(1,Ajmin(p+1)-klu), lda, grj(klu+1,1,p), 2*klu)

    do j=1,nrhs,1
!      call T_DEC(AXPY)(klu,ONE_PREC,grj(klu+1,j,p),1,f(Ajmin(p)+mysize-klu,j),1)
      call T_DEC(AXPY)(klu,ONE_PREC,grj(klu+1,j,p),1,f(Ajmin(p+1)-klu,j),1)
    enddo

    call T_DEC(TBSM)('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, f(1,1), ldf)
  endif

  if(p > 1 + spm(23) .and. p < nbpart) then 
    mystart = nrhs*(p-(2+spm(23)))+1
    g(:,mystart:mystart+nrhs-1) = ZERO_PREC
    call T_DEC(LACPY)('F', klu, nrhs, grj(1,1,p), 2*klu, g(1,mystart), mysize)
    call T_DEC(LACPY)('F', klu, nrhs, grj(klu+1,1,p), 2*klu, g(mysize-klu+1,mystart), mysize)
    call T_DEC(TBSM)('U',trans,'N', mysize, nrhs, ku, A(1,Ajmin(p)),    lda, g(1,mystart), mysize)

!DTBSMPY combines the solve and sum, to save us trips to memory.
    call T_DEC(TBSMPY)('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, g(1,mystart), mysize, f(Ajmin(p),1), ldf)
! Keeping this around just in case...
!    call T_DEC(TBSM)('L',trans,'U', mysize, nrhs, kl, A(ku+1,Ajmin(p)), lda, g(1,mystart), mysize)
!    do j=1,nrhs,1
!      call T_DEC(AXPY)(mysize,ONE_PREC,g(1,mystart+j-1),1,f(Ajmin(p),j),1)
!    enddo
!    deallocate(g)
  endif

  if(p==nbpart) then 
    call T_DEC(TBSM)('L',trans,'N',klu,nrhs,kl,A(ku+1,Ajmin(p)),lda,grj(1,1,p),2*klu)
    do j=1,nrhs,1
      call T_DEC(AXPY)(klu,ONE_PREC,grj(1,j,p),1,f(Ajmin(p),j),1)
    enddo
    call T_DEC(TBSM)('U',trans,'U',mysize,nrhs,ku,A(1,Ajmin(p)),   lda,f(Ajmin(p),1),ldf)
  endif

  if(p > 1 .and. p <= 1 + spm(23)) then 
    mystart = nrhs*(p-2)+1

    if(keys(2,p-1) == omp_get_thread_num()) then
      g2(1:mysize/2,mystart:mystart+nrhs-1) = ZERO_PREC
      call T_DEC(LACPY)('F', klu, nrhs, grj(1,1,p), 2*klu, g2(1,mystart), mysize)
    endif

    if(keys(3,p-1) == omp_get_thread_num()) then
      g2(mysize/2+1:mysize,mystart:mystart+nrhs-1) = ZERO_PREC
      call T_DEC(LACPY)('F', klu, nrhs, grj(klu+1,1,p), 2*klu, g2(mysize-klu+1,mystart), mysize)
    endif

    call T_DEC(spike_GBTRSk2)(.false.,trans,Ajmin(p+1)-Ajmin(p), nrhs, kl, ku, A(1,Ajmin(p)), lda, g2(1,mystart), mysize, ipiv_dummy, keys(1,p-1), spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1),.false.)
!    call T_DEC(spike_simple_solve_transpose)(Ajmin(p+1)-Ajmin(p),trans,nrhs,kl,ku,A(1,Ajmin(p)),lda,f(Ajmin(p),1), ldf,info,keys(1,p-1),spike2_red(1,(p-2)*2*klu+1),spike2_grj(1,(p-2)*nrhs+1))

    if(keys(2,p-1) == omp_get_thread_num()) then
      do j=1,nrhs,1
        call T_DEC(AXPY)(mysize/2,ONE_PREC,g2(1,mystart-1+j),1,f(Ajmin(p),j),1)
      enddo
    endif
    if(keys(3,p-1) == omp_get_thread_num()) then
      do j=1,nrhs,1
        call T_DEC(AXPY)(mysize/2,ONE_PREC,g2(mysize/2+1,mystart-1+j),1,f(Ajmin(p)+mysize/2,j),1)
      enddo
    endif

  endif

  timing_array(i,2) = OMP_GET_WTIME() - timing_array(i,2)
enddo
parallel_end_time_2 = OMP_GET_WTIME()


if(spm(23) > 0) then
  call T_ALLOC(wdeallocate_2)(g2)
endif

if(nbpart-(spm(23)+2) > 0) then
  call T_ALLOC(wdeallocate_2)(g)
endif

call T_ALLOC(wdeallocate_2)(spike2_red)
call T_ALLOC(wdeallocate_2)(spike2_grj)
call T_ALLOC(wdeallocate_3)(grj)
call wdeallocate_2i(keys)
call wdeallocate_1i(Ajmin)
endif

end_time = OMP_GET_WTIME()

if(spm(1) .eq. 1) then
  call wwrite_s(fout,'------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Solve Time (Transpose)|')
  call wwrite_n(fout)
  call wwrite_s(fout,'------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|  Solve Sweeps on blocks    |----------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'| Thread | Partition |         Time 1        |          Time 2       |')
  call wwrite_n(fout)
  do i=1,spm(22)
    call threadmap(i,p,nbpart,spm(23))
    call wwrite_s(fout,'|    ')
    call wwrite_i(fout,i)
    call wwrite_s(fout,'   |     ')
    call wwrite_i(fout,p)
    call wwrite_s(fout,'     | ')
    call wwrite_d(fout,timing_array(i,1))
    call wwrite_s(fout,' | ')
    call wwrite_d(fout,timing_array(i,2))
    call wwrite_s(fout,' |')
    call wwrite_n(fout)
  enddo 
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Blocks Solve        | ')
  call wwrite_d(fout,(parallel_end_time_1-parallel_start_time_1))
  call wwrite_s(fout,' | ')
  call wwrite_d(fout,(parallel_end_time_2-parallel_start_time_2))
  call wwrite_s(fout,' |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Reduced System Solve| ')
  call wwrite_d(fout,(redsys_end_time-redsys_start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'|Overall Solve       | ')
  call wwrite_d(fout,(end_time-start_time))
  call wwrite_s(fout,'                         |')
  call wwrite_n(fout)
  call wwrite_s(fout,'----------------------------------------------------------------------')
  call wwrite_n(fout)
endif

call wdeallocate_2d(timing_array)

end subroutine T_DEC(spike_gbtrst2)



!Iterative refinment
subroutine T_DEC(spike_itrefinement)(spm,trans,n,kl,ku,nrhs,A,work,lda,C,ldc,B,oB,ldb)
!=====================================================================
!  Braegan Spring  2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: maxiter,n,kl,ku,lda,ldb,nrhs,ldc
MAIN_TYPE, dimension(lda,*) :: A
MAIN_TYPE, dimension(ldb,*) :: B
MAIN_TYPE, dimension(ldb,*) :: oB
MAIN_TYPE, dimension(ldc,*) :: C
character :: trans
REAL_TYPE,dimension(:),  pointer :: normB
REAL_TYPE,dimension(:),  pointer :: normRes
integer,  dimension(:),pointer :: Ajmin
!double precision, dimension(ldb,nrhs) :: res
MAIN_TYPE, dimension(:,:), pointer :: res
integer :: i,j,p
integer :: work_complete
MAIN_TYPE, dimension(*) :: work
REAL_TYPE :: res_max
REAL_TYPE, dimension(:),pointer  :: res_max_temp
integer :: infoloc
integer(8),parameter :: fout=6
integer :: norm
!double precision, external :: OMP_GET_WTIME
double precision :: t1,t2


norm=spm(14)


j=0
work_complete = 0
maxiter = spm(11)

call wallocate_1i(Ajmin,spm(20)+1,infoloc)
call T_ALLOC(wallocate_2)(res,ldb,nrhs,infoloc)
call T_REAL_ALLOC(wallocate_1)(normB      ,nrhs          ,infoloc)
call T_REAL_ALLOC(wallocate_1)(normRes    ,nrhs          ,infoloc)

call spikerl_calc_size_partitions(spm,Ajmin,n)

!Each right hand side should have an individual norm. 

call T_DEC(spike_vector_norm)(spm,norm,oB,n,nrhs,normB)


call T_REAL_ALLOC(wallocate_1)(res_max_temp,spm(20),infoloc)

do while ((j<maxiter) .and. (work_complete .eq. 0)) 
  j=j+1

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    do i=1,nrhs
      res(Ajmin(p):Ajmin(p+1)-1,1:nrhs)=oB(Ajmin(p):Ajmin(p+1)-1,1:nrhs)
    enddo
  enddo
  res_max_temp = ONE_REAL_PREC

!How far are we from the desired value? 
! Find RelativeResidual = norm(oB-AB)/norm(oB) 
! (since B contains the result) 

! Set res = oB-AB
t1 = OMP_GET_WTIME()
  call T_DEC(spike_matmul)(spm,Ajmin,trans,kl,ku,nrhs,C,lda,B,ldb,res)
t2 = OMP_GET_WTIME()
! RRES = norm(res)/norm(oB)

  call T_DEC(spike_vector_norm)(spm,norm,res,n,nrhs,normRes)
  res_max = ZERO_REAL_PREC
  do i=1,nrhs
    res_max = max(res_max,normRes(i)/normB(i))
  enddo

  if(spm(1) .eq. 1) then
    call wwrite_s(fout,'Num Refinements so far: ')
    call wwrite_i(fout,j-1)
    call wwrite_n(fout)
    call wwrite_s(fout,'Current relative residual norm: ')
    call T_REAL_PRINT(wwrite_)(fout,res_max)
    call wwrite_n(fout)
  endif

  if(-log10(abs(res_max)) >= spm(SPM_PREC_ENTRY)) then 
   work_complete=1
    if(spm(1) .eq. 1) then
      call wwrite_s(fout,'Desired Residual Reached')
      call wwrite_n(fout)
    endif
  endif

 
  if(work_complete .eq. 0) then
    call T_DEC(spike_gbtrs)(spm,trans,n,kl,ku,nrhs,A,lda,work,res,ldb)
  
    !$OMP PARALLEL DO default(shared) private(i)
    do p=1,spm(20)
      do i=1,nrhs
        call T_DEC(AXPY)(Ajmin(p+1)-Ajmin(p), ONE_PREC, res(Ajmin(p),i), 1, B(Ajmin(p),i), 1)
      enddo
    enddo
  endif

enddo
  call T_REAL_ALLOC(wdeallocate_1)(normB)
  call T_REAL_ALLOC(wdeallocate_1)(normRes) 
  call T_REAL_ALLOC(wdeallocate_1)(res_max_temp)
  call T_ALLOC(wdeallocate_2)(res)
  call wdeallocate_1i(Ajmin)

end subroutine T_DEC(spike_itrefinement)




subroutine T_DEC(spike_matmul)(spm,Ajmin,trans,kl,ku,nrhs,C,ldc,B,ldb,res)
!=====================================================================
!  Braegan Spring - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer, dimension(64) :: spm
integer :: kl,ku,ldb,nrhs,ldc
MAIN_TYPE, dimension(ldb,*) :: B
MAIN_TYPE, dimension(ldc,*) :: C
character :: trans
integer, dimension(spm(20)+1) :: Ajmin
MAIN_TYPE, dimension(ldb,*) :: res
MAIN_TYPE,  dimension(:,:),pointer :: temp_upper
MAIN_TYPE,  dimension(:,:),pointer :: temp_lower
integer :: i,p
integer :: infoloc

if(trans .eq. 'n' .or. trans .eq. 'N') then 

  call T_ALLOC(wallocate_2)(temp_upper,ku,spm(20)*nrhs,infoloc)
  call T_ALLOC(wallocate_2)(temp_lower,kl,spm(20)*nrhs,infoloc)
  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    if(p<spm(20)) then
      call T_DEC(LACPY)('F',ku,nrhs,B(Ajmin(p+1), 1),ldb,temp_upper(1,(p-1)*nrhs+1),ku)
    endif
    if(p>1) then
      call T_DEC(LACPY)('F',kl,nrhs,B(Ajmin(p)-kl,1),ldb,temp_lower(1,(p-1)*nrhs+1),kl)
    endif
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    !Use the fancy structurally symmetric multiply if we can, otherwise use the general multiply. 
    if(kl == ku) then
      call T_DEC(SBMM)('F', Ajmin(p+1)-Ajmin(p), nrhs, KL, -ONE_PREC, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,ONE_PREC, res(Ajmin(p),1), LDB)
    else
      call T_DEC(GBMM)(trans,'N', Ajmin(p+1)-Ajmin(p), nrhs, KL, KU, -ONE_PREC, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,ONE_PREC, res(Ajmin(p),1), LDB)
    endif
!The triangular bits
    if(p<spm(20)) then
      call T_DEC(TRMM)('L','L',trans,'N', ku, nrhs, ONE_PREC, C(1,Ajmin(p+1)), ldc-1, temp_upper(1,(p-1)*nrhs+1),ku)
      do i=1,nrhs 
        call T_DEC(AXPY)(ku, -ONE_PREC, temp_upper(1,(p-1)*nrhs+i), 1, res(Ajmin(p+1)-ku,i), 1)
      enddo
    endif
    if(p>1) then
      call T_DEC(TRMM)('L','U',trans,'N', kl, nrhs, ONE_PREC, C(ku+kl+1,Ajmin(p)-kl), ldc-1, temp_lower(1,(p-1)*nrhs+1),kl)
      do i=1,nrhs 
        call T_DEC(AXPY)(kl, -ONE_PREC, temp_lower(1,(p-1)*nrhs+i), 1, res(Ajmin(p),i), 1)
      enddo
    endif
  enddo

else
  call T_ALLOC(wallocate_2)(temp_upper,kl,spm(20)*nrhs,infoloc)
  call T_ALLOC(wallocate_2)(temp_lower,ku,spm(20)*nrhs,infoloc)

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    if(p<spm(20)) then
      call T_DEC(LACPY)('F',kl,nrhs,B(Ajmin(p+1), 1),ldb,temp_upper(1,(p-1)*nrhs+1),kl)
    endif
    if(p>1) then
      call T_DEC(LACPY)('F',ku,nrhs,B(Ajmin(p)-ku,1),ldb,temp_lower(1,(p-1)*nrhs+1),ku)
    endif
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
    call T_DEC(GBMM)(trans,'N', Ajmin(p+1)-Ajmin(p), nrhs, KU, KL, -ONE_PREC, C(1,Ajmin(p)), ldc, B(Ajmin(p),1),LDB,ONE_PREC, res(Ajmin(p),1), LDB)
  enddo

  !$OMP PARALLEL DO default(shared) private(i)
  do p=1,spm(20)
!The triangular bits
    if(p>1) then
    !Multiply by the upper band-- so this should be the stuff that was pulled from the current partition (because transpose)
      call T_DEC(TRMM)('L','L',trans,'N', ku, nrhs, ONE_PREC, C(1,Ajmin(p)), ldc-1, temp_lower(1,(p-1)*nrhs+1),ku)
      do i=1,nrhs 
        call T_DEC(AXPY)(ku, -ONE_PREC, temp_lower(1,(p-1)*nrhs+i), 1, res(Ajmin(p),i), 1)
      enddo
    endif

    if(p<spm(20)) then
      call T_DEC(TRMM)('L','U',trans,'N', kl, nrhs, ONE_PREC, C(ku+kl+1,Ajmin(p+1)-kl), ldc-1, temp_upper(1,(p-1)*nrhs+1),kl)
      do i=1,nrhs 
        call T_DEC(AXPY)(kl, -ONE_PREC, temp_upper(1,(p-1)*nrhs+i), 1, res(Ajmin(p+1)-kl,i), 1)
      enddo
    endif
  enddo
endif

call T_ALLOC(wdeallocate_2)(temp_upper)
call T_ALLOC(wdeallocate_2)(temp_lower)

end subroutine T_DEC(spike_matmul)





subroutine T_DEC(spike_multi)(klu,nrhs,A,xb,ldb,xt,ldt)
! Purpose:
! Solve of the reduced system A (different format possible) by x (xb,xt), result in x
! This is an internally used function, unlikely to be useful for the user.
! --------
! Arguments: 
!       A (input) real(kind=(kind(1.0d0))) array, input matrix 
!       xb (input/output) real(kind=(kind(1.0d0))) array, the multiplier and bottom half
!       of the result 
!       xt (input/output) real(kind=(kind(1.0d0))) array, the multiplier and the top
!       half of the result  
!       ldb,ldt leading dimensions
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer:: klu,nrhs,ldb,ldt
MAIN_TYPE, dimension(2*klu,*):: A
MAIN_TYPE, dimension(ldb,*) :: xb
MAIN_TYPE, dimension(ldt,*) :: xt
integer:: Info_lap
MAIN_TYPE, dimension(:,:),pointer :: A_aux, f
integer, dimension(:),pointer :: IPIV
character(len=1) :: TRANSA,TRANSB,SIDE
integer :: infoloc

TRANSA='N'
TRANSB = 'N'
SIDE='L'



!TODO Fix this 
if(.true.) then
  ! I-VW storage
  
  call wallocate_1i(IPIV,klu,infoloc)
  call T_ALLOC(wallocate_2)(A_aux,klu,klu,infoloc)
!!! step (1) solve (I-WV)xt=xt-W*xb
  call T_DEC(GEMM)(TRANSA, TRANSB, klu, nrhs, klu, -ONE_PREC, A(klu+1,1), 2*klu, xb, ldb,ONE_PREC, xt, ldt) 
!!! step (1) solve (I-WV)xt=xt-W*xb
  call T_DEC(LACPY)('A',klu,klu,A,2*klu,A_aux,klu)

  call T_DEC(GESV)( klu, nrhs, A_aux, klu, IPIV, xt, ldt, INFO_lap) 
!!! step (2) compute xb=xb-Vxt
  call T_DEC(GEMM )(TRANSA, TRANSB, klu, nrhs, klu, -ONE_PREC, A(1,klu+1), 2*klu, xt, ldt,ONE_PREC, xb, ldb) 

  call wdeallocate_1i(IPIV)
  call T_ALLOC(wdeallocate_2)(A_aux)
else
  ! Explicit storage
  call wallocate_1i(IPIV,2*klu,infoloc)
  call T_ALLOC(wallocate_2)(A_aux,2*klu,2*klu,infoloc)
  call T_ALLOC(wallocate_2)(f,2*klu,nrhs,infoloc)
!!! step (1) solve (I-WV)xt=xt-W*xb
!!! step (1) solve (I-WV)xt=xt-W*xb
  call T_DEC(LACPY)('A',2*klu,2*klu,A,2*klu,A_aux,2*klu)
  call T_DEC(LACPY)('A',  klu, nrhs,xt,ldt,f(    1,1),2*klu)
  call T_DEC(LACPY)('A',  klu, nrhs,xb,ldb,f(klu+1,1),2*klu)

  call T_DEC(GESV)( 2*klu, nrhs, A_aux, 2*klu, IPIV, f, 2*klu, INFO_lap) 

!!! step (2) compute xb=xb-Vxt

  call T_DEC(LACPY)('A', klu, nrhs,f(    1,1),2*klu,xt,ldt)
  call T_DEC(LACPY)('A', klu, nrhs,f(klu+1,1),2*klu,xb,ldb)
 
  call T_ALLOC(wdeallocate_2)(A_aux)
  call T_ALLOC(wdeallocate_2)(f)
  call wdeallocate_1i(IPIV)
  
endif

end subroutine  T_DEC(spike_multi)





subroutine T_DEC(spike_invred)(klu,V,W,redA)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Purpose: Invert  or store the truncated reduced system by blocs
! matrix to invert as the following form
!
!                |I V|^(-1) or |I-WV       V|
!                |W I|         |W          I|
! Arguments: 
!       V (input) MAIN_TYPE array, 1x2 block of input matrix 
!       W (input) MAIN_TYPE array, 2x1 block of the input matrix 
!       readA (output) MAIN_TYPE array, result matrix  
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
integer :: klu,i
MAIN_TYPE, dimension(2*klu,*):: redA
MAIN_TYPE, dimension(2*klu,*):: V
MAIN_TYPE, dimension(2*klu,*):: W
character(len=1) :: TRANSA,TRANSB
!character(len=100) :: format
!format = '(F6.2, F6.2, F6.2, F6.2, F6.2, F6.2, F6.2, F6.2, F6.2, F6.2)'

!TODO Fix this 
if(.true.) then
redA(1:klu,1:klu)             = ZERO_PREC
redA(klu+1:2*klu,klu+1:2*klu) = ZERO_PREC
redA(klu+1:2*klu,1:klu)       = W(1:klu,1:klu)
redA(1:klu,klu+1:2*klu)       = V(klu+1:2*klu,1:klu)


do i=1,klu
  redA(i,i)=ONE_PREC
  redA(klu+i,klu+i)=ONE_PREC
end do


!! Step I-WV
TRANSA='N'
TRANSB='N'
call T_DEC(GEMM )( TRANSA, TRANSB, klu, klu, klu, -ONE_PREC, redA(klu+1,1), 2*klu, redA(1,klu+1), 2*klu,ONE_PREC, redA(1,1), 2*klu ) 

else

!W(1:2*klu,1:2*klu)=-1.0d0
!V(1:2*klu,1:2*klu)=1.0d0

!redA(klu+1:2*klu,1:klu)=W(1:klu,1:klu)
!redA(1:klu,klu+1:2*klu)=V(klu+1:2*klu,1:klu)

redA(1:klu,1:klu)             = ZERO_PREC
redA(klu+1:2*klu,klu+1:2*klu) = ZERO_PREC
redA(klu+1:2*klu,1:klu)       = V(klu+1:2*klu,1:klu)
redA(1:klu,klu+1:2*klu)       = W(1:klu,1:klu)

do i=1,2*klu
  redA(i,i)=ONE_PREC
end do


endif

end subroutine T_DEC(spike_invred)




subroutine T_DEC(spike_prep_recn)(kl,ku,rV,rW,redA,nbpart,nblevel) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Prepares the reduced system. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer :: nbpart,nblevel 
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
integer :: j,sub,nbsub,l,i,kl,ku,klu
integer:: p
MAIN_TYPE,dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rV,rW
MAIN_TYPE,dimension(2*max(kl,ku),2*max(kl,ku),*) :: redA
MAIN_TYPE, dimension(:,:,:), pointer :: aux1,aux2
character(len=1) :: TRANSA,TRANSB
integer :: t1,t2,tim,t3,t4
integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
klu=max(kl,ku)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
call T_ALLOC(wallocate_3)(aux1,klu,klu,nbpart,infoloc)
call T_ALLOC(wallocate_3)(aux2,klu,klu,nbpart,infoloc)
!!!!!
!!! We associate to each CPU 1.0d0 particular reduce system (which is required for the recursion)
!!! Example(nbpart=8) for level 1 we get the reduce systems for CPU=1,3,5,7  
    !!                    for level 2 we get the reduce systems for CPU=2,6    
    !!                    for level 3 we get the reduce systems for CPU=4 
!!!!!

TRANSA='N'
TRANSB ='N'


do j=1,nblevel-1
!!! compute the spikes level j
  nbsub=2**(nblevel-j) !! Number of subsystem to solve for the level j
  sub=2**(j-1) !! Number of spikes by subsystem

!$omp parallel do private(l,p,t1,t2,t3,t4,tim)
  do i=1,nbpart,1
    p=i
    if ((mod(p-sub,nbpart/nbsub) == 0) .and. p < nbpart-nbpart/nbsub) then
      aux1(:,:,p)=ZERO_PREC
      call T_DEC(LACPY)('A',klu,klu,rV(1,1,j,p+1),2*klu,aux2(1,1,p),klu)
      call T_DEC(spike_multi)(klu,klu,redA(1,1,p),aux1(1,1,p),klu,aux2(1,1,p),klu) !solution of the  truncated reduce system 

      if (j>1) then
        do l=p-sub+1,p-1  !! CPU(/=i) involved in the calculation of the new V spikes (Vtop)
            !send i -> l  aux2 
            !recv l <- i  aux2
          call T_DEC(GEMM)( TRANSA, TRANSB, 2*klu, klu, klu, -ONE_PREC, rV(1,1,j,l), 2*klu, aux2(1,1,p), klu,ZERO_PREC, rV(1,1,j+1,l), 2*klu )
        enddo
      endif
      call T_DEC(GEMM)( TRANSA, TRANSB, klu, klu, klu, -ONE_PREC, rV(1,1,j,p), 2*klu, aux2(1,1,p), klu,ZERO_PREC, rV(1,1,j+1,p), 2*klu )
      call T_DEC(LACPY)('A',klu,klu,aux1(1,1,p),klu,rV(klu+1,1,j+1,p),2*klu)

      do l=p+1,p+sub !! CPU(/=i) involved in the calculation of the new V spikes (Vbottom)
        !send i -> l  aux1
        !recv l <- i  aux1
        if (l/=p+1) then
          call T_DEC(LACPY)('A',2*klu,klu,rV(1,1,j,l),2*klu,rV(1,1,j+1,l),2*klu)
          call T_DEC(GEMM)( TRANSA, TRANSB, 2*klu, klu, klu, -ONE_PREC, rW(1,1,j,l), 2*klu, aux1(1,1,p), klu,ONE_PREC,rV(1,1,j+1,l), 2*klu )
        else
          !send i -> i+1  aux2 
          !recv i+1 <- i  aux2
          call T_DEC(LACPY)('A',klu,klu,aux2(1,1,p),klu,rV(1,1,j+1,l),2*klu)
          call T_DEC(LACPY)('A',klu,klu,rV(klu+1,1,j,l),2*klu,rV(klu+1,1,j+1,l),2*klu)
          call T_DEC(GEMM)( TRANSA, TRANSB, klu, klu, klu, -ONE_PREC, rW(klu+1,1,j,l), 2*klu, aux1(1,1,p), klu,ONE_PREC, rV(klu+1,1,j+1,l), 2*klu )
        endif
      enddo
    endif

    p=i-1
    if ((p>= sub+nbpart/nbsub) .and. (mod(p-sub,nbpart/nbsub)==0)) then
      aux2(:,:,p-1)=ZERO_PREC
      call T_DEC(LACPY)('A',klu,klu,rW(klu+1,1,j,p),2*klu,aux1(1,1,p-1),klu)
      call T_DEC(spike_multi)(klu,klu,redA(1,1,p),aux1(1,1,p-1),klu,aux2(1,1,p-1),klu) !! solution of the  truncated reduce system 

      if(j>1) then
        do l=p-sub+1,p-1  !! CPU(/=i) involved in the calculation of the new W spikes (Wtop)
          !send i -> l  aux4 
          !recv l <- i  aux4
          call T_DEC(LACPY)('A',2*klu,klu,rW(1,1,j,l),2*klu,rW(1,1,j+1,l),2*klu)
          call T_DEC(GEMM )( TRANSA, TRANSB, 2*klu, klu, klu, -ONE_PREC, rV(1,1,j,l), 2*klu, aux2(1,1,p-1), klu,ONE_PREC, rW(1,1,j+1,l), 2*klu )
        enddo
      endif

      call T_DEC(LACPY)('A',klu,klu,rW(1,1,j,p),2*klu,rW(1,1,j+1,p),2*klu)
      call T_DEC(GEMM)( TRANSA, TRANSB, klu, klu, klu, -ONE_PREC, rV(1,1,j,p), 2*klu, aux2(1,1,p-1), klu,ONE_PREC, rW(1,1,j+1,p), 2*klu )
      call T_DEC(LACPY)('A',klu,klu,aux1(1,1,p-1),klu,rW(klu+1,1,j+1,p),2*klu)

      do l=p+1,p+sub !! CPU(/=i) involved in the calculation of the new W spikes (Wbottom)
        !send i -> l  aux3
        !recv l <- i  aux3
        if (l/=p+1) then
          call T_DEC(GEMM )( TRANSA, TRANSB, 2*klu, klu, klu, -ONE_PREC, rW(1,1,j,l), 2*klu, aux1(1,1,p-1), klu,ZERO_PREC, rW(1,1,j+1,l), 2*klu ) 
        else
          !send i -> i+1  aux4 
          !recv i+1 <- i  aux4
          call T_DEC(LACPY)('A',klu,klu,aux2(1,1,p-1),klu,rW(1,1,j+1,l),2*klu)
          call T_DEC(GEMM )( TRANSA, TRANSB, klu, klu, klu, -ONE_PREC, rW(klu+1,1,j,l), 2*klu, aux1(1,1,p-1), klu,ZERO_PREC, rW(klu+1,1,j+1,l), 2*klu )
        endif
      enddo
    endif
  enddo

!!!<<< all threads synchronize before computing reduce system for level j+1
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! New reduce systems for level j+1 !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  nbsub=2**(nblevel-j-1) !! Number of subsystem to solve for the level j
  sub=2**j !! Number of spikes by subsystem

!$omp parallel do default(shared) private(p)
  do i=1,nbpart
    p=i
    if ((mod(p-sub,nbpart/nbsub) == 0)) then
!!!!! Inverse of the truncated reduced system  
      call T_DEC(spike_invred)(klu,rV(1,1,j+1,p),rW(1,1,j+1,p+1),redA(1,1,p)) !! give the inverse of the reduce system in the new bloc numbering
    endif
  enddo
enddo

! Remove Fortran runtime dependency
call T_ALLOC(wdeallocate_3)(aux1)
call T_ALLOC(wdeallocate_3)(aux2)
! End of removal

end subroutine T_DEC(spike_prep_recn)





subroutine T_DEC(spike_solve_recn)(kl,ku,nrhs,rV,rW,redA,nbpart,nblevel,grj)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! Solves the reduced system
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer :: nbpart,nblevel
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads
integer :: j,sub,nbsub,l,i,kl,ku,klu,nrhs
character(len=1) :: TRANSA,TRANSB
MAIN_TYPE, dimension(:,:,:),pointer :: aux1,aux2
MAIN_TYPE, dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rV,rW
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku),*) :: redA
MAIN_TYPE, dimension(2*max(kl,ku),nrhs,*) :: grj
integer :: infoloc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
klu=max(kl,ku)

! Remove Fortran runtime dependency
call T_ALLOC(wallocate_3)(aux1,klu,nrhs,nbpart,infoloc)
call T_ALLOC(wallocate_3)(aux2,klu,nrhs,nbpart,infoloc)
        
! End of removal
TRANSA='N'
TRANSB='N'

do j=1,nblevel

!!! compute the spikes level j
  nbsub=2**(nblevel-j) !! Number of subsystem to solve for the level j
  sub=2**(j-1) !! Number of spikes by subsystem
  aux1=ZERO_PREC
  aux2=ZERO_PREC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! CALCULATION INVOLVING  V and W SPIKES at level j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !$omp parallel do
  do i=sub,nbpart,nbpart/nbsub

    call T_DEC(LACPY)('A',klu,nrhs,grj(1,1,i+1),  2*klu,aux2(1,1,i),klu)
    call T_DEC(spike_multi)(klu,nrhs,redA(1,1,i),grj(klu+1,1,i),2*klu,aux2(1,1,i),klu) !! solution of the  truncated reduce system 
    call T_DEC(LACPY)('A',klu,nrhs,grj(klu+1,1,i),2*klu,aux1(1,1,i),klu)

  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel do
  do i=sub,nbpart,nbpart/nbsub
    do l=i-sub+1,i-1  !! CPU(/=i) involved in the calculation with the V spikes 
             !send i -> l  aux2 
             !recv l <- i  aux2
      call T_DEC(GEMM )( TRANSA, TRANSB, 2*klu, nrhs, klu, -ONE_PREC, rV(1,1,j,l), 2*klu, aux2(1,1,i), klu,ONE_PREC, grj(1,1,l), 2*klu )
    end do
    call T_DEC(GEMM )( TRANSA, TRANSB, klu, nrhs, klu, -ONE_PREC, rV(1,1,j,i), 2*klu, aux2(1,1,i), klu,ONE_PREC, grj(1,1,i), 2*klu )
    do l=i+1,i+sub  !! CPU(/=i) involved in the calculation of the W spikes 
             !send i -> l  aux1
             !recv l <- i  aux1
      if (l/=i+1) then
                   !         grj(:,:)=grj(:,:)-matmul(sW(:,:,j),aux1) 
        call T_DEC(GEMM )( TRANSA, TRANSB, 2*klu, nrhs, klu, -ONE_PREC, rW(1,1,j,l), 2*klu, aux1(1,1,i), klu,ONE_PREC, grj(1,1,l), 2*klu )
      else
                   !send i -> i+1  aux2 
                   !recv i+1 <- i  aux2
        call T_DEC(LACPY)('A',klu,nrhs,aux2(1,1,i),klu,grj(1,1,l),2*klu)
        call T_DEC(GEMM )( TRANSA, TRANSB, klu, nrhs, klu, -ONE_PREC, rW(klu+1,1,j,l), 2*klu, aux1(1,1,i), klu,ONE_PREC, grj(klu+1,1,l), 2*klu )
      end if
    end do
  end do
enddo

! Remove Fortran runtime dependency 
call T_ALLOC(wdeallocate_3)(aux1)
call T_ALLOC(wdeallocate_3)(aux2)
! End of removal
end subroutine T_DEC(spike_solve_recn)



subroutine T_DEC(spike_solve_recn_transpose)(trans,kl,ku,nrhs,rV,rW,redA,nbpart,nblevel,grj) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Transpose solve for the reduced system 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - 2018
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
character :: trans
integer   :: klu,kl,ku,nrhs
integer   :: nblevel, nbpart, nbsub, p
integer   :: i,j,k
MAIN_TYPE, dimension(2*max(kl,ku),max(kl,ku),nblevel,*) :: rW,rV
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku),*)       :: redA
MAIN_TYPE, dimension(2*max(kl,ku),nrhs,*)               :: grj

klu=max(kl,ku)
! grj holds the right hand side
! It is broken up by the number of partitons, because each partition works on a section at a time.
! redA holds the reduced systems 
! These are all the same size, and there is one per partition. 

! rV and rW hold the spike tips.
! The storage here is a bit weird
! The spikes for each level are part of the D matrix for the next level.
! They double in size as we move up the levels, starting at a size of 2xklu
! But these tips are stored in a weird way, because the dimensions of rV and rW are 2*klu x klu...
! Because this is the transpose case, we're going 'backwards' across the levels of the recursive systems. 

do i=nblevel,1,-1
! nsub is the number of subsystems that we'll be working on for this level. 
! it starts at 1, because the first matrix we work on is s_3 (because everything has been reversed in the transpose case), which has only one reduced system (in the middle)

! p = (2*j-1)*2**(i-1) maps j and i to the partiton (subsystem), p 
! The idea:
! Take the number of partitions total, and divide that by the number of partitions we'll be working on for this level to get the number of
! partitions we step over to get from one "partition we are working on in this level" to the next. 
! The starting partition is gotten by dividing the step size by two. 
! But since everything is a power of two, it can be done as below. 
! ... it makes sense if you draw it on graph paper!
  nbsub = 2**(nblevel-i)

! Build the modified right hand side 
  !$OMP PARALLEL DO default(none) private(p,k) shared(trans,i,klu,nrhs,rV,grj,nbpart,nbsub,rW)
  do j=1,nbsub

    p = (2*j-1)*2**(i-1)
    ! Now we have to build the modified right hand side.
    call T_DEC(GEMM)(trans,'N',klu,nrhs,  klu,-ONE_PREC,rV(1,1,i,p),2*klu,grj(1,1,p),2*klu,ONE_PREC,grj(1,1,p+1),2*klu)
    do k=1,nbpart/(2*nbsub)-1,1
      call T_DEC(GEMM)(trans,'N',klu,nrhs,2*klu,-ONE_PREC,rV(1,1,i,p-k),2*klu,grj(1,1,p-k),2*klu,ONE_PREC,grj(1,1,p+1),2*klu)
    enddo
  enddo

  !$OMP PARALLEL DO default(none) private(p,k) shared(trans,i,klu,nrhs,rV,grj,nbpart,nbsub,rW)
   do j=1,nbsub
    p = (2*j-1)*2**(i-1)
    call T_DEC(GEMM)(trans,'N',klu,nrhs,  klu,-ONE_PREC,rW(klu+1,1,i,p+1),2*klu,grj(klu+1,1,p+1),2*klu,ONE_PREC,grj(klu+1,1,p),2*klu)
    do k=1,nbpart/(2*nbsub)-1,1
      call T_DEC(GEMM)(trans,'N',klu,nrhs,2*klu,-ONE_PREC,rW(1    ,1,i,p+k+1),2*klu,grj(1,1  ,p+k+1),2*klu,ONE_PREC,grj(klu+1,1,p),2*klu)
    enddo
  enddo
  !$OMP PARALLEL DO default(shared) private(p)
  do j=1,nbsub
    p = (2*j-1)*2**(i-1)
    call T_DEC(spike_multi_transpose)(trans,klu,nrhs,redA(1,1,p),grj(1,1,p+1),2*klu,grj(klu+1,1,p),2*klu)
  enddo
enddo

end subroutine T_DEC(spike_solve_recn_transpose)




subroutine T_DEC(spike_multi_transpose)(trans,klu,nrhs,red,xb,ldxb,xt,ldxt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Recursive spike calls SPIKE on the reduced system for the main system.
! This solves the reduced system for the SPIKE calls on the main reduced system.
!!!!!!!!!!!!!!!!!!!!!!!!!!!
!=====================================================================
!  Braegan Spring - 2018
!=====================================================================
implicit none
include 'f90_noruntime_interface.fi'
integer   :: klu, nrhs,ldxb,ldxt
integer   :: info
character :: trans
integer, dimension(klu)                   :: IPIV
MAIN_TYPE, dimension(2*klu,*)    :: xb,xt
MAIN_TYPE, dimension(2*klu,*) :: red
!MAIN_TYPE, dimension(klu,*)     :: red_aux 
MAIN_TYPE, pointer, dimension(:,:)     :: red_aux 
integer :: infoloc

call T_ALLOC(wallocate_2)(red_aux,klu,klu,infoloc)
! The mathematical format of the reduced system after block LU is 
!(where w and v are actually the tips of their respective spikes):
! __   __    __       __
! |I  0 |    | I    V  |
! |W  I |    | 0  I-WV |
! --   --    --       --
! Which transposes to 
! __         __    __     __
! |I  0       |    | I  WT | = LU
! |VT (I-WV)T |    | 0  I  |
! --         --    --     --
! But it is actually stored with (I-WV) in the top-left corner
! Because of the identities in the top left and bottom right corners, 
! most of this solve can actually be done with DGEMM's

call T_DEC(LACPY)('A',klu,klu,red,2*klu,red_aux,klu)
! Set up modified right hand side (to take advantage of the fact that we already know the top
! part of the vector(s) that we are solving the reduced system for)
call T_DEC(GEMM)(trans,'N',klu,nrhs,klu,-ONE_PREC,red(1,klu+1),2*klu,xt,ldxt,ONE_PREC,xb,ldxb)

call T_DEC(GETRF)(klu,klu,red_aux,klu,ipiv,info)
call T_DEC(GETRS)(trans,klu,nrhs,red_aux,klu,ipiv,xb,ldxb,info)
! Now for U we already know the bottom part of the solution, so we just need to do a multiply and 
! subtract for the top
call T_DEC(GEMM)(trans,'N',klu,nrhs,klu,-ONE_PREC,red(klu+1,1),2*klu,xb,ldxb,ONE_PREC,xt,ldxt)
call T_ALLOC(wdeallocate_2)(red_aux)
end subroutine T_DEC(spike_multi_transpose)




!The simple two partition case.
subroutine T_DEC(spike_GBTRFk2)(pivot,n,kl,ku,A,lda,info,ipiv,keys)
!=====================================================================
!  Braegan Spring  - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer :: kl,ku,klu,lda,n,info,n1,n2
MAIN_TYPE, dimension(lda,*) :: A
!MAIN_TYPE, dimension(max(kl,ku),2*max(kl,ku)) :: vw
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel,omp_get_ancestor_thread_num
!MAIN_TYPE, dimension(kl,kl) :: w
!MAIN_TYPE, dimension(ku,ku) :: v
MAIN_TYPE, pointer, dimension(:,:) :: w
MAIN_TYPE, pointer, dimension(:,:) :: v
REAL_TYPE :: norma,nzero
integer, dimension(4) :: keys
integer :: infoloc
logical :: pivot
integer, dimension(*) :: ipiv
!double precision,external :: OMP_GET_WTIME
double precision :: start_time_fac, end_time_fac
!character(len=100) :: format
!write(format,'(A1,i2,A6)') '(', 1, 'F6.2)'

klu=max(kl,ku)
norma = ZERO_REAL_PREC
nzero = ZERO_REAL_PREC


if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1
endif
! Conceptually the thing we're pulling out of the top is a lower triangular
! and the thing we're pulling out of the bottom is upper triangular.
! But they look opposite because banded is weird.

! i=omp_get_thread_num()+1
! Top partition
! Solving A[n/2-ku:n/2,n/2-(kl+ku):n/2]*V[n/2-ku:n/2]=B=A[n/2-ku:n/2,n/2:n/2+ku]
! But in banded.

!lu factorize

! Here is the logic that will make us happy even if this is only called with one thread:
! If there is no block, block other threads, and claim the key for the first section. 
! Then, unblock other threads and go to factorize the first partition. 
! If the key has been claimed for the first partition by some other thread, jump over that code. 
! Place a block after the end of the first partition factorization section. 
! Claim the key for the second partition factorization if it is unclaimed. 
! And so on. 

! Using this scheme we will be able to avoid talking to all other threads/using critical sections if it is unnecessary...

! Test the key. If the key is not claimed, enter a critical region. Test again (because the key may have been claimed between the time you checked and the time you entered the critical region)

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_1) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_1)
endif

if(keys(2) == omp_get_thread_num()) then
  if(pivot) then
    START_TIME_FAC = OMP_GET_WTIME()
    call T_DEC(GBTRF)(n1,n1,kl,ku,A(1+klu-kl,1),lda,ipiv(1),info)
    END_TIME_FAC = OMP_GET_WTIME()
!    write (*,format), END_TIME_FAC-START_TIME_FAC
  else
    call T_ALLOC(wallocate_2)(v,ku,ku,infoloc)
    v=ZERO_PREC
    call T_DEC(GBALU)(n1,kl,ku,A(1,1),lda,nzero,norma,info)
    call T_DEC(LACPY)('L',ku,ku,A(1,n1+1),lda-1,v(1,1),ku) 
    call T_DEC(TBSM)('L','N','U', ku, ku, kl, A(ku+1,n1+1-ku), lda, v(1,1), ku)
    call T_DEC(LACPY)('L',ku,ku,v(1,1),ku,A(1,n1+1),lda-1)
    call T_ALLOC(wdeallocate_2)(v)
  endif
endif

if(keys(3) == -1) then
!$omp critical (keyclaim_region_2) 
  !$omp flush
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    endif
!$omp end critical (keyclaim_region_2)
endif

if(keys(3) == omp_get_thread_num()) then
  if(pivot) then
  !  call T_ALLOC(wallocate_2)(w,2*klu,kl,infoloc)

    call T_DEC(GBTRFUL)(n2,n2,kl,ku,A(1,n1+1),lda,ipiv(n1+1),info)

  else
    call T_ALLOC(wallocate_2)(w,kl,kl,infoloc)
    w=ZERO_PREC
    call T_DEC(GBAUL)(n2,kl,ku,A(1,n1+1),lda,nzero,norma,info)
    call T_DEC(LACPY)('U',kl,kl,A(kl+ku+1, n1+1-kl),kl+ku,w(1,1),kl)
    call T_DEC(TBSM)('U','N','U',kl, kl, ku, A(1, n1+1), lda, w(1,1), kl)
    call T_DEC(LACPY)('U',kl,kl,w(1,1),kl,A(kl+ku+1, n1+1-kl),kl+ku)
    call T_ALLOC(wdeallocate_2)(w)
  endif
endif

end subroutine T_DEC(spike_GBTRFk2)








subroutine T_DEC(spike_GBTRSk2)(pivot,trans,n,nrhs,kl,ku,A,lda,B,ldb,ipiv,keys,red,grj,zeroskip)
!=====================================================================
!  Braegan Spring  2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer                                    :: n,nrhs,kl,ku,lda,ldb,info
MAIN_TYPE, dimension(lda,*)                :: A
MAIN_TYPE, dimension(ldb,*)                :: B
integer, dimension(*)                      :: ipiv 
logical                                    :: zeroskip,pivot
character                                  :: trans
integer, dimension(4)                      :: keys
MAIN_TYPE, dimension(2*max(kl,ku),*)       :: red
MAIN_TYPE, dimension(2*max(kl,ku),*)       :: grj 
!MAIN_TYPE, target, dimension(2*max(kl,ku),2*max(kl,ku)+nrhs) :: work
!MAIN_TYPE, dimension(:,:), pointer :: red
!MAIN_TYPE, dimension(:,:), pointer :: grj
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel

if(trans == 'N' .or. trans=='n') then
!  red => work(1:2*max(kl,ku),1:2*max(kl,ku))
!  grj => work(1:2*max(kl,ku),2*max(kl,ku)+1:2*max(kl,ku)+nrhs)
  call T_DEC(spike_simple_solve)(pivot,n,nrhs,kl,ku,A,lda,B,ldb,info,ipiv,keys,red,grj,zeroskip)
else
!if(omp_get_thread_num() == 1) then
!endif
!  red => work(1:2*max(kl,ku),1:2*max(kl,ku))
!  grj => work(1:2*max(kl,ku),2*max(kl,ku)+1:2*max(kl,ku)+nrhs)
  call T_DEC(spike_simple_solve_transpose)(n,trans,nrhs,kl,ku,A,lda,B,ldb,info,keys,red,grj)
endif

end subroutine T_DEC(spike_GBTRSk2)




subroutine T_DEC(spike_simple_solve)(pivot,n,nrhs,kl,ku,A,lda,f,ldf,info,ipiv,keys,red,grj,zeroskip)
!=====================================================================
!  Braegan Spring - Eric Polizzi - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer :: j,i,kl,kd,ku,klu,nrhs,local_start,local_nrhs,n,lda,info,n1,n2,ldf,iastartu
logical :: invred
logical :: zeroskip,pivot
integer, dimension(*) :: ipiv
MAIN_TYPE, dimension(lda,*):: A
!integer,external :: omp_get_thread_num,omp_get_num_threads,omp_get_max_threads,omp_in_parallel
MAIN_TYPE, dimension(ldf,*):: f
MAIN_TYPE, dimension(2*max(kl,ku),*) :: red
MAIN_TYPE, dimension(2*max(kl,ku),*) :: grj
!integer, dimension(2*max(kl,ku)) :: ipiv
integer, dimension(2*max(kl,ku)) :: ipiv_red
integer , dimension(4) :: keys
integer :: infoloc
MAIN_TYPE, pointer, dimension(:,:) :: w
MAIN_TYPE, pointer, dimension(:,:) :: v
MAIN_TYPE, pointer, dimension(:,:) :: aux
integer, pointer, dimension(:) :: ipiv_aux

invred=.false.
klu=max(kl,ku)

!Note : In general, keys(2) belongs to the one that works one the top part of the f vector, and keys(3) belongs to the one that works on the bottom part of the f vector

!call wallocate_1i(ipiv,2*max(kl,ku),infoloc)
!allocate(ipiv(2*max(kl,ku)))

!allocate(red(2*klu,2*klu))
!allocate(grj(2*klu,nrhs))


if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1 
endif

! Automatic keyclaiming region -- if you haven't set the keys the first thread to get here will grab a key.

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(2) == omp_get_thread_num()) then   

  red(1:klu,1:2*klu)=ZERO_PREC
  do i=1,klu
    red(i,i)=ONE_PREC
  enddo

  if(zeroskip) then
    local_start = ku+1
    local_nrhs  = kl
  else
    local_start = 1
    local_nrhs  = nrhs
  endif
 

  if(pivot) then
    iastartu = klu-kl+1
    kd = kl+ku

    call T_ALLOC(wallocate_2)(v,2*klu,ku,infoloc)
    call wallocate_1i(ipiv_aux,2*klu,infoloc)
    call T_DEC(GBTRSL)('N',n1,KL,KU,local_nrhs,A(klu-kl+1,1),LDA,ipiv(1),f(1,local_start),ldf,INFO)
    v=ZERO_PREC

    call T_DEC(LACPY)( 'L',ku , ku, A(1,1) , lda-1,v(1+2*klu-ku,1),2*klu)

    ipiv_aux(1:2*klu) = ipiv(n1+1-2*klu : n1) - (n1 - 2*klu)
    call T_DEC(GBTRSL)('N',2*KLU,KL,KU,KU,A(klu-kl+1,n1+1-2*klu),LDA,ipiv_aux,v(1,1),2*klu,INFO)
    call T_DEC(TBSM)('U','N','N', klu, ku, kl+ku, A(klu-kl+1, n1+1-klu), lda, v(klu+1,1), 2*klu)

    call T_DEC(LACPY)('A', klu, ku, v(klu+1,1), 2*klu, red(1,1+klu), 2*klu)

    call T_ALLOC(wdeallocate_2)(v)
  else
    !! Solving Ajgj=fj
    ! Unpacking the semi-solved v from A and doing the U sweep to form the tip of v.
    iastartu = 1
    kd = ku

    ! grj(1:klu,1:nrhs) = 0.0d0
    ! Now begin the real work of the solve stage
    call T_DEC(TBSM)('L','N','U', n1, local_nrhs,kl,A(ku+1,1),lda,f(1,local_start),ldf)
 
    call T_DEC(LACPY)('L',ku,ku,A(1,n1+1),kl+ku,red(klu-ku+1,1+klu),2*klu)
    call T_DEC(TBSM)('U','N','N', klu, ku, ku, A(1,n1+1-klu), lda, red(1,1+klu), 2*klu)
  endif

  call T_DEC(LACPY)('F', klu, nrhs,f(n1+1-klu,1),ldf, grj(1,1),2*klu)
  call T_DEC(TBSM)('U','N','N',klu,nrhs,kd,A(iastartu,n1+1-klu),lda,grj(1,1),2*klu)

  !$omp atomic
  keys(1)= keys(1)+ 1
endif


! Automatic keyclaiming region -- if you haven't set the keys the first thread to get here will grab a key.
if(keys(3) == -1) then
  !$omp critical (keyclaim_region_4) 
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_4)
endif

if(keys(3) == omp_get_thread_num()) then
!! solving Anbpart gnbpart =f nbpart
! Unpack the semi-solved w from A and do the L sweep.

  red(klu+1:2*klu,1:2*klu) = ZERO_PREC
  do i=1,klu
    red(klu+i,klu+i)=ONE_PREC
  enddo

  if(zeroskip) then
    local_nrhs  = ku
  else
    local_nrhs  = nrhs
  endif

  if(pivot) then
    if(zeroskip) then 
      local_start=kl+1
    else
      local_start=1
    endif
    iastartu = klu-ku+1
    kd = kl+ku

    call T_ALLOC(wallocate_2)(w,2*klu,kl,infoloc)
    call wallocate_1i(ipiv_aux,2*klu,infoloc)

    call T_DEC(GBTRSL)('N',n2,KU,KL,local_nrhs,A(klu-ku+1,n1+1),LDA,ipiv(n1+1),f(n1+1,local_start),ldf,INFO)

!    do j=1,nrhs
!      call T_DEC(COPY)(klu,f(n+1-klu,j),1,grj(klu+1,nrhs+1-j),1)
!    enddo
    call T_DEC(LACPY)('F', klu, nrhs, f(n+1-klu,1),ldf, grj(klu+1,1), 2*klu )

    call T_DEC(TBSM)('U','N','N', klu, nrhs, kd,A(klu-ku+1, n+1-klu), lda, grj(klu+1,1),2*klu)

    do j=1,nrhs/2
      call T_DEC(SWAP)(klu,grj(klu+1,j),1,grj(klu+1,nrhs+1-j),-1)
    enddo
    if(mod(nrhs,2) .eq. 1) then
      call T_DEC(SWAP)(klu/2,grj(klu+1,nrhs/2+1),1,grj(klu+1+(klu+1)/2,nrhs/2+1),-1)
    endif

    w=ZERO_PREC
    call T_DEC(LACPY)('L', kl, kl, A(klu+1,n1+1), lda-1, w(1+2*klu-kl,1), 2*klu)

    ipiv_aux(1:2*klu) = ipiv(n+1-2*klu : n) - (n2 - 2*klu)

    call T_DEC(GBTRSL)('N',2*KLU,KU,KL,KL,A(klu-ku+1,n+1-2*klu),LDA,ipiv_aux,w(1,1),2*klu,INFO)
    call T_DEC(TBSM)('U','N','N', klu, kl, kd, A(klu-ku+1, n+1-klu), lda, w(klu+1,1), 2*klu)

    !Need to swap the order of w somewhere... might as well do it in the w array...

    do j=1,kl/2
      call T_DEC(SWAP)(klu,w(klu+1,j),1,w(klu+1,kl+1-j),-1)
    enddo
    if(mod(kl,2) .eq. 1) then
      call T_DEC(SWAP)(klu/2,w(klu+1,kl/2+1),1,w(klu+1+(klu+1)/2,kl/2+1),-1)
    endif

    call T_DEC(LACPY)('A', klu, kl, w(klu+1,1), 2*klu, red(klu+1,1+klu-kl), 2*klu)
    
    call T_ALLOC(wdeallocate_2)(w)
  else

    call T_DEC(TBSM)('U','N','U', n2, local_nrhs, ku,A(1,n1+1),lda,f(n1+1,1),ldf)

    call T_DEC(LACPY)('F', klu, nrhs, f(n1+1,1),ldf, grj(klu+1,1), 2*klu )
    call T_DEC(TBSM)('L','N','N',klu,nrhs,kl,A(ku+1,n1+1),lda,grj(klu+1,1),2*klu)

    call T_DEC(LACPY)('U',kl,kl,A(kl+ku+1,n1+1-kl),lda-1,red(klu+1,klu-kl+1),2*klu)
    call T_DEC(TBSM)('L','N','N',klu, kl, kl, A(ku+1,n1+1), lda, red(klu+1,klu-kl+1),2*klu)
  endif


  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

!do j=1,2*klu
!enddo

do while(keys(1)< 2)
!$omp flush
enddo

if(keys(2) == omp_get_thread_num()) then
  !$omp flush
  call T_DEC(GESV)( 2*klu, nrhs, red, 2*klu, ipiv_red, grj, 2*klu, info) 
  !$omp atomic
  keys(1)= keys(1)+ 1
endif


do while(keys(1)< 3)
!$omp flush
enddo

!<<<<<<<<<<<<<<<<<<<<<<<< Retrieval  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>!

if(keys(2) == omp_get_thread_num()) then
  if(pivot) then

    call T_ALLOC(wallocate_2)(aux,2*klu,nrhs,info)
    aux = ZERO_PREC
    call T_DEC(TRMM)('L','L','N','N',ku,nrhs,-ONE_PREC,A(1,1),lda-1,grj(klu+1,1),2*klu)
    call T_DEC(LACPY)('A', ku, nrhs, grj(klu+1,1), 2*klu, aux(2*klu-ku+1,1), 2*klu)
    call T_DEC(GBTRSL)('N',2*KLU,KL,KU,NRHS,A(klu-kl+1,n1+1-2*klu),LDA,ipiv_aux,aux(1,1),2*klu,INFO)

    do i=1,nrhs
      call T_DEC(AXPY)(2*klu,ONE_PREC,aux(1,i),1,f(n/2+1-2*klu,i),1) 
    enddo

    call T_DEC(TBSM)('U','N','N', n1, nrhs, kl+ku, A(klu-kl+1, 1), lda, f(1,1), ldf)

    call wdeallocate_1i(ipiv_aux)
    call T_ALLOC(wdeallocate_2)(aux)
  else
    call T_DEC(TRMM)('L','L','N','N',ku,nrhs,-ONE_PREC,A(1 ,n1+1),lda-1,grj(klu+1,1),2*klu)

    do i=1,nrhs
      call T_DEC(AXPY)(ku,ONE_PREC,grj(klu+1,i),1,f(n/2+1-ku,i),1) 
    enddo
    call T_DEC(TBSM)('U','N','N', n1, nrhs, ku,A(1,1),lda,f(1,1),ldf)
  endif
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

if(keys(3) == omp_get_thread_num()) then


  if(pivot) then
    
    call T_ALLOC(wallocate_2)(aux,2*klu,nrhs,info)
    aux = ZERO_PREC

    do j=1,nrhs/2
      call T_DEC(SWAP)(kl,grj(klu-kl+1,j),1,grj(klu-kl+1,nrhs+1-j),-1)
    enddo
    if(mod(nrhs,2) .eq. 1) then
      call T_DEC(SWAP)(kl/2,grj(klu-kl+1,nrhs/2+1),1,grj((klu-kl+1)+(kl+1)/2,nrhs/2+1),-1)
    endif

    call T_DEC(TRMM)('L','L','N','N',kl,nrhs,-ONE_PREC,A(klu+1,n1+1),lda-1,grj(klu+1-kl,1),2*klu)
    call T_DEC(LACPY)('A', kl, nrhs, grj(klu+1-kl,1), 2*klu, aux(2*klu-kl+1,1), 2*klu)

    call T_DEC(GBTRSL)('N',2*KLU,KU,KL,NRHS,A(klu-ku+1,n+1-2*klu),LDA,ipiv_aux,aux(1,1),2*klu,INFO)

    do i=1,nrhs
      call T_DEC(AXPY)(2*klu,ONE_PREC,aux(1,i),1,f(n+1-2*klu,i),1) 
    enddo

    call T_DEC(TBSM)('U','N','N', n2, nrhs, kl+ku, A(klu-ku+1, n1+1), lda, f(n1+1,1), ldf)

    call wdeallocate_1i(ipiv_aux)
    call T_ALLOC(wdeallocate_2)(aux)
  else

    call T_DEC(TRMM)('L','U','N','N',kl,nrhs,-ONE_PREC,A(kl+ku+1,n1+1-kl),lda-1,grj(klu-kl+1,1),2*klu)

    do i=1,nrhs
      call T_DEC(AXPY)(kl,ONE_PREC,grj(klu-kl+1,i),1,f(n1+1,i),1) 
    enddo

    call T_DEC(TBSM)('L','N','N', n2, nrhs, kl,A(ku+1,n1+1),lda,f(n1+1,1),ldf)

  endif

  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

do while(keys(1)< 5)
!$omp flush
enddo
!deallocate(grj,red)

end subroutine T_DEC(spike_simple_solve)


subroutine T_DEC(spike_simple_solve_transpose) &
(n,trans,nrhs,kl,ku,A,lda,f,ldf,info,keys,red,grj)
!=====================================================================
!  Braegan Spring  - 2018
!=====================================================================
use omp_lib
implicit none
include 'f90_noruntime_interface.fi'
integer :: j,kl,ku,klu,nrhs,n,lda,info,n1,n2,ldf
MAIN_TYPE, dimension(lda,*):: A
MAIN_TYPE, dimension(ldf,*):: f
!real(kind=(kind(ONE_PREC))), dimension(2*max(kl,ku),nrhs) :: aux
MAIN_TYPE, pointer, dimension(:,:) :: aux
MAIN_TYPE, dimension(2*max(kl,ku),2*max(kl,ku)) :: red
MAIN_TYPE, dimension(2*max(kl,ku),nrhs) :: grj
integer, dimension(2*max(kl,ku)) :: ipiv_red
integer :: infoloc
!character(len=30) :: format
integer, dimension(4) :: keys
!integer, external :: omp_get_thread_num
character :: trans


call T_ALLOC(wallocate_2)(aux,2*max(kl,ku),nrhs,infoloc)
klu=max(kl,ku)
!allocate(red(2*klu,2*klu))
!allocate(grj(2*klu,nrhs))
!allocate(aux(2*klu,nrhs))

if(mod(n,2)==0) then
  n1 = n/2
  n2 = n/2
else
  n1 = n/2
  n2 = n/2+1 
endif

if(keys(3) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(3) == -1) then
    keys(3) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(3) == omp_get_thread_num()) then
  call T_DEC(LACPY)('F', klu, nrhs, f(n1+1,1), ldf, grj(klu+1,1), 2*klu)
  red(1:2*klu,1:klu) = ZERO_PREC
!    red(1:klu,1:2*klu) = 0.0d0
  do j=1,klu
    red(j,j)=ONE_PREC
  enddo
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

if(keys(2) == -1) then
  !$omp critical (keyclaim_region_3) 
  if(keys(2) == -1) then
    keys(2) = omp_get_thread_num()
    !$omp flush
  endif
  !$omp end critical (keyclaim_region_3)
endif

if(keys(2) == omp_get_thread_num()) then
  call T_DEC(LACPY)('F', klu, nrhs, f(n1+1-klu,1), ldf, grj(1,1), 2*klu)
  red(1:2*klu,1+klu:2*klu) = ZERO_PREC
!    red(klu+1:2*klu,1:2*klu) = 0.0d0
  do j=1,klu
    red(klu+j,klu+j)=ONE_PREC
  enddo
  !$omp atomic
  keys(1)= keys(1)+ 1
endif 

do while(keys(1)< 2)
!$omp flush
enddo

!Set up f-vector - save tip and augment it with zeroes
if(keys(2) == omp_get_thread_num()) then

  call T_DEC(LACPY)('L',ku,ku,A(1,n1+1),lda-1,red(klu-ku+1,1+klu),2*klu)
  call T_DEC(TBSM)('U','N','N', klu, ku, ku, A(1,n1+1-klu), lda, red(1,1+klu), 2*klu)

  f(n1-klu+1:n1,1:nrhs) = ZERO_PREC
! Perform the up sweep on f, then multiply it by the already L solved v (the L solve was done in the factorization stage)
  call T_DEC(TBSM)('U',trans,'N',n1, nrhs, ku, A(1,1), lda, f(1,1), ldf)
  call T_DEC(LACPY)('F', klu, nrhs, f(n1+1-klu,1), ldf,aux(1,1), 2*klu)

! This multiples by the v~ 
!    call T_DEC(TRMM)('L','L',trans,'N',ku, nrhs, ONE_PREC, A(1,n1+1), kl+ku, f(n1-ku+1,1), n) ! "upper" band 
  call T_DEC(TRMM)('L','L',trans,'N',ku, nrhs, ONE_PREC, A(1,n1+1), lda-1, aux(klu-ku+1,1), 2*klu) ! "upper" band 
! Later we'll assume that vw contains the tips of the actual v and w spikes (instead of the semi-solved version from the factorizaiton stage), so do their sweep here.
  do j=1,nrhs
!      call T_DEC(AXPY)(ku,-ONE_PREC,f(n1-ku+1,j),1, grj(klu+1,j),1)
    call T_DEC(AXPY)(ku,-ONE_PREC,aux(klu-ku+1,j),1, grj(klu+1,j),1)
  enddo
!  call T_DEC(LACPY)('F', klu, nrhs, aux(1,1),2*klu,f(n1+1-klu,1), 2*klu)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif  

if(keys(3) == omp_get_thread_num()) then
  call T_DEC(LACPY)('U',kl,kl,A(kl+ku+1,n1+1-kl),lda-1,red(klu+1,klu-kl+1),2*klu)
  call T_DEC(TBSM)('L','N','N',klu, kl, kl, A(ku+1,n1+1), lda, red(klu+1,klu-kl+1),2*klu)

  f(n1+1:n1+klu,1:nrhs) = ZERO_PREC
  call T_DEC(TBSM)('L',trans,'N',n2, nrhs, kl, A(ku+1,n1+1), lda, f(n1+1,1), ldf)
  call T_DEC(LACPY)('F', klu, nrhs, f(n1+1,1), ldf, aux(klu+1,1), 2*klu)
!    call T_DEC(TRMM)('L','U',trans,'N',kl, nrhs,  1.0d0, A(kl+ku+1,n1+1-kl), kl+ku, f(n1+1,1),n) ! "lower" band 
  call T_DEC(TRMM)('L','U',trans,'N',kl, nrhs, ONE_PREC, A(kl+ku+1,n1+1-kl), lda-1, aux(klu+1,1),2*klu) ! "lower" band 
  do j=1,nrhs
!      call T_DEC(AXPY)(kl,-1.0d0,f(n1+1,j),1, grj(klu-kl+1,j),1)
  call T_DEC(AXPY)(kl,-ONE_PREC,aux(klu+1,j),1, grj(klu-kl+1,j),1)
  enddo
!    call T_DEC(LACPY)('F', klu, nrhs, aux(klu+1,1), 2*klu, f(n1+1,1), n)
  !$omp atomic
  keys(1)= keys(1)+ 1
endif

do while(keys(1)< 4)
!$omp flush
enddo

!Communication -- only effects red, which is shared
if(keys(2) == omp_get_thread_num()) then
do j=1,2*klu 
enddo

  call T_DEC(GETRF)(2*klu,2*klu,red,2*klu,ipiv_red,info)
  call T_DEC(GETRS)(trans,2*klu,nrhs,red,2*klu,ipiv_red,grj,2*klu,info)

do j=1,2*klu 
enddo

  !$omp atomic
  keys(1)= keys(1)+ 1
endif

do while(keys(1)< 5)
!$omp flush
enddo

if(keys(2) == omp_get_thread_num()) then
  call T_DEC(TBSM)('U',trans,'N',klu, nrhs, ku, A(1,n1-klu+1), lda, grj(1,1), 2*klu)
  do j=1,nrhs
    call T_DEC(AXPY)(klu,ONE_PREC,grj(1,j),1, f(n1-klu+1,j),1)
  enddo
  call T_DEC(TBSM)('L',trans,'U',n1, nrhs,kl,A(1+ku,1),lda,f(1,1),ldf)
endif


if(keys(3) == omp_get_thread_num()) then
  call T_DEC(TBSM)('L',trans,'N',klu,  nrhs, kl, A(ku+1,n1+1),   lda, grj(klu+1,1),2*klu)
  do j=1,nrhs
    call T_DEC(AXPY)(klu,ONE_PREC,grj(klu+1,j),1, f(n1+1,j),1)
  enddo
  call T_DEC(TBSM)('U',trans,'U', n2, nrhs, ku, A(1,n1+1), lda, f(n1+1,1), ldf)
endif
 
!deallocate(grj,red)
call T_ALLOC(wdeallocate_2)(aux)

end subroutine T_DEC(spike_simple_solve_transpose)



subroutine T_DEC(spike_vector_norm)(spm,norm,x,n,nrhs,dnorm)
! Purpose
! -------
! Computes the norm of a given vector or set of vectors. 
! options are norm one, norm 2 and norm infinity.
! Retains the SPIKE partitioning scheme. 
! -------
!
! Arguments  
! ---------
! spm   (in) : The spike paramater array. Described in spikeinit, in the file spike_smp_utilities.f90. Contains information about the size of each partition.
!
! norm  (in) : They type of norm computed. 0 for infinorm, 1 for norm one, and 2 for infinorm.
! 
! x     (in) : The vector or set of vectors. Dimension (N by nrhs)
!
! n     (in) : The number of rows in x. This is also the leading dimensions of x. 
!
! nrhs  (in) : The number of vectors in x.
!
! dnorm (out): Array to hold the norms of the x vectors. Dimension nrhs. The norm for the i'th column of x is stores in dnorm(i)
!
!---------
!=====================================================================
!  Braegan Spring - Eric Polizzi 2018
!=====================================================================

implicit none
include 'f90_noruntime_interface.fi'
MAIN_TYPE, dimension(n,*)          :: x
integer, dimension(64)             :: spm
integer                            :: norm
integer                            :: n,nrhs
REAL_TYPE, dimension(nrhs)         :: dnorm
REAL_TYPE, dimension(:,:), pointer :: dnormTemp
integer                            :: imax
integer                            :: p,i,infoloc
REAL_TYPE, external                :: T_REAL_DEC(asum),T_REAL_DEC(nrm2)
integer, external                  :: T_INORM_DEC(amax)
#ifdef SPIKECOMPLEX 
REAL_TYPE, external                :: T_NORM_DEC(asum),T_NORM_DEC(nrm2)
integer, external                  :: T_INORM_REAL_DEC(amax)
#endif
integer, dimension(spm(20)+1)      :: Ajmin

call spikerl_calc_size_partitions(spm,Ajmin,n)
call T_REAL_ALLOC(wallocate_2)(dnormTemp,spm(20),nrhs,infoloc)

!One norm
if(norm==1) then
! Take the norm of each partition, then take the norm of those norms to get the overall norm. 
  !$OMP PARALLEL DO default(shared) private(i)  
  do p=1,spm(20)  
    do i=1,nrhs 
      dnormTemp(p,i) = T_NORM_DEC(asum)(Ajmin(p+1)-Ajmin(p), x(Ajmin(p),i),1)
    enddo
  enddo

  do i=1,nrhs 
    dnorm(i) = T_REAL_DEC(asum)(spm(20),dnormTemp(1,i),1)
  enddo
endif


!Infinorm
if(norm==0) then
  !$OMP PARALLEL DO default(shared) private(i,imax)  
  do p=1,spm(20)  
    do i=1,nrhs 
      imax = T_INORM_DEC(amax)(Ajmin(p+1)-Ajmin(p),x(Ajmin(p),i),1)
      dnormTemp(p,i) = T_NORM_DEC(asum)(1, x(Ajmin(p)+imax,i), 1) ! Use Xasum to get absolute value from LAPACK rather than fortran.
    enddo
  enddo
  do i=1,nrhs
    imax = T_INORM_REAL_DEC(amax)(spm(20),dnormTemp(1,i),1)
    dnorm(i) = dnormTemp(imax,i)
  enddo
endif


!2 norm (or euclidean)
if(norm==2) then
  !$OMP PARALLEL DO default(shared) private(i)  
  do p=1,spm(20)  
    do i=1,nrhs 
      dnormTemp(p,i) = T_NORM_DEC(nrm2)(Ajmin(p+1)-Ajmin(p) ,x(Ajmin(p),i),1)
    enddo
  enddo

  do i=1,nrhs 
    dnorm(i) = T_REAL_DEC(nrm2)(spm(20),dnormTemp(1,i),1)
  enddo
endif

call T_REAL_ALLOC(wdeallocate_2)(dnormTemp)

end subroutine T_DEC(spike_vector_norm)

