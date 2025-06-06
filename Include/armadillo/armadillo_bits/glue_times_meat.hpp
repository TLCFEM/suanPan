// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2008-2016 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2008-2016 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


//! \addtogroup glue_times
//! @{



template<bool do_inv_detect>
template<typename T1, typename T2>
inline
void
glue_times_redirect2_helper<do_inv_detect>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> U1(X.A);
  const partial_unwrap<T2> U2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = U1.M;
  const typename partial_unwrap<T2>::stored_type& B = U2.M;
  
  constexpr bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const     eT       alpha = use_alpha ? (U1.get_val() * U2.get_val()) : eT(0);
  
  if( (is_cx<eT>::no) && (resolves_to_rowvector<T1>::value && resolves_to_colvector<T2>::value) )
    {
    arma_debug_print("glue_times: dot product optimisation");
    
    arma_conform_assert_mul_size(A, B, U1.do_trans, U2.do_trans, "matrix multiplication");
    
    const eT val = op_dot::direct_dot(A.n_elem, A.memptr(), B.memptr());
    
    out.set_size(1,1);
    
    out[0] = (use_alpha) ? (val * alpha) : (val);
    
    return;
    }
  
  const bool alias = U1.is_alias(out) || U2.is_alias(out);
  
  if(alias == false)
    {
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      use_alpha
      >
      (out, A, B, alpha);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      use_alpha
      >
      (tmp, A, B, alpha);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_redirect2_helper<true>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(arma_config::optimise_invexpr && (strip_inv<T1>::do_inv_gen || strip_inv<T1>::do_inv_spd))
    {
    // replace inv(A)*B with solve(A,B)
    
    arma_debug_print("glue_times_redirect<2>::apply(): detected inv(A)*B");
    
    const strip_inv<T1> A_strip(X.A);
    
    typedef typename strip_inv<T1>::stored_type T1_stripped;
    
    if( (is_cx<eT>::no) && (strip_inv<T1>::do_inv_gen) && (is_Mat<T1_stripped>::value) && (is_Mat<T2>::value) )
      {
      const unwrap<T1_stripped> UA(A_strip.M);
      const unwrap<T2         > UB(X.B);
      
      const typename unwrap<T1_stripped>::stored_type& A = UA.M;
      const typename unwrap<T2         >::stored_type& B = UB.M;
      
      const uword N = A.n_rows;
      
      if( (N > 0) && (N <= uword(3)) && (N == A.n_cols) && (N == B.n_rows) && (void_ptr(&out) != void_ptr(&B)) )
        {
        arma_debug_print("glue_times_redirect<2>::apply(): inv tiny matrix optimisation");
        
        Mat<eT> AA(N, N, arma_nozeros_indicator());
        
        arrayops::copy(AA.memptr(), A.memptr(), AA.n_elem);
        
        bool inv_status = false;
        
        if(N == 1)  { const eT a = AA[0]; AA[0] = eT(1) / a; inv_status = (a != eT(0)); }
        if(N == 2)  { inv_status = op_inv_gen_full::apply_tiny_2x2(AA); }
        if(N == 3)  { inv_status = op_inv_gen_full::apply_tiny_3x3(AA); }
        
        if(inv_status)  { glue_times::apply<eT,false,false,false>(out, AA, B, eT(0)); return; }
        
        arma_debug_print("glue_times_redirect<2>::apply(): inv tiny matrix optimisation failed");
        
        // fallthrough if optimisation failed
        }
      }
    
    Mat<eT> A = A_strip.M;
    
    arma_conform_check( (A.is_square() == false), "inv(): given matrix must be square sized" );
    
    if( (strip_inv<T1>::do_inv_spd) && (arma_config::check_conform) && (auxlib::rudimentary_sym_check(A) == false) )
      {
      if(is_cx<eT>::no )  { arma_warn(1, "inv_sympd(): given matrix is not symmetric"); }
      if(is_cx<eT>::yes)  { arma_warn(1, "inv_sympd(): given matrix is not hermitian"); }
      }
    
    const unwrap_check<T2> B_tmp(X.B, out);
    const Mat<eT>& B = B_tmp.M;
    
    arma_conform_assert_mul_size(A, B, "matrix multiplication");
    
    const bool is_sym = (strip_inv<T1>::do_inv_spd) ? false : ( arma_config::optimise_sym && (auxlib::crippled_lapack(A) == false) && (is_sym_expr<T1>::eval(X.A) || sym_helper::is_approx_sym(A, uword(100))) );
    
    const bool status = (strip_inv<T1>::do_inv_spd) ? auxlib::solve_sympd_fast(out, A, B) : ( (is_sym) ? auxlib::solve_sym_fast(out, A, B) : auxlib::solve_square_fast(out, A, B) );
    
    if(status == false)
      {
      out.soft_reset();
      arma_stop_runtime_error("matrix multiplication: problem with matrix inverse; suggest to use solve() instead");
      }
    
    return;
    }
  
  if(arma_config::optimise_invexpr && strip_inv<T2>::do_inv_spd)
    {
    // replace A*inv_sympd(B) with trans( solve(trans(B),trans(A)) )
    // transpose of B is avoided as B is explicitly marked as symmetric
    
    arma_debug_print("glue_times_redirect<2>::apply(): detected A*inv_sympd(B)");
    
    const Mat<eT> At = trans(X.A);
    
    const strip_inv<T2> B_strip(X.B);
    
    Mat<eT> B = B_strip.M;
    
    arma_conform_check( (B.is_square() == false), "inv_sympd(): given matrix must be square sized" );
    
    if( (arma_config::check_conform) && (auxlib::rudimentary_sym_check(B) == false) )
      {
      if(is_cx<eT>::no )  { arma_warn(1, "inv_sympd(): given matrix is not symmetric"); }
      if(is_cx<eT>::yes)  { arma_warn(1, "inv_sympd(): given matrix is not hermitian"); }
      }
    
    arma_conform_assert_mul_size(At.n_cols, At.n_rows, B.n_rows, B.n_cols, "matrix multiplication");
    
    const bool status = auxlib::solve_sympd_fast(out, B, At);
    
    if(status == false)
      {
      out.soft_reset();
      arma_stop_runtime_error("matrix multiplication: problem with matrix inverse; suggest to use solve() instead");
      }
    
    out = trans(out);
    
    return;
    }
  
  glue_times_redirect2_helper<false>::apply(out, X);
  }



template<bool do_inv_detect>
template<typename T1, typename T2, typename T3>
inline
void
glue_times_redirect3_helper<do_inv_detect>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue<T1,T2,glue_times>, T3, glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // we have exactly 3 objects
  // hence we can safely expand X as X.A.A, X.A.B and X.B
  
  const partial_unwrap<T1> U1(X.A.A);
  const partial_unwrap<T2> U2(X.A.B);
  const partial_unwrap<T3> U3(X.B  );
  
  const typename partial_unwrap<T1>::stored_type& A = U1.M;
  const typename partial_unwrap<T2>::stored_type& B = U2.M;
  const typename partial_unwrap<T3>::stored_type& C = U3.M;
  
  constexpr bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times || partial_unwrap<T3>::do_times;
  const     eT       alpha = use_alpha ? (U1.get_val() * U2.get_val() * U3.get_val()) : eT(0);
  
  const bool alias = U1.is_alias(out) || U2.is_alias(out) || U3.is_alias(out);
  
  if(alias == false)
    {
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      partial_unwrap<T3>::do_trans,
      use_alpha
      >
      (out, A, B, C, alpha);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      partial_unwrap<T3>::do_trans,
      use_alpha
      >
      (tmp, A, B, C, alpha);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2, typename T3>
inline
void
glue_times_redirect3_helper<true>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue<T1,T2,glue_times>, T3, glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  if(arma_config::optimise_invexpr && (strip_inv<T1>::do_inv_gen || strip_inv<T1>::do_inv_spd))
    {
    // replace inv(A)*B*C with solve(A,B*C);
    
    arma_debug_print("glue_times_redirect<3>::apply(): detected inv(A)*B*C");
    
    const strip_inv<T1> A_strip(X.A.A);
    
    Mat<eT> A = A_strip.M;
    
    arma_conform_check( (A.is_square() == false), "inv(): given matrix must be square sized" );
    
    const partial_unwrap<T2> U2(X.A.B);
    const partial_unwrap<T3> U3(X.B  );
    
    const typename partial_unwrap<T2>::stored_type& B = U2.M;
    const typename partial_unwrap<T3>::stored_type& C = U3.M;
    
    constexpr bool use_alpha = partial_unwrap<T2>::do_times || partial_unwrap<T3>::do_times;
    const     eT       alpha = use_alpha ? (U2.get_val() * U3.get_val()) : eT(0);
    
    Mat<eT> BC;
    
    glue_times::apply
      <
      eT,
      partial_unwrap<T2>::do_trans,
      partial_unwrap<T3>::do_trans,
      use_alpha
      >
      (BC, B, C, alpha);
    
    arma_conform_assert_mul_size(A, BC, "matrix multiplication");
    
    if( (strip_inv<T1>::do_inv_spd) && (arma_config::check_conform) && (auxlib::rudimentary_sym_check(A) == false)  )
      {
      if(is_cx<eT>::no )  { arma_warn(1, "inv_sympd(): given matrix is not symmetric"); }
      if(is_cx<eT>::yes)  { arma_warn(1, "inv_sympd(): given matrix is not hermitian"); }
      }
    
    const bool is_sym = (strip_inv<T1>::do_inv_spd) ? false : ( arma_config::optimise_sym && (auxlib::crippled_lapack(A) == false) && (is_sym_expr<T1>::eval(X.A.A) || sym_helper::is_approx_sym(A, uword(100))) );
    
    const bool status = (strip_inv<T1>::do_inv_spd) ? auxlib::solve_sympd_fast(out, A, BC) : ( (is_sym) ? auxlib::solve_sym_fast(out, A, BC) : auxlib::solve_square_fast(out, A, BC) );
    
    if(status == false)
      {
      out.soft_reset();
      arma_stop_runtime_error("matrix multiplication: problem with matrix inverse; suggest to use solve() instead");
      }
    
    return;
    }
  
  
  if(arma_config::optimise_invexpr && (strip_inv<T2>::do_inv_gen || strip_inv<T2>::do_inv_spd))
    {
    // replace A*inv(B)*C with A*solve(B,C)
    
    arma_debug_print("glue_times_redirect<3>::apply(): detected A*inv(B)*C");
    
    const strip_inv<T2> B_strip(X.A.B);
    
    Mat<eT> B = B_strip.M;
    
    arma_conform_check( (B.is_square() == false), "inv(): given matrix must be square sized" );
    
    const quasi_unwrap<T3> U3(X.B);
    const Mat<eT>& C =     U3.M;
    
    arma_conform_assert_mul_size(B, C, "matrix multiplication");
    
    if( (strip_inv<T2>::do_inv_spd) && (arma_config::check_conform) && (auxlib::rudimentary_sym_check(B) == false)  )
      {
      if(is_cx<eT>::no )  { arma_warn(1, "inv_sympd(): given matrix is not symmetric"); }
      if(is_cx<eT>::yes)  { arma_warn(1, "inv_sympd(): given matrix is not hermitian"); }
      }
    
    Mat<eT> solve_result;
    
    const bool is_sym = (strip_inv<T1>::do_inv_spd) ? false : ( arma_config::optimise_sym && (auxlib::crippled_lapack(B) == false) && (is_sym_expr<T2>::eval(X.A.B) || sym_helper::is_approx_sym(B, uword(100))) );
    
    const bool status = (strip_inv<T2>::do_inv_spd) ? auxlib::solve_sympd_fast(solve_result, B, C) : ( (is_sym) ? auxlib::solve_sym_fast(solve_result, B, C) : auxlib::solve_square_fast(solve_result, B, C) );
    
    if(status == false)
      {
      out.soft_reset();
      arma_stop_runtime_error("matrix multiplication: problem with matrix inverse; suggest to use solve() instead");
      return;
      }
    
    const partial_unwrap<T1> U1(X.A.A);
    
    const typename partial_unwrap<T1>::stored_type& A = U1.M;
    
    constexpr bool use_alpha = partial_unwrap<T1>::do_times;
    const     eT       alpha = use_alpha ? U1.get_val() : eT(0);
    
    if(U1.is_alias(out))
      {
      Mat<eT> tmp;
      
      glue_times::apply<eT, partial_unwrap<T1>::do_trans, false, use_alpha>(tmp, A, solve_result, alpha);
      
      out.steal_mem(tmp);
      }
    else
      {
      glue_times::apply<eT, partial_unwrap<T1>::do_trans, false, use_alpha>(out, A, solve_result, alpha);
      }
    
    return;
    }
  
  
  glue_times_redirect3_helper<false>::apply(out, X);
  }



template<uword N>
template<typename T1, typename T2>
inline
void
glue_times_redirect<N>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const partial_unwrap<T1> U1(X.A);
  const partial_unwrap<T2> U2(X.B);
  
  const typename partial_unwrap<T1>::stored_type& A = U1.M;
  const typename partial_unwrap<T2>::stored_type& B = U2.M;
  
  constexpr bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times;
  const     eT       alpha = use_alpha ? (U1.get_val() * U2.get_val()) : eT(0);
  
  const bool alias = U1.is_alias(out) || U2.is_alias(out);
  
  if(alias == false)
    {
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      use_alpha
      >
      (out, A, B, alpha);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      use_alpha
      >
      (tmp, A, B, alpha);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
glue_times_redirect<2>::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  glue_times_redirect2_helper< is_supported_blas_type<eT>::value >::apply(out, X);
  }



template<typename T1, typename T2, typename T3>
inline
void
glue_times_redirect<3>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue<T1,T2,glue_times>, T3, glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  glue_times_redirect3_helper< is_supported_blas_type<eT>::value >::apply(out, X);
  }



template<typename T1, typename T2, typename T3, typename T4>
inline
void
glue_times_redirect<4>::apply(Mat<typename T1::elem_type>& out, const Glue< Glue< Glue<T1,T2,glue_times>, T3, glue_times>, T4, glue_times>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  // there is exactly 4 objects
  // hence we can safely expand X as X.A.A.A, X.A.A.B, X.A.B and X.B
  
  const partial_unwrap<T1> U1(X.A.A.A);
  const partial_unwrap<T2> U2(X.A.A.B);
  const partial_unwrap<T3> U3(X.A.B  );
  const partial_unwrap<T4> U4(X.B    );
  
  const typename partial_unwrap<T1>::stored_type& A = U1.M;
  const typename partial_unwrap<T2>::stored_type& B = U2.M;
  const typename partial_unwrap<T3>::stored_type& C = U3.M;
  const typename partial_unwrap<T4>::stored_type& D = U4.M;
  
  constexpr bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times || partial_unwrap<T3>::do_times || partial_unwrap<T4>::do_times;
  const     eT       alpha = use_alpha ? (U1.get_val() * U2.get_val() * U3.get_val() * U4.get_val()) : eT(0);
  
  const bool alias = U1.is_alias(out) || U2.is_alias(out) || U3.is_alias(out) || U4.is_alias(out);
  
  if(alias == false)
    {
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      partial_unwrap<T3>::do_trans,
      partial_unwrap<T4>::do_trans,
      use_alpha
      >
      (out, A, B, C, D, alpha);
    }
  else
    {
    Mat<eT> tmp;
    
    glue_times::apply
      <
      eT,
      partial_unwrap<T1>::do_trans,
      partial_unwrap<T2>::do_trans,
      partial_unwrap<T3>::do_trans,
      partial_unwrap<T4>::do_trans,
      use_alpha
      >
      (tmp, A, B, C, D, alpha);
    
    out.steal_mem(tmp);
    }
  }



template<typename T1, typename T2>
inline
void
glue_times::apply(Mat<typename T1::elem_type>& out, const Glue<T1,T2,glue_times>& X)
  {
  arma_debug_sigprint();
  
  constexpr uword N_mat = 1 + depth_lhs< glue_times, Glue<T1,T2,glue_times> >::num;
  
  arma_debug_print(arma_str::format("glue_times::apply(): N_mat: %u") % N_mat);
  
  glue_times_redirect<N_mat>::apply(out, X);
  }



template<typename T1>
inline
void
glue_times::apply_inplace(Mat<typename T1::elem_type>& out, const T1& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  Mat<eT> tmp = out * X;
  
  out.steal_mem(tmp);
  }



template<typename T1, typename T2>
inline
void
glue_times::apply_inplace_plus(Mat<typename T1::elem_type>& out, const Glue<T1, T2, glue_times>& X, const sword sign)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type            eT;
  typedef typename get_pod_type<eT>::result  T;
  
  if( X.is_alias(out) || (is_outer_product<T1>::value) || (has_op_inv_any<T1>::value) || (has_op_inv_any<T2>::value) )
    {
    // handle aliasing and partial workaround for corner cases
    
    const Mat<eT> tmp(X);
    
    if(sign > sword(0))  { out += tmp; }  else  { out -= tmp; }
    
    return;
    }
  
  const partial_unwrap<T1> U1(X.A);
  const partial_unwrap<T2> U2(X.B);
  
  typedef typename partial_unwrap<T1>::stored_type TA;
  typedef typename partial_unwrap<T2>::stored_type TB;
  
  const TA& A = U1.M;
  const TB& B = U2.M;
  
  constexpr bool do_trans_A = partial_unwrap<T1>::do_trans;
  constexpr bool do_trans_B = partial_unwrap<T2>::do_trans;
  
  const bool use_alpha = partial_unwrap<T1>::do_times || partial_unwrap<T2>::do_times || (sign < sword(0));
  
  const eT       alpha = use_alpha ? ( U1.get_val() * U2.get_val() * ( (sign > sword(0)) ? eT(1) : eT(-1) ) ) : eT(0);
  
  arma_conform_assert_mul_size(A, B, do_trans_A, do_trans_B, "matrix multiplication");
  
  const uword result_n_rows = (do_trans_A == false) ? (TA::is_row ? 1 : A.n_rows) : (TA::is_col ? 1 : A.n_cols);
  const uword result_n_cols = (do_trans_B == false) ? (TB::is_col ? 1 : B.n_cols) : (TB::is_row ? 1 : B.n_rows);
  
  arma_conform_assert_same_size(out.n_rows, out.n_cols, result_n_rows, result_n_cols, ( (sign > sword(0)) ? "addition" : "subtraction" ) );
  
  if(out.n_elem == 0)  { return; }
  
  if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == false) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) )  { gemv<true,         false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<false,        false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else                                                             { gemm<false, false, false, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == true) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) )  { gemv<true,         true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<false,        true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else                                                             { gemm<false, false, true, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == false) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no)  )  { gemv<true,        false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if(  (B.n_cols == 1) || (TB::is_col)                      )  { gemv<true,        false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no)  )  { syrk<true,        false, true>::apply(out,          A,             alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::yes) )  { herk<true,        false, true>::apply(out,          A,              T(0),  T(1)); }
    else                                                              { gemm<true, false, false, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == true) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<true,        true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<true,        true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no) )  { syrk<true,        true, true>::apply(out,          A,             alpha, eT(1)); }
    else                                                             { gemm<true, false, true, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == false) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no)  )  { gemv<false,       false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no)  )  { gemv<false,       false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no)  )  { syrk<false,       false, true>::apply(out,          A,             alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::yes) )  { herk<false,       false, true>::apply(out,          A,              T(0),  T(1)); }
    else                                                              { gemm<false, true, false, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == true) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) )  { gemv<false,       true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) )  { gemv<false,       true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no) )  { syrk<false,       true, true>::apply(out,          A,             alpha, eT(1)); }
    else                                                             { gemm<false, true, true, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == false) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<false,      false, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) )  { gemv<true,       false, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else                                                             { gemm<true, true, false, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == true) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<false,      true, true>::apply(out.memptr(), B, A.memptr(), alpha, eT(1)); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) )  { gemv<true,       true, true>::apply(out.memptr(), A, B.memptr(), alpha, eT(1)); }
    else                                                             { gemm<true, true, true, true>::apply(out,          A, B,          alpha, eT(1)); }
    }
  }



template<typename eT, const bool do_trans_A, const bool do_trans_B, typename TA, typename TB>
arma_inline
uword
glue_times::mul_storage_cost(const TA& A, const TB& B)
  {
  const uword final_A_n_rows = (do_trans_A == false) ? ( TA::is_row ? 1 : A.n_rows ) : ( TA::is_col ? 1 : A.n_cols );
  const uword final_B_n_cols = (do_trans_B == false) ? ( TB::is_col ? 1 : B.n_cols ) : ( TB::is_row ? 1 : B.n_rows );
  
  return final_A_n_rows * final_B_n_cols;
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool use_alpha,
  typename   TA,
  typename   TB
  >
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const TA&      A,
  const TB&      B,
  const eT       alpha
  )
  {
  arma_debug_sigprint();
  
  //arma_conform_assert_mul_size(A, B, do_trans_A, do_trans_B, "matrix multiplication");
  arma_conform_assert_trans_mul_size<do_trans_A, do_trans_B>(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
  
  const uword final_n_rows = (do_trans_A == false) ? (TA::is_row ? 1 : A.n_rows) : (TA::is_col ? 1 : A.n_cols);
  const uword final_n_cols = (do_trans_B == false) ? (TB::is_col ? 1 : B.n_cols) : (TB::is_row ? 1 : B.n_rows);
  
  out.set_size(final_n_rows, final_n_cols);
  
  if( (A.n_elem == 0) || (B.n_elem == 0) )  { out.zeros(); return; }
  
  if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == false) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) )  { gemv<true,         false, false>::apply(out.memptr(), B, A.memptr()); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<false,        false, false>::apply(out.memptr(), A, B.memptr()); }
    else                                                             { gemm<false, false, false, false>::apply(out,          A, B         ); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == false) && (use_alpha == true) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) )  { gemv<true,         true, false>::apply(out.memptr(), B, A.memptr(), alpha); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<false,        true, false>::apply(out.memptr(), A, B.memptr(), alpha); }
    else                                                             { gemm<false, false, true, false>::apply(out,          A, B,          alpha); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == false) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no)  )  { gemv<true,        false, false>::apply(out.memptr(), B, A.memptr()); }
    else if(  (B.n_cols == 1) || (TB::is_col)                      )  { gemv<true,        false, false>::apply(out.memptr(), A, B.memptr()); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no)  )  { syrk<true,        false, false>::apply(out,          A            ); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::yes) )  { herk<true,        false, false>::apply(out,          A            ); }
    else                                                              { gemm<true, false, false, false>::apply(out,          A, B         ); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == false) && (use_alpha == true) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<true,        true, false>::apply(out.memptr(), B, A.memptr(), alpha); }
    else if(  (B.n_cols == 1) || (TB::is_col)                     )  { gemv<true,        true, false>::apply(out.memptr(), A, B.memptr(), alpha); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no) )  { syrk<true,        true, false>::apply(out,          A,             alpha); }
    else                                                             { gemm<true, false, true, false>::apply(out,          A, B,          alpha); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == false) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no)  )  { gemv<false,       false, false>::apply(out.memptr(), B, A.memptr()); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no)  )  { gemv<false,       false, false>::apply(out.memptr(), A, B.memptr()); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no)  )  { syrk<false,       false, false>::apply(out,          A            ); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::yes) )  { herk<false,       false, false>::apply(out,          A            ); }
    else                                                              { gemm<false, true, false, false>::apply(out,          A, B         ); }
    }
  else
  if( (do_trans_A == false) && (do_trans_B == true) && (use_alpha == true) )
    {
         if( ((A.n_rows == 1) || (TA::is_row)) && (is_cx<eT>::no) ) { gemv<false,       true, false>::apply(out.memptr(), B, A.memptr(), alpha); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) ) { gemv<false,       true, false>::apply(out.memptr(), A, B.memptr(), alpha); }
    else if( (void_ptr(&A) == void_ptr(&B))    && (is_cx<eT>::no) ) { syrk<false,       true, false>::apply(out,          A,             alpha); }
    else                                                            { gemm<false, true, true, false>::apply(out,          A, B,          alpha); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == false) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<false,      false, false>::apply(out.memptr(), B, A.memptr()); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) )  { gemv<true,       false, false>::apply(out.memptr(), A, B.memptr()); }
    else                                                             { gemm<true, true, false, false>::apply(out,          A, B         ); }
    }
  else
  if( (do_trans_A == true) && (do_trans_B == true) && (use_alpha == true) )
    {
         if( ((A.n_cols == 1) || (TA::is_col)) && (is_cx<eT>::no) )  { gemv<false,      true, false>::apply(out.memptr(), B, A.memptr(), alpha); }
    else if( ((B.n_rows == 1) || (TB::is_row)) && (is_cx<eT>::no) )  { gemv<true,       true, false>::apply(out.memptr(), A, B.memptr(), alpha); }
    else                                                             { gemm<true, true, true, false>::apply(out,          A, B,          alpha); }
    }
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool do_trans_C,
  const bool use_alpha,
  typename   TA,
  typename   TB,
  typename   TC
  >
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const TA&      A,
  const TB&      B,
  const TC&      C,
  const eT       alpha
  )
  {
  arma_debug_sigprint();
  
  Mat<eT> tmp;
  
  const uword storage_cost_AB = glue_times::mul_storage_cost<eT, do_trans_A, do_trans_B>(A, B);
  const uword storage_cost_BC = glue_times::mul_storage_cost<eT, do_trans_B, do_trans_C>(B, C);
  
  if(storage_cost_AB <= storage_cost_BC)
    {
    // out = (A*B)*C
    
    glue_times::apply<eT, do_trans_A, do_trans_B, use_alpha>(tmp, A,   B, alpha);
    glue_times::apply<eT, false,      do_trans_C, false    >(out, tmp, C, eT(0));
    }
  else
    {
    // out = A*(B*C)
    
    glue_times::apply<eT, do_trans_B, do_trans_C, use_alpha>(tmp, B, C,   alpha);
    glue_times::apply<eT, do_trans_A, false,      false    >(out, A, tmp, eT(0));
    }
  }



template
  <
  typename   eT,
  const bool do_trans_A,
  const bool do_trans_B,
  const bool do_trans_C,
  const bool do_trans_D,
  const bool use_alpha,
  typename   TA,
  typename   TB,
  typename   TC,
  typename   TD
  >
inline
void
glue_times::apply
  (
        Mat<eT>& out,
  const TA&      A,
  const TB&      B,
  const TC&      C,
  const TD&      D,
  const eT       alpha
  )
  {
  arma_debug_sigprint();
  
  Mat<eT> tmp;
  
  const uword storage_cost_AC = glue_times::mul_storage_cost<eT, do_trans_A, do_trans_C>(A, C);
  const uword storage_cost_BD = glue_times::mul_storage_cost<eT, do_trans_B, do_trans_D>(B, D);
  
  if(storage_cost_AC <= storage_cost_BD)
    {
    // out = (A*B*C)*D
    
    glue_times::apply<eT, do_trans_A, do_trans_B, do_trans_C, use_alpha>(tmp, A, B, C, alpha);
    
    glue_times::apply<eT, false, do_trans_D, false>(out, tmp, D, eT(0));
    }
  else
    {
    // out = A*(B*C*D)
    
    glue_times::apply<eT, do_trans_B, do_trans_C, do_trans_D, use_alpha>(tmp, B, C, D, alpha);
    
    glue_times::apply<eT, do_trans_A, false, false>(out, A, tmp, eT(0));
    }
  }



//
// glue_times_diag


template<typename T1, typename T2>
inline
void
glue_times_diag::apply(Mat<typename T1::elem_type>& actual_out, const Glue<T1, T2, glue_times_diag>& X)
  {
  arma_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const strip_diagmat<T1> S1(X.A);
  const strip_diagmat<T2> S2(X.B);
  
  typedef typename strip_diagmat<T1>::stored_type T1_stripped;
  typedef typename strip_diagmat<T2>::stored_type T2_stripped;
  
  if( (strip_diagmat<T1>::do_diagmat == true) && (strip_diagmat<T2>::do_diagmat == false) )
    {
    arma_debug_print("glue_times_diag::apply(): diagmat(A) * B");
    
    const diagmat_proxy<T1_stripped> A(S1.M);
    
    const quasi_unwrap<T2> UB(X.B);
    const Mat<eT>& B     = UB.M;
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    const uword A_length = (std::min)(A_n_rows, A_n_cols);
    
    const uword B_n_rows = B.n_rows;
    const uword B_n_cols = B.n_cols;
    
    arma_conform_assert_mul_size(A_n_rows, A_n_cols, B_n_rows, B_n_cols, "matrix multiplication");
    
    const bool is_alias = (A.is_alias(actual_out) || UB.is_alias(actual_out));
    
    if(is_alias)  { arma_debug_print("glue_times_diag::apply(): aliasing detected"); }
    
    Mat<eT>  tmp;
    Mat<eT>& out = (is_alias) ? tmp : actual_out;
    
    out.zeros(A_n_rows, B_n_cols);
    
    for(uword col=0; col < B_n_cols; ++col)
      {
            eT* out_coldata = out.colptr(col);
      const eT*   B_coldata =   B.colptr(col);
      
      for(uword i=0; i < A_length; ++i)  { out_coldata[i] = A[i] * B_coldata[i]; }
      }
    
    if(is_alias)  { actual_out.steal_mem(tmp); }
    }
  else
  if( (strip_diagmat<T1>::do_diagmat == false) && (strip_diagmat<T2>::do_diagmat == true) )
    {
    arma_debug_print("glue_times_diag::apply(): A * diagmat(B)");
    
    const quasi_unwrap<T1> UA(X.A);
    const Mat<eT>& A     = UA.M;
    
    const diagmat_proxy<T2_stripped> B(S2.M);
    
    const uword A_n_rows = A.n_rows;
    const uword A_n_cols = A.n_cols;
    
    const uword B_n_rows = B.n_rows;
    const uword B_n_cols = B.n_cols;
    const uword B_length = (std::min)(B_n_rows, B_n_cols);
    
    arma_conform_assert_mul_size(A_n_rows, A_n_cols, B_n_rows, B_n_cols, "matrix multiplication");
    
    const bool is_alias = (UA.is_alias(actual_out) || B.is_alias(actual_out));
    
    if(is_alias)  { arma_debug_print("glue_times_diag::apply(): aliasing detected"); }
    
    Mat<eT>  tmp;
    Mat<eT>& out = (is_alias) ? tmp : actual_out;
    
    out.zeros(A_n_rows, B_n_cols);
    
    for(uword col=0; col < B_length; ++col)
      {
      const eT  val = B[col];
      
            eT* out_coldata = out.colptr(col);
      const eT*   A_coldata =   A.colptr(col);
      
      for(uword i=0; i < A_n_rows; ++i)  { out_coldata[i] = A_coldata[i] * val; }
      }
    
    if(is_alias)  { actual_out.steal_mem(tmp); }
    }
  else
  if( (strip_diagmat<T1>::do_diagmat == true) && (strip_diagmat<T2>::do_diagmat == true) )
    {
    arma_debug_print("glue_times_diag::apply(): diagmat(A) * diagmat(B)");
    
    const diagmat_proxy<T1_stripped> A(S1.M);
    const diagmat_proxy<T2_stripped> B(S2.M);
    
    arma_conform_assert_mul_size(A.n_rows, A.n_cols, B.n_rows, B.n_cols, "matrix multiplication");
    
    const bool is_alias = (A.is_alias(actual_out) || B.is_alias(actual_out));
    
    if(is_alias)  { arma_debug_print("glue_times_diag::apply(): aliasing detected"); }
    
    Mat<eT>  tmp;
    Mat<eT>& out = (is_alias) ? tmp : actual_out;
    
    out.zeros(A.n_rows, B.n_cols);
    
    const uword A_length = (std::min)(A.n_rows, A.n_cols);
    const uword B_length = (std::min)(B.n_rows, B.n_cols);
    
    const uword N = (std::min)(A_length, B_length);
    
    for(uword i=0; i < N; ++i)  { out.at(i,i) = A[i] * B[i]; }
    
    if(is_alias)  { actual_out.steal_mem(tmp); }
    }
  }



//! @}
