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


//! \addtogroup Proxy
//! @{


// within each specialisation of the Proxy class:
// 
// elem_type        = type of the elements obtained from object Q
// pod_type         = underlying type of elements if elem_type is std::complex
// stored_type      = type of Q object
// ea_type          = type of object that provides access to elements via operator[i]
// aligned_ea_type  = type of object that provides access to elements via at_alt(i)
// 
// use_at           = boolean to indicate at(row,col) must be used to get elements
// use_mp           = boolean to indicate OpenMP can be used while processing elements
// has_subview      = boolean to indicate Q object has a subview
// 
// is_row           = boolean to indicate Q object can be treated a row vector
// is_col           = boolean to indicate Q object can be treated a column vector
// is_xvec          = boolean to indicate Q object is a vector with unknown orientation
// 
// Q                = object that can be unwrapped via unwrap family of classes (ie. Q must be convertible to Mat)
// 
// get_n_rows()     = return number of rows in Q
// get_n_cols()     = return number of columns in Q
// get_n_elem()     = return number of elements in Q
// 
// operator[i]      = linear element accessor; valid only if 'use_at' boolean is false
// at(row,col)      = access elements via (row,col); valid only if 'use_at' boolean is true
// at_alt(i)        = aligned linear element accessor; valid only if 'use_at' boolean is false and is_aligned() returns true
// 
// get_ea()         = return object that provides linear access to elements via operator[i]
// get_aligned_ea() = return object that provides linear access to elements via at_alt(i); valid only if is_aligned() returns true
// 
// is_alias(X)      = return true/false to indicate Q object aliases matrix X
// has_overlap(X)   = return true/false to indicate Q object has overlap with subview X
// is_aligned()     = return true/false to indicate Q object has aligned memory



template<typename T1>
struct Proxy_default
  {
  inline Proxy_default(const T1&)
    {
    arma_type_check(( is_arma_type<T1>::value == false ));
    }
  };



template<typename T1>
struct Proxy_fixed
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef T1                                       stored_type;
  typedef const elem_type*                         ea_type;
  typedef const T1&                                aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = T1::is_row;
  static constexpr bool is_col  = T1::is_col;
  static constexpr bool is_xvec = T1::is_xvec;
  
  arma_aligned const T1& Q;
  
  inline explicit Proxy_fixed(const T1& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  //// this may require T1::n_elem etc to be declared as static constexpr inline variables (C++17)
  //// see also the notes in Mat::fixed
  //// https://en.cppreference.com/w/cpp/language/static
  //// https://en.cppreference.com/w/cpp/language/inline
  // 
  // static constexpr uword get_n_rows() { return T1::n_rows; }
  // static constexpr uword get_n_cols() { return T1::n_cols; }
  // static constexpr uword get_n_elem() { return T1::n_elem; }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : T1::n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : T1::n_cols; }
  arma_inline uword get_n_elem() const { return              T1::n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&Q) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const
    {
    #if defined(ARMA_HAVE_ALIGNED_ATTRIBUTE)
      return true;
    #else
      return memory::is_aligned(Q.memptr());
    #endif
    }
  };



template<typename T1, bool condition>
struct Proxy_redirect {};

template<typename T1>
struct Proxy_redirect<T1, false> { typedef Proxy_default<T1> result; };

template<typename T1>
struct Proxy_redirect<T1, true>  { typedef Proxy_fixed<T1>   result; };



template<typename T1>
struct Proxy : public Proxy_redirect<T1, is_Mat_fixed<T1>::value>::result
  {
  inline Proxy(const T1& A)
    : Proxy_redirect<T1, is_Mat_fixed<T1>::value>::result(A)
    {
    }
  };



template<typename eT>
struct Proxy< Mat<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Mat<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const Mat<eT>& Q;
  
  inline explicit Proxy(const Mat<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&Q) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< Col<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Col<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Col<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const Col<eT>& Q;
  
  inline explicit Proxy(const Col<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&Q) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< Row<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Row<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Row<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const Row<eT>& Q;
  
  inline explicit Proxy(const Row<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  constexpr   uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword, const uword c) const { return Q[c];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&Q) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename gen_type>
struct Proxy< Gen<T1, gen_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Gen<T1, gen_type>                        stored_type;
  typedef const Gen<T1, gen_type>&                 ea_type;
  typedef const Gen<T1, gen_type>&                 aligned_ea_type;
  
  static constexpr bool use_at      = Gen<T1, gen_type>::use_at;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = Gen<T1, gen_type>::is_row;
  static constexpr bool is_col  = Gen<T1, gen_type>::is_col;
  static constexpr bool is_xvec = Gen<T1, gen_type>::is_xvec;
  
  const Gen<T1, gen_type>& Q;
  
  inline explicit Proxy(const Gen<T1, gen_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return (is_row ? 1 : Q.n_rows);                           }
  arma_inline uword get_n_cols() const { return (is_col ? 1 : Q.n_cols);                           }
  arma_inline uword get_n_elem() const { return (is_row ? 1 : Q.n_rows) * (is_col ? 1 : Q.n_cols); }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];       }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c); }
  arma_inline elem_type at_alt     (const uword i)                const { return Q[i];       }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  constexpr bool is_aligned() const { return Gen<T1, gen_type>::is_simple; }
  };



template<typename T1, typename eop_type>
struct Proxy< eOp<T1, eop_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eOp<T1, eop_type>                        stored_type;
  typedef const eOp<T1, eop_type>&                 ea_type;
  typedef const eOp<T1, eop_type>&                 aligned_ea_type;
  
  static constexpr bool use_at      = eOp<T1, eop_type>::use_at;
  static constexpr bool use_mp      = eOp<T1, eop_type>::use_mp;
  static constexpr bool has_subview = eOp<T1, eop_type>::has_subview;
  
  static constexpr bool is_row  = eOp<T1, eop_type>::is_row;
  static constexpr bool is_col  = eOp<T1, eop_type>::is_col;
  static constexpr bool is_xvec = eOp<T1, eop_type>::is_xvec;
  
  const eOp<T1, eop_type>& Q;
  
  inline explicit Proxy(const eOp<T1, eop_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.get_n_cols(); }
  arma_inline uword get_n_elem() const { return Q.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Q.P.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Q.P.has_overlap(X); }
  
  arma_inline bool is_aligned() const { return Q.P.is_aligned(); }
  };



template<typename T1, typename T2, typename eglue_type>
struct Proxy< eGlue<T1, T2, eglue_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef eGlue<T1, T2, eglue_type>                stored_type;
  typedef const eGlue<T1, T2, eglue_type>&         ea_type;
  typedef const eGlue<T1, T2, eglue_type>&         aligned_ea_type;
  
  static constexpr bool use_at      = eGlue<T1, T2, eglue_type>::use_at;
  static constexpr bool use_mp      = eGlue<T1, T2, eglue_type>::use_mp;
  static constexpr bool has_subview = eGlue<T1, T2, eglue_type>::has_subview;
  
  static constexpr bool is_row  = eGlue<T1, T2, eglue_type>::is_row;
  static constexpr bool is_col  = eGlue<T1, T2, eglue_type>::is_col;
  static constexpr bool is_xvec = eGlue<T1, T2, eglue_type>::is_xvec;
  
  const eGlue<T1, T2, eglue_type>& Q;
  
  inline explicit Proxy(const eGlue<T1, T2, eglue_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.get_n_cols(); }
  arma_inline uword get_n_elem() const { return Q.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (Q.P1.is_alias(X) || Q.P2.is_alias(X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return (Q.P1.has_overlap(X) || Q.P2.has_overlap(X)); }
  
  arma_inline bool is_aligned() const { return (Q.P1.is_aligned() && Q.P2.is_aligned()); }
  };



template<typename T1, typename op_type>
struct Proxy< Op<T1, op_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = Op<T1, op_type>::is_row;
  static constexpr bool is_col  = Op<T1, op_type>::is_col;
  static constexpr bool is_xvec = Op<T1, op_type>::is_xvec;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const Op<T1, op_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename T2, typename glue_type>
struct Proxy< Glue<T1, T2, glue_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = Glue<T1, T2, glue_type>::is_row;
  static constexpr bool is_col  = Glue<T1, T2, glue_type>::is_col;
  static constexpr bool is_xvec = Glue<T1, T2, glue_type>::is_xvec;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const Glue<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }

  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename T2>
struct Proxy< Glue<T1, T2, glue_min> >
  {
  typedef        Glue<T1, T2, glue_min>   this_Glue_type;
  typedef Proxy< Glue<T1, T2, glue_min> > this_Proxy_type;
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef this_Glue_type                           stored_type;
  typedef const this_Proxy_type&                   ea_type;
  typedef const this_Proxy_type&                   aligned_ea_type;
  
  static constexpr bool use_at      = (Proxy<T1>::use_at      || Proxy<T2>::use_at     );
  static constexpr bool use_mp      = (Proxy<T1>::use_mp      || Proxy<T2>::use_mp     );
  static constexpr bool has_subview = (Proxy<T1>::has_subview || Proxy<T2>::has_subview);
  
  static constexpr bool is_row  = this_Glue_type::is_row;
  static constexpr bool is_col  = this_Glue_type::is_col;
  static constexpr bool is_xvec = this_Glue_type::is_xvec;
  
  const this_Glue_type& Q;
  const Proxy<T1>       P1;
  const Proxy<T2>       P2;
  
  arma_lt_comparator<elem_type> comparator;
  
  inline explicit Proxy(const this_Glue_type& X)
    : Q (X  )
    , P1(X.A)
    , P2(X.B)
    {
    arma_debug_sigprint();
    
    arma_conform_assert_same_size(P1, P2, "element-wise min()");
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : P1.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : P1.get_n_cols(); }
  arma_inline uword get_n_elem() const { return P1.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                const { const elem_type A = P1[i];        const elem_type B = P2[i];        return comparator(A,B) ? A : B; }
  arma_inline elem_type at         (const uword r, const uword c) const { const elem_type A = P1.at(r,c);   const elem_type B = P2.at(r,c);   return comparator(A,B) ? A : B; }
  arma_inline elem_type at_alt     (const uword i)                const { const elem_type A = P1.at_alt(i); const elem_type B = P2.at_alt(i); return comparator(A,B) ? A : B; }
  
  arma_inline         ea_type         get_ea() const { return *this; }
  arma_inline aligned_ea_type get_aligned_ea() const { return *this; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (P1.is_alias(X) || P2.is_alias(X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return (P1.has_overlap(X) || P2.has_overlap(X)); }
  
  arma_inline bool is_aligned() const { return (P1.is_aligned() && P2.is_aligned()); }
  };



template<typename T1, typename T2>
struct Proxy< Glue<T1, T2, glue_max> >
  {
  typedef        Glue<T1, T2, glue_max>   this_Glue_type;
  typedef Proxy< Glue<T1, T2, glue_max> > this_Proxy_type;
  
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef this_Glue_type                           stored_type;
  typedef const this_Proxy_type&                   ea_type;
  typedef const this_Proxy_type&                   aligned_ea_type;
  
  static constexpr bool use_at      = (Proxy<T1>::use_at      || Proxy<T2>::use_at     );
  static constexpr bool use_mp      = (Proxy<T1>::use_mp      || Proxy<T2>::use_mp     );
  static constexpr bool has_subview = (Proxy<T1>::has_subview || Proxy<T2>::has_subview);
  
  static constexpr bool is_row  = this_Glue_type::is_row;
  static constexpr bool is_col  = this_Glue_type::is_col;
  static constexpr bool is_xvec = this_Glue_type::is_xvec;
  
  const this_Glue_type& Q;
  const Proxy<T1>       P1;
  const Proxy<T2>       P2;
  
  arma_gt_comparator<elem_type> comparator;
  
  inline explicit Proxy(const this_Glue_type& X)
    : Q (X  )
    , P1(X.A)
    , P2(X.B)
    {
    arma_debug_sigprint();
    
    arma_conform_assert_same_size(P1, P2, "element-wise max()");
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : P1.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : P1.get_n_cols(); }
  arma_inline uword get_n_elem() const { return P1.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                const { const elem_type A = P1[i];        const elem_type B = P2[i];        return comparator(A,B) ? A : B; }
  arma_inline elem_type at         (const uword r, const uword c) const { const elem_type A = P1.at(r,c);   const elem_type B = P2.at(r,c);   return comparator(A,B) ? A : B; }
  arma_inline elem_type at_alt     (const uword i)                const { const elem_type A = P1.at_alt(i); const elem_type B = P2.at_alt(i); return comparator(A,B) ? A : B; }
  
  arma_inline         ea_type         get_ea() const { return *this; }
  arma_inline aligned_ea_type get_aligned_ea() const { return *this; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (P1.is_alias(X) || P2.is_alias(X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return (P1.has_overlap(X) || P2.has_overlap(X)); }
  
  arma_inline bool is_aligned() const { return (P1.is_aligned() && P2.is_aligned()); }
  };



template<typename out_eT, typename T1, typename op_type>
struct Proxy< mtOp<out_eT, T1, op_type> >
  {
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  typedef          const Mat<out_eT>&           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = mtOp<out_eT, T1, op_type>::is_row;
  static constexpr bool is_col  = mtOp<out_eT, T1, op_type>::is_col;
  static constexpr bool is_xvec = mtOp<out_eT, T1, op_type>::is_xvec;
  
  const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtOp<out_eT, T1, op_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r,c);   }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename out_eT, typename T1, typename T2, typename glue_type>
struct Proxy< mtGlue<out_eT, T1, T2, glue_type> >
  {
  typedef          out_eT                       elem_type;
  typedef typename get_pod_type<out_eT>::result pod_type;
  typedef          Mat<out_eT>                  stored_type;
  typedef          const elem_type*             ea_type;
  typedef          const Mat<out_eT>&           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = mtGlue<out_eT, T1, T2, glue_type>::is_row;
  static constexpr bool is_col  = mtGlue<out_eT, T1, T2, glue_type>::is_col;
  static constexpr bool is_xvec = mtGlue<out_eT, T1, T2, glue_type>::is_xvec;
  
  const Mat<out_eT> Q;
  
  inline explicit Proxy(const mtGlue<out_eT, T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r,c);   }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename op_type>
struct Proxy< CubeToMatOp<T1, op_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = CubeToMatOp<T1, op_type>::is_row;
  static constexpr bool is_col  = CubeToMatOp<T1, op_type>::is_col;
  static constexpr bool is_xvec = CubeToMatOp<T1, op_type>::is_xvec;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const CubeToMatOp<T1, op_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1>
struct Proxy< CubeToMatOp<T1, op_vectorise_cube_col> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const unwrap_cube<T1> U;
  const Mat<elem_type>  Q;
  
  inline explicit Proxy(const CubeToMatOp<T1, op_vectorise_cube_col>& A)
    : U(A.m)
    , Q(const_cast<elem_type*>(U.M.memptr()), U.M.n_elem, 1, false, true)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename op_type>
struct Proxy< SpToDOp<T1, op_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = SpToDOp<T1, op_type>::is_row;
  static constexpr bool is_col  = SpToDOp<T1, op_type>::is_col;
  static constexpr bool is_xvec = SpToDOp<T1, op_type>::is_xvec;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const SpToDOp<T1, op_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< SpToDOp<SpMat<eT>, op_sp_nonzeros> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const SpMat<elem_type>& R;
  const Mat<elem_type>    Q;
  
  inline explicit Proxy(const SpToDOp<SpMat<eT>, op_sp_nonzeros>& A)
    : R(A.m)
    , Q(const_cast<elem_type*>(R.values), R.n_nonzero, 1, false, true)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, typename T2, typename glue_type>
struct Proxy< SpToDGlue<T1, T2, glue_type> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = SpToDGlue<T1, T2, glue_type>::is_row;
  static constexpr bool is_col  = SpToDGlue<T1, T2, glue_type>::is_col;
  static constexpr bool is_xvec = SpToDGlue<T1, T2, glue_type>::is_xvec;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const SpToDGlue<T1, T2, glue_type>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }

  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< subview<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview<eT>                              stored_type;
  typedef const subview<eT>&                       ea_type;
  typedef const subview<eT>&                       aligned_ea_type;
  
  static constexpr bool use_at      = true;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const subview<eT>& Q;
  
  inline explicit Proxy(const subview<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];       }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c); }
  arma_inline elem_type at_alt     (const uword i)                const { return Q[i];       }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&(Q.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Q.check_overlap(X); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT>
struct Proxy< subview_col<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_col<eT>                          stored_type;
  typedef const eT*                                ea_type;
  typedef const subview_col<eT>&                   aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_col<eT>& Q;
  
  inline explicit Proxy(const subview_col<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.colmem; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;        }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&(Q.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Q.check_overlap(X); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.colmem); }
  };



template<typename eT>
struct Proxy< subview_cols<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Mat<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const subview_cols<eT>& sv;
  const Mat<eT>           Q;
  
  inline explicit Proxy(const subview_cols<eT>& A)
    : sv(A)
    , Q ( const_cast<eT*>( A.colptr(0) ), A.n_rows, A.n_cols, false, false )
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r,c);   }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&(sv.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return sv.check_overlap(X); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< subview_row<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_row<eT>                          stored_type;
  typedef const subview_row<eT>&                   ea_type;
  typedef const subview_row<eT>&                   aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const subview_row<eT>& Q;
  
  inline explicit Proxy(const subview_row<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  constexpr   uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword, const uword c) const { return Q[c]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&(Q.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Q.check_overlap(X); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT, typename T1>
struct Proxy< subview_elem1<eT,T1> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_elem1<eT,T1>                     stored_type;
  typedef const Proxy< subview_elem1<eT,T1> >&     ea_type;
  typedef const Proxy< subview_elem1<eT,T1> >&     aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_elem1<eT,T1>& Q;
  const Proxy<T1>             R;
  
  inline explicit Proxy(const subview_elem1<eT,T1>& A)
    : Q(A)
    , R(A.a.get_ref())
    {
    arma_debug_sigprint();
    
    const bool R_is_vec   = ((R.get_n_rows() == 1) || (R.get_n_cols() == 1));
    const bool R_is_empty = (R.get_n_elem() == 0);
    
    arma_conform_check( ((R_is_vec == false) && (R_is_empty == false)), "Mat::elem(): given object must be a vector" );
    }
  
  arma_inline uword get_n_rows() const { return R.get_n_elem(); }
  constexpr   uword get_n_cols() const { return 1;              }
  arma_inline uword get_n_elem() const { return R.get_n_elem(); }
  
  arma_inline elem_type operator[] (const uword i)              const { const uword ii = (Proxy<T1>::use_at) ? R.at(i,0) : R[i]; arma_conform_check_bounds( (ii >= Q.m.n_elem), "Mat::elem(): index out of bounds" ); return Q.m[ii]; }
  arma_inline elem_type at         (const uword r, const uword) const { const uword ii = (Proxy<T1>::use_at) ? R.at(r,0) : R[r]; arma_conform_check_bounds( (ii >= Q.m.n_elem), "Mat::elem(): index out of bounds" ); return Q.m[ii]; }
  arma_inline elem_type at_alt     (const uword i)              const { const uword ii = (Proxy<T1>::use_at) ? R.at(i,0) : R[i]; arma_conform_check_bounds( (ii >= Q.m.n_elem), "Mat::elem(): index out of bounds" ); return Q.m[ii]; }
  
  arma_inline         ea_type         get_ea() const { return (*this); }
  arma_inline aligned_ea_type get_aligned_ea() const { return (*this); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return ( (void_ptr(&X) == void_ptr(&(Q.m))) || (R.is_alias(X)) ); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT, typename T1, typename T2>
struct Proxy< subview_elem2<eT,T1,T2> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Mat<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const Mat<eT> Q;
  
  inline explicit Proxy(const subview_elem2<eT,T1,T2>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< diagview<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef diagview<eT>                             stored_type;
  typedef const diagview<eT>&                      ea_type;
  typedef const diagview<eT>&                      aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const diagview<eT>& Q;
  
  inline explicit Proxy(const diagview<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];       }
  arma_inline elem_type at         (const uword r, const uword) const { return Q.at(r, 0); }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i];       }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (is_same_type<eT,eT2>::yes) && (void_ptr(&(Q.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T1>
struct Proxy_diagvec_mat
  {
  inline Proxy_diagvec_mat(const T1&) {}
  };



template<typename T1>
struct Proxy_diagvec_mat< Op<T1, op_diagvec> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef diagview<elem_type>                      stored_type;
  typedef const diagview<elem_type>&               ea_type;
  typedef const diagview<elem_type>&               aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const Mat<elem_type>&     R;
  const diagview<elem_type> Q;
  
  inline explicit Proxy_diagvec_mat(const Op<T1, op_diagvec>& A)
    : R(A.m), Q( R.diag() )
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];       }
  arma_inline elem_type at         (const uword r, const uword) const { return Q.at(r, 0); }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i];       }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&R) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T1>
struct Proxy_diagvec_expr
  {
  inline Proxy_diagvec_expr(const T1&) {}
  };



template<typename T1>
struct Proxy_diagvec_expr< Op<T1, op_diagvec> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy_diagvec_expr(const Op<T1, op_diagvec>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q.at(r, 0);  }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, bool condition>
struct Proxy_diagvec_redirect {};

template<typename T1>
struct Proxy_diagvec_redirect< Op<T1, op_diagvec>, true > { typedef Proxy_diagvec_mat < Op<T1, op_diagvec> > result; };

template<typename T1>
struct Proxy_diagvec_redirect< Op<T1, op_diagvec>, false> { typedef Proxy_diagvec_expr< Op<T1, op_diagvec> > result; };



template<typename T1>
struct Proxy< Op<T1, op_diagvec> >
  : public Proxy_diagvec_redirect< Op<T1, op_diagvec>, is_Mat<T1>::value >::result
  {
  typedef typename Proxy_diagvec_redirect< Op<T1, op_diagvec>, is_Mat<T1>::value >::result Proxy_diagvec;
  
  inline explicit Proxy(const Op<T1, op_diagvec>& A)
    : Proxy_diagvec(A)
    {
    arma_debug_sigprint();
    }
  };



template<typename T1>
struct Proxy< Op<T1, op_diagvec2> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const Mat<elem_type> Q;
  
  inline explicit Proxy(const Op<T1, op_diagvec2>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q.at(r, 0);  }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1>
struct Proxy_xtrans_default
  {
  inline Proxy_xtrans_default(const T1&) {}
  };



template<typename T1>
struct Proxy_xtrans_default< Op<T1, op_htrans> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef xtrans_mat<elem_type,true>               stored_type;
  typedef const xtrans_mat<elem_type,true>&        ea_type;
  typedef const xtrans_mat<elem_type,true>&        aligned_ea_type;
  
  static constexpr bool use_at      = true;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const unwrap<T1>                 U;
  const xtrans_mat<elem_type,true> Q;
  
  inline explicit Proxy_xtrans_default(const Op<T1, op_htrans>& A)
    : U(A.m)
    , Q(U.M)
    {
    arma_debug_sigprint();
    }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return void_ptr(&(U.M)) == void_ptr(&X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T1>
struct Proxy_xtrans_default< Op<T1, op_strans> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef xtrans_mat<elem_type,false>              stored_type;
  typedef const xtrans_mat<elem_type,false>&       ea_type;
  typedef const xtrans_mat<elem_type,false>&       aligned_ea_type;
  
  static constexpr bool use_at      = true;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const unwrap<T1>                  U;
  const xtrans_mat<elem_type,false> Q;
  
  inline explicit Proxy_xtrans_default(const Op<T1, op_strans>& A)
    : U(A.m)
    , Q(U.M)
    {
    arma_debug_sigprint();
    }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return void_ptr(&(U.M)) == void_ptr(&X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T1>
struct Proxy_xtrans_vector
  {
  inline Proxy_xtrans_vector(const T1&) {}
  };



template<typename T1>
struct Proxy_xtrans_vector< Op<T1, op_htrans> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = quasi_unwrap<T1>::has_subview;
  
  // NOTE: the Op class takes care of swapping row and col for op_htrans
  static constexpr bool is_row  = Op<T1, op_htrans>::is_row;
  static constexpr bool is_col  = Op<T1, op_htrans>::is_col;
  static constexpr bool is_xvec = Op<T1, op_htrans>::is_xvec;
  
  const quasi_unwrap<T1> U; // avoid copy if T1 is a Row, Col or subview_col
  const Mat<elem_type>   Q;
  
  inline Proxy_xtrans_vector(const Op<T1, op_htrans>& A)
    : U(A.m)
    , Q(const_cast<elem_type*>(U.M.memptr()), U.M.n_cols, U.M.n_rows, false, false)
    {
    arma_debug_sigprint();
    }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return U.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1>
struct Proxy_xtrans_vector< Op<T1, op_strans> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = quasi_unwrap<T1>::has_subview;
  
  // NOTE: the Op class takes care of swapping row and col for op_strans
  static constexpr bool is_row  = Op<T1, op_strans>::is_row;
  static constexpr bool is_col  = Op<T1, op_strans>::is_col;
  static constexpr bool is_xvec = Op<T1, op_strans>::is_xvec;
  
  const quasi_unwrap<T1> U; // avoid copy if T1 is a Row, Col or subview_col
  const Mat<elem_type>   Q;
  
  inline Proxy_xtrans_vector(const Op<T1, op_strans>& A)
    : U(A.m)
    , Q(const_cast<elem_type*>(U.M.memptr()), U.M.n_cols, U.M.n_rows, false, false)
    {
    arma_debug_sigprint();
    }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return U.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1, bool condition>
struct Proxy_xtrans_redirect {};

template<typename T1>
struct Proxy_xtrans_redirect<T1, false> { typedef Proxy_xtrans_default<T1> result; };

template<typename T1>
struct Proxy_xtrans_redirect<T1, true>  { typedef Proxy_xtrans_vector<T1>  result; };



template<typename T1>
struct Proxy< Op<T1, op_htrans> >
  : public
    Proxy_xtrans_redirect
      <
      Op<T1, op_htrans>,
      ((is_cx<typename T1::elem_type>::no) && ((Op<T1, op_htrans>::is_row) || (Op<T1, op_htrans>::is_col)) )
      >::result
  {
  typedef
  typename
  Proxy_xtrans_redirect
    <
    Op<T1, op_htrans>,
    ((is_cx<typename T1::elem_type>::no) && ((Op<T1, op_htrans>::is_row) || (Op<T1, op_htrans>::is_col)) )
    >::result
  Proxy_xtrans;
  
  typedef typename Proxy_xtrans::elem_type       elem_type;
  typedef typename Proxy_xtrans::pod_type        pod_type;
  typedef typename Proxy_xtrans::stored_type     stored_type;
  typedef typename Proxy_xtrans::ea_type         ea_type;
  typedef typename Proxy_xtrans::aligned_ea_type aligned_ea_type;
  
  static constexpr bool use_at      = Proxy_xtrans::use_at;
  static constexpr bool use_mp      = Proxy_xtrans::use_mp;
  static constexpr bool has_subview = Proxy_xtrans::has_subview;
  
  static constexpr bool is_row  = Proxy_xtrans::is_row;
  static constexpr bool is_col  = Proxy_xtrans::is_col;
  static constexpr bool is_xvec = Proxy_xtrans::is_xvec;
  
  using Proxy_xtrans::Q;
  
  inline explicit Proxy(const Op<T1, op_htrans>& A)
    : Proxy_xtrans(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Proxy_xtrans::get_ea();         }
  arma_inline aligned_ea_type get_aligned_ea() const { return Proxy_xtrans::get_aligned_ea(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Proxy_xtrans::is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Proxy_xtrans::has_overlap(X); }
  
  arma_inline bool is_aligned() const { return Proxy_xtrans::is_aligned(); }
  };



template<typename T1>
struct Proxy< Op<T1, op_strans> >
  : public
    Proxy_xtrans_redirect
      <
      Op<T1, op_strans>,
      ( (Op<T1, op_strans>::is_row) || (Op<T1, op_strans>::is_col) )
      >::result
  {
  typedef
  typename
  Proxy_xtrans_redirect
    <
    Op<T1, op_strans>,
    ( (Op<T1, op_strans>::is_row) || (Op<T1, op_strans>::is_col) )
    >::result
  Proxy_xtrans;
  
  typedef typename Proxy_xtrans::elem_type       elem_type;
  typedef typename Proxy_xtrans::pod_type        pod_type;
  typedef typename Proxy_xtrans::stored_type     stored_type;
  typedef typename Proxy_xtrans::ea_type         ea_type;
  typedef typename Proxy_xtrans::aligned_ea_type aligned_ea_type;
  
  static constexpr bool use_at      = Proxy_xtrans::use_at;
  static constexpr bool use_mp      = Proxy_xtrans::use_mp;
  static constexpr bool has_subview = Proxy_xtrans::has_subview;
  
  static constexpr bool is_row  = Proxy_xtrans::is_row;
  static constexpr bool is_col  = Proxy_xtrans::is_col;
  static constexpr bool is_xvec = Proxy_xtrans::is_xvec;
  
  using Proxy_xtrans::Q;
  
  inline explicit Proxy(const Op<T1, op_strans>& A)
    : Proxy_xtrans(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.n_rows; }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem;              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Proxy_xtrans::get_ea();         }
  arma_inline aligned_ea_type get_aligned_ea() const { return Proxy_xtrans::get_aligned_ea(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Proxy_xtrans::is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Proxy_xtrans::has_overlap(X); }
  
  arma_inline bool is_aligned() const { return Proxy_xtrans::is_aligned(); }
  };



template<typename eT>
struct Proxy_subview_row_htrans_cx
  {
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  typedef subview_row_htrans<eT>            stored_type;
  typedef const subview_row_htrans<eT>&     ea_type;
  typedef const subview_row_htrans<eT>&     aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_row_htrans<eT> Q;
  
  inline explicit Proxy_subview_row_htrans_cx(const Op<subview_row<eT>, op_htrans>& A)
    : Q(A.m)
    {
    arma_debug_sigprint();
    }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.sv_row.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  };



template<typename eT>
struct Proxy_subview_row_htrans_non_cx
  {
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  typedef subview_row_strans<eT>            stored_type;
  typedef const subview_row_strans<eT>&     ea_type;
  typedef const subview_row_strans<eT>&     aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_row_strans<eT> Q;
  
  inline explicit Proxy_subview_row_htrans_non_cx(const Op<subview_row<eT>, op_htrans>& A)
    : Q(A.m)
    {
    arma_debug_sigprint();
    }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.sv_row.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  };



template<typename eT, bool condition>
struct Proxy_subview_row_htrans_redirect {};

template<typename eT>
struct Proxy_subview_row_htrans_redirect<eT, true>  { typedef Proxy_subview_row_htrans_cx<eT>      result; };

template<typename eT>
struct Proxy_subview_row_htrans_redirect<eT, false> { typedef Proxy_subview_row_htrans_non_cx<eT>  result; };



template<typename eT>
struct Proxy< Op<subview_row<eT>, op_htrans> >
  : public
    Proxy_subview_row_htrans_redirect
      <
      eT,
      is_cx<eT>::yes
      >::result
  {
  typedef
  typename
  Proxy_subview_row_htrans_redirect
      <
      eT,
      is_cx<eT>::yes
      >::result
  Proxy_sv_row_ht;
  
  typedef typename Proxy_sv_row_ht::elem_type   elem_type;
  typedef typename Proxy_sv_row_ht::pod_type    pod_type;
  typedef typename Proxy_sv_row_ht::stored_type stored_type;
  typedef typename Proxy_sv_row_ht::ea_type     ea_type;
  typedef typename Proxy_sv_row_ht::ea_type     aligned_ea_type;
  
  static constexpr bool use_at      = Proxy_sv_row_ht::use_at;
  static constexpr bool use_mp      = Proxy_sv_row_ht::use_mp;
  static constexpr bool has_subview = Proxy_sv_row_ht::has_subview;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  using Proxy_sv_row_ht::Q;
  
  inline explicit Proxy(const Op<subview_row<eT>, op_htrans>& A)
    : Proxy_sv_row_ht(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Proxy_sv_row_ht::is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return Proxy_sv_row_ht::has_overlap(X); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT>
struct Proxy< Op<subview_row<eT>, op_strans> >
  {
  typedef eT                                elem_type;
  typedef typename get_pod_type<eT>::result pod_type;
  typedef subview_row_strans<eT>            stored_type;
  typedef const subview_row_strans<eT>&     ea_type;
  typedef const subview_row_strans<eT>&     aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_row_strans<eT> Q;
  
  inline explicit Proxy(const Op<subview_row<eT>, op_strans>& A)
    : Q(A.m)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.sv_row.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T>
struct Proxy< Op< Row< std::complex<T> >, op_htrans> >
  {
  typedef typename std::complex<T>  eT;
  
  typedef typename std::complex<T>  elem_type;
  typedef T                         pod_type;
  typedef xvec_htrans<eT>           stored_type;
  typedef const xvec_htrans<eT>&    ea_type;
  typedef const xvec_htrans<eT>&    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const xvec_htrans<eT> Q;
  const Row<eT>&        src;
  
  inline explicit Proxy(const Op< Row< std::complex<T> >, op_htrans>& A)
    : Q  (A.m.memptr(), A.m.n_rows, A.m.n_cols)
    , src(A.m)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return void_ptr(&src) == void_ptr(&X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T>
struct Proxy< Op< Col< std::complex<T> >, op_htrans> >
  {
  typedef typename std::complex<T>  eT;

  typedef typename std::complex<T>  elem_type;
  typedef T                         pod_type;
  typedef xvec_htrans<eT>           stored_type;
  typedef const xvec_htrans<eT>&    ea_type;
  typedef const xvec_htrans<eT>&    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const xvec_htrans<eT> Q;
  const Col<eT>&        src;
  
  inline explicit Proxy(const Op< Col< std::complex<T> >, op_htrans>& A)
    : Q  (A.m.memptr(), A.m.n_rows, A.m.n_cols)
    , src(A.m)
    {
    arma_debug_sigprint();
    }
  
  constexpr   uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword, const uword c) const { return Q[c]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return void_ptr(&src) == void_ptr(&X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T>
struct Proxy< Op< subview_col< std::complex<T> >, op_htrans> >
  {
  typedef typename std::complex<T>  eT;
  
  typedef typename std::complex<T>  elem_type;
  typedef T                         pod_type;
  typedef xvec_htrans<eT>           stored_type;
  typedef const xvec_htrans<eT>&    ea_type;
  typedef const xvec_htrans<eT>&    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = true;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const xvec_htrans<eT>  Q;
  const subview_col<eT>& src;
  
  inline explicit Proxy(const Op< subview_col< std::complex<T> >, op_htrans>& A)
    : Q  (A.m.colptr(0), A.m.n_rows, A.m.n_cols)
    , src(A.m)
    {
    arma_debug_sigprint();
    }
  
  constexpr   uword get_n_rows() const { return 1;        }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword, const uword c) const { return Q[c]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return void_ptr(&src.m) == void_ptr(&X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename T1>
struct Proxy< Op<T1, op_htrans2> >
  {
  typedef typename T1::elem_type                           elem_type;
  typedef typename get_pod_type<elem_type>::result         pod_type;
  typedef       eOp< Op<T1, op_htrans>, eop_scalar_times>  stored_type;
  typedef const eOp< Op<T1, op_htrans>, eop_scalar_times>& ea_type;
  typedef const eOp< Op<T1, op_htrans>, eop_scalar_times>& aligned_ea_type;
  
  static constexpr bool use_at      = eOp< Op<T1, op_htrans>, eop_scalar_times>::use_at;
  static constexpr bool use_mp      = eOp< Op<T1, op_htrans>, eop_scalar_times>::use_mp;
  static constexpr bool has_subview = eOp< Op<T1, op_htrans>, eop_scalar_times>::has_subview;
  
  // NOTE: the Op class takes care of swapping row and col for op_htrans
  static constexpr bool is_row  = eOp< Op<T1, op_htrans>, eop_scalar_times>::is_row;
  static constexpr bool is_col  = eOp< Op<T1, op_htrans>, eop_scalar_times>::is_col;
  static constexpr bool is_xvec = eOp< Op<T1, op_htrans>, eop_scalar_times>::is_xvec;
  
  const      Op<T1, op_htrans>                     R;
  const eOp< Op<T1, op_htrans>, eop_scalar_times > Q;
  
  inline explicit Proxy(const Op<T1, op_htrans2>& A)
    : R(A.m)
    , Q(R, A.aux)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return is_row ? 1 : Q.get_n_rows(); }
  arma_inline uword get_n_cols() const { return is_col ? 1 : Q.get_n_cols(); }
  arma_inline uword get_n_elem() const { return Q.get_n_elem();              }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r, c);  }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return Q.P.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return Q.P.is_aligned(); }
  };



template<typename eT>
struct Proxy< subview_row_strans<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_row_strans<eT>                   stored_type;
  typedef const subview_row_strans<eT>&            ea_type;
  typedef const subview_row_strans<eT>&            aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_row_strans<eT>& Q;
  
  inline explicit Proxy(const subview_row_strans<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.sv_row.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT>
struct Proxy< subview_row_htrans<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef subview_row_htrans<eT>                   stored_type;
  typedef const subview_row_htrans<eT>&            ea_type;
  typedef const subview_row_htrans<eT>&            aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const subview_row_htrans<eT>& Q;
  
  inline explicit Proxy(const subview_row_htrans<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i]; }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r]; }
  arma_inline elem_type at_alt     (const uword i)              const { return Q[i]; }
  
  arma_inline         ea_type         get_ea() const { return Q; }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q; }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return (void_ptr(&(Q.sv_row.m)) == void_ptr(&X)); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  constexpr bool is_aligned() const { return false; }
  };



template<typename eT, bool do_conj>
struct Proxy< xtrans_mat<eT, do_conj> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Mat<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = false;
  
  const Mat<eT> Q;
  
  inline explicit Proxy(const xtrans_mat<eT, do_conj>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r,c);   }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename eT>
struct Proxy< xvec_htrans<eT> >
  {
  typedef eT                                       elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<eT>                                  stored_type;
  typedef const eT*                                ea_type;
  typedef const Mat<eT>&                           aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = false;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = false;
  static constexpr bool is_xvec = true;
  
  const Mat<eT> Q;
  
  inline explicit Proxy(const xvec_htrans<eT>& A)
    : Q(A)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  arma_inline uword get_n_cols() const { return Q.n_cols; }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)                const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword c) const { return Q.at(r,c);   }
  arma_inline elem_type at_alt     (const uword i)                const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  constexpr bool is_alias(const Mat<eT2>&) const { return false; }
  
  template<typename eT2>
  constexpr bool has_overlap(const subview<eT2>&) const { return false; }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1>
struct Proxy_vectorise_col_mat
  {
  inline Proxy_vectorise_col_mat(const T1&) {}
  };



template<typename T1>
struct Proxy_vectorise_col_mat< Op<T1, op_vectorise_col> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef Mat<elem_type>                           stored_type;
  typedef const elem_type*                         ea_type;
  typedef const Mat<elem_type>&                    aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = false;
  static constexpr bool has_subview = true;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const unwrap<T1>     U;
  const Mat<elem_type> Q;
  
  inline explicit Proxy_vectorise_col_mat(const Op<T1, op_vectorise_col>& A)
    : U(A.m)
    , Q(const_cast<elem_type*>(U.M.memptr()), U.M.n_elem, 1, false, false)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return Q.n_rows; }
  constexpr   uword get_n_cols() const { return 1;        }
  arma_inline uword get_n_elem() const { return Q.n_elem; }
  
  arma_inline elem_type operator[] (const uword i)              const { return Q[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return Q[r];        }
  arma_inline elem_type at_alt     (const uword i)              const { return Q.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return Q.memptr(); }
  arma_inline aligned_ea_type get_aligned_ea() const { return Q;          }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return ( void_ptr(&X) == void_ptr(&(U.M)) ); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return memory::is_aligned(Q.memptr()); }
  };



template<typename T1>
struct Proxy_vectorise_col_expr
  {
  inline Proxy_vectorise_col_expr(const T1&) {}
  };



template<typename T1>
struct Proxy_vectorise_col_expr< Op<T1, op_vectorise_col> >
  {
  typedef typename T1::elem_type                   elem_type;
  typedef typename get_pod_type<elem_type>::result pod_type;
  typedef          Op<T1, op_vectorise_col>        stored_type;
  typedef typename Proxy<T1>::ea_type              ea_type;
  typedef typename Proxy<T1>::aligned_ea_type      aligned_ea_type;
  
  static constexpr bool use_at      = false;
  static constexpr bool use_mp      = Proxy<T1>::use_mp;
  static constexpr bool has_subview = Proxy<T1>::has_subview;
  
  static constexpr bool is_row  = false;
  static constexpr bool is_col  = true;
  static constexpr bool is_xvec = false;
  
  const Op<T1, op_vectorise_col>& Q;
  const Proxy<T1>                 R;
  
  inline explicit Proxy_vectorise_col_expr(const Op<T1, op_vectorise_col>& A)
    : Q(A)
    , R(A.m)
    {
    arma_debug_sigprint();
    }
  
  arma_inline uword get_n_rows() const { return R.get_n_elem(); }
  constexpr   uword get_n_cols() const { return 1;              }
  arma_inline uword get_n_elem() const { return R.get_n_elem(); }
  
  arma_inline elem_type operator[] (const uword i)              const { return R[i];        }
  arma_inline elem_type at         (const uword r, const uword) const { return R.at(r, 0);  }
  arma_inline elem_type at_alt     (const uword i)              const { return R.at_alt(i); }
  
  arma_inline         ea_type         get_ea() const { return R.get_ea();         }
  arma_inline aligned_ea_type get_aligned_ea() const { return R.get_aligned_ea(); }
  
  template<typename eT2>
  arma_inline bool is_alias(const Mat<eT2>& X) const { return R.is_alias(X); }
  
  template<typename eT2>
  arma_inline bool has_overlap(const subview<eT2>& X) const { return is_alias(X.m); }
  
  arma_inline bool is_aligned() const { return R.is_aligned(); }
  };



template<typename T1, bool condition>
struct Proxy_vectorise_col_redirect {};

template<typename T1>
struct Proxy_vectorise_col_redirect< Op<T1, op_vectorise_col>, true > { typedef Proxy_vectorise_col_mat < Op<T1, op_vectorise_col> > result; };

template<typename T1>
struct Proxy_vectorise_col_redirect< Op<T1, op_vectorise_col>, false> { typedef Proxy_vectorise_col_expr< Op<T1, op_vectorise_col> > result; };



template<typename T1>
struct Proxy< Op<T1, op_vectorise_col> >
  : public Proxy_vectorise_col_redirect< Op<T1, op_vectorise_col>, (Proxy<T1>::use_at) >::result
  {
  typedef typename Proxy_vectorise_col_redirect< Op<T1, op_vectorise_col>, (Proxy<T1>::use_at) >::result Proxy_vectorise_col;
  
  inline explicit Proxy(const Op<T1, op_vectorise_col>& A)
    : Proxy_vectorise_col(A)
    {
    arma_debug_sigprint();
    }
  };



//! @}
