// Determinant from an existing CPU SPIKE-SMP factorization.
//
// Mirrors the argument list of the original Fortran factorization:
//
//   call dspike_gbtrf(spm, n, kl, ku, A, lda, work, info)
//   call dspike_gbtrfp(spm, n, kl, ku, A, lda, work, ipiv, info)
//
// This header consumes the *outputs* of those calls (A, work, ipiv, info)
// and returns det(A).  It performs no factorization of A; the only extra
// work is an LU of the small interface system (2*klu*(nbpart-1) unknowns),
// which is rebuilt from the level-1 spikes stored in work.
//
// Identity used (single-level SPIKE, valid for any number of partitions):
//   det(A) = prod_p det(A_pp) * det(R),
//   R = block tridiagonal interface matrix built from the saved spikes.
// Each LU factor contributes (-1)^swaps * prod U_kk.
//
// Layout assumptions (from Xspike_smp.f90 / spike_smp_utilities.f90):
//   spm(20)=nbpart, spm(21)=nblevel, spm(23)=# two-thread partitions,
//   spm(3)==1 selects the pivoted path.  Fortran indices are 1-based; all
//   arrays below are passed as plain C++ pointers, index 0 == Fortran 1.
//   A      : band storage, lda = kl+ku+1 (unpivoted) or
//            klu+kl+ku+1 (pivoted), each partition factored in place.
//   work   : [rV][rW][red], each rV/rW block has 2*klu*klu*nblevel entries
//            per partition, stored Fortran-column-major.
//   ipiv   : LAPACK-style 1-based pivot indices, one entry per global row,
//            required by the pivoted path.
//   info   : 1 = diagonal boosting was used (result approximate),
//            2 = illegal value (determinant reported as zero).
//
// Note: det() for near-singular matrices is ill-conditioned; the sign is
// accumulated independently of the magnitude and is more robust.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace gspike_cpu {

    namespace detail {

        // Mantissa in [0.5,1) plus a base-2 exponent; combining never
        // overflows or underflows before the final ldexp.
        template<class T>
        struct scaled_det {
            T mant = T(1);
            long long exp = 0;
            int sign = 1;
            bool zero = false;

            void mul(const T v) {
                if(v == T(0)) {
                    zero = true;
                    return;
                }
                int e = 0;
                const T m = std::frexp(std::fabs(v), &e);
                mant *= m;
                exp += e;
                if(mant < T(0.5)) {
                    mant *= T(2);
                    --exp;
                }
                if(v < T(0)) sign = -sign;
            }

            void mul(const scaled_det& o) {
                mant *= o.mant;
                exp += o.exp;
                sign *= o.sign;
                zero = zero || o.zero;
                if(mant < T(0.5)) {
                    mant *= T(2);
                    --exp;
                }
            }

            [[nodiscard]] T value() const {
                if(zero) return T(0);
                const long long e = exp;
                const int ec = (int)(e > 2000 ? 2000 : (e < -2000 ? -2000 : e));
                const T v = std::ldexp(mant, ec);
                return sign < 0 ? -v : v;
            }
        };

        // Replicates spikerl_calc_size_partitions (spike_smp_utilities.f90).
        inline void partition_starts(const int* spm, int n, std::vector<int>& ajmin) {
            const int nbpart = spm[19]; // spm(20)
            const int y = spm[22];      // spm(23), two-thread partitions
            const int x = nbpart - (2 + y);
            const double r = (double)(spm[1] == 1 ? 10 : spm[3]); // spm(2)/spm(4)
            const double r_large = 0.1 * r, r_small = 0.2 * r;
            const double denom = 2.0 + x / r_small + y / r_large;

            int size_large = (int)(n / (r_large * denom));
            size_large += size_large % 2;
            const int size_small = (int)(n / (r_small * denom));
            const int size_huge = (n - (x * size_small + y * size_large)) / 2;

            ajmin.assign((size_t)nbpart + 1, 0);
            ajmin[1] = size_huge;
            for(int i = 2; i <= 1 + y; ++i) ajmin[i] = ajmin[i - 1] + size_large;
            for(int i = 2 + y; i <= nbpart - 1; ++i) ajmin[i] = ajmin[i - 1] + size_small;
            ajmin[nbpart] = n;
        }

        // Dense LU of the interface system; returns its scaled determinant.
        template<class T>
        scaled_det<T> interface_det(const T* M0, int m) {
            scaled_det<T> r;
            if(m == 0) return r;
            std::vector<T> M(M0, M0 + (size_t)m * m);
            const auto at = [&](const int i, const int j) -> T& { return M[(size_t)j * m + i]; };
            for(int j = 0; j < m; ++j) {
                int pv = j;
                for(int i = j + 1; i < m; ++i)
                    if(std::fabs(at(i, j)) > std::fabs(at(pv, j))) pv = i;
                if(pv != j) {
                    for(int c = 0; c < m; ++c) std::swap(at(j, c), at(pv, c));
                    r.sign = -r.sign;
                }
                const T d = at(j, j);
                if(d == T(0)) {
                    r.zero = true;
                    return r;
                }
                r.mul(d);
                for(int i = j + 1; i < m; ++i) {
                    at(i, j) /= d;
                    for(int c = j + 1; c < m; ++c) at(i, c) -= at(i, j) * at(j, c);
                }
            }
            return r;
        }

    } // namespace detail

    // Return det(A) from the factorized A / work / ipiv produced by
    // spike_gbtrf or spike_gbtrfp.  If sign is non-null it receives -1/0/+1.
    template<class T>
    T spike_det(const int* spm, int n, int kl, int ku,
                const T* A, int lda, const T* work,
                int info, const int* ipiv = nullptr, int* sign = nullptr) {
        const int klu0 = std::max(kl, ku);
        const bool needs_work = spm && spm[19] > 1 && klu0 > 0;
        if(!spm || !A || (needs_work && !work) || n < 1 || kl < 0 || ku < 0 || lda < kl + ku + 1) {
            if(sign) *sign = 0;
            return T(0);
        }
        // Two-thread partitions (spm(23)>0) factor the two halves separately
        // and keep the determinant of their internal interface in a temporary
        // array that is not returned in work.  That contribution cannot be
        // recovered here, so report 0 instead of a wrong value.
        if(spm[22] > 0) {
            if(sign) *sign = 0;
            return T(0);
        }
        if(info == 2) { // illegal value reported by the factorization
            if(sign) *sign = 0;
            return T(0);
        }

        const int klu = std::max(kl, ku);
        const int nbpart = spm[19]; // spm(20)
        const int nblevel = spm[20]; // spm(21)
        const int ntwo = spm[22]; // spm(23)
        const bool pivot = spm[2] == 1; // spm(3)

        std::vector<int> ajmin;
        detail::partition_starts(spm, n, ajmin);

        detail::scaled_det<T> total;

        // ---- determinants of the factored diagonal blocks -------------------
        // For GBALU/GBAUL (no pivoting) the stored diagonal is at row ku+1
        // (C index ku).  For GBTRF it is at C index klu+ku.  GBTRFUL swaps
        // band columns with a -1 stride, i.e. it factorizes J*A*J; det is
        // unchanged, but its diagonal sits at C index klu+kl.
        const int diag_f = pivot ? klu + ku : ku;
        const int diag_b = pivot ? klu + kl : ku;
        for(int p = 0; p < nbpart; ++p) {
            const int start = ajmin[p], size = ajmin[p + 1] - ajmin[p];
            const bool two_thread = (p >= 1 && p < 1 + ntwo && size > 1);
            const bool reversed = pivot && (two_thread || p == nbpart - 1);

            auto block = [&](const int from, const int count, const int diag_off) {
                for(int j = 0; j < count; ++j)
                    total.mul(A[(size_t)(from + j) * lda + diag_off]);
                if(pivot && ipiv) {
                    for(int j = 0; j < count; ++j)
                        if(ipiv[from + j] != j + 1) total.sign = -total.sign;
                }
            };

            if(two_thread) {
                const int n1 = size / 2;
                block(start, n1, diag_f);
                block(start + n1, size - n1, diag_b);
            }
            else {
                block(start, size, reversed ? diag_b : diag_f);
            }
        }

        // ---- rebuild the level-0 interface system from the saved spikes -----
        // rV/rW: (2*klu, klu, nblevel, nbpart), Fortran column major.
        // Level 1 is untouched by spike_prep_recn, so the original spikes
        // are still available there after factorization.
        if(nbpart > 1 && klu > 0) {
            const size_t plane = (size_t)2 * klu * klu * (size_t)nblevel;
            const size_t rw0 = plane * (size_t)nbpart;
            const int m = 2 * klu * (nbpart - 1);
            std::vector<T> R((size_t)m * m, T(0));
            const auto r = [&](const int i, const int j) -> T& { return R[(size_t)j * m + i]; };
            // rV(i,j,l,p) / rW(i,j,l,p): 0-based i,j,l,p.
            const auto v = [&](const int i, const int j, const int p) -> T {
                return work[(size_t)p * plane + (size_t)j * 2 * klu + i];
            };
            const auto w = [&](const int i, const int j, const int p) -> T {
                return work[rw0 + (size_t)p * plane + (size_t)j * 2 * klu + i];
            };

            for(int i = 0; i < nbpart - 1; ++i) {
                const int I0 = 2 * klu * i;
                for(int a = 0; a < klu; ++a) {
                    r(I0 + a, I0 + a) = T(1);
                    r(I0 + klu + a, I0 + klu + a) = T(1);
                    for(int b = 0; b < klu; ++b) {
                        r(I0 + a, I0 + klu + b) = v(klu + a, b, i);       // Vp_i^bot
                        r(I0 + klu + a, I0 + b) = w(a, b, i + 1);          // Vm_{i+1}^top
                        if(i > 0)
                            r(I0 + a, I0 - 2 * klu + b) = w(klu + a, b, i); // Vm_i^bot
                        if(i + 1 < nbpart - 1)
                            r(I0 + klu + a, I0 + 3 * klu + b) = v(a, b, i + 1); // Vp_{i+1}^top
                    }
                }
            }
            total.mul(detail::interface_det<T>(R.data(), m));
        }

        total.sign = total.zero ? 0 : total.sign;
        if(sign) *sign = total.sign;
        return total.value();
    }

} // namespace gspike_cpu
